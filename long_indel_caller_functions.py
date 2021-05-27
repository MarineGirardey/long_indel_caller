#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr 14 08:49:55 2021
@author: Marine Girardey

This program gives the long INDEL in a bam file after NGS analysis.
"""

# Import modules
import pysam
import os
import pyfastx
import pandas as pd
from collections import Counter
from difflib import SequenceMatcher
import heapq
from operator import itemgetter


def check_index(bam_file):
    """
    Check if a bam file have an index file if not create one

    Returns
    ----------
    Create a indexed bam file or do nothing
    """

    pysam.index(bam_file)


# FILTERING PHASE
def filter_bam(bam_file, span, position, mapping_quality, filtered_bam):
    """
    Filter a bam file on a specific position too keep reads that cover a region
    and keep reads with CIGAR containing 'S', 'D' or 'I' or with low mapping quality

    Parameters
    ----------
    bam_file : complete path to the original bam file
    position : specified position to filter the bam file
    span : span around the specified position
    mapping_quality : mapping quality
    filtered_bam : path of the output bam

    Returns
    ----------
    Create a filtered filtered bam file
    """
    # Extract position to map according to the input position
    position = position.split(':')
    nuc_pos = int(position[1])
    # Chromosome, position start and stop
    chr_pos = position[0]
    nuc_start_pos = nuc_pos - span
    nuc_end_pos = nuc_pos + span

    # Letters corresponding respectively to Deletion, Insertion and Split reads
    del_ins_split = ['D', 'I', 'S']

    # Open the file and read in it in binary
    sam_file = pysam.AlignmentFile(bam_file, "rb")

    # Open the filtered bam file to write inside
    filtered_bam = pysam.AlignmentFile(filtered_bam, 'wbu', template=sam_file)

    # For lines (reads) in the bam file, sort depending on the specified positions
    for read in sam_file.fetch(chr_pos, nuc_start_pos, nuc_end_pos):

        # Grep the CIGAR string
        cigar = read.cigarstring
        read_mapping_quality = read.mapping_quality

        # If for NoneType cigar (* = unmapped reads)
        if cigar is None:
            pass

        else:
            # If one or more letter from the del_ins_split list is in cigar string
            if any(x in del_ins_split for x in cigar) or read_mapping_quality < mapping_quality:
                # Then write the associated read in the new bam file
                filtered_bam.write(read)

    # Close original bam file and filtered bam file
    filtered_bam.close()
    sam_file.close()


# CIGAR STATISTICS PHASE
def stats_cigar_indel(filtered_bam):
    """
    Statistics on indel positions in cigarstrings

    Parameters
    ----------
    filtered_bam : the complete path to the filtered bam

    Returns
    ----------
    Dictionary with position and read id associated
    """

    # Each CIGAR string is associated to a value (based on pysam doc)
    cigar_translate = {'M': 0, 'D': 2, 'I': 1, 'S': 4}
    # Final dictionary
    read_id_pos_dict = {}
    res_dict = {}

    # Open the file and read in it in binary
    sam_file = pysam.AlignmentFile(filtered_bam, "rb")

    # For lines (reads) in the bam file
    for read in sam_file.fetch():
        # Get the read name/id
        read_id = str(read.query_name)
        flag = read.flag
        read_id_flag = read_id, flag

        # Get the chromosome name
        chromosome = read.reference_name

        # Get the cigar of the read as a tuple with a value associated to every possible cigar letter.
        # Ex cigar tuple : (0, 104), (1, 1), (0, 5)
        cigar_tuples = read.cigartuples

        # Get the start position of the variation in the chromosome
        start_pos_read = read.reference_start

        # Create a new dictionary which exclude the 'M' cigar letter only to detect indel in the cigar
        cigar_translate_except_match = {key: val for key, val in cigar_translate.items() if key != 'M'}

        # Browse each key/value of the previous dictionary
        for cigar_event, cigar_event_value in cigar_translate_except_match.items():

            # Specific case to excluse : if cigar is a Nontype (*), don't pay attention
            if cigar_tuples is None:
                pass

            else:
                res_dict = create_pos_and_read_dict(cigar_event_value, cigar_tuples, chromosome, start_pos_read,
                                                    read_id_flag, read_id_pos_dict)

    sam_file.close()
    return res_dict


def create_pos_and_read_dict(cigar_event_value, cigar_tuples, chromosome, start_pos_read,
                             read_id_flag, read_id_pos_dict):
    """
    Create two dictionary, one with indel start position & number of reads concerned and one with indel start position
    and associated reads_id.

    Parameters
    ----------
    cigar_event_value : dictionary with cigar event and associated pysam value
    cigar_tuples : list of tuples which associated value from  a cigar event
    chromosome : number of the chromosome
    start_pos_read : start position of the read
    read_id_flag : flag of the read to know if it is R1 or R2 strand
    read_id_pos_dict : dictionary containing indel position

    Returns
    ----------
    Dictionary with position and number of reads associated
    """

    # for each tuple in cigar_tuple. Ex tuple in cigar tuple : (0, 104)
    for tup in cigar_tuples:

        # If a tuple value correspond to a value in the dictionary
        # Ex : 'D' (Deletion) correspond to value 2 and the cigar tuple is (2, 5)
        if tup[0] == cigar_event_value:

            # Apply function to return position of indel
            indel_position = len_before_indel(chromosome, start_pos_read, cigar_tuples, tup)

            # If position not in dictionary
            if indel_position not in read_id_pos_dict.keys():
                read_id_pos_dict[indel_position] = read_id_flag

            # If position is in dictionary
            else:

                if not isinstance(read_id_pos_dict[indel_position], list):
                    # Then make it list
                    read_id_pos_dict[indel_position] = [read_id_pos_dict[indel_position]]

                # Append the value in list
                read_id_pos_dict[indel_position].append(read_id_flag)

    return read_id_pos_dict


def len_before_indel(chromosome, start_pos_read, cigar_tuples, tup):
    """
    Compute the fragment length before the INDEL thanks to cigartuple and return position of the INDEL

    Parameters
    ----------
    chromosome : chromosome that matches the read
    start_pos_read : start position of the alignment of the read
    cigar_tuples : the dictionary to update
    tup : Indel tuple

    Returns
    ----------
    indel_chr_position : chromosomal position of the concerned INDEL
    """

    # Variable to add nucleotide to the start position of the read until the I or D cigar letter is found
    fragment_len_before_indel = 0

    # Then takes the index of this tuple
    index_tuple = cigar_tuples.index(tup)
    # Use the index to store all the tuples before this one (the one with the indel)
    tuples_of_interest_before_indel = cigar_tuples[:index_tuple]

    if len(tuples_of_interest_before_indel) == 0:

        fragment_len_before_indel = 0
        indel_chr_position = chromosome + ':' + str(start_pos_read + fragment_len_before_indel)

        return indel_chr_position

    else:
        # For each tuple before the tuple that contain the indel
        for tuple_before_indel in tuples_of_interest_before_indel:
            # Store the length of the match until the indel is reached
            fragment_len_before_indel += tuple_before_indel[1] +1
            # Add this length to the alignment start position of the read (position on the reference) (O based)
            indel_chr_position = chromosome + ':' + str(start_pos_read + fragment_len_before_indel)

    return indel_chr_position


# TODO: miss docstring
def create_data_frame_in_file(read_id_pos_dict, general_file):

    enriched_pos_dict = {}
    for key, val in read_id_pos_dict.items():
        enriched_pos_dict[key] = len(val)

    data_frame_pos = general_file + '_df.txt'
    data_frame = pd.DataFrame(enriched_pos_dict.items())
    data_frame.to_csv(data_frame_pos, '\t')
    return data_frame


# ENRICHED POSITION PHASE
def create_enriched_position_read_id_dict(pos_read_id, nb_of_max_val):
    """
    Create a list with all read_id of the enriched position

    Parameters
    ----------
    pos_read_id : dictionary containing reads of each positions
    nb_of_max_val : number of enriched position to analysed (specified by user)

    Returns
    ----------
    read_enriched_pos_dict: A dictionary of read id of enriched pos to analyse (number of position precised by user)
    """
    enriched_position_list = []
    read_enriched_pos_dict = {}
    pos_read_number_dict = {}

    for key, val in pos_read_id.items():
        pos_read_number_dict[key] = len(val)

    max_val = heapq.nlargest(nb_of_max_val, pos_read_number_dict.items(), key=itemgetter(1))  # Use .iteritems() on Py2
    max_val_dict = dict(max_val)

    for pos_key, nb_read_val in max_val_dict.items():
        enriched_position_list.append(pos_key)

    for enriched_pos in enriched_position_list:

        # For positions and all associated read id in dictionary
        for position, read_id in pos_read_id.items():

            # If enriched position is a key
            if position == enriched_pos:
                # Then get all the associated read id
                read_enriched_pos_dict[position] = pos_read_id[position]

    return read_enriched_pos_dict


def create_fasta_on_enriched_position(pos_read_id, nb_of_max_val, filtered_bam, general_file):
    """
    Run a blat on the most represented position.

    Parameters
    ----------

    pos_read_id : dictionary containing reads of each positions
    nb_of_max_val
    filtered_bam : filtered bam file
    general_file : path to the fasta file to create

    Returns
    ----------
    Create fasta files and return a list of generated file paths
    """
    all_fasta_created = []
    # Run the function to create a list with read_id for the most represented indel position
    read_id_dict = create_enriched_position_read_id_dict(pos_read_id, nb_of_max_val)

    for position_key in read_id_dict.keys():
        # Open a file that will contain the read to blat (fasta format)
        pos_fasta_path = general_file + '_' + position_key.replace(':', '_') + '.fa'

        all_fasta_created.append(pos_fasta_path)
        read_to_blat_fasta = open(pos_fasta_path, 'w')

        # Open the file and read in it in binary
        filtered_sam_file = pysam.AlignmentFile(filtered_bam, "rb")

        for read in filtered_sam_file:

            # Get the read name/id
            read_id_bam = str(read.query_name).strip(' ')
            flag = read.flag
            read_id_flag_tup = (read_id_bam, flag)

            if read_id_flag_tup in read_id_dict[position_key]:

                # Get the read complete sequence
                read_seq = str(read.query_sequence)
                # Write in the fasta file the read id and the read sequence
                read_to_blat_fasta.write('>' + read_id_bam + '_' + str(flag) + '\n' + read_seq + '\n')

        read_to_blat_fasta.close()
        filtered_sam_file.close()

    return all_fasta_created


# BLAT ANALYSIS PHASE
def blat(reference, fasta, blat_output):
    """
    Run a BLAT for each read in the fasta file.

    Parameters
    ----------
    reference : reference database for the blat
    fasta : the complete path to the fasta file
    blat_output : the complete path to the output blat file

    Returns
    ----------
    Create a multiple blat file
    """
    # Prepare the command all_line to run a blat for the read
    blat_command = 'blat -stepSize=5 -repMatch=2253 -minScore=20 -minIdentity=0 ' + reference + ' ' + fasta + ' ' + blat_output
    # Run the command all_line
    os.system(blat_command)


def blat_analysis(blat_file, chromosome):
    """
    Analyse BLAT for each blat all_line and only return positions of interesting one

    Parameters
    ----------
    blat_file : the blat file to analyse
    chromosome
    Returns
    ----------
    A dictionary with general raw positions
    """
    # Create a dictionary to store all blat positions
    positions_dict = {}

    # Open a file that will contain all blat lines (psl format)
    blat_file = open(blat_file, 'r')

    # For each all_line of the blat header
    for header_line in blat_file:

        # When the last all_line of the blat header is reached
        if '-----' in header_line:
            break

    # Then begin to read the blat lines
    for blat_line in blat_file:

        # Remove the \n and ',' of in blat all_line
        blat_line = blat_line.strip('\n')
        blat_line = blat_line.strip(',')
        # Put the blat all_line in a list with each column as an element
        split_blat_line = blat_line.split('\t')

        # Store read name
        read_id = split_blat_line[9].split('_')[0]
        flag = split_blat_line[9].split('_')[1]
        # Store start and end of the target to compute the span
        t_start = split_blat_line[15]
        t_end = split_blat_line[16]
        # Store position all_line
        read_positions = flag, split_blat_line[13:]

        # Compute span to chose the blat output
        span = int(t_end) - int(t_start) + 1

        if chromosome == split_blat_line[13]:

            # If the dictionary is not empty (the first blat all_line has been treated)
            if len(positions_dict) != 0:

                # And if the read has not been seen yet
                if read_id not in positions_dict.keys():
                    # Store the read and the associated blat all_line containing all positions
                    positions_dict[read_id] = span, read_positions

                # If the read is already in the dictionary
                else:
                    # And if the span of this other blat output for the same read is highest then before
                    if span > positions_dict[read_id][0]:
                        # Replace the read in the dictionary with the new one
                        positions_dict[read_id] = span, read_positions

            # If the dictionary is empty
            else:
                positions_dict[read_id] = span, read_positions

    # Close the blat file
    blat_file.close()

    # Loop to remove the span of the dictionary
    for key, val in positions_dict.items():
        positions_dict[key] = val[1]

    return positions_dict


def create_alt_ref_block_positions_dict(positions_dict, span_ref_alt, min_block_size):
    """
    Compute dictionary positions (chromosome pos and read pos) to grep the REF sequence and ALT sequence.

    Parameters
    ----------
    positions_dict : dictionary containing target position and query position corresponding to a blat all_line
    span_ref_alt : span to define around the break point
    min_block_size : minimum length for block size

    Returns
    ----------
    ref_alt_pos_dict : dictionary with read id and associated positions of alternative and reference blocks.
    """
    # Create final dictionary to return
    ref_alt_pos_dict = {}

    # For each blat output selected in positions_dict dictionary (after the blat analysis function)
    for key, val in positions_dict.items():

        # Store read flag, blat data (from T name to tStarts), number of aligned blocks
        flag = val[0]
        val = val[1]
        number_aligned_block = val[4]

        # Only analyse blat all_line if number of aligned blocks is 2
        if number_aligned_block == '2':

            # TARGET BASED POSITIONS : Store start and end of alignment (cf. blat file)
            t_start = int(val[2])
            t_end = int(val[3])

            # Store align blocks length
            block_sizes = val[5].split(',')
            q_len_block1 = int(block_sizes[0])
            q_len_block2 = int(block_sizes[1])

            # Store break point position (target positions)
            t_break_point_block1 = t_start + q_len_block1
            t_break_point_block2 = t_end - q_len_block2

            # QUERY BASED POSITIONS : Store starts of the two aligned blocks
            q_starts = val[6].split(',')
            q_start_block1 = int(q_starts[0])
            q_start_block2 = int(q_starts[1])
            q_end_block2 = q_start_block2 + q_len_block2

            # Store break point position (query positions)
            q_break_point1 = q_start_block1 + q_len_block1
            q_break_point2 = q_start_block2

            # Store DELETION positions to grep the sequence
            deletion_length = t_break_point_block2 - t_break_point_block1
            deletion_pos = [t_break_point_block1, deletion_length]
            # Store INSERTION positions to grep the sequence
            insertion_length = q_break_point2 - q_break_point1
            insertion_pos = [q_break_point1, insertion_length]
            # Store DELINS positions to grep the sequence
            delins_pos = [insertion_pos, deletion_pos]

            # Pass blat lines with too small aligned blocks
            if q_len_block1 < min_block_size or q_len_block2 < min_block_size or q_len_block1 + q_len_block2 < span_ref_alt*2:
                pass

            # DELETION
            if q_start_block1 + q_len_block1 == q_start_block2:

                # Define the two REF blocks (around the breaking point)
                target_block1_pos = [t_break_point_block1 - span_ref_alt, t_break_point_block1 + span_ref_alt]
                target_block2_pos = [t_break_point_block2 - span_ref_alt, t_break_point_block2 + span_ref_alt]

                # REF(2 blocks) and ALT(1 block) positions
                ref_pos = [target_block1_pos, target_block2_pos]
                alt_pos = [q_break_point1 - span_ref_alt, q_break_point1 + span_ref_alt]

                # Adjust ALT positions if one block is smallest than the specified span
                alt_pos = adjust_position_first_block(alt_pos)
                alt_pos = adjust_position_second_block(alt_pos, q_end_block2)

                # Update the dictionary for the concerned read (key)
                ref_alt_pos_dict[key] = flag, alt_pos, ref_pos, deletion_pos

            # INSERTION OR SUBSTITUTION
            else:
                # INSERTION
                if t_start + q_len_block1 == t_end - q_len_block2:

                    # Define the two ALT blocks (around the breaking point)
                    query_block1_pos = [q_break_point1 - span_ref_alt, q_break_point1 + span_ref_alt]
                    query_block2_pos = [q_break_point2 - span_ref_alt, q_break_point2 + span_ref_alt]

                    # TODO : Code the case of overlap between alt1 and alt2 (if indel size > 15)
                    # Adjust ALT positions if one block is smallest than the specified span
                    query_block1_pos = adjust_position_first_block(query_block1_pos)
                    query_block2_pos = adjust_position_first_block(query_block2_pos)
                    query_block1_pos = adjust_position_second_block(query_block1_pos, q_end_block2)
                    query_block2_pos = adjust_position_second_block(query_block2_pos, q_end_block2)

                    # REF(1 blocks) and ALT(2 block) positions
                    alt_pos = [query_block1_pos, query_block2_pos]
                    ref_pos = [t_break_point_block1 - span_ref_alt, t_break_point_block1 + span_ref_alt]

                    # Update the dictionary for the concerned read (key)
                    ref_alt_pos_dict[key] = flag, alt_pos, ref_pos, insertion_pos

                # DELETION-INSERTION
                else:

                    # Define the two ALT blocks (around the breaking point)
                    query_block1_pos = [q_break_point1 - span_ref_alt, q_break_point1 + span_ref_alt]
                    query_block2_pos = [q_break_point2 - span_ref_alt, q_break_point2 + span_ref_alt]

                    # Define the two REF blocks (around the breaking point)
                    target_block1_pos = [t_break_point_block1 - span_ref_alt, t_break_point_block1 + span_ref_alt]
                    target_block2_pos = [t_break_point_block2 - span_ref_alt, t_break_point_block2 + span_ref_alt]

                    # Adjust ALT positions if one block is smallest than the specified span
                    query_block1_pos = adjust_position_first_block(query_block1_pos)
                    query_block2_pos = adjust_position_first_block(query_block2_pos)
                    query_block1_pos = adjust_position_second_block(query_block1_pos, q_end_block2)
                    query_block2_pos = adjust_position_second_block(query_block2_pos, q_end_block2)

                    # REF(2 blocks) and ALT(2 block) positions
                    ref_pos = [target_block1_pos, target_block2_pos]
                    alt_pos = [query_block1_pos, query_block2_pos]

                    # Update the dictionary for the concerned read (key)
                    ref_alt_pos_dict[key] = flag, alt_pos, ref_pos, delins_pos

    return ref_alt_pos_dict


def adjust_position_first_block(position_to_adjust):
    """
    Adjust the position of the left block (before the INDEL) is smaller than the defined span

    Parameters
    ----------
    position_to_adjust : entire block to move

    Returns
    ----------
    position_to_adjust : block with new alt positions (more span for the right block and less for the left)
    """
    if position_to_adjust[0] < 0:
        position_to_adjust = [position_to_adjust[0] + 1, position_to_adjust[1] + 1]
        return adjust_position_first_block(position_to_adjust)
    else:
        return position_to_adjust


def adjust_position_second_block(position_to_adjust, limit_position):
    """
    Adjust the position of the right block (after the INDEL) is smaller than the defined span

    Parameters
    ----------
    position_to_adjust : entire block to move
    limit_position : the end position of the right block
    Returns
    ----------
    position_to_adjust : block with new alt positions (more span for the left block and less for the right)
    """
    if position_to_adjust[1] > limit_position:
        position_to_adjust = [position_to_adjust[0] - 1, position_to_adjust[1] - 1]
        return adjust_position_second_block(position_to_adjust, limit_position)
    else:
        return position_to_adjust


def insertion_or_deletion(pos_dict):
    """
    Analyse the alteration type of each read

    Parameters
    ----------
    pos_dict : dictionary with read id as key and as value the flag and blocks positions of associated REF and ALT

    Returns
    ----------
    pos_dict : new dictionary with only reads that correspond to criteria of the most represented alteration
    """
    # INSERTION, DELETION and DELETION-INSERTION counter
    ins_counter = 0
    del_counter = 0
    delins_counter = 0

    # Read each positions blocks (one by read)
    for positions_val in pos_dict.values():
        # If the REF block positions contain only 1 block
        if isinstance(positions_val[2][0], int):
            # Then it is an INSERTION
            ins_counter += 1

        # If the ALT block positions contain only 1 block
        if isinstance(positions_val[1][0], int):
            # Then it is a DELETION
            del_counter += 1

        # If the ALT and the REF blocks positions contain both 2 blocks
        if isinstance(positions_val[1][0], list) and isinstance(positions_val[2][0], list):
            # Then it is a DELINS
            delins_counter += 1

    # Check which alteration is the most represented in the read set and update the dictionary only with
    # reads that correspond to this alteration
    if max(ins_counter, del_counter, delins_counter) == ins_counter:
        pos_dict = {key: val for key, val in pos_dict.items() if isinstance(val[2][0], int)}
        pos_dict['alteration_type'] = 'insertion'

    if max(ins_counter, del_counter, delins_counter) == del_counter:
        pos_dict = {key: val for key, val in pos_dict.items() if isinstance(val[1][0], int)}
        pos_dict['alteration_type'] = 'deletion'

    if max(ins_counter, del_counter, delins_counter) == delins_counter:
        pos_dict = {key: val for key, val in pos_dict.items() if isinstance(val[1][0], list) and isinstance(val[2][0], list)}
        pos_dict['alteration_type'] = 'deletion_insertion'

    return pos_dict


def alteration_pathway(pos_dict, filtered_bam, chromosome, general_file, fasta_ref):
    """
    Depending on the alteration type of the set, run associated functions

    Parameters
    ----------
    pos_dict : sorted dictionary with reads and associated block positions
    filtered_bam : filtered bam to grep ALT sequences in it
    chromosome : the number of the concerned chromosome to search in the REF fasta (bed file creation)
    general_file : output path an name associated to the position analyzed
    fasta_ref : reference genome to search REF sequences

    Returns
    ----------
    Tuple of two lists
    alt_seq_list : list with all ALT sequences found
    fasta_file_list : list with all REF sequences found
    """
    # INSERTION
    if pos_dict['alteration_type'] == 'insertion':
        del pos_dict['alteration_type']
        # Grep ALT sequences
        alt_seq_list = grep_alt(pos_dict, filtered_bam)
        # Create a bed file
        bed_file_list = insertion_create_bed_file(pos_dict, chromosome, general_file)
        bed = bed_file_list[0]
        # Grep REF sequences using the bed file
        fasta_file_list = insertion_grep_sequence(fasta_ref, bed, general_file)

        return alt_seq_list, fasta_file_list

    # DELETION
    if pos_dict['alteration_type'] == 'deletion':
        del pos_dict['alteration_type']

        alt_seq_list = grep_alt(pos_dict, filtered_bam)
        bed_file_list = deletion_and_delins_create_bed_file(pos_dict, chromosome, general_file)
        # Two bed file because two REF sequences for a deletion
        bed1 = bed_file_list[0]
        bed2 = bed_file_list[1]
        bed_indel = bed_file_list[2]
        fasta_file_list = deletion_and_delins_grep_sequence(fasta_ref, bed1, bed2, bed_indel, general_file)

        return alt_seq_list, fasta_file_list

    # DELINS
    if pos_dict['alteration_type'] == 'deletion_insertion':
        del pos_dict['alteration_type']

        alt_seq_list = grep_alt(pos_dict, filtered_bam)
        bed_file_list = deletion_and_delins_create_bed_file(pos_dict, chromosome, general_file)
        bed1 = bed_file_list[0]
        bed2 = bed_file_list[1]
        bed_indel = bed_file_list[2]
        fasta_file_list = deletion_and_delins_grep_sequence(fasta_ref, bed1, bed2, bed_indel, general_file)

        return alt_seq_list, fasta_file_list


# GREP SEQUENCES PHASE
# ins
def insertion_create_bed_file(pos_dict, chromosome, general_file):
    """
    Create a BED file to grep REF sequences in the reference fasta

    Parameters
    ----------
    pos_dict : sorted dictionary with reads and associated block positions
    chromosome : chromosome concerned by the alteration in the reference
    general_file : general path to create bed file in the good position

    Returns
    ----------
    String
    bed : string which is the path of the newly created bed file
    """
    # Define the name of the bed file and open it in write mode
    bed = general_file + '_ref1.bed'
    write_bed = open(bed, 'w')

    # for each element of the REF ALT positions dictionary
    for read_id, position in pos_dict.items():
        ref_block = position[2]
        # Write positions to search of the sequences around indel in the bed file
        write_bed.write(chromosome + '\t' + str(ref_block[0]) + '\t' + str(ref_block[1]) + '\t' + read_id + '\n')

    write_bed.close()

    return [bed]


def deletion_and_delins_create_bed_file(pos_dict, chromosome, general_file):
    """
    Create a BED file to grep REF sequences in the reference fasta

    Parameters
    ----------
    pos_dict : sorted dictionary with reads and associated block positions
    chromosome : chromosome concerned by the alteration in the reference
    general_file : general path to create bed file in the good position

    Returns
    ----------
    List
    bed1 : string which is the path of the newly created block 1 BED file
    bed2 : string which is the path of the newly created block 2 BED file
    """
    # Create and open the two BED files
    bed1, bed2, bed_indel = general_file + '_ref1.bed', general_file + '_ref2.bed', general_file + '_ref_indel.bed'
    write_bed1, write_bed2, write_bed_indel = open(bed1, 'w'), open(bed2, 'w'), open(bed_indel, 'w')

    # for each element of the REF ALT positions dictionary
    for read_id, position in pos_dict.items():
        # Define REF block 1 and block 2
        block1 = position[2][0]
        block2 = position[2][1]
        indel = position[3]
        # DELINS
        if isinstance(indel[0], list):
            del_ref = position[3][1]
        else:
            del_ref = position[3]

        # Write positions to search REF sequences in BED file
        # REF block 1
        write_bed1.write(chromosome + '\t' + str(block1[0]) + '\t' + str(block1[1]) + '\t' + read_id + '\n')
        # REF block 2
        write_bed2.write(chromosome + '\t' + str(block2[0]) + '\t' + str(block2[1]) + '\t' + read_id + '\n')
        # REF INDEL
        # Minus one for the annotation with snpEff
        write_bed_indel.write(chromosome + '\t' + str(del_ref[0] - 1) + '\t' + str(del_ref[0] + del_ref[1]) + '\t' + read_id + '\n')

    write_bed1.close()
    write_bed2.close()
    write_bed_indel.close()
    return [bed1, bed2, bed_indel]


def grep_alt(pos_dict, filtered_bam):
    """
    Grep ALT sequences from read containing an insertion or a delins

    Parameters
    ----------
    pos_dict : sorted dictionary with reads and associated block positions
    filtered_bam : filtered bam to grep ALT sequences in it

    Returns
    ----------
    A tuple of two list or just one list
    """
    # Initialisation of the list and open the BAM file in read mode
    alt_list = []
    indel_list = []
    filtered_bam = pysam.AlignmentFile(filtered_bam, "rb")

    # Browse all reads in the filtered BAM
    for read in filtered_bam:
        # For each read only store read id and associated flag
        read_id_bam = str(read.query_name)
        flag_bam = int(read.flag)

        # When the read id in the dictionary is found in BAM
        if read_id_bam in pos_dict.keys():
            # Get the read id associated value and flag in the dictionary (which is query based sequences positions)
            positions = pos_dict[read_id_bam]
            flag_dict = int(positions[0])

            # When the read id and same flag (same strand) is found
            if flag_bam == flag_dict:
                # Store the sequence of the read
                read_seq = str(read.query_sequence)

                # INS OR DELINS
                if isinstance(positions[1][0], list):
                    # Extract block 1 and block 2 lists positions
                    block1 = positions[1][0]
                    block2 = positions[1][1]
                    # Only grep sequences that correspond to positions of block 1 and block 2
                    alt_seq_block1 = read_seq[block1[0]: block1[1]]
                    alt_seq_block2 = read_seq[block2[0]: block2[1]]
                    alt_seq = [alt_seq_block1, alt_seq_block2]

                    # DELINS
                    if isinstance(positions[3][0], list):
                        # Grep DELINS sequence
                        delins_alt = positions[3][0]
                        indel_alt_seq = read_seq[delins_alt[0]:delins_alt[0] + delins_alt[1]]

                    # INS
                    else:
                        # Grep INS sequence
                        insertion = positions[3]
                        indel_alt_seq = read_seq[insertion[0]:insertion[0] + insertion[1]]

                    indel_list.append(indel_alt_seq)

                # DEL
                else:
                    # Extract the only one block of ALT  positions
                    block = positions[1]
                    # Only grep sequences that correspond to positions of ALT block
                    alt_seq = read_seq[block[0]: block[1]]

                    # # Grep DEL sequence
                    # deletion = positions[3]
                    # indel_alt_seq = read_seq[deletion[0]:deletion[0] + deletion[1]]

                # Add this sequence to the ALT list
                alt_list.append(alt_seq)

    filtered_bam.close()

    # INS or DELINS
    if len(indel_list) != 0:
        return alt_list, indel_list
    # DEL
    else:
        return alt_list


def insertion_grep_sequence(fasta_ref, bed, general_file):
    """
    Grep REF sequences for an INSERTION in reference fasta file using a BED file

    Parameters
    ----------
    fasta_ref : reference genome in fasta format
    bed : BED with positions of the REF block
    general_file : general path to create
    Returns
    ----------
    List with the path file of the fasta created
    """
    fasta1 = general_file + '_ref.fa'

    # Grep the indel sequence and the sequence around indel in fasta (reference and alignment fasta)
    bedtools_getfasta_ref_block1 = 'bedtools getfasta -fi ' + fasta_ref + ' -bed ' + bed + ' > ' + fasta1

    # Run both command
    os.system(bedtools_getfasta_ref_block1)

    return [fasta1]


def deletion_and_delins_grep_sequence(fasta_ref, bed1, bed2, bed_indel, general_file):
    """
    Grep REF sequences for a DELETION or a DELINS in reference fasta file using a BED file

    Parameters
    ----------
    fasta_ref : reference genome in fasta format
    bed1 : BED with positions of block 1
    bed2 : BED with positions of block 2
    general_file : general path to create
    Returns
    ----------
    Tuple with a list and a file path
    """
    fasta1 = general_file + '_ref_block1.fa'
    fasta2 = general_file + '_ref_block2.fa'
    fasta3 = general_file + '_indel_ref.fa'

    # Grep the indel sequence and the sequence around indel in fasta (reference and alignment fasta)
    bedtools_getfasta_ref_block1 = 'bedtools getfasta -fi ' + fasta_ref + ' -bed ' + bed1 + ' > ' + fasta1
    bedtools_getfasta_ref_block2 = 'bedtools getfasta -fi ' + fasta_ref + ' -bed ' + bed2 + ' > ' + fasta2
    bedtools_getfasta_ref_delins = 'bedtools getfasta -fi ' + fasta_ref + ' -bed ' + bed_indel + ' > ' + fasta3

    # Run both command
    os.system(bedtools_getfasta_ref_block1)
    os.system(bedtools_getfasta_ref_block2)
    os.system(bedtools_getfasta_ref_delins)

    return [fasta1, fasta2], fasta3


# COUNT PHASE
def sequence_count(sequences):
    """
    Statistics on all ALT and REF sequences found
    sequences : list of ALT sequences or string with the fasta path to REF sequences
    Parameters
    ----------
    Returns
    ----------
    Dictionaries which associate to a sequence the number of time it is found
    Two dictionary, first one the indel sequences and values, second one the around indel sequences and values
    """

    # ALT sequences
    if isinstance(sequences, list):
        # Stats on ALT sequences
        seq_stats = Counter(sequences)

        return seq_stats

    # REF sequences
    else:
        seq_list = []
        # Read the fasta file
        fasta = pyfastx.Fasta(sequences)
        for seq in fasta:
            seq_list.append(str(seq))
        # Stats on REF sequences
        seq_stats = Counter(seq_list)

        return seq_stats


def grep_seq_on_samples(seq_stats, fastq_r1, fastq_r2):
    """
    Count for the most common sequence (indel and around indel) the number of match with samples

    Parameters
    ----------
    seq_stats
    fastq_r1 : sample 1 which is a fastq
    fastq_r2 : sample 2 which is a fastq

    Returns
    ----------
    Dictionary with percentage for R1 and R2
    """
    # For the sequence found the most
    sequence_found = max(seq_stats, key=seq_stats.get)

    # Grep this sequence in both samples
    sequence_grep_r1 = os.popen('grep ' + sequence_found + ' ' + fastq_r1 + ' | wc -l').readlines()
    sequence_grep_r2 = os.popen('grep ' + sequence_found + ' ' + fastq_r2 + ' | wc -l').readlines()
    # Store the result
    sequence_grep_r1 = int(sequence_grep_r1[0].strip('\n'))
    sequence_grep_r2 = int(sequence_grep_r2[0].strip('\n'))

    return {'R1': round(sequence_grep_r1), 'R2': round(sequence_grep_r2)}


def compute_percentage(alt_count_dict, ref_count_dict):
    """
    Compute the percentage of the INDEL in samples

    Returns
    ----------
    Return VAF percentage of the alteration found in both samples
    """
    # Store ALT count found in R1 and R2
    alt_r1, alt_r2 = list(alt_count_dict.values())[0], list(alt_count_dict.values())[1]
    # Store REF count found in R1 and R2
    ref_r1, ref_r2 = list(ref_count_dict.values())[0], list(ref_count_dict.values())[1]
    # Compute sum for REF and ALT in R1 and R2
    sum_r1, sum_r2 = alt_r1 + ref_r1, alt_r2 + ref_r2

    try:
        # Compute percentage in R1 and R2
        alt_percentage_r1, alt_percentage_r2 = (alt_r1 / sum_r1) * 100, (alt_r2 / sum_r2) * 100

        return {'R1': str(round(alt_percentage_r1, 1)), 'R2': str(round(alt_percentage_r2, 1))}

    except sum_r1 == 0 or sum_r2 == 0:
        return 'ERROR'


def ins_dup(seq_stats_alt1, seq_stats_alt2):
    """
    Determine if the insertion is duplicated sequence or not by searching matches in block 1 and block 2 of the ALT seq

    Parameters
    ----------
    seq_stats_alt1 :
    seq_stats_alt2 :

    Returns
    ----------
    Return number of match found
    """
    # For mostly found ALT sequences
    max_alt1_seq = max(seq_stats_alt1, key=seq_stats_alt1.get)
    max_alt2_seq = max(seq_stats_alt2, key=seq_stats_alt2.get)

    # Use SequenceMatcher to search matches between 2 strings
    match = SequenceMatcher(None, max_alt1_seq, max_alt2_seq).find_longest_match(0, len(max_alt1_seq), 0, len(max_alt2_seq))
    seq_matches = max_alt1_seq[match.a: match.a + match.size]

    return seq_matches


def vcf_file(general_file, position_enriched, alt, nt_before_del=None):
    """
    Create a VCF format file for the annotation phase

    Parameters
    ----------

    Returns
    ----------
    IDK
    """
    chrom = position_enriched.split('_')[0]
    pos = position_enriched.split('_')[1]

    # DATAFRAME FILE CREATION
    file_vcf = general_file + '_' + position_enriched + '.vcf'
    vcf_data = {'#CHROM': chrom, 'POS': pos, 'ID': '.'}
    if nt_before_del:
        alt_in_file(alt, vcf_data, nt_before_del)
    else:
        alt_in_file(alt, vcf_data)

    vcf_data['QUAL'] = '.'
    vcf_data['FILTER'] = '.'
    vcf_data['INFO'] = '.'

    data_frame = pd.DataFrame(vcf_data, columns=vcf_data.keys(), index=[0])
    data_frame.to_csv(file_vcf, sep='\t', na_rep='NA', index=False)

    return file_vcf


def alt_in_file(alt, data, nt_before_del=None):
    """
    Store associated information in the dictionary depending on the alteration type

    Parameters
    ----------
    alt : dictionary containing alteration type
    data : dictionary to create the final result file

    Returns
    ----------
    Nothing
    """
    if len(alt.keys()) == 1:
        if 'insertion' in alt.keys():
            data['REF'] = '.'
            data['ALT'] = alt['insertion']
        if 'deletion' in alt.keys():
            if nt_before_del:
                data['REF'] = nt_before_del + alt['deletion']
                data['ALT'] = nt_before_del
            else:
                data['REF'] = alt['deletion']
                data['ALT'] = '.'
    else:
        data['REF'] = alt['deletion']
        data['ALT'] = alt['insertion']


def create_annotation(vcf_input, snpeff_path, ram=None):
    """
    Create an annotation for the concerning variant

    Parameters
    ----------
    vcf_input :
    snpeff_path :
    ram :

    Returns
    ----------
    IDK
    """
    vcf_output = vcf_input.split('.')[0] + '_ann' + '.vcf'
    snpeff_cmd = 'java ' + '-Xmx' + ram + ' -jar ' + snpeff_path + ' eff -noStats hg19 ' + vcf_input + ' > ' + vcf_output
    os.system(snpeff_cmd)

    return vcf_output


def output_file(output_folder, main_output):
    """
    Create a tabulate output file

    Parameters
    ----------
    output_folder :
    main_output :

    Returns
    ----------
    Nothing
    """
    data = {'chromosome': main_output[1].split('_')[0], 'position': main_output[1].split('_')[1]}

    alt_in_file(main_output[9], data)

    data['gene_name'] = main_output[10][0]
    data['NM'] = main_output[10][1]
    data['C.'] = main_output[10][2]
    data['p.'] = main_output[10][3]
    data['TYPE'] = main_output[2]

    # COMPUTE VAF
    for key, val in main_output[5].items():
        data['alt_count_' + str(key)] = val

    for key, val in main_output[6].items():
        data['ref_count_' + str(key)] = val

    mean_vaf = 0
    for key, val in main_output[7].items():
        data['VAF_' + str(key)] = val
        mean_vaf += float(val)

    mean_vaf = mean_vaf / 2
    data['VAF_mean'] = mean_vaf

    # SEQ BLOCK ALT
    for key, val in main_output[3].items():
        data['alt_block' + str(key)] = val

    # SEQ BLOCK REF
    for key, val in main_output[4].items():
        data['ref_block' + str(key)] = val

    # ANALYSIS STATE
    print(main_output[8])

    # DATAFRAME FILE CREATION
    dataframe_file = output_folder + str(main_output[0]) + '_' + str(main_output[1]) + '_results_df.txt'
    data_frame = pd.DataFrame(data, columns=data.keys(), index=[0])
    data_frame.to_csv(dataframe_file, sep='\t', na_rep='NA')


def bwa(reference, fastq1, fastq2, bam_name, optional_folder):
    """
    Store associated information in the dictionary depending on the alteration type

    Parameters
    ----------
    reference : reference genome
    fastq1 : patient sample
    fastq2 : patient sample
    bam_name : name of the bam file
    optional_folder : folder to delete by default at the end of the analysis

    Returns
    ----------
    The bwa bam file path
    """

    path_to_barcode = optional_folder + bam_name

    command_bwa = 'bwa mem ' + reference + ' ' + fastq1 + ' ' + fastq2 + ' > ' + path_to_barcode + '_bwa.sam'
    sam_to_bam_command = 'samtools view -S -b ' + path_to_barcode + '_bwa.sam' + ' > ' + path_to_barcode + '_bwa.bam'
    sorted_bam_command = 'samtools sort -o ' + path_to_barcode + '_bwa.bam' + ' ' + path_to_barcode + '_bwa.bam'

    os.system(command_bwa)
    os.system(sam_to_bam_command)
    os.system(sorted_bam_command)

    return path_to_barcode + '_bwa.bam'


def parse_vcf_ann_file(vcf_ann_file, ann_nm_file=None):
    """
    Compute the fragment length before the INDEL thanks to cigartuple and return position of the INDEL

    Parameters
    ----------
    vcf_ann_file : a vcf file for the snpeff annotation
    ann_nm_file : a file with NM to get the corresponding annotation if it is possible

    Returns
    ----------
    list with the c. and the p. annotation
    """

    new_column_name_list = []
    # Open the VCF file
    with open(vcf_ann_file, 'r') as vcf_ann:

        # Browse lines of the VCF file
        for all_line in vcf_ann:

            # Only store lines before header
            if '##INFO=<ID=ANN' in all_line:
                all_line = all_line.split(',')
                description = all_line[3]
                column_name_list = description.split("'")[1].split('|')

                for elem in column_name_list:
                    new_elem = elem.replace(' ', '')
                    new_column_name_list.append(new_elem)

                index_num_col_gene_name = new_column_name_list.index('Gene_Name')
                index_nm = new_column_name_list.index('Feature_ID')
                index_c_ann = new_column_name_list.index('HGVS.c')
                index_p_ann = new_column_name_list.index('HGVS.p')

            if '#CHROM' in all_line:
                break

        for line in vcf_ann:
            vcf_line = line.split('\t')
            ann_line = vcf_line[7].strip('ANN=')
            ann_line = ann_line.split(',')

            for annotation in ann_line:
                ann_col = annotation.split('|')
                gene_name = ann_col[index_num_col_gene_name]
                nm = ann_col[index_nm].split('.')[0]
                c_ann = ann_col[index_c_ann]
                p_ann = ann_col[index_p_ann]

                if ann_nm_file:
                    with open(ann_nm_file, 'r') as nm_file:
                        for file_nm_line in nm_file:
                            file_nm = file_nm_line.strip('\n').replace(' ', '').split('\t')
                            if file_nm[-1] == nm:
                                return [gene_name, nm, c_ann, p_ann]

        ann_line = ann_line[0]
        ann_col = ann_line.split('|')
        gene_name = ann_col[index_num_col_gene_name]
        nm = ann_col[index_nm]
        c_ann = ann_col[index_c_ann]
        p_ann = ann_col[index_p_ann]

        return [gene_name, nm, c_ann, p_ann]
