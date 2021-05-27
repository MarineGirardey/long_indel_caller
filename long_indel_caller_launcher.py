#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr 14 08:49:55 2021
@author: Marine Girardey
"""

# Import modules
from long_indel_caller_functions import *
import argparse
import os
import sys
import shutil

program_description = '''
Detect long INDELs in a BAM file an return a informations about the INDEL and
the Variant Allele Frequency (VAF).
'''

parser = argparse.ArgumentParser(add_help=True, description=program_description)

parser.add_argument('-fp', '--filter_pos', help='Position to filter the BAM file', default=None, metavar='filter_position', type=str, required=True)
parser.add_argument('-o', '--output_folder', help='PPath to the output folder', default=None, metavar='output', type=str, required=True)
parser.add_argument('-b', '--bam', help='Path to the bam file', default=None, metavar='BAM', type=str, required=True)
parser.add_argument('-r', '--reference', help='Path to the reference genome', default=None, metavar='ref', type=str, required=True)
parser.add_argument('-fq1', '--fastq1', help='Path to the fastq sample', default=None, metavar='fastq1', type=str, required=True)
parser.add_argument('-fq2', '--fastq2', help='Path to the fastq sample', default=None, metavar='fastq2', type=str, required=True)
parser.add_argument('-snpeff', '--snpeff', help='Path to snpeff', default=None, metavar='snpeff', type=str, required=True)
parser.add_argument('-ram', '--ram', help='Ram for snpeff', default='8g', metavar='ram', type=str, required=False)
parser.add_argument('-nm', '--nm', help='NM for annotation', metavar='nm', type=str, required=False)
parser.add_argument('-bwa', '--bwa', help='Create bwa bam and realize the analysis on it', default=None, action='store_true', required=False)

parser.add_argument('-sp', '--span', help='Position span to filter BAM (default: %(default)s)', default=200, metavar='span', type=int, required=False)
parser.add_argument('-mq', '--map_quality', help='Minimum mapping quality for the filtering part (default: %(default)s)', default=10, metavar='map_qua', type=int, required=False)
parser.add_argument('-bsp', '--blocks_span', help='Block span to map the break point (default: %(default)s)', default=15, metavar='min_block_size', type=int, required=False)
parser.add_argument('-mbs', '--min_len_blocks_size', help='Min blocks size length acceptable (default: %(default)s)', default=5, metavar='span', type=int, required=False)
parser.add_argument('-nbp', '--nb_pos_analyse', help='Number of position to analyse (default: %(default)s)', default=1, metavar='nb_pos', type=int, required=False)
parser.add_argument('-kf', '--keep_folder', help='Keep the folder with other files', action='store_true', required=False)

args = parser.parse_args()


def main():
    """
    Main function that run the application 'long_indel_detection_functions.py'

    Returns
    ----------
    A result vcf file
    """

    # Preparation
    bam_name = args.bam.split('/')[-1].split('.')[0]
    optional_folder = args.output_folder + bam_name + '_others/'
    fasta_folder = args.output_folder + bam_name + '_fasta/'

    try:
        os.mkdir(optional_folder)
        os.mkdir(fasta_folder)
    except OSError:
        pass

    if args.bwa:
        bam = bwa(args.reference, args.fastq1, args.fastq2, bam_name, optional_folder)
        bam_name = bam.split('/')[-1].split('.')[0]
        check_index(bam)

    else:
        bam = args.bam
        check_index(args.bam)

    fasta_folder_bam = fasta_folder + bam_name
    optional_folder_bam_name = optional_folder + bam_name

    filtered_bam = optional_folder_bam_name + '_filtered.bam'
    # FILTERING PHASE
    filter_bam(bam, args.span, args.filter_pos, args.map_quality, filtered_bam)
    check_index(filtered_bam)

    # CIGAR STATISTICS PHASE
    read_id_by_position = stats_cigar_indel(filtered_bam)
    create_data_frame_in_file(read_id_by_position, optional_folder_bam_name)

    # Creation of FASTA
    fasta_created = create_fasta_on_enriched_position(read_id_by_position,
                                                      args.nb_pos_analyse, filtered_bam, fasta_folder_bam)
    # RUN BLAT PHASE
    for fasta_f in fasta_created:

        # blat preparation
        chromosome = fasta_f.split('/')[-1].split('_')[-2]
        pos = fasta_f.split('/')[-1].split('_')[-1].split('.')[0]
        position_enriched = chromosome + '_' + pos

        path_position_enriched = optional_folder_bam_name + '_' + position_enriched.replace(':', '_')
        blat_output = path_position_enriched + '.psl'
        # Run blat tool
        blat(args.reference, fasta_f, blat_output)

        # BLAT ANALYSIS PHASE
        blat_output_dict = blat_analysis(blat_output, chromosome)

        pos_dict = create_alt_ref_block_positions_dict(blat_output_dict, args.blocks_span, args.min_len_blocks_size)

        # TODO : Remove this
        if len(pos_dict) < 10:
            pass

        else:
            # ANALYSE ALTERATION TYPE PHASE
            pos_dict_sorted_on_alteration_type = insertion_or_deletion(pos_dict)
            alt_type = pos_dict_sorted_on_alteration_type['alteration_type']
            alt_ref_seq = alteration_pathway(pos_dict_sorted_on_alteration_type, filtered_bam, chromosome,
                                             path_position_enriched, args.reference)

            # VARIABLE PREPARATION DEPENDING ON ALTERATION TYPE
            if alt_type == 'deletion':
                alt_seq_list = alt_ref_seq[0]
                ref_fasta = alt_ref_seq[1][0]
                indel_seq_list = alt_ref_seq[1][1]

            if alt_type == 'insertion':
                # Output analyse preparation
                alt_seq_list = alt_ref_seq[0][0]
                indel_seq_list = alt_ref_seq[0][1]
                ref_fasta = alt_ref_seq[1][0]

            if alt_type == 'deletion_insertion':
                # Output analyse preparation
                alt_seq_list = alt_ref_seq[0][0]
                indel_seq_list = alt_ref_seq[0][1]
                ref_fasta = alt_ref_seq[1][0]

            ############ ALT SEQUENCE ############

            alt = {}
            seq_alt = {}
            seq_ref = {}
            seq_alt1 = []
            seq_alt2 = []
            ref_count_dict = {}
            alt_count_dict = {}
            # INS or DELINS
            # INS
            if isinstance(alt_seq_list[0], list):
                # Insert in list
                for alt_seq in alt_seq_list:
                    seq_alt1.append(alt_seq[0])  # BLOCK 1
                    seq_alt2.append(alt_seq[1])  # BLOCK 2

                # COUNT
                dict_stats_alt1 = sequence_count(seq_alt1)
                dict_stats_alt2 = sequence_count(seq_alt2)
                insertion = sequence_count(indel_seq_list)
                alt['insertion'] = max(insertion, key=insertion.get)
                seq_alt[1] = max(dict_stats_alt1, key=dict_stats_alt1.get)
                seq_alt[2] = max(dict_stats_alt2, key=dict_stats_alt2.get)
                # GREP

                alt1_count_dict = grep_seq_on_samples(dict_stats_alt1, args.fastq1, args.fastq2)
                alt2_count_dict = grep_seq_on_samples(dict_stats_alt2, args.fastq1, args.fastq2)
                seq_match = ins_dup(dict_stats_alt1, dict_stats_alt2)

                # DUPLICATED INSERTION
                if len(seq_match) == args.blocks_span:
                    alt_count_dict = alt1_count_dict
                else:
                    alt_count_mean_r1 = list(alt1_count_dict.values())[0] + list(alt2_count_dict.values())[0] / 2
                    alt_count_mean_r2 = list(alt1_count_dict.values())[1] + list(alt2_count_dict.values())[1] / 2
                    alt_count_dict['R1'] = round(alt_count_mean_r1)
                    alt_count_dict['R2'] = round(alt_count_mean_r2)

                # DELINS
                if isinstance(ref_fasta, list):
                    indel_fasta = alt_ref_seq[1][1]
                    deletion = sequence_count(indel_fasta)
                    alt['deletion'] = max(deletion, key=deletion.get)

            # DEL
            else:
                # COUNT
                deletion = sequence_count(indel_seq_list)
                alt['deletion'] = max(deletion, key=deletion.get)
                nt_before_del = alt['deletion'][0]
                alt['deletion'] = alt['deletion'][1:]
                dict_stats_alt = sequence_count(alt_seq_list)
                seq_alt[1] = max(dict_stats_alt, key=dict_stats_alt.get)
                # GREP
                alt_count_dict = grep_seq_on_samples(dict_stats_alt, args.fastq1, args.fastq2)

            ############ REF SEQUENCE ############

            # DEL or DELINS
            if len(ref_fasta) == 2:
                ref_fasta1 = ref_fasta[0]
                ref_fasta2 = ref_fasta[1]
                # COUNT
                dict_stats_ref1 = sequence_count(ref_fasta1)
                dict_stats_ref2 = sequence_count(ref_fasta2)
                seq_ref[1] = max(dict_stats_ref1, key=dict_stats_ref1.get)
                seq_ref[2] = max(dict_stats_ref2, key=dict_stats_ref2.get)
                # GREP
                ref1_count_dict = grep_seq_on_samples(dict_stats_ref1, args.fastq1, args.fastq2)
                ref2_count_tuple = grep_seq_on_samples(dict_stats_ref2, args.fastq1, args.fastq2)
                ref_count_mean_r1 = list(ref1_count_dict.values())[0] + list(ref2_count_tuple.values())[0] / 2
                ref_count_mean_r2 = list(ref1_count_dict.values())[1] + list(ref2_count_tuple.values())[1] / 2
                ref_count_dict['R1'] = round(ref_count_mean_r1)
                ref_count_dict['R2'] = round(ref_count_mean_r2)
            # INS
            else:

                # COUNT
                dict_stats_ref = sequence_count(ref_fasta)
                seq_ref[1] = max(dict_stats_ref, key=dict_stats_ref.get)
                # GREP
                ref_count_dict = grep_seq_on_samples(dict_stats_ref, args.fastq1, args.fastq2)

            # COMPUTE PERCENTAGE
            dict_percentage = compute_percentage(alt_count_dict, ref_count_dict)

        # ANALYSIS SUCCEEDED OR NOT
        try:
            if dict_percentage:
                vcf_input = vcf_file(optional_folder_bam_name, position_enriched, alt, nt_before_del)
                vcf_output = create_annotation(vcf_input, args.snpeff, args.ram)
                if args.nm:
                    annotations = parse_vcf_ann_file(vcf_output, args.nm)
                else:
                    annotations = parse_vcf_ann_file(vcf_output)

                analyse_state = 'Analysis succeeded - result file created in folder'

                if args.keep_folder:
                    pass
                else:
                    shutil.rmtree(optional_folder)

                return bam_name, position_enriched, alt_type, seq_alt, seq_ref, alt_count_dict, ref_count_dict, dict_percentage, analyse_state, alt, annotations
                break

        except UnboundLocalError:
            analyse_state = '\nThe analysis failed for these position. You can try to : \n' \
                            '- Check your path and/or adjust your parameters\n' \
                            '- Analyse another position\n'

            print('\n----------------------')
            print('POSITION : ' + position_enriched)
            print(analyse_state)

    if args.keep_folder:
        pass
    else:
        shutil.rmtree(optional_folder)
    return '\nAnalysis failed for position(s) precised. If not already, try BWA or a manual review\n'


if __name__ == '__main__':

    main_output = main()

    # IF ANALYSIS SUCCEEDED
    if isinstance(main_output, tuple):
        output_file(args.output_folder, main_output)

    # IF ANALYSIS FAILED
    else:
        print(main_output)
