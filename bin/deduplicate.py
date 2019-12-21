#!/usr/bin/env python3

'''
BabrahamLinkON: Analysis pipeline for VDJ-seq
Copyright (C) 2017  Peter Chovanec

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <https://www.gnu.org/licenses/>.
'''

import os
import pysam
import glob
from pathlib import Path
import subprocess
import shutil
import re
from collections import defaultdict, Counter
import fileinput
import shlex
import argparse
import pandas as pd
import numpy as np
#Do not hard code matplotlib backend, use export MPLBACKEND=pdf instead if running on headless node
# import matplotlib
# matplotlib.use('pdf')
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.pyplot as plt

from babrahamlinkon import general, presets, umi_correction, UMI_seqlogo, deduplication_general
import multiprocessing
import math

import itertools
from joblib import Parallel, delayed
import logging
# import pickle
# import warnings

import Levenshtein
from tqdm import tqdm
import pyximport
from babrahamlinkon._dedup_umi import edit_distance
from babrahamlinkon.version import __version__
# from memory_profiler import profile




def parse_args():
    parser = argparse.ArgumentParser(description='BabrahamLinkON Deduplicate')

    parser.add_argument('--version', action='version', version='%(prog)s ' + __version__)

    sub = parser.add_subparsers(dest='action', description='Choose pipeline')

    sp1 = sub.add_parser('umi')
    sp2 = sub.add_parser('short')
    sp3 = sub.add_parser('short_anchor')
    sp4 = sub.add_parser('no_anchor')
    sp5 = sub.add_parser('umi_seq_logo')
    sp6 = sub.add_parser('reverse_complement')

    for sp in [sp1,sp2,sp3,sp4,sp5]: #common to all 3

        sp.add_argument('--input_dir', dest='input', type=str, required=True, help='Input directory (created for/by preclean) or can specify a file')
        sp.add_argument('--out', dest='out_dir', type=str, help='Output directory, default: creates Deduplicated in main directory')

    for sp in [sp1,sp2,sp3,sp4]:
        sp.add_argument('--threads', dest='nthreads', default=1, type=int, help='Number of threads to use (if aligning), default: 1')

        sp.add_argument('--mismatch', dest='mismatch', type=int, default=5, help='Number of mismatches allowed in consensus sequence comparison [5]')
        sp.add_argument('--min_reads', dest='minreads', type=int, default=2, help='Minimum number of reads in UMI group, if less than or equal to [2] then discard')
        sp.add_argument('--gt_ratio', dest='gtratio', type=float, default=1, help='Ratio of good to total reads to mark UMI group as erroneous 0-1 [1]')

        sp.add_argument('--stats', action='store_true', help='Output stats from UMI deduplication [False]')

        sp.add_argument('--umi_correction', action='store_true', help='Perform correction of errors that might be present in UMI')
        sp.add_argument('--threshold', dest='threshold', type=int, default=1, help='Number of mismatches allowed in UMI when doing UMI correction [1]')

        sp.add_argument('--skip_unclear', action='store_true', help='Skip unclear J reads [False]')
        sp.add_argument('--keep_mh', action='store_true', help='Keep multiple hit unclear J reads [False]')

        sp.add_argument('--use_j', action='store_true', help='Deduplicate using J identity')
        sp.add_argument('--ignore_umi', action='store_true', help='Deduplicate without using UMI')

        sp.add_argument('--j_trim', dest='j_trim', default=25, type=int, help='Trim J primer when comparing to consensus, default: 25')
  
        sp.add_argument('--no_msa', dest='no_msa', action='store_true', help='Don\'t use msa to derive consensus sequence [False]')
        sp.add_argument('--fq', dest='fq', action='store_true', help='Output fastq instead of fasta')
        sp.add_argument('--cons_no_qual', dest='cons_no_qual', action='store_true', help='Make consensus without using quality scores')
        sp.add_argument('--with_N', dest='with_N', action='store_true', help='Output consensus sequences with ambigious N bases')

    for sp in [sp1, sp3, sp5, sp6]:


        sp.add_argument('--an1', dest='an1', default='GACTCGT', type=str, help='Default: GACTCGT')
        sp.add_argument('--an2', dest='an2', default='CTGCTCCT', type=str, help='Default: CTGCTCCT')



    sp6.add_argument('--input', dest='input', type=str, required=True, help='Input file/directory with files')
    sp6.add_argument('--fq', dest='fq', action='store_true', help='Convert fastq')


    sp1.set_defaults(short=False, no_anchor=False, seqlogo=False, rev_comp=False)
    sp2.set_defaults(short=True, no_anchor=True, fq=False, seqlogo=False, rev_comp=False)
    sp3.set_defaults(short=True, no_anchor=False, seqlogo=False, rev_comp=False)
    sp4.set_defaults(short=False, no_anchor=True, seqlogo=False, rev_comp=False)
    sp5.set_defaults(seqlogo=True, rev_comp=False, cons_no_qual=False, no_msa=False, no_anchor=False)
    sp6.set_defaults(seqlogo=False, rev_comp=True, cons_no_qual=False, no_msa=False, no_anchor=False, in_dir=None)


    opts = parser.parse_args()

    return opts



def main():

    #argparse
    opts = parse_args()

    #combinations which can't be used together
    if opts.cons_no_qual and opts.fq:
        raise Exception('Can\'t output fastq without producing consensus quality')
    if opts.no_msa and opts.fq:
        raise Exception('Can\'t output fastq without producing consensus quality which is only used with msa at the moment')

    if opts.no_anchor:
        an1 = ''
        an2 = ''
    # elif opts.short:
    #     an1 = ''
    #     an2 = ''
    else:
        an1 = opts.an1
        an2 = opts.an2

    #directory or file submitted
    if os.path.isfile(opts.input):
        in_dir = os.path.basename(opts.input)
        in_file = True
    #if only directory provided
    elif os.path.isdir(opts.input):
        in_dir = opts.input
        in_file = False


    if not opts.rev_comp:
        #initiate deduplication object
        dedup = deduplicate(file_directory=opts.input, an1=an1, an2=an2)
        dedup.create_dirs_assembled(out_dir=opts.out_dir)


    if opts.seqlogo:
        if opts.no_anchor:
            if in_file:
                jv_fastq = opts.input
            else:
                jv_fastq = glob.glob(in_dir + '/*all_j')[0]
            UMI_seqlogo.umi_seq_logo(jv_fastq, dedup.out_dir + '/' + dedup.jv_prefix + '.eps')
        else:
            if in_file:
                raise ValueError('Cannot submit multiple file paths, use directory instead')
            else:
                jv_fastq_an1 = glob.glob(in_dir + '/*all_j*' + an1)[0]
                jv_fastq_an2 = glob.glob(in_dir + '/*all_j*' + an2)[0]
            UMI_seqlogo.umi_seq_logo(jv_fastq_an1, dedup.out_dir + '/' + dedup.jv_prefix + '_' + an1 + '.eps')
            UMI_seqlogo.umi_seq_logo(jv_fastq_an2, dedup.out_dir + '/' + dedup.jv_prefix + '_' + an2 + '.eps')

    #reverse complement for partis (needs to be VDJ, default output is JDV)
    elif opts.rev_comp:
        #if file path is provided
        if opts.input.endswith('fastq'):
            print('Reverse complementing:', opts.input)
            rev_comp_fq(os.path.abspath(opts.input), fq=True)
        elif opts.input.endswith('fasta'):
            print('Reverse complementing:', opts.input)
            rev_comp_fq(os.path.abspath(opts.input), fq=False)
        #if only directory provided
        else:
            if opts.fq:
                jv_fnames = glob.glob(opts.input + '/*.fastq')
            else:
                jv_fnames = glob.glob(opts.input + '/*.fasta')
            for jv_fname in jv_fnames:
                print('Reverse complementing:', jv_fname)
                rev_comp_fq(os.path.abspath(jv_fname), fq=opts.fq)

    else:

        logging.basicConfig(level=logging.DEBUG, filename=dedup.out_dir +'/' + dedup.jv_prefix + '.log', filemode='a+',
                            format='%(asctime)-15s %(levelname)-8s %(message)s')

        logging.info(opts)
        print('Starting deduplication')
        dedup.deduplicate_reads(threshold=opts.threshold, min_reads=opts.minreads, threads=opts.nthreads,
                                            mismatch=opts.mismatch, gt_threshold=opts.gtratio,
                                            j_trim=opts.j_trim, stats=opts.stats,
                                            ignore_umi=opts.ignore_umi, use_j=opts.use_j,
                                            skip_unclear=opts.skip_unclear, keep_mh=opts.keep_mh, no_msa=opts.no_msa,
                                            umi_cor=opts.umi_correction,
                                            no_anchor=opts.no_anchor, short=opts.short, fq=opts.fq,
                                            cons_no_qual=opts.cons_no_qual, with_N=opts.with_N)




################################################################################
#Aggregate reads and bundle with read name seq count
################################################################################


def make_bundle(fastq, ignore_umi, use_j, skip_unclear, keep_mh, no_anchor=False, short=False):
    '''bundle reads
    '''
    unclear_skip = 0
    #Deduplication without alignment
    reads_dict = defaultdict(lambda: defaultdict(dict))
    qual_dict = defaultdict(list)

    with general.file_open(fastq) as jv_in:
        # lines = jv_in.read().splitlines()
        for item in general.fastq_parse(jv_in):
            qname = item[0]
            seq = item[1]
            thrd = item[2]
            qual = item[3]

            #Get UMI from qname
            if ignore_umi:
                umi = ''
            else:
                umi = qname.split(' ')[0].split('_')[-1]

            if no_anchor:
                anchor = ''
            else:
                anchor = qname.split(' ')[0].split('_')[-2]


            if skip_unclear:
                if 'unclear' in qname:
                    unclear_skip += 1
                    continue
            elif not keep_mh:
                if 'unclear' in qname and '-mh_' in qname:
                    unclear_skip += 1
                    continue

            if short and anchor:
                j_idn = qname.split('_')[-4]
            else:
                j_idn = qname.split('_')[-3]


            if use_j:
                key = j_idn + anchor
            else:
                key = anchor


            #create dictionary of sequence and quality...
            qual_dict[seq].append(qual)

            try:
                reads_dict[key][umi]['count'] += 1
                reads_dict[key][umi]['seq'].update([seq]) #add all the seqs for consensus

            except KeyError:
                reads_dict[key][umi]['count'] = 1
                reads_dict[key][umi]['read'] = qname.split(' ')[0] #+ ' ' + v_seq #put v_seq into qname for stats
                reads_dict[key][umi]['seq'] = Counter([seq]) #add all the seqs for consensus

        #if same sequence has 2 quals take highest

        for k,v in qual_dict.items():
            if len(v) > 1:
                qual_lofls = [list(item) for item in v]
                qual_dict[k] = [deduplication_general.qual_highest(qual_lofls)]

    return (reads_dict, unclear_skip, qual_dict)





def reduce_clusters_worker(bundle, clusters, counts, j_trim, no_msa, qual_dict,
                           gt_threshold, cons_no_qual, short, mismtch, with_N):

    reads = []
    consensus_seqs = []
    consensus_quals = []
    final_umis = []
    umi_counts = []
    low_gt = 0
    corrected = 0
    low_gt_corrected = 0

    cons_diffs = defaultdict()
    cons_algn = defaultdict()

    for cluster in clusters:
        umi_in_cluster = len(cluster)
        #Consensus for read loss filter
        out_dict = {umi:bundle[umi]['seq'] for umi in cluster}

        assert len(out_dict) > 0, 'No sequence from umi'

        if umi_in_cluster > 1: #contains umi with 1 error
            corrected += 1

        if no_msa:
            alignment = ''
            gt_ratio, consensus_seq, consensus_qual, diffs_from_cons = deduplication_general.read_loss(out_dict, qual_dict, differences=mismtch,
                                                                                short=short, no_msa=no_msa, cons_no_qual=cons_no_qual,
                                                                                j_trim=j_trim, with_N=with_N) #umi=umi_cons
        else:
            alignment, new_qual_dict = deduplication_general.kalign_msa(out_dict, qual_dict)
            gt_ratio, consensus_seq, consensus_qual, diffs_from_cons = deduplication_general.read_loss(alignment, new_qual_dict, differences=mismtch,
                                                                                 no_msa=no_msa, cons_no_qual=cons_no_qual,
                                                                                 j_trim=j_trim, with_N=with_N) #umi=umi_cons
        #keep record of the distance between sequence and consensus per UMI bases (clustered UMIs seperated by ,)
        if not isinstance(diffs_from_cons, int):
            cons_algn[','.join(cluster)] = ','.join(x for x in alignment)
            cons_diffs[','.join(cluster)] = ','.join(x for x in diffs_from_cons)
        else:
            cons_algn[','.join(cluster)] = diffs_from_cons
            cons_diffs[','.join(cluster)] = diffs_from_cons


        if gt_ratio >= gt_threshold:
            #Parent umi = highest count umi which account for the cluster
            parent_umi = deduplication_general.get_best_higher_counts(cluster, counts)
            reads.append(bundle[parent_umi]['read'])
            consensus_seqs.append(consensus_seq.replace('-', '')) #remove padding or indels from msa
            consensus_quals.append(consensus_qual.replace('#', '')) #should not have any # qual as preclean removed them

            final_umis.append(parent_umi)
            #Number of UMI's in the cluster (how many have been collapsed)
            umi_counts.append(sum([counts[x] for x in cluster]))
        else:
            low_gt += 1

            if umi_in_cluster > 1: #contains umi with 1 error and low ratio
                low_gt_corrected += 1

    return [reads, consensus_seqs, consensus_quals, final_umis, umi_counts, low_gt,
            corrected, low_gt_corrected, cons_diffs, cons_algn]




################################################################################
#Deduplication in parallel
################################################################################


class results():
    '''Holder for results
    '''
    def __init__(self):

                self.reads = []
                self.consensus_seqs = []
                self.final_umis = []
                self.umi_counts = []
                self.low_gt = 0
                self.corrected = 0
                self.low_gt_corrected = 0
                self.num_input = 0
                self.stats_pre_df_dict = {'UMI': [], 'counts': []}
                self.cons_diffs = defaultdict()
                self.cons_algn = defaultdict()
                self.cons_all = defaultdict(lambda: defaultdict())
                self.consensus_quals = []



#TODO need to update
def deduplication_worker_umi(bundle, threshold, stats, mismatch, gt_threshold, j_trim,
                             no_msa, short, qual_dict, cons_no_qual, with_N):
    ''' worker for deduplicate_bundle_parallel '''


    reads, consensus_seqs, consensus_quals, final_umis, umi_counts, low_gt, corrected, low_gt_corrected, cons_diffs, cons_algn =\
    umi_correction.run_dir_adj(bundle, threshold=threshold, mismatches=mismatch,
                gt_threshold=gt_threshold, qual_dict=qual_dict, no_msa=no_msa, short=short,
                cons_no_qual=cons_no_qual, with_N=with_N, j_trim=j_trim)


    num_input = sum([bundle[umi]['count'] for umi in bundle])
    # collect pre-dudupe stats
    stats_pre_df_dict = {'UMI': [], 'counts': []}
    # pre_average_distance = ''
    if stats:
        stats_pre_df_dict['UMI'].extend(bundle) #umi + read
        stats_pre_df_dict['counts'].extend([bundle[UMI]['count'] for UMI in bundle]) #umi counts


    return [reads, consensus_seqs, final_umis, umi_counts, low_gt, corrected, low_gt_corrected,
            num_input, stats_pre_df_dict, cons_diffs, cons_algn, consensus_quals] #, pre_average_distance]



def deduplication_worker(bundle, threshold, stats, mismatch, gt_threshold, j_trim,
                         no_msa, short, qual_dict, cons_no_qual, with_N):
    ''' worker for deduplicate_bundle_parallel without umi correction'''

    rc_results = results()

    #reduce bundle
    umis = bundle.keys()

    #clusters need to be a list of sets
    clusters = []
    for key in bundle.keys():
        clusters.append(set([key]))

    counts = {umi: bundle[umi]['count'] for umi in umis}
    # bundle, clusters, counts, stats, mismtch, gt_threshold, no_msa=False

    reads, consensus_seqs, consensus_quals, final_umis, umi_counts, low_gt, corrected, low_gt_corrected, cons_diffs, cons_algn =\
    deduplication_general.reduce_clusters_single(bundle, clusters, counts, mismatch, gt_threshold,
                           qual_dict, j_trim, no_msa, short, cons_no_qual, with_N)


    num_input = sum([bundle[umi]['count'] for umi in bundle])
    # collect pre-dudupe stats
    stats_pre_df_dict = {'UMI': [], 'counts': []}
    # pre_average_distance = ''
    if stats:
        stats_pre_df_dict['UMI'].extend(bundle) #umi + read
        stats_pre_df_dict['counts'].extend([bundle[UMI]['count'] for UMI in bundle]) #umi counts

        # pre_average_distance = get_average_umi_distance(bundle.keys()) #v_seq + umi
    rc_results.reads = reads
    rc_results.consensus_seqs = consensus_seqs
    rc_results.final_umis = final_umis
    rc_results.umi_counts = umi_counts
    rc_results.low_gt = low_gt
    rc_results.corrected = corrected
    rc_results.low_gt_corrected = low_gt_corrected
    rc_results.num_input = num_input
    rc_results.stats_pre_df_dict = stats_pre_df_dict
    rc_results.cons_diffs = cons_diffs
    rc_results.cons_algn = cons_algn
    rc_results.consensus_quals = consensus_quals

    return [rc_results] #, pre_average_distance]
    # return [reads, consensus_seqs, final_umis, umi_counts, low_gt, corrected,
    #low_gt_corrected, num_input, stats_pre_df_dict, cons_diffs, consensus_quals] #, pre_average_distance]


def chunk_it(seq, num):
  avg = len(seq) / float(num)
  out = []
  last = 0.0

  while last < len(seq):
    out.append(seq[int(last):int(last + avg)])
    last += avg

  return out



def deduplicate_bundle_parallel(reads_dict, qual_dict, threshold,
                    mismatch, stats, threads, gt_threshold, j_trim,
                    no_msa, umi_correction, no_anchor=False, short=False,
                    cons_no_qual=False, use_j=False, with_N=False):
    '''
    Deduplicated reads in parralel with either UMI correction or without
    '''
    #restrict number of threads to number of bundles + subprocess that will be spawned
    num_bundles = len(reads_dict.values())
    if num_bundles > threads:
        use_threads = threads
    elif num_bundles <= threads:
        use_threads = num_bundles

    if umi_correction:

        dir_adj_results = []

        #will return a list of results (two lists that need to be changed into results format)
        dir_adj_bundles = Parallel(n_jobs=use_threads)(delayed(deduplication_worker_umi)(bundle, threshold, stats,
        mismatch, gt_threshold, j_trim, no_msa, short, qual_dict, cons_no_qual, with_N) for bundle in reads_dict.values())

        for i in range(num_bundles):
            reads, consensus_seqs, final_umis, umi_counts, low_gt, corrected, low_gt_corrected, \
            num_input, stats_pre_df_dict, cons_diffs, cons_algn, consensus_quals = dir_adj_bundles[i]

            rc_results = results()

            rc_results.reads = reads
            rc_results.consensus_seqs = consensus_seqs
            rc_results.final_umis = final_umis
            rc_results.umi_counts = umi_counts
            rc_results.low_gt = low_gt
            rc_results.corrected = corrected
            rc_results.low_gt_corrected = low_gt_corrected
            rc_results.num_input = num_input
            rc_results.stats_pre_df_dict = stats_pre_df_dict
            rc_results.cons_diffs = cons_diffs
            rc_results.cons_algn = cons_algn
            rc_results.consensus_quals = consensus_quals
            # print(dict(itertools.islice(cons_algn.items(), 0, 5)))
            #merge UMI cluster sequence differences and MSA alignments into one object a dict of dicts
            for k,v in cons_diffs.items():
                try:
                    rc_results.cons_all[k]['diffs_from_cons'] += '_' + v

                except KeyError:
                    rc_results.cons_all[k]['diffs_from_cons'] = v


            for k,v in cons_algn.items():
                try:
                    rc_results.cons_all[k]['alignments'] += '_' + v
                except KeyError:
                    rc_results.cons_all[k]['alignments'] = v

            dir_adj_results.append(rc_results)


    else:

        if no_anchor:
            if not no_msa:
                use_threads = use_threads/2


        dir_adj_results = []

        #dict within a dict
        for bundle in tqdm(reads_dict.values()): #only a single bundle present
            #do in parallel
            #reduce bundle
            umis = bundle.keys()

            bundle_results = results()
            #clusters need to be a list of sets
            clusters = []
            for key in bundle.keys():
                clusters.append(set([key]))

            list_of_clusters = chunk_it(clusters, threads)

            counts = {umi: bundle[umi]['count'] for umi in umis}

            #num_input
            bundle_results.num_input += sum([bundle[umi]['count'] for umi in bundle])

            # pre_average_distance = ''
            if stats:
                bundle_results.stats_pre_df_dict['UMI'].extend(bundle) #umi + read
                bundle_results.stats_pre_df_dict['counts'].extend([bundle[UMI]['count'] for UMI in bundle]) #umi counts
            # if stats:
            #     dir_adj_results[0][8]['UMI'].extend(bundle) #umi + read
            #     dir_adj_results[0][8]['counts'].extend([bundle[UMI]['count'] for UMI in bundle]) #umi counts

            dir_adj_results_lists = \
            Parallel(n_jobs=threads)(delayed(reduce_clusters_worker)(bundle, clusters, counts, j_trim,
            no_msa, qual_dict, gt_threshold, cons_no_qual, short, mismatch, with_N) for clusters in list_of_clusters)

            for reads_s, consensus_seqs_s, consensus_quals_s, final_umis_s, umi_counts_s, low_gt_s, \
            corrected_s, low_gt_corrected_s, cons_diffs_s, cons_algn_s in dir_adj_results_lists:

                bundle_results.reads.extend(reads_s)
                bundle_results.consensus_seqs.extend(consensus_seqs_s)
                bundle_results.consensus_quals.extend(consensus_quals_s)
                bundle_results.final_umis.extend(final_umis_s)
                bundle_results.umi_counts.extend(umi_counts_s)
                bundle_results.low_gt += low_gt_s
                bundle_results.corrected += corrected_s
                bundle_results.low_gt += low_gt_corrected_s
                for k,v in cons_diffs_s.items():
                    try:
                        bundle_results.cons_all[k]['diffs_from_cons'] += '_' + v

                    except KeyError:
                        bundle_results.cons_all[k]['diffs_from_cons'] = v


                for k,v in cons_algn_s.items():
                    try:
                        bundle_results.cons_all[k]['alignments'] += '_' + v
                    except KeyError:
                        bundle_results.cons_all[k]['alignments'] = v

            dir_adj_results.append(bundle_results)

    return dir_adj_results




def write_out_deduplicated(dir_adj_results, low_umi_out, out, stats, min_reads,
                           no_anchor, fq, pdf_out):
    '''
    Write out reads from result object produced from deduplication and
    return stats
    '''
    stats_pre_df_dict_all = {'UMI': [], 'counts': []}
    stats_post_df_dict = {'UMI': [], 'counts': []}
    pre_cluster_stats = []
    post_cluster_stats = []
    stats_cons_diffs = defaultdict(lambda: defaultdict())


    num_input_all, num_output = 0, 0

    low_gt_reads = 0
    corrected_reads = 0
    low_gt_corrected_reads = 0
    low_umi_count = 0


    print('Writing out')

    for bundle in range(len(dir_adj_results)):
        # reads, consensus, final_umis, umi_counts, low_gt, corrected, low_gt_corrected,\
        # topologies, nodes, num_input, stats_pre_df_dict, pre_average_distance = results
        # print('Unique:', len(set(dir_adj_results[bundle][0])), 'All:', len(dir_adj_results[bundle][0]))
        #umi_counts
        #in case no reads are output
        if dir_adj_results[bundle].umi_counts:
            #*-operator to unpack the arguments out of a list or tuple #dir_adj_results[bundle][3]
            labels, values = zip(*Counter(dir_adj_results[bundle].umi_counts).items())
            non_dedup_values = tuple(l*v for l, v in zip(labels, values))

            if min_reads != None:
                # cut_point = count_change(labels, non_dedup_values)

                cut_point = min_reads

                plt.figure()
                plt.bar(labels, non_dedup_values)
                #extract title from read
                # j_nam = re.search('J\d+', dir_adj_results[bundle][0][0].split(' ')[0]).group(0)
                j_nam = dir_adj_results[bundle].reads[0].split(' ')[0].split('_')[-3] #dir_adj_results[bundle][0][0]
                if no_anchor:
                    anchor = ''
                else:
                    anchor = dir_adj_results[bundle].reads[0].split(' ')[0].split('_')[-2] #dir_adj_results[bundle][0][0]

                plt.title(j_nam + ' ' + anchor + ' Cut point: ' + str(cut_point), ha='center') #need to put name to know which bundle J is being processed
                my_plot = plt.axvline(cut_point, linestyle='dashed', linewidth=2).get_figure()
                pdf_out.savefig(my_plot)
                plt.close('all')

        num_input_all += dir_adj_results[bundle].num_input #dir_adj_results[bundle][7] #num_input

        #remove low umi counts 1-5
        indx = 0
        for count in dir_adj_results[bundle].umi_counts: #dir_adj_results[bundle][3]: #umi_counts
            if count <= cut_point:
                if fq:
                    assert len(dir_adj_results[bundle].consensus_seqs[indx]) == len(dir_adj_results[bundle].consensus_quals[indx]), \
                    'Consensus sequence not same length as consensus quality'
                    # dir_adj_results[bundle][1]    dir_adj_results[bundle][10]

                    low_umi_out.write(dir_adj_results[bundle].reads[indx].split(' ')[0] + '_' + str(count) +'\n' +
                                      dir_adj_results[bundle].consensus_seqs[indx] + '\n' +
                                      '+' + '\n' +
                                      dir_adj_results[bundle].consensus_quals[indx] + '\n')
                else:
                    #write out fasta
                    low_umi_out.write(dir_adj_results[bundle].reads[indx].split(' ')[0].replace('@', '>') + '_' + str(count) +'\n' +
                    dir_adj_results[bundle].consensus_seqs[indx] + '\n')
                low_umi_count += 1
            else:
                if fq:
                    assert len(dir_adj_results[bundle].consensus_seqs[indx]) == len(dir_adj_results[bundle].consensus_quals[indx]), \
                    'Consensus sequence not same length as consensus quality'

                    out.write(dir_adj_results[bundle].reads[indx].split(' ')[0] + '_' + str(count) + '\n' +
                              dir_adj_results[bundle].consensus_seqs[indx] + '\n' +
                              '+' + '\n' +
                              dir_adj_results[bundle].consensus_quals[indx] + '\n')
                else:
                    #write out fasta
                    out.write(dir_adj_results[bundle].reads[indx].split(' ')[0].replace('@', '>') + '_' + str(count) + '\n' +
                    dir_adj_results[bundle].consensus_seqs[indx] + '\n')
                num_output += 1
            indx += 1


        if stats:
            # assert len(reads)  == len(consensus) == len(umi_counts), 'Reads, consensus and counts differ in length'

            ##################

            low_gt_reads += dir_adj_results[bundle].low_gt #dir_adj_results[bundle][4] #low_gt
            corrected_reads +=  dir_adj_results[bundle].corrected #dir_adj_results[bundle][5] #corrected
            low_gt_corrected_reads += dir_adj_results[bundle].low_gt_corrected #dir_adj_results[bundle][6] #low_gt_corrected

            stats_pre_df_dict_all['UMI'].extend(dir_adj_results[bundle].stats_pre_df_dict['UMI']) #stats_pre_df_dict umi + read #dir_adj_results[bundle][8]
            stats_pre_df_dict_all['counts'].extend(dir_adj_results[bundle].stats_pre_df_dict['counts']) #stats_pre_df_dict umi counts

            # pre_cluster_stats.append(pre_average_distance)

            #aggregate errors per cluster for each bundle
            cons_diffs = dir_adj_results[bundle].cons_all #dir_adj_results[bundle][9]
            for k,v in cons_diffs.items():
                try:
                    # print(stats_cons_diffs[k], '_', v)
                    stats_cons_diffs[k]['diffs_from_cons'] += '_' + v['diffs_from_cons']
                    stats_cons_diffs[k]['alignments'] += '_' + v['alignments']
                except KeyError:
                    stats_cons_diffs[k]['diffs_from_cons'] = v['diffs_from_cons']
                    stats_cons_diffs[k]['alignments'] = v['alignments']

            # collect post-dudupe stats
            #v_seq + umi
            # post_cluster_umis = [qname.split(' ')[-1] for qname in dir_adj_results[bundle][0]] #reads are just qnames
            stats_post_df_dict['UMI'].extend(dir_adj_results[bundle].final_umis) #final_umis #dir_adj_results[bundle][2]
            stats_post_df_dict['counts'].extend(dir_adj_results[bundle].umi_counts) #umi_counts #dir_adj_results[bundle][3]


    return [stats_pre_df_dict_all, stats_post_df_dict, pre_cluster_stats, post_cluster_stats,
    num_input_all, num_output, low_gt_reads, corrected_reads, low_gt_corrected_reads, low_umi_count,
    stats_cons_diffs]




def aggregate_stats_df(stats_df):
    ''' return a data frame with aggregated counts per UMI
    '''

    total_counts = stats_df.pivot_table(
        columns="UMI", values="counts", aggfunc=np.sum).T
    total_counts.rename(columns={'counts':'total_counts'}, inplace=True)
    median_counts = stats_df.pivot_table(
        columns="UMI", values="counts", aggfunc=np.median).T
    median_counts.rename(columns={'counts':'median_counts'}, inplace=True)
    times_observed = stats_df.pivot_table(
        columns="UMI", values="counts", aggfunc=len).T
    times_observed.rename(columns={'counts':'times_observed'}, inplace=True)

    agg_df = pd.concat([total_counts, median_counts, times_observed], axis=1)

    return agg_df


################################################################################
#Deduplication run class
################################################################################

class deduplicate:
    '''Deduplicate using J, V start and UMI
    anchor 1 GACTCGT  anchor 2 CTGCTCCT
    '''

    def __init__(self, file_directory, an1, an2):
        '''
        :param file_directory: where files are
        :param single: seperate V and J (not assembled)
        '''
        self.an1 = an1
        self.an2 = an2

        if os.path.isfile(file_directory):
            in_dir = os.path.basename(file_directory)
            self.file_directory = in_dir
            self.jv_fastq_an1 = file_directory
        #if only directory provided
        elif os.path.isdir(file_directory):
            in_dir = file_directory
            self.file_directory = in_dir

            try:
                self.jv_fastq_an1 = glob.glob(self.file_directory + '/*all_j*' + an1)[0]
                self.jv_fastq_an2 = glob.glob(self.file_directory + '/*all_j*' + an2)[0]
            except IndexError:
                raise Exception('No input files found')


        self.out_dir = ''

        self.jv_prefix = ''
        self.header_postfix_jv = ''


    def create_dirs_assembled(self, out_dir=None):
        '''Create directories
        '''

        dir_main = Path(os.path.abspath(self.jv_fastq_an1)).parents[1] #1 dir up, create outside of preclean directory
        self.jv_prefix = re.split('(_all)', os.path.basename(self.jv_fastq_an1))[0]


        #final output dir
        if out_dir == None:
            self.out_dir = str(dir_main) + '/' + str(self.jv_prefix) + '_Deduplicated'
        else:
            self.out_dir = os.path.abspath(out_dir)

        try:
            os.mkdir(self.out_dir)
        except FileExistsError:
            print('Directory', self.out_dir, 'already exists')
            pass

        #Get 1:N:0:GCCAAT from header
        with general.file_open(self.jv_fastq_an1) as f:
                first_line = f.readline()
                self.header_postfix_jv = first_line.decode('utf-8').split(' ')[1].rstrip()




    def deduplicate_reads(self, threshold, min_reads, threads, mismatch, gt_threshold, j_trim,
                                      stats=False, ignore_umi=False, use_j=False, skip_unclear=False,
                                      keep_mh=False, no_msa=False, umi_cor=False, no_anchor=False,
                                      short=False, fq=False, cons_no_qual=False,
                                      with_N=False):
        '''Main deduplication function
        First bundles are created from input reads_s
        Second deduplication function is run
        Throughout logs are prepared
        '''
        if no_anchor:
            reads_dict, unclear_skip_an1, qual_dict = make_bundle(self.jv_fastq_an1,
                                                            ignore_umi=ignore_umi, use_j=use_j,
                                                            skip_unclear=skip_unclear, keep_mh=keep_mh,
                                                            no_anchor=no_anchor, short=short)
            unclear_skip_an2 = 0

        else:

            reads_dict_an1, unclear_skip_an1, qual_dict_an1 = make_bundle(self.jv_fastq_an1,
                                                            ignore_umi=ignore_umi, use_j=use_j,
                                                            skip_unclear=skip_unclear, keep_mh=keep_mh, short=short)
            reads_dict_an2, unclear_skip_an2, qual_dict_an2 = make_bundle(self.jv_fastq_an2,
                                                            ignore_umi=ignore_umi, use_j=use_j,
                                                            skip_unclear=skip_unclear, keep_mh=keep_mh, short=short)

            # merge the two dict_keys
            reads_dict = {**reads_dict_an1, **reads_dict_an2}

            #second dict values will overwrite those from the first
            # qual_dict = {**qual_dict_an1, **qual_dict_an2}

            qual_dict = defaultdict(list)
            for k,v in itertools.chain(qual_dict_an1.items(), qual_dict_an2.items()):
                qual_dict[k].extend(v)

            #keep only highest quality for duplicates
            for k,v in qual_dict.items():
                if len(v) > 1:
                    qual_lofls = [list(item) for item in v]
                    qual_dict[k] = [deduplication_general.qual_highest(qual_lofls)]



        ########################

        if stats:
            print('Unclear skipped:', unclear_skip_an1 + unclear_skip_an2)
            all_unclear =  unclear_skip_an1 + unclear_skip_an2
            logging.info('Unclear skipped:' +  str(all_unclear))

            # set up arrays to hold stats data
            stats_pre_df_dict_all = {'UMI': [], 'counts': []}
            stats_post_df_dict_all = {'UMI': [], 'counts': []}
            pre_cluster_stats_all = []
            post_cluster_stats_all = []


            num_input_an1, num_output_an1 = 0, 0
            num_input_an2, num_output_an2 = 0, 0
            # line_fmt = "@{0!s}\n{1!s}\n+\n{2!s}\n"
            low_gt_reads_an1, low_gt_reads_an2 = 0, 0
            corrected_reads_an1, corrected_reads_an2 = 0, 0
            low_gt_corrected_reads_an1, low_gt_corrected_reads_an2 = 0, 0
            low_umi_count_an1, low_umi_count_an2 = 0, 0

        #Write both anchors into same file
        #Can't split into DJ and V

        if fq:
            out_name = self.out_dir + '/' + self.jv_prefix + '_dedup.fastq'
            low_name = self.out_dir + '/' + self.jv_prefix + '_low_umi.fastq'
        else:
            out_name = self.out_dir + '/' + self.jv_prefix + '_dedup.fasta'
            low_name = self.out_dir + '/' + self.jv_prefix + '_low_umi.fasta'

        with open(out_name, 'w') as jv_out, open(low_name, 'w') as low_umi_out, \
        PdfPages(self.out_dir + '/' + self.jv_prefix + '_histogram.pdf') as pdf:

            #run an1 and an2 side by side (not tested!)
            # if len(reads_dict) >= threads:
            #     nprocs = 1
            # else:
            #     nprocs = int(threads/len(reads_dict)) #how many unused cores are available?

            #an1+2

            deduplication_results =\
            deduplicate_bundle_parallel(reads_dict, qual_dict, threshold=threshold,
                            mismatch=mismatch, gt_threshold=gt_threshold,
                            stats=stats, threads=threads, j_trim=j_trim, no_msa=no_msa,
                            umi_correction=umi_cor, no_anchor=no_anchor, short=short, cons_no_qual=cons_no_qual,
                            with_N=with_N)

            stats_pre_df_dict, stats_post_df_dict, pre_cluster_stats, post_cluster_stats, \
            num_input, num_output, low_gt_reads, corrected_reads, \
            low_gt_corrected_reads, low_umi_count, stats_cons_diffs=\
            write_out_deduplicated(deduplication_results, low_umi_out, jv_out, stats=stats, min_reads=min_reads,
                                   no_anchor=no_anchor, fq=fq, pdf_out=pdf)

            #stats
            if stats:
                num_input_an1 += num_input
                num_output_an1 += num_output
                low_gt_reads_an1 += low_gt_reads
                corrected_reads_an1 += corrected_reads
                low_gt_corrected_reads_an1 += low_gt_corrected_reads
                low_umi_count_an1 += low_umi_count

                stats_pre_df_dict_all.update(stats_pre_df_dict)
                stats_post_df_dict_all.update(stats_post_df_dict)

                pre_cluster_stats_all.extend(pre_cluster_stats)
                post_cluster_stats_all.extend(post_cluster_stats)




        if stats:
            print('Number of input reads:', num_input_an1)
            print('Number of output reads:', num_output_an1)
            logging.info('Number of input reads:' + str(num_input_an1))
            logging.info('Number of output reads:' + str(num_output_an1))

            print('Number of clusters with low ratio discarded:' + str(low_gt_reads_an1))
            logging.info('Number of clusters with low ratio discarded:' + str(low_gt_reads_an1))

            if umi_cor:
                print('Number of directional-adjacency corrected clusters:', corrected_reads_an1)
                logging.info('Number of directional-adjacency corrected clusters:' + str(corrected_reads_an1))
                print('Number of corrected clusters with low ratio discarded:', low_gt_corrected_reads_an1)
                logging.info('Number of corrected clusters with low ratio discarded:' +  str(low_gt_corrected_reads_an1))

            print('Number of low UMI count groups:', low_umi_count_an1)
            logging.info('Number of low UMI count groups:' + str(low_umi_count_an1))


            ##########################################

            #From UMI_tools
            stats_pre_df = pd.DataFrame(stats_pre_df_dict_all)
            stats_post_df = pd.DataFrame(stats_post_df_dict_all)

            # print(pd.DataFrame.from_dict(Counter(stats_post_df_dict['counts']), orient='index').reset_index())

            # generate histograms of counts per UMI at each position
            UMI_counts_df_pre = pd.DataFrame(stats_pre_df.pivot_table(
                columns=stats_pre_df['counts'], values='counts', aggfunc=len)).T #not sure why it need to be transposed now when it didnt before!

            UMI_counts_df_post = pd.DataFrame(stats_post_df.pivot_table(
                columns=stats_post_df['counts'], values='counts', aggfunc=len)).T

            # UMI_counts_df_pre.columns = ['instances']
            # UMI_counts_df_post.columns = ['instances']
            UMI_counts_df_pre.rename(columns={'counts':'instances'}, inplace=True)
            UMI_counts_df_post.rename(columns={'counts':'instances'}, inplace=True)

            # print(stats_pre_df.pivot_table(columns="UMI", values="counts", aggfunc=np.sum))

            UMI_counts_df = pd.merge(UMI_counts_df_pre, UMI_counts_df_post,
                                     how='outer', left_index=True, right_index=True,
                                     sort=True, suffixes=['_pre', '_post'])

            UMI_counts_df = UMI_counts_df.fillna(0).astype(int)

            UMI_counts_df.to_csv(self.out_dir + '/' + self.jv_prefix + '_per_umi_per_position.tsv', sep='\t')

            ##########################

             # aggregate stats pre/post per UMI
            agg_pre_df = aggregate_stats_df(stats_pre_df)
            agg_post_df = aggregate_stats_df(stats_post_df)

            agg_df = pd.merge(agg_pre_df, agg_post_df, how='left',
                              left_index=True, right_index=True,
                              sort=True, suffixes=['_pre', '_post'])

            agg_df = agg_df.fillna(0).astype(int)

            # stats_consensus_difference = pd.DataFrame(stats_cons_diffs, index=[0])
            stats_consensus_difference = pd.DataFrame(stats_cons_diffs)
            stats_consensus_difference = stats_consensus_difference.T
            stats_consensus_difference.columns = ['Alignments', 'Consensus_differences']

            # stats_consensus_difference['UMI'] = stats_consensus_difference.index

            agg_df = pd.merge(agg_df, stats_consensus_difference, how='left',
                              left_index=True, right_index=True,
                              sort=True,)

            agg_df.to_csv(self.out_dir + '/' + self.jv_prefix + '_per_umi.tsv', sep="\t")



            # stats_consensus_difference.to_csv(self.out_dir + '/' + self.jv_prefix + '_consensus_difference.tsv', sep="\t")


def rev_comp_fq(path, fq):
    '''Reverse complement a fastx file

    Args:
        path (str): path to fastx file
        fq (logical): is the supplied file a fastq file?

    Example:
        path = '/tset/sava/af.aef.fastq.gz'

    Return:
        Writes out a fastx file with _rv.fastx appended in the same directory as the
        input path.
    '''
    if path.endswith('.gz'):
        *fname_path, fname_end, fname_gz = path.split('.')
    else:
        *fname_path, fname_end = path.split('.')
    with general.file_open(path) as in_fname, open('.'.join(fname_path) + '_rv.' + fname_end, 'w') as out_fname:
        if fq:
            for qname, seq, thrd, qual in general.fastq_parse(in_fname):
                rv_seq = general.reverse_complement(seq)
                out_fname.write(qname + '\n' + rv_seq + '\n' + thrd + '\n' + qual + '\n')
        else:
            for qname, seq in general.fasta_iter(path):
                rv_seq = general.reverse_complement(seq)
                out_fname.write('>' + qname + '\n' + rv_seq + '\n')




if __name__ == "__main__":
    main()
