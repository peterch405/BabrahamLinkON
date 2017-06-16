#!/usr/bin/env python3

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

from babrahamlinkon import general, presets, umi_correction, UMI_seqlogo
# from skbio import DNA, TabularMSA
import multiprocessing
import math

import itertools
from joblib import Parallel, delayed
import logging
# import pickle
# import warnings
import operator
import Levenshtein

import pyximport
from babrahamlinkon._dedup_umi import edit_distance



################################################################################
#Aggregate reads and bundle with read name seq count
################################################################################


def make_bundle(fastq, v_len, j_len, ignore_umi, ignore_j, ignore_v, skip_unclear, skip_mh, no_anchor=False, short=False):
    '''bundle reads
    '''
    unclear_skip = 0
    #Deduplication without alignment
    reads_dict = defaultdict(lambda: defaultdict(dict))


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

            if short:
                #Only J seq present, trim all to 50bp and take 8bp from there
                if v_len > 0:
                    v_seq = seq[:50][-v_len:]
            else:
                if v_len > 0:
                    v_seq = seq[-v_len:]

            if skip_unclear:
                if 'unclear' in qname:
                    unclear_skip += 1
                    continue
            elif skip_mh:
                if 'unclear' in qname and '-mh_' in qname:
                    unclear_skip += 1
                    continue

            if short and anchor:
                j_idn = qname.split('_')[-4]
            else:
                j_idn = qname.split('_')[-3]


            if ignore_j:
                key = anchor
            else:
                key = j_idn + anchor


            if ignore_v:
                dedup_seq = umi
            elif int(j_len) > 0: #if more than 0
                j_seq = seq[:50][-j_len:]
                dedup_seq = j_seq + umi
            else:
                dedup_seq = v_seq + umi

            #TODO: include quality for consensus

            try:
                reads_dict[key][dedup_seq]['count'] += 1
                reads_dict[key][dedup_seq]['seq'].update([seq]) #add all the seqs for consensus

            except KeyError:
                reads_dict[key][dedup_seq]['count'] = 1
                reads_dict[key][dedup_seq]['read'] = qname.split(' ')[0] #+ ' ' + v_seq #put v_seq into qname for stats
                reads_dict[key][dedup_seq]['seq'] = Counter([seq]) #add all the seqs for consensus


    return (reads_dict, unclear_skip)



################################################################################
#consensus and read loss
################################################################################

#TODO: use quality as an additional metric to resolve ambigious bases (highest quality base prevail)?
def consensus(list_of_lists):
    '''
    :param list_of_lists: string split into individual letters
                          [['C', 'A', 'C', 'A', 'T', 'A', 'T', 'A'],
                          ['G', 'A', 'T', 'A', 'T', 'A', 'T', 'A'],
                          ['G', 'A', 'T', 'A', 'T', 'A', 'T', 'A']]
    If base ambigous return N
    Else return most frequent base
    '''

    consensus_seq = ''
    # count = 0
    for pos in zip(*list_of_lists):
        base_counts = Counter(pos)
        most_freq = base_counts.most_common(2)
        if len(most_freq) > 1:
            if most_freq[0][1] == most_freq[1][1]:
                consensus_seq += 'N' #ambigous base
            else:
                consensus_seq += most_freq[0][0]
        else:
            consensus_seq += most_freq[0][0]

    return consensus_seq



def consensus_unequal(list_of_lists):
    '''
    :param list_of_lists: string split into individual letters
                          [['C', 'A', 'C', 'A', 'T', 'A', 'T', 'A'],
                          ['G', 'A', 'T', 'A', 'T', 'A', 'T', 'A'],
                          ['G', 'A', 'T', 'A', 'T', 'A', 'T', 'A']]
    If base ambigous return N
    Else return most frequent base
    '''

    consensus_seq = ''
    # count = 0
    for pos in itertools.zip_longest(*list_of_lists):
        base_counts = Counter(pos)
        most_freq = base_counts.most_common(2)
        if len(most_freq) > 1:
            if most_freq[0][1] == most_freq[1][1]:
                consensus_seq += 'N' #ambigous base
            elif most_freq[0][0] == None:
                #skip position
                continue
            else:
                consensus_seq += most_freq[0][0]
        elif most_freq[0][0] == None:
            continue
        else:
            consensus_seq += most_freq[0][0]

    return consensus_seq



def kalign_msa(seq_counter_dict): #, umi=None
    '''Read loss analysis
    :param dict seq_counter: dict of Counter object with sequences
    :param differences: number of differences from consensus allowed
    '''
    #Only aligning single copy of a duplicate sequence and then subsequently
    #multiplying the output alignment by Counter

    #Convert Counter into fasta input for kalign
    #Need to sort beacuse of Instability in progressive multiple sequence alignment algorithms?
    seq_fasta = ''
    count = 0
    reads = 0
    for umi, cntr in sorted(seq_counter_dict.items(), key=operator.itemgetter(0)):
        for seq, freq in sorted(cntr.items(), key=lambda x:x[0]):
            seq_fasta = seq_fasta + '>' + str(count) + '_' + umi + '_' + str(freq) + '\n' + seq + '\n'
            count += 1  #how many different reads
            reads += freq #how many total reads


    #Can't get consensus from single sequence, return single sequence
    if count == 1:
        assert len(list(seq_counter_dict.values())) == 1, 'Not a single sequence'
        seq_out= re.sub('\n', '', seq_fasta)
        seq_out = list(filter(None, re.split('(>\d+_[A-Z]+_\d+)', seq_out)))
        return seq_out

    #Multiple sequence alignment
    #http://msa.sbc.su.se/cgi-bin/msa.cgi
    kalign_cmd = ['kalign', '-f', 'fasta']


    p = subprocess.Popen(kalign_cmd, stdout=subprocess.PIPE, stdin=subprocess.PIPE, stderr=subprocess.STDOUT)
    kalign_stdout = p.communicate(input=seq_fasta.encode('utf-8'))[0]

    #parse kalign output and derive consensus seq
    head, sep, tail = kalign_stdout.partition(b'>') #remove head kalign intro

    alignment = sep.decode() + tail.decode()
    alignment = re.sub('\n', '', alignment)
    alignment = list(filter(None, re.split('(>\d+_[A-Z]+_\d+)', alignment)))

    # assert len(alignment) == count*2, 'Alignment output fewer reads'
    return alignment



def consensus_difference(seq_counter_dict, msa=False, short=False):
    '''
    :param alignment: output from msa
    :return: number of differences between two umi group consensus sequences
    '''

    if msa:
        seq_dict = defaultdict(list)

        for item in general.fasta_parse(seq_counter_dict):
            qname = item[0]
            seq = item[1]

            freq = int(qname.split('_')[-1]) #freq saved in name
            #split fasta into umi groups (head and child)
            umi = qname.split('_')[-2]
            for i in range(freq):
                seq_dict[umi].append(seq)

        assert len(seq_dict.values()) == 2, 'More than two groups'

        lst_1, lst_2 = list(seq_dict.values())

        lst_lists_1 = [list(item) for item in lst_1]
        lst_lists_2 = [list(item) for item in lst_2]

        cons_seq_1 = consensus(lst_lists_1)
        cons_seq_2 = consensus(lst_lists_2)

    else:

        cntr_1, cntr_2 = list(seq_counter_dict.values())

        lst_lists_1 = [list(item) for item in cntr_1.elements()]
        lst_lists_2 = [list(item) for item in cntr_2.elements()]

        cons_seq_1 = consensus_unequal(lst_lists_1)
        cons_seq_2 = consensus_unequal(lst_lists_2)

        if short: #need to trim sequence
            min_len = min(len(cons_seq_1),len(cons_seq_2))
            cons_seq_1 = cons_seq_1[:min_len]
            cons_seq_2 = cons_seq_2[:min_len]
        else:
            #Need to pad seq if length unequal!
            len_diff = len(cons_seq_1) - len(cons_seq_2)
            if len_diff < 0: #pad cons_seq_1
                cons_seq_1 = cons_seq_1 + '-'*abs(len_diff)
            elif len_diff > 0:
                cons_seq_2 = cons_seq_2 + '-'*len_diff

        assert len(cons_seq_1) == len(cons_seq_2), 'Sequences for hamming distance not same length!'

    num_diffs = edit_distance(cons_seq_1.encode('utf-8'), cons_seq_2.encode('utf-8'))

    return num_diffs





def read_loss(seq_counter_dict, differences=5, msa=False, short=False): #, umi=None
    '''Read loss analysis
    :param alignment: fasta from msa function
    :param differences: number of differences from consensus allowed
    '''
    if msa:
        seq_dict = defaultdict(list)

        for item in general.fasta_parse(seq_counter_dict):
            qname = item[0]
            seq = item[1]

            freq = int(qname.split('_')[-1]) #freq saved in name
            for i in range(freq):
                seq_dict[qname+str(i)] = seq #each seq has unique key

        #all values in list of lists
        seq_list_lsts = [list(item) for item in seq_dict.values()]

        cons_seq = consensus(seq_list_lsts)

        #Skip poor consensus (i.e. if building consensus from only two different seqs)
        if cons_seq.count('N') > differences: #5
            return(0, cons_seq, '0')
        if len(cons_seq) < 25:
            return(0, cons_seq, '0')

        good = 0
        total = 0
        diffs = len(cons_seq)
        best_seq = ''
        diffs_from_cons = []

        for item in general.fasta_parse(seq_counter_dict):
            qname = item[0]
            seq = item[1]

            freq = int(qname.split('_')[-1])
            total += 1

            assert len(seq) == len(cons_seq), 'Length of sequences into hamming distance unequal'
            #TODO: make the trimming more robust?
            consensus_diff = edit_distance(seq[25:].encode('utf-8'), cons_seq[25:].encode('utf-8')) #how many mismatches present
            diffs_from_cons.append(','.join([str(consensus_diff)]*freq))

            if consensus_diff < diffs: #need only if many N's in cons_seq
                diffs = consensus_diff
                best_seq = str(seq)
            if consensus_diff <= differences:
                good += 1


    else:

        list_of_lists = []
        for umi, seqs in seq_counter_dict.items():
            for seq in seqs.elements():
                list_of_lists.append(list(seq))

        # cons_seq = consensus(seq_list_lsts)
        cons_seq = consensus_unequal(list_of_lists)

        #Skip poor consensus (i.e. if building consensus from only two different seqs)
        if cons_seq.count('N') > differences: #5
            return(0, cons_seq, '0')
        if len(cons_seq) < 25:
            return(0, cons_seq, '0')

        good = 0
        total = 0
        diffs = cons_seq.count('N')+2 #some wiggle room
        best_seq = ''
        diffs_from_cons = []
        for umi, seqs in seq_counter_dict.items():
            for seq in seqs.elements():
                total += 1

                if short: #trim sequence to shortest one
                    min_len = min(len(seq),len(cons_seq))
                    cons_seq_el = cons_seq[:min_len]
                    seq = seq[:min_len]

                #Need to pad seq if length unequal!
                elif len(seq) != len(cons_seq):
                    len_diff = len(cons_seq) - len(seq)
                    if len_diff < 0: #pad cons_seq_1
                        cons_seq_el = cons_seq + '-'*abs(len_diff)
                    elif len_diff > 0:
                        seq = seq + '-'*len_diff
                        cons_seq_el = cons_seq
                else:
                    cons_seq_el = cons_seq

                assert len(seq) == len(cons_seq_el), 'Length of sequences into hamming distance unequal'

                consensus_diff = edit_distance(seq[25:].encode('utf-8'), cons_seq_el[25:].encode('utf-8')) #how many mismatches present
                diffs_from_cons.append(str(consensus_diff))
                # if consensus_diff > 5:
                    # print('seq', seq, 'cons', cons_seq)
                # consensus_diff = Levenshtein.distance(seq, cons_seq)
                #replace possible N's in consensus sequence with best seq
                if consensus_diff < diffs: #need only if many N's in cons_seq
                    diffs = consensus_diff
                    best_seq = str(seq)
                if consensus_diff <= differences:
                    good += 1

    #replace N's in cons_seq by values in best_seq
    if 'N' in cons_seq:
        new_str = ''
        for best_str, cons_str in zip(list(best_seq), list(cons_seq)):
            if best_str == cons_str:
                new_str += best_str
            elif cons_str == 'N':
                new_str += best_str
            else:
                new_str += cons_str
        #0 bad 1 good
        return (good/total, new_str, diffs_from_cons)
    else:
        return (good/total, cons_seq, diffs_from_cons)


# import pickle

#
# def rcs_check(func):
#     def func_wrapper(bundle, clusters, counts, stats, mismtch, gt_threshold, msa, short):
#         count = len(glob.glob('/media/chovanec/My_Passport/test/ratio_file*'))
#         with open('/media/chovanec/My_Passport/test/ratio_file' + str(count) + '.pkl', 'wb') as out_pkl:
#             reads, consensus, final_umis, umi_counts, low_gt, corrected, low_gt_corrected, gt_list = func(bundle, clusters, counts, stats, mismtch, gt_threshold, msa=False, short=False)
#             pickle.dump(gt_list, out_pkl)
#         return reads, consensus, final_umis, umi_counts, low_gt, corrected, low_gt_corrected
#     return func_wrapper
#
# @rcs_check
def reduce_clusters_single(bundle, clusters, counts, stats, mismtch, gt_threshold, msa=False, short=False):
    ''' collapse clusters down to the UMI which accounts for the cluster
    and return the list of final UMIs using consensus sequence'''

    reads = []
    consensus = []
    final_umis = []
    umi_counts = []
    low_gt = 0
    corrected = 0
    low_gt_corrected = 0
    total = 0

    gt_list = []

    # parent_umi_dict = defaultdict()
    cons_diffs = defaultdict()

    for cluster in clusters:
        total += 1
        umi_in_cluster = len(cluster)
        #Consensus for read loss filter
        out_dict = {umi:bundle[umi]['seq'] for umi in cluster}


        assert len(out_dict) != 0, 'No sequence from umi'

        if umi_in_cluster > 1: #contains umi with 1 error
            corrected += 1

        if msa:
            alignment = kalign_msa(out_dict)
            gt_ratio, consensus_seq, diffs_from_cons = read_loss(alignment, differences=mismtch, msa=True) #umi=umi_cons
        else:
            gt_ratio, consensus_seq, diffs_from_cons = read_loss(out_dict, differences=mismtch, short=short) #umi=umi_cons
        #keep record of the distance between sequence and consensus per UMI bases (clustered UMIs seperated by ,)
        if not isinstance(diffs_from_cons, int):
            cons_diffs[','.join(cluster)] = ','.join(x for x in diffs_from_cons)
        else:
            cons_diffs[','.join(cluster)] = diffs_from_cons

        gt_list.append(gt_ratio)
        if gt_ratio >= gt_threshold:
            #Parent umi = highest count umi which account for the cluster
            parent_umi = umi_correction.get_best_higher_counts(cluster, counts)
            reads.append(bundle[parent_umi]['read'])
            consensus.append(consensus_seq.replace('-', '')) #remove padding
            final_umis.append(parent_umi)
            #Number of UMI's in the cluster (how many have been collapsed)
            umi_counts.append(sum([counts[x] for x in cluster]))
        else:
            low_gt += 1

            if umi_in_cluster > 1: #contains umi with 1 error and low ratio
                low_gt_corrected += 1

    assert len(set(reads)) == len(reads), 'Not all reads unique!'

    #list of reads, final umi's used, list of umi counts within clusters
    return reads, consensus, final_umis, umi_counts, low_gt, corrected, low_gt_corrected, cons_diffs#, gt_list



#
# def reduce_clusters_single_parallel(bundle, clusters, counts, nprocs, stats, mismtch, gt_threshold):
#     ''' collapse clusters down to the UMI which accounts for the cluster
#     using the adjacency dictionary and return the list of final UMIs'''
#
#     def worker(bundle, clusters, counts, out_q, stats, gt_threshold):
#
#         inter_results = {'reads':[], 'consensus':[], 'final_umis':[], 'umi_counts':[], 'low_gt':0, 'corrected':0, 'low_gt_corrected':0}
#         group = 0
#
#         for cluster in clusters:
#
#             umis_in_cluster = 0
#
#             #Consensus for read loss filter
#             out = [] #This is what want a consensus of
#             for umi in cluster:
#                 umis_in_cluster += 1
#                 try:
#                     out.append(bundle[umi]['seq'])
#                 except KeyError:
#                     print('UMI not in bundle')
#
#             assert len(out) != 0, 'No sequence from umi'
#
#             if umis_in_cluster > 1:
#                 inter_results['corrected'] += 1
#
#             #Get consensus sequence and good to total ratio
#             # alignment =  msa(out)
#             gt_ratio, consensus_seq = read_loss(alignment, differences=mismtch)
#
#
#             if gt_ratio >= gt_threshold:
#                 #Parent umi = highest count umi which account for the cluster
#                 parent_umi = umi_correction.get_best_higher_counts(cluster, counts)
#                 #get name from highest count but use consensus seq
#                 inter_results['reads'].append(bundle[parent_umi]['read'])
#                 inter_results['consensus'].append(consensus_seq.replace('-', ''))
#
#             else:
#                 inter_results['low_gt'] += 1
#
#                 if umis_in_cluster > 1:
#                     inter_results['low_gt_corrected'] += 1
#
#                 continue
#
#             if stats:
#                 inter_results['final_umis'].append(parent_umi)
#                 #Number of UMI's in the cluster (how many have been collapsed)
#                 inter_results['umi_counts'].append(sum([counts[x] for x in cluster]))
#
#         out_q.put(inter_results)
#
#
#     # Each process will get 'chunksize' nums and a queue to put its out dict into
#     out_q = multiprocessing.Queue()
#     cluster_chunk = int(math.ceil(len(clusters)/float(nprocs)))
#     procs = []
#
#     for i in range(nprocs):
#         p = multiprocessing.Process(target=worker,
#         args=(bundle, clusters[cluster_chunk * i:cluster_chunk * (i+1)], counts, out_q, stats, gt_threshold))
#         procs.append(p)
#         p.start()
#
#     # Collect all results into a single result dict
#     reads = []
#     consensus = []
#     final_umis = []
#     umi_counts = []
#     low_gt = 0
#     corrected = 0
#     low_gt_corrected = 0
#     for i in range(nprocs):
#         out_dict = out_q.get()
#         reads.extend(out_dict['reads'])
#         consensus.extend(out_dict['consensus'])
#         final_umis.extend(out_dict['final_umis'])
#         umi_counts.extend(out_dict['umi_counts'])
#         low_gt += out_dict['low_gt']
#         corrected += out_dict['corrected']
#         low_gt_corrected += out_dict['low_gt_corrected']
#
#     # Wait for all worker processes to finish
#     for p in procs:
#         p.join()
#
#
#     #list of reads, final umi's used, list of umi counts within clusters
#     return reads, consensus, final_umis, umi_counts, low_gt, corrected, low_gt_corrected
#



################################################################################
#Deduplication in parallel
################################################################################

def deduplication_worker_umi(bundle, threshold, stats, mismatch, min_reads, nprocs, gt_threshold, msa, short):
    ''' worker for deduplicate_bundle_parallel '''

    reads, consensus, final_umis, umi_counts, low_gt, corrected, low_gt_corrected =\
    umi_correction.run_dir_adj(bundle, threshold=threshold, stats=stats, mismatches=mismatch,
                nprocs=nprocs, gt_threshold=gt_threshold, msa=msa, short=short)


    num_input = sum([bundle[umi]['count'] for umi in bundle])
    # collect pre-dudupe stats
    stats_pre_df_dict = {'UMI': [], 'counts': []}
    # pre_average_distance = ''
    if stats:
        stats_pre_df_dict['UMI'].extend(bundle) #umi + read
        stats_pre_df_dict['counts'].extend([bundle[UMI]['count'] for UMI in bundle]) #umi counts

        # pre_average_distance = get_average_umi_distance(bundle.keys()) #v_seq + umi

    return [reads, consensus, final_umis, umi_counts, low_gt, corrected, low_gt_corrected, num_input, stats_pre_df_dict, cons_diffs] #, pre_average_distance]


def deduplication_worker(bundle, threshold, stats, mismatch, min_reads, nprocs, gt_threshold, msa, short):
    ''' worker for deduplicate_bundle_parallel without umi correction'''

    #reduce bundle
    umis = bundle.keys()

    #clusters need to be a list of sets
    clusters = []
    for key in bundle.keys():
        clusters.append(set([key]))

    counts = {umi: bundle[umi]['count'] for umi in umis}
    # bundle, clusters, counts, stats, mismtch, gt_threshold, msa=False

    reads, consensus, final_umis, umi_counts, low_gt, corrected, low_gt_corrected, cons_diffs =\
    reduce_clusters_single(bundle, clusters, counts, stats, mismatch, gt_threshold, msa, short)


    num_input = sum([bundle[umi]['count'] for umi in bundle])
    # collect pre-dudupe stats
    stats_pre_df_dict = {'UMI': [], 'counts': []}
    # pre_average_distance = ''
    if stats:
        stats_pre_df_dict['UMI'].extend(bundle) #umi + read
        stats_pre_df_dict['counts'].extend([bundle[UMI]['count'] for UMI in bundle]) #umi counts

        # pre_average_distance = get_average_umi_distance(bundle.keys()) #v_seq + umi

    return [reads, consensus, final_umis, umi_counts, low_gt, corrected, low_gt_corrected, num_input, stats_pre_df_dict, cons_diffs] #, pre_average_distance]



def deduplicate_bundle_parallel(reads_dict, low_umi_out, out, threshold, min_reads, mismatch,
                   stats, threads, pdf_out, gt_threshold, nprocs, msa,
                   skip_umi_correction, no_anchor=False, short=False):
    '''
    :param reads_dict:
    :param min_reads: minimun number of reads required in a umi group [5]
    '''

    if skip_umi_correction:
        #list of result stats
        dir_adj_results = Parallel(n_jobs=threads)(delayed(deduplication_worker)(bundle, threshold, stats,
        mismatch, min_reads, nprocs, gt_threshold, msa, short) for bundle in reads_dict.values())
    else:
        dir_adj_results = Parallel(n_jobs=threads)(delayed(deduplication_worker_umi)(bundle, threshold, stats,
        mismatch, min_reads, nprocs, gt_threshold, msa, short) for bundle in reads_dict.values())


    stats_pre_df_dict_all = {'UMI': [], 'counts': []}
    stats_post_df_dict = {'UMI': [], 'counts': []}
    pre_cluster_stats = []
    post_cluster_stats = []
    stats_cons_diffs = defaultdict()


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
        labels, values = zip(*Counter(dir_adj_results[bundle][3]).items()) #*-operator to unpack the arguments out of a list or tuple
        non_dedup_values = tuple(l*v for l, v in zip(labels, values))

        if min_reads != None:
            # cut_point = count_change(labels, non_dedup_values)

            cut_point = min_reads

            plt.figure()
            plt.bar(labels, non_dedup_values)
            #extract title from read
            # j_nam = re.search('J\d+', dir_adj_results[bundle][0][0].split(' ')[0]).group(0)
            j_nam = dir_adj_results[bundle][0][0].split(' ')[0].split('_')[-3]
            if no_anchor:
                anchor = ''
            else:
                anchor = dir_adj_results[bundle][0][0].split(' ')[0].split('_')[-2]

            plt.title(j_nam + ' ' + anchor + ' Cut point: ' + str(cut_point), ha='center') #need to put name to know which bundle J is being processed
            my_plot = plt.axvline(cut_point, linestyle='dashed', linewidth=2).get_figure()
            pdf_out.savefig(my_plot)
            plt.close('all')

        num_input_all += dir_adj_results[bundle][7] #num_input

        #remove low umi counts 1-5
        indx = 0
        for count in dir_adj_results[bundle][3]: #umi_counts
            if count <= cut_point:
                #write out fasta
                low_umi_out.write(dir_adj_results[bundle][0][indx].split(' ')[0].replace('@', '>') + '_' + str(count) +'\n' + dir_adj_results[bundle][1][indx] + '\n')
                low_umi_count += 1
            else:
                #write out fasta
                out.write(dir_adj_results[bundle][0][indx].split(' ')[0].replace('@', '>') + '_' + str(count) + '\n' + dir_adj_results[bundle][1][indx] + '\n')
                num_output += 1
            indx += 1


        if stats:
            # assert len(reads)  == len(consensus) == len(umi_counts), 'Reads, consensus and counts differ in length'

            ##################

            low_gt_reads += dir_adj_results[bundle][4] #low_gt
            corrected_reads +=  dir_adj_results[bundle][5] #corrected
            low_gt_corrected_reads += dir_adj_results[bundle][6] #low_gt_corrected

            # # collect pre-dudupe stats
            # stats_pre_df_dict['UMI'].extend(bundle) #umi + read
            # stats_pre_df_dict['counts'].extend([bundle[UMI]['count'] for UMI in bundle]) #umi counts
            #
            # pre_average_distance = get_average_umi_distance(bundle.keys()) #v_seq + umi
            stats_pre_df_dict_all['UMI'].extend(dir_adj_results[bundle][8]['UMI']) #stats_pre_df_dict umi + read
            stats_pre_df_dict_all['counts'].extend(dir_adj_results[bundle][8]['counts']) #stats_pre_df_dict umi counts

            # pre_cluster_stats.append(pre_average_distance)

            #aggregate errors per cluster for each bundle
            cons_diffs = dir_adj_results[bundle][9]
            for k,v in cons_diffs.items():
                try:
                    # print(stats_cons_diffs[k], '_', v)
                    stats_cons_diffs[k] += '_' + v
                except KeyError:
                    stats_cons_diffs[k] = v

            # collect post-dudupe stats
            #v_seq + umi
            # post_cluster_umis = [qname.split(' ')[-1] for qname in dir_adj_results[bundle][0]] #reads are just qnames
            stats_post_df_dict['UMI'].extend(dir_adj_results[bundle][2]) #final_umis
            stats_post_df_dict['counts'].extend(dir_adj_results[bundle][3]) #umi_counts



    return [stats_pre_df_dict_all, stats_post_df_dict, pre_cluster_stats, post_cluster_stats,
    num_input_all, num_output, low_gt_reads, corrected_reads, low_gt_corrected_reads, low_umi_count,
    stats_cons_diffs]



######### Stats ################


def aggregate_Stats_df(stats_df):
    ''' return a data from with aggregated counts per UMI'''

    agg_df_dict = {}

    agg_df_dict['total_counts'] = stats_df.pivot_table(
        columns="UMI", values="counts", aggfunc=np.sum)

    # agg_df_dict['median_counts'] = stats_df.pivot_table(
    #     columns="UMI", values="counts", aggfunc=np.median)

    agg_df_dict['times_observed'] = stats_df.pivot_table(
        columns="UMI", values="counts", aggfunc=len)

    return pd.DataFrame(agg_df_dict)


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

        self.file_directory = file_directory
        self.an1 = an1
        self.an2 = an2
        self.out_dir = ''

        try:
            self.jv_fastq_an1 = glob.glob(self.file_directory + '/*all_jv*' + an1)[0]
            self.jv_fastq_an2 = glob.glob(self.file_directory + '/*all_jv*' + an2)[0]
        except IndexError:
            raise Exception('No input files found')

        self.jv_prefix = ''
        self.header_postfix_jv = ''


    def create_dirs_assembled(self, out_dir=None):
        '''#Create directories
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
        with open(self.jv_fastq_an1, 'r') as f:
                first_line = f.readline()
                self.header_postfix_jv = first_line.split(' ')[1].rstrip()




    def v_start_j_umi_dedup_assembled(self, threshold, min_reads, threads, mismatch, gt_threshold, v_len, j_len,
                                      stats=False, ignore_umi=False, ignore_j=False, ignore_v=False, skip_unclear=False,
                                      skip_mh=False, msa=False, skip_umi_cor=False, no_anchor=False, short=False):
        '''Determine start position of v reads
        some elements inspired by umi_tools by tom smith cagt

        '''
        if no_anchor:
            reads_dict, unclear_skip_an1 = make_bundle(self.jv_fastq_an1, v_len=v_len, j_len=j_len, ignore_umi=ignore_umi, ignore_j=ignore_j, ignore_v=ignore_v,
                                                           skip_unclear=skip_unclear, skip_mh=skip_mh, no_anchor=no_anchor, short=short)
            unclear_skip_an2 = 0
        else:

            reads_dict_an1, unclear_skip_an1 = make_bundle(self.jv_fastq_an1, v_len=v_len, j_len=j_len, ignore_umi=ignore_umi,
                                                           ignore_j=ignore_j, ignore_v=ignore_v,
                                                           skip_unclear=skip_unclear, skip_mh=skip_mh, short=short)
            reads_dict_an2, unclear_skip_an2 = make_bundle(self.jv_fastq_an2, v_len=v_len, j_len=j_len, ignore_umi=ignore_umi,
                                                           ignore_j=ignore_j, ignore_v=ignore_v,
                                                           skip_unclear=skip_unclear, skip_mh=skip_mh, short=short)

            # merge the two dict_keys
            reads_dict = {**reads_dict_an1, **reads_dict_an2}

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

        # print(reads_dict.keys())
        # with pysam.AlignmentFile(self.tmp_dir + '/' + self.v_prefix_an1 + '_' + self.an1 + '_dedup.bam', "wb", template=sam_algn_v_an1) as out_file:
        with open(self.out_dir + '/' + self.jv_prefix + '_dedup.fasta', 'w') as jv_out, \
        open(self.out_dir + '/' + self.jv_prefix + '_low_umi.fasta', 'w') as low_umi_out, \
        PdfPages(self.out_dir + '/' + self.jv_prefix + '_histogram.pdf') as pdf:

            #run an1 and an2 side by side (not tested!)
            if len(reads_dict) >= threads:
                nprocs = 1
            else:
                nprocs = int(threads/len(reads_dict)) #how many unused cores are available?

            #an1+2
            stats_pre_df_dict, stats_post_df_dict, pre_cluster_stats, post_cluster_stats, \
            num_input, num_output, low_gt_reads, corrected_reads, \
            low_gt_corrected_reads, low_umi_count, stats_cons_diffs=\
            deduplicate_bundle_parallel(reads_dict, low_umi_out, jv_out, threshold=threshold, min_reads=min_reads,
                           mismatch=mismatch, gt_threshold=gt_threshold,
                           stats=stats, threads=threads, pdf_out=pdf, nprocs=nprocs, msa=msa,
                           skip_umi_correction=skip_umi_cor, no_anchor=no_anchor, short=short)

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
                columns=stats_pre_df['counts'], values='counts', aggfunc=len))
            UMI_counts_df_post = pd.DataFrame(stats_post_df.pivot_table(
                columns=stats_post_df['counts'], values='counts', aggfunc=len))

            UMI_counts_df_pre.columns = ['instances']
            UMI_counts_df_post.columns = ['instances']

            UMI_counts_df = pd.merge(UMI_counts_df_pre, UMI_counts_df_post,
                                     how='outer', left_index=True, right_index=True,
                                     sort=True, suffixes=['_pre', '_post'])

            UMI_counts_df = UMI_counts_df.fillna(0).astype(int)

            UMI_counts_df.to_csv(self.out_dir + '/' + self.jv_prefix + '_per_umi_per_position.tsv', sep='\t')

            ##########################

             # aggregate stats pre/post per UMI
            agg_pre_df = aggregate_Stats_df(stats_pre_df)
            agg_post_df = aggregate_Stats_df(stats_post_df)

            agg_df = pd.merge(agg_pre_df, agg_post_df, how='left',
                              left_index=True, right_index=True,
                              sort=True, suffixes=['_pre', '_post'])

            agg_df = agg_df.fillna(0).astype(int)

            stats_consensus_difference = pd.DataFrame(stats_cons_diffs, index=[0])
            stats_consensus_difference = stats_consensus_difference.T
            stats_consensus_difference.columns = ['Consensus_differences']
            # stats_consensus_difference['UMI'] = stats_consensus_difference.index

            agg_df = pd.merge(agg_df, stats_consensus_difference, how='left',
                              left_index=True, right_index=True,
                              sort=True,)
            agg_df.to_csv(self.out_dir + '/' + self.jv_prefix + '_per_umi.tsv', sep="\t")



            # stats_consensus_difference.to_csv(self.out_dir + '/' + self.jv_prefix + '_consensus_difference.tsv', sep="\t")


def parse_args():
    parser = argparse.ArgumentParser(description='BabrahamLinkON Deduplicate')

    parser.add_argument('--input_dir', dest='in_dir', type=str, required=True, help='Input directory (created for/by preclean)')
    parser.add_argument('--cores', dest='nthreads', default=1, type=int, help='Number of cores to use (if aligning), default: 1')
    # parser.add_argument('--species', dest='species', default='mmu', type=str, help='Which species (mmu, hsa), default: mmu')
    parser.add_argument('--an1', dest='an1', default='GACTCGT', type=str, help='Default: GACTCGT')
    parser.add_argument('--an2', dest='an2', default='CTGCTCCT', type=str, help='Default: CTGCTCCT')
    # parser.add_argument('--verbose', action='store_true', help='Print detailed progress')
    parser.add_argument('--out', dest='out_dir', type=str, help='Output directory, default: creates Deduplicated in main directory')
    parser.add_argument('--no_anchor', action='store_true', help='No anchor sequence present')
    parser.add_argument('--short', action='store_true', help='Short sequences present')

    parser.add_argument('--ignore_umi', action='store_true', help='Deduplicate without using UMI')
    #TODO: make using these default
    parser.add_argument('--ignore_j', action='store_true', help='Deduplicate without using J end identity')
    parser.add_argument('--ignore_v', action='store_true', help='Deduplicate without using 8bp V end')
    parser.add_argument('--v_len', dest='v_len', type=int, default=0, help='Length V end sequence to add to the UMI [0]')
    parser.add_argument('--j_len', dest='j_len', type=int, default=0, help='Length of J end sequence, 50bp into read (-), to add to the UMI (for short UMIs) [0]')

    parser.add_argument('--skip_umi_correction', action='store_true', help='Skip correcting errors that might be present in UMI')


    parser.add_argument('--stats', action='store_true', help='Output stats from UMI deduplication [False]')
    parser.add_argument('--skip_unclear', action='store_true', help='Skip unclear J reads [False]')
    parser.add_argument('--skip_mh', action='store_true', help='Skip multiple hits unclear J reads [False]')
    # parser.add_argument('--assembled', action='store_true', help='Assembled reads are being provided as input (from PEAR) [False]')
    parser.add_argument('--mismatch', dest='mismatch', type=int, default=5, help='Number of mismatches allowed in consensus sequence comparison [5]')
    parser.add_argument('--threshold', dest='threshold', type=int, default=1, help='Number of mismatches allowed in UMI [1]')
    parser.add_argument('--min_reads', dest='minreads', type=int, default=5, help='Minimum number of reads in UMI group, if less than or equal to [5] then discard')
    parser.add_argument('--gt_ratio', dest='gtratio', type=float, default=1, help='Ratio of good to total reads to mark UMI group as early PCR error 0-1 [1]')
    parser.add_argument('--msa', dest='msa', action='store_true', help='Use msa to derive consensus sequence (Slow) [False]')

    parser.add_argument('--umi_seq_logo', dest='seqlogo', action='store_true', help='Make seqlogo from UMIs')

    opts = parser.parse_args()

    return opts



def main():

    #argparse
    opts = parse_args()

    if opts.no_anchor:
        an1 = ''
        an2 = ''
    elif opts.short:
        an1 = ''
        an2 = ''
    else:
        an1 = opts.an1
        an2 = opts.an2

    dedup = deduplicate(file_directory=opts.in_dir, an1=an1, an2=an2)

    dedup.create_dirs_assembled(out_dir=opts.out_dir)

    if opts.seqlogo:
        if opts.no_anchor:
            jv_fastq = glob.glob(opts.in_dir + '/*all_jv')[0]
            UMI_seqlogo.umi_seq_logo(jv_fastq, dedup.out_dir + '/' + dedup.jv_prefix + '.eps')
        else:
            jv_fastq_an1 = glob.glob(opts.in_dir + '/*all_jv*' + an1)[0]
            jv_fastq_an2 = glob.glob(opts.in_dir + '/*all_jv*' + an2)[0]
            UMI_seqlogo.umi_seq_logo(jv_fastq_an1, dedup.out_dir + '/' + dedup.jv_prefix + '_' + an1 + '.eps')
            UMI_seqlogo.umi_seq_logo(jv_fastq_an2, dedup.out_dir + '/' + dedup.jv_prefix + '_' + an2 + '.eps')
    else:


        logging.basicConfig(level=logging.DEBUG, filename=dedup.out_dir +'/' + dedup.jv_prefix + '.log', filemode='a+',
                            format='%(asctime)-15s %(levelname)-8s %(message)s')

        logging.info(opts)
        print('Starting deduplication')
        dedup.v_start_j_umi_dedup_assembled(threshold=opts.threshold, min_reads=opts.minreads, threads=opts.nthreads,
                                            mismatch=opts.mismatch, gt_threshold=opts.gtratio, v_len=opts.v_len, j_len=opts.j_len,
                                            stats=opts.stats,
                                            ignore_umi=opts.ignore_umi, ignore_j=opts.ignore_j, ignore_v=opts.ignore_v,
                                            skip_unclear=opts.skip_unclear, skip_mh=opts.skip_mh, msa=opts.msa,
                                            skip_umi_cor=opts.skip_umi_correction,
                                            no_anchor=opts.no_anchor, short=opts.short)



#TODO: do a collision measure like in presto, could perhaps remove collision UMIs distance-to-nearest plot
#TODO: could speeed up UMI correction if first doing hierarchical clustering?

if __name__ == "__main__":
    main()
