#!/usr/bin/env python3

import os
import pysam
import glob
from pathlib import Path
import subprocess
import shutil
import re
import collections
import fileinput
import shlex
import argparse
import itertools
import pandas as pd
import numpy as np
from babrahamlinkon import general, presets
from skbio import DNA, TabularMSA
import multiprocessing
import math
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import itertools
from joblib import Parallel, delayed
import logging
# import pickle
import warnings
import operator
import Levenshtein
from copy import deepcopy

import pyximport

### Parts modified from UMI tools: ###

from babrahamlinkon._dedup_umi import edit_distance


'''A functor that clusters a bundle of reads,
indentifies the parent UMIs and returns the selected reads, umis and counts
The initiation of the functor defines the methods:
  ** get_adj_list ** - returns the edges connecting the UMIs
  ** connected_components ** - returns clusters of connected components
                               using the edges in the adjacency list
  ** get_best ** - returns the parent UMI(s) in the connected_components
  ** reduce_clusters ** - loops through the connected components in a
                          cluster and returns the unique reads. Optionally
                          returns lists of umis and counts per umi also
Note: The get_adj_list and connected_components methods are not required by
all custering methods. Where there are not required, the methods return
None or the input parameters.

directional-adjacency
from:
https://github.com/CGATOxford/UMI-tools/blob/master/umi_tools/dedup.py

The MIT License (MIT)

Copyright (c) 2015 CGAT

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.


'''

# Use Cpython edit_distance
# def edit_dist(first, second):
#     ''' returns the edit distance/hamming distances between
#     its two arguements '''
#
#     dist = sum([not a == b for a, b in zip(first, second)])
#     return dist


def breadth_first_search(node, adj_list):
    searched = set()
    found = set()
    queue = set()
    queue.update((node,))
    found.update((node,))

    while len(queue)>0:
        node=(list(queue))[0]

        found.update(adj_list[node])
        queue.update(adj_list[node])

        searched.update((node,))
        queue.difference_update(searched)

    return found


######### read loss #########

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
        base_counts = collections.Counter(pos)
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
        base_counts = collections.Counter(pos)
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



def msa(seq_counter_dict): #, umi=None
    '''Read loss analysis
    :param dict seq_counter: dict of Counter object with sequences
    :param differences: number of differences from consensus allowed
    '''
    #Only aligning single copy of a duplicate sequence and then subsequently
    #multiplying the output alignment by Counter

    #Convert Counter into fasta input for kalign
    #Need to sort beacuse of Instability in progressive multiple sequence alignment algorithms
    # seq_fasta = ''
    # count = 0
    # reads = 0
    # for umi, cntr in sorted(seq_counter_dict.items(), key=operator.itemgetter(0)):
    #     for seq, freq in sorted(cntr.items(), key=lambda x:x[0]):
    #         seq_fasta = seq_fasta + '>' + str(count) + '_' + umi + '_' + str(freq) + '\n' + seq + '\n'
    #         count += 1  #how many different reads
    #         reads += freq #how many total reads
    #
    #
    # #Can't get consensus from single sequence, return single sequence
    # if count == 1:
    #     assert len(list(seq_counter_dict.values())) == 1, 'Not a single sequence'
    #     seq_out= re.sub('\n', '', seq_fasta)
    #     seq_out = list(filter(None, re.split('(>\d+_[A-Z]+_\d+)', seq_out)))
    #     return seq_out
    #
    #
    # #Multiple sequence alignment
    # #http://msa.sbc.su.se/cgi-bin/msa.cgi
    # kalign_cmd = ['kalign', '-f', 'fasta']
    #
    #
    # p = subprocess.Popen(kalign_cmd, stdout=subprocess.PIPE, stdin=subprocess.PIPE, stderr=subprocess.STDOUT)
    # kalign_stdout = p.communicate(input=seq_fasta.encode('utf-8'))[0]
    #
    # #parse kalign output and derive consensus seq
    # head, sep, tail = kalign_stdout.partition(b'>') #remove head kalign intro
    #
    # alignment = sep.decode() + tail.decode()
    # alignment = re.sub('\n', '', alignment)
    # alignment = list(filter(None, re.split('(>\d+_[A-Z]+_\d+)', alignment)))
    #
    # # assert len(alignment) == count*2, 'Alignment output fewer reads'
    # return alignment

    # msa_algn = general.mafft()
    #
    # single_seq = msa_algn.write_fasta(seq_counter_dict)
    #
    # if single_seq is None:
    #     alignment = msa_algn.align(threads=4)
    #
    #     msa_algn.del_tmp()
    #     # print(alignment)
    #     return alignment
    # else:
    #     return single_seq

# def consensus_difference(alignment):
def consensus_difference(seq_counter_dict):
    '''
    :param alignment: output from msa
    :return: number of differences between two umi group consensus sequences
    '''

    # seq_dict = collections.defaultdict(list)
    #
    # for item in general.fasta_parse(alignment):
    #     qname = item[0]
    #     seq = item[1]
    #
    #     freq = int(qname.split('_')[-1]) #freq saved in name
    #     #split fasta into umi groups (head and child)
    #     umi = qname.split('_')[-2]
    #     for i in range(freq):
    #         seq_dict[umi].append(seq)



    # assert len(seq_dict.values()) == 2, 'More than two UMI groups, only two allowed'

    # lst_1, lst_2 = list(seq_dict.values())
    #
    # lst_lists_1 = [list(item) for item in lst_1]
    # lst_lists_2 = [list(item) for item in lst_2]

    cntr_1, cntr_2 = list(seq_counter_dict.values())

    lst_lists_1 = [list(item) for item in cntr_1.elements()]
    lst_lists_2 = [list(item) for item in cntr_2.elements()]

    cons_seq_1 = consensus_unequal(lst_lists_1)
    cons_seq_2 = consensus_unequal(lst_lists_2)

    # num_diffs = edit_distance(cons_seq_1.encode('utf-8'), cons_seq_2.encode('utf-8'))

    #Need to pad seq if length unequal!
    len_diff = len(cons_seq_1) - len(cons_seq_2)
    if len_diff < 0: #pad cons_seq_1
        cons_seq_1 = cons_seq_1 + '-'*abs(len_diff)
    elif len_diff > 0:
        cons_seq_2 = cons_seq_2 + '-'*len_diff

    assert len(cons_seq_1) == len(cons_seq_2), 'Sequences for hamming distance not same length!'

    num_diffs = edit_distance(cons_seq_1.encode('utf-8'), cons_seq_2.encode('utf-8'))

    return num_diffs



# def read_loss(alignment, differences=5):
def read_loss(seq_counter_dict, differences=5): #, umi=None
    '''Read loss analysis
    :param alignment: fasta from msa function
    :param differences: number of differences from consensus allowed
    '''
    # seq_dict = collections.defaultdict(list)
    #
    # for item in general.fasta_parse(alignment):
    #     qname = item[0]
    #     seq = item[1]
    #
    #     freq = int(qname.split('_')[-1]) #freq saved in name
    #     for i in range(freq):
    #         seq_dict[qname+str(i)] = seq #each seq has unique key
    #
    # #all values in list of lists
    # seq_list_lsts = [list(item) for item in seq_dict.values()]

    list_of_lists = []
    for umi, seqs in seq_counter_dict.items():
        for seq in seqs.elements():
            list_of_lists.append(list(seq))

    # cons_seq = consensus(seq_list_lsts)
    cons_seq = consensus_unequal(list_of_lists)

    #Skip poor consensus (i.e. if building consensus from only two different seqs)
    if cons_seq.count('N') > differences: #5
        return(0, cons_seq)

    good = 0
    total = 0
    diffs = len(cons_seq)
    best_seq = ''
    for umi, seqs in seq_counter_dict.items():
        for seq in seqs.keys():
            total += 1

            #Need to pad seq if length unequal!
            len_diff = len(cons_seq) - len(seq)
            if len_diff < 0: #pad cons_seq_1
                cons_seq = cons_seq + '-'*abs(len_diff)
            elif len_diff > 0:
                seq = seq + '-'*len_diff

            assert len(seq) == len(cons_seq), 'Length of sequences into hamming distance unequal'

            consensus_diff = edit_distance(seq.encode('utf-8'), cons_seq.encode('utf-8')) #how many mismatches present
            # if consensus_diff > 5:
                # print('seq', seq, 'cons', cons_seq)
            # consensus_diff = Levenshtein.distance(seq, cons_seq)
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

        return (good/total, new_str)
    else:
        return (good/total, cons_seq)


######## "get_best" methods ##########

def get_best_higher_counts(cluster, counts):
    ''' return the UMI with the highest counts'''
    count = 0

    if len(cluster) == 1:
        return list(cluster)[0]
    else:

        sorted_nodes = sorted(cluster, key=lambda x: counts[x],
                              reverse=True)
        return sorted_nodes[0]



######## "get_adj_list" methods ##########


def ham_adj_list_directional_adjacency(umis, counts, threshold=1):
    ''' identify all umis within the hamming distance threshold (1 is best)
    and where the counts of the first umi is > (2 * second umi counts)-1
    will have duplicates'''

    # return {umi: [umi2 for umi2 in umis if
    #               edit_distance(umi.encode('utf-8'),
    #                             umi2.encode('utf-8')) == threshold and
    #               counts[umi] >= (counts[umi2]*2)-1] for umi in umis}

    #should be 50% faster
    adj_list = collections.defaultdict(list)
    for i,umi in enumerate(umis):
        a1 = adj_list[umi]
        c1 = counts[umi]
        for j in range(i+1,len(umis)):
            umi2 = umis[j] #dict_keys object doesn't support indexing
            if edit_distance(umi.encode('utf-8'), umi2.encode('utf-8')) == threshold:
                c2 = counts[umi2]
                if c1 >= (c2*2)-1:
                    adj_list[umi].append(umi2)
                if c2 >= (c1*2)-1:
                    adj_list[umi2].append(umi)
    return adj_list

######## "get_connected_components" methods ##########

def get_connected_components_adjacency(umis, graph, counts):
    ''' find the connected UMIs within an adjacency dictionary'''

    found = list()
    components = list()

    for node in sorted(graph, key=lambda x: counts[x], reverse=True):
        if node not in found:
            component = breadth_first_search(node, graph)
            found.extend(component)
            components.append(component)

    if 'GGAGCAAGCTTTCA' in components:
        print('components', components)
    return components


# def merge(sets):
#     '''
#     Merge shared sets in a list
#     '''
#     merged = 1
#     while merged:
#         merged = 0
#         results = []
#         while sets:
#             common, rest = sets[0], sets[1:]
#             sets = []
#             for x in rest:
#                 if x.isdisjoint(common): #sets are disjoint only if their intersection is the empty set
#                     sets.append(x) #don't share values with other sets
#                 else:
#                     merged = 1
#                     common.update(x)
#             results.append(common)
#         sets = results
#     return sets
#
#
# #TODO: merge get_shared and merge functions?
# def get_shared(network_dict):
#     shared_list = [] #which head nodes share values
#     found = set()
#     for node in network_dict: #head_node
#         if node not in found:
#             shared_values = set()
#             for key, values in network_dict.items(): #loop through all networks to check if they share values with node
#                 if len(network_dict[node].intersection(values)) > 0: #setA network_dict[node] share values with setB (values)
#                     shared_values.add(key) #will also add initial head node itself
#                     found.add(key)
#             if len(shared_values) != 0:
#                 shared_list.append(shared_values)
#     return shared_list


def merge_dict(in_dict):
    #Messes up in_dict!
    k_v = [(set([k]), v) for k, v in in_dict.items()]
    merged = 1
    while merged:
        merged = 0
        results = []
        while k_v:
            common, rest = k_v[0], k_v[1:]
            k_v = []
            for x in rest:
                if x[1].isdisjoint(common[1]): #sets are disjoint only if their intersection is the empty set
                    k_v.append((x[0], x[1])) #don't share values with other sets
                else:
                    merged = 1
                    common[1].update(x[1])
                    common[0].update(x[0])
            results.append(common)
        k_v = results
    #return only keys (shared keys)
    return [tp[0] for tp in k_v]



def resolve_clusters(bundle, clusters, counts, differences, gt_threshold):
    '''
    Which shared nodes belong to which head node
    '''
    single_components = []
    network_dict = collections.defaultdict(set)
    # network_dict_2 = collections.defaultdict()
    cont_comp = 0

    for comp in clusters:
        pr_comp = False
        if 'CCTAGGGACGCTGG' in set(comp):
            print('comp_CCTAGGGACGCTGG', comp)
            pr_comp = True

        if 'CCTAGGGACGATGC' in set(comp):
            print('comp_CCTAGGGACGATGC', comp)
            pr_comp = True
        if len(comp) == 1: #not a cluster
            single_components.append(comp)
        else:
            cont_comp += len(comp)
            ordered_network = sorted(comp, key=lambda x: counts[x], reverse=True)
            if pr_comp:
                print(ordered_network)
            network_dict[ordered_network[0]] = set(ordered_network[1:])
            # network_dict[ordered_network[0]].update(ordered_network[1:])


    #which clusters share components
    print('single_component', len(single_components), 'network_dict', len(network_dict))#, 'network_dict_2', len(network_dict_2))

    # my_list = ['TGAGGTCCCAGCGG', 'TGAGGTCCACCCGA', 'TGAGGTCCACGCGA', 'TGAGGTCCACGCGG', 'TGAGGTCCCCGCGG']
    # for item in my_list:
    #     try:
    #         print(item, bundle[item]['seq'])
    #     except KeyError:
    #         pass

    shared_list = merge_dict(deepcopy(network_dict))
    # network_dict_2 = {item[0]:item[1] for item in sorted(network_dict.items(), key=operator.itemgetter(1))}
    # shared_list_2 = merge_dict(network_dict_2)
    #
    # single_components = []
    # network_dict = collections.defaultdict(set)
    # # network_dict_2 = collections.defaultdict()
    # cont_comp = 0
    #
    # for comp in clusters:
    #     pr_comp = False
    #     if 'CCTAGGGACGCTGG' in set(comp):
    #         print('comp_CCTAGGGACGCTGG', comp)
    #         pr_comp = True
    #
    #     if 'CCTAGGGACGATGC' in set(comp):
    #         print('comp_CCTAGGGACGATGC', comp)
    #         pr_comp = True
    #     if len(comp) == 1: #not a cluster
    #         single_components.append(comp)
    #     else:
    #         cont_comp += len(comp)
    #         ordered_network = sorted(comp, key=lambda x: counts[x], reverse=True)
    #         if pr_comp:
    #             print(ordered_network)
    #         network_dict[ordered_network[0]] = set(ordered_network[1:])
    #         # network_dict[ordered_network[0]].update(ordered_network[1:])




    for item in shared_list:
        if item == {'CCTAGGGACGATGC', 'CCTAGGGACGCTGG'}:
            print('present', {'CCTAGGGACGATGC', 'CCTAGGGACGCTGG'})
        if item == {'GGAGCAAGCTTTCA', 'GGAGCAAGGTTTCG'}:
            print('present', {'GGAGCAAGCTTTCA', 'GGAGCAAGGTTTCG'})


    my_count = collections.Counter()
    for item in shared_list:
        my_count.update([len(item)])

    print('counter', my_count)


    connect_comp = []
    added_key = set()
    duplicates_removed = 0
    head_node_removed = 0
    single_share = 0
    umbigous = 0



    shared_num = 0
    combinations = 0
    # print(shared_list)
    shared_list.sort(key=len, reverse=True)
    # shared_list_lst = []
    # for item in shared_list:
    #     shared_list_lst.append(sorted(item))
    # print(shared_list_lst)
    remove_dict = collections.defaultdict(set)
    #Process shared one at a time (create subset of network_dict)
    for shared in shared_list:
        shr = False
        if shared == {'CCTAGGGACGATGC', 'CCTAGGGACGCTGG'}:
            print('present_shared', {'CCTAGGGACGATGC', 'CCTAGGGACGCTGG'})
            shr = True
        if shared == {'GGAGCAAGCTTTCA', 'GGAGCAAGGTTTCG'}:
            print('present_shared', {'GGAGCAAGCTTTCA', 'GGAGCAAGGTTTCG'})
            shr = True

        shared_num += 1

        #get all the netwroks that share components
        # shared_network = {key:network_dict[key] for key in shared}
        # shared_network = {key:value for key,value in network_dict.items() if key in shared}
        appended_df = []
        #skip next part if head node not sharing item with other head nodes
        if len(shared) > 1:
            # print('shared', sorted(shared))

            remove_dict = collections.defaultdict(set)
            rm = False
            for head_1, head_2 in itertools.combinations(list(shared),2):
                combinations += 1
                # shared_values = shared_network[head_1].intersection(shared_network[head_2])
                h_1 = network_dict[head_1]
                h_2 = network_dict[head_2]
                shared_values = h_1.intersection(h_2)
                if shr:
                    print('shared_values',shared_values)
                    print('hd1', head_1, 'hd2', head_2)
                    # print('sh_h1', shared_network[head_1], 'sh_h2', shared_network[head_2])
                    print('n_h1', network_dict[head_1], 'n_h2', network_dict[head_2])
                    if network_dict[head_1] == network_dict[head_2]:
                        print('same')
                    # print('shared_values_2', shared_values_2)
                # print(sorted(shared_values), len(shared), shared_values==shared_values_2)
                columns = [head_1, head_2]
                df = pd.DataFrame(index = shared_values, columns=columns)

                if len(shared_values) > 0: #head nodes need to share values
                    # head_1_gt, head_1_cons = head_cons_dict[head_1]
                    # head_2_gt, head_2_cons = head_cons_dict[head_2]

                    to_1 = 0
                    to_2 = 0
                    for item in sorted(shared_values):

                        # if item in {'TGAGGTCCACGCGA', 'TGAGGTCCACGCGG', 'TGAGGTCCCCGCGG'}:
                        #     print(item, bundle[item]['seq'])
                        #     rm = True
                        #
                        # if head_1 in {'TGAGGTCCCAGCGG', 'TGAGGTCCACCCGA'}:
                        #     print(head_1, bundle[head_1]['seq'])
                        #     resolve_dict_1 = {head_1:bundle[head_1]['seq'], item:bundle[item]['seq']}
                        #     print('dict_1', resolve_dict_1)
                        #     rm = True
                        # if head_2 in {'TGAGGTCCCAGCGG', 'TGAGGTCCACCCGA'}:
                        #     print(head_2, bundle[head_2]['seq'])
                        #     resolve_dict_2 = {head_2:bundle[head_2]['seq'], item:bundle[item]['seq']}
                        #     print('dict_2', resolve_dict_2)
                        #     rm = True

                        resolve_dict_1 = {head_1:bundle[head_1]['seq'], item:bundle[item]['seq']}
                        resolve_dict_2 = {head_2:bundle[head_2]['seq'], item:bundle[item]['seq']}
                        #Align head 1 and value
                        # algn_1 = msa(resolve_dict_1)
                        # diff_1 = consensus_difference(algn_1)
                        diff_1 = consensus_difference(resolve_dict_1)
                        #Align head 2 and value
                        # algn_2 = msa(resolve_dict_2)
                        # diff_2 = consensus_difference(algn_2)
                        diff_2 = consensus_difference(resolve_dict_2)
                        if math.isnan(diff_1):
                            print(resolve_dict_1, 'not a number')
                        if math.isnan(diff_2):
                            print(resolve_dict_2, 'not a number')


                        df.set_value(item, head_1, diff_1)
                        df.set_value(item, head_2, diff_2)

                        #which ever is lower asign to that
                        # print('diff_1', diff_1, 'diff_2', diff_2)
                        if diff_1 < diff_2:
                            remove_dict[head_2].update([item])
                        elif diff_1 > diff_2:
                            remove_dict[head_1].update([item])
                        elif diff_1 == diff_2: #remove from both
                            umbigous += 1
                            # print('rs1', diff_1, resolve_dict_1)
                            # print('rs2', diff_2, resolve_dict_2)
                            # print('shared', shared)
                            rm = True
                            # print(item, head_1, head_2, 'shared_values', shared_network)
                            # print('resolve_dict_1', resolve_dict_1)
                            # print('resolve_dict_2', resolve_dict_2)
                            remove_dict[head_2].update([item])
                            remove_dict[head_1].update([item])

                        else:
                            print('Something else')
                        #
                        # if rm:
                        #     print('remove_dict', remove_dict, 'diff_1', diff_1, 'diff_2', diff_2, 'rd_1', resolve_dict_1, 'rd_2', resolve_dict_2)

                appended_df.append(df)
            # print(pd.concat(appended_df, axis=1))
            # if rm:
            #     print('remove_dict', remove_dict)
            # print('sn before', shared_network)
            for key, value in remove_dict.items():
                # shared_network[key].difference_update(value)
                network_dict[key].difference_update(value)
            # print('sn after', shared_network)
        else:
            single_share += 1


        # for key, value in shared_network.items():
        #     if key not in added_key:
        #
        #         node = set([key]).union(value)
        #         connect_comp.append(node)
        #     else:
        #         duplicates_removed +=1

    for key, value in network_dict.items():
        if key not in added_key:

            node = set([key]).union(value)
            connect_comp.append(node)
        else:
            duplicates_removed +=1

    # 'shared_list_2', sum(len(item) for item in shared_list_2),
    print('shared_list', sum(len(item) for item in shared_list),  'connect_comp', len(connect_comp), 'shared_num', shared_num, 'combinations', combinations)

    out_components = single_components + connect_comp

    all_connected = 0
    for item in connect_comp:
        all_connected += len(item) #difference in how clusters are resolved

    print('connect_comp:', len(connect_comp), 'cont_comp', cont_comp, 'all_connected', all_connected, 'out_component', len(out_components),\
    'single_share', single_share, 'dup_removed', duplicates_removed, 'umbigous', umbigous)
    # print('connect_comp set:', len(set(frozenset(cluster) for cluster in connect_comp)))
    # if len(connect_comp) != len(set(frozenset(cluster) for cluster in connect_comp)):
    #     # print(connect_comp)
    #     with open('/media/chovanec/My_Passport/Sync/BabrahamLinkON/notebooks/network_dict.json', 'wb') as f:
    #         pickle.dump(network_dict, f)
    #     with open('/media/chovanec/My_Passport/Sync/BabrahamLinkON/notebooks/bundle.json', 'wb') as f:
    #         pickle.dump(bundle, f)
    #
    # # assert  len(connect_comp) == len(set(frozenset(cluster) for cluster in connect_comp)), 'End'
    # print('single_components:', len(single_components))
    # print('single_components set:', len(set(frozenset(cluster) for cluster in single_components)))
    # assert duplicates_removed == 0, 'Duplicates present'

    print('Head nodes removed:', head_node_removed)
    # print('Final out:',len(out_components))
    return out_components


######## "reduce_clusters" methods ##########


def reduce_clusters_single(bundle, clusters, counts, stats, mismtch, gt_threshold):
    ''' collapse clusters down to the UMI which accounts for the cluster
    using the adjacency dictionary and return the list of final UMIs
    using consensus sequence'''

    reads = []
    consensus = []
    final_umis = []
    umi_counts = []
    low_gt = 0
    corrected = 0
    low_gt_corrected = 0
    total = 0


    parent_umi_dict = collections.defaultdict()
    # cons_not_match = 0
    # print('all clusters:',len(clusters))
    # print('unique clusters:', len(set(frozenset(cluster) for cluster in clusters)))

    for cluster in clusters:
        total += 1
        umi_in_cluster = len(cluster)
        #Consensus for read loss filter
        # out = [] #This is what want a consensus of
        # umi_cons = []

        out_dict = {umi:bundle[umi]['seq'] for umi in cluster}

        # for umi in cluster: #if not corrected should have single UMI in cluster
        #     umi_in_cluster += 1
        #     try:
        #         out.append(bundle[umi]['seq'])
        #         # umi_cons.append(umi)
        #     except KeyError:
        #         print('UMI not in bundle')


        assert len(out_dict) != 0, 'No sequence from umi'

        if umi_in_cluster > 1: #contains umi with 1 error
            corrected += 1
        # print(umi_cons)
        #TODO: Make sure there aren't any long streches of - which would lead to low gt_ratio
        # alignment = msa(out_dict)
        # gt_ratio, consensus_seq = read_loss(alignment, differences=mismtch) #umi=umi_cons
        gt_ratio, consensus_seq = read_loss(out_dict, differences=mismtch) #umi=umi_cons

        if gt_ratio >= gt_threshold:
            #Parent umi = highest count umi which account for the cluster
            parent_umi = get_best_higher_counts(cluster, counts)
            # if consensus_seq != bundle[parent_umi]['seq'].most_common()[0][0]:
            #     cons_not_match += 1
            #     if cons_not_match < 100:
            #         print(consensus_seq, bundle[parent_umi]['seq'].most_common()[0][0])
            reads.append(bundle[parent_umi]['read'])
            #make sure all bundles are length 1
            # if len(bundle[parent_umi]['read']) != 1:
            #     print('parent umi:', parent_umi)
            #     print('bundle', bundle[parent_umi]['read'])
            # assert len(bundle[parent_umi]['read']) == 1, 'More reads than 1'
            # parent_umi_dict[bundle[parent_umi]['read']] = cluster
            consensus.append(consensus_seq.replace('-', ''))
            final_umis.append(parent_umi)
            #Number of UMI's in the cluster (how many have been collapsed)
            umi_counts.append(sum([counts[x] for x in cluster]))
        else:
            #REVIEW: if the gt ratio low seperate seq into unique seqs? UMI not unique enought and count too low?
            low_gt += 1
            # print('Cluster:', cluster, 'seqs:', gt_ratio)
            if umi_in_cluster > 1: #contains umi with 1 error and low ratio
                low_gt_corrected += 1
    #
    # if total == 14176:
    #     with open('/media/chovanec/My_Passport/Dan_paper_datasets/mu/test_cluster', 'w') as out_file:
    #         out_file.write(out_lst)

    print('Total:', total, 'Output:', len(reads), 'low_gt', low_gt, 'umi_count', len(umi_counts), 'corrected', corrected)
    # print(gt_ratio_lst)
    assert len(set(reads)) == len(reads), 'Not all reads unique!'
    # if len(set(reads)) != len(reads):
    #     print('bad bundle')
        # print(parent_umi_dict)
    # print('Not match:', cons_not_match)
    #list of reads, final umi's used, list of umi counts within clusters
    return reads, consensus, final_umis, umi_counts, low_gt, corrected, low_gt_corrected




def reduce_clusters_single_parallel(bundle, clusters, counts, nprocs, stats, mismtch, gt_threshold):
    ''' collapse clusters down to the UMI which accounts for the cluster
    using the adjacency dictionary and return the list of final UMIs'''

    def worker(bundle, clusters, counts, out_q, stats, gt_threshold):

        inter_results = {'reads':[], 'consensus':[], 'final_umis':[], 'umi_counts':[], 'low_gt':0, 'corrected':0, 'low_gt_corrected':0}
        group = 0

        for cluster in clusters:

            umis_in_cluster = 0

            #Consensus for read loss filter
            out = [] #This is what want a consensus of
            for umi in cluster:
                umis_in_cluster += 1
                try:
                    out.append(bundle[umi]['seq'])
                except KeyError:
                    print('UMI not in bundle')

            assert len(out) != 0, 'No sequence from umi'

            if umis_in_cluster > 1:
                inter_results['corrected'] += 1

            #Get consensus sequence and good to total ratio
            # alignment =  msa(out)
            gt_ratio, consensus_seq = read_loss(alignment, differences=mismtch)


            if gt_ratio >= gt_threshold:
                #Parent umi = highest count umi which account for the cluster
                parent_umi = get_best_higher_counts(cluster, counts)
                #get name from highest count but use consensus seq
                inter_results['reads'].append(bundle[parent_umi]['read'])
                inter_results['consensus'].append(consensus_seq.replace('-', ''))

            else:
                inter_results['low_gt'] += 1

                if umis_in_cluster > 1:
                    inter_results['low_gt_corrected'] += 1

                continue

            if stats:
                inter_results['final_umis'].append(parent_umi)
                #Number of UMI's in the cluster (how many have been collapsed)
                inter_results['umi_counts'].append(sum([counts[x] for x in cluster]))

        out_q.put(inter_results)


    # Each process will get 'chunksize' nums and a queue to put its out dict into
    out_q = multiprocessing.Queue()
    cluster_chunk = int(math.ceil(len(clusters)/float(nprocs)))
    procs = []

    for i in range(nprocs):
        p = multiprocessing.Process(target=worker,
        args=(bundle, clusters[cluster_chunk * i:cluster_chunk * (i+1)], counts, out_q, stats, gt_threshold))
        procs.append(p)
        p.start()

    # Collect all results into a single result dict
    reads = []
    consensus = []
    final_umis = []
    umi_counts = []
    low_gt = 0
    corrected = 0
    low_gt_corrected = 0
    for i in range(nprocs):
        out_dict = out_q.get()
        reads.extend(out_dict['reads'])
        consensus.extend(out_dict['consensus'])
        final_umis.extend(out_dict['final_umis'])
        umi_counts.extend(out_dict['umi_counts'])
        low_gt += out_dict['low_gt']
        corrected += out_dict['corrected']
        low_gt_corrected += out_dict['low_gt_corrected']

    # Wait for all worker processes to finish
    for p in procs:
        p.join()


    #list of reads, final umi's used, list of umi counts within clusters
    return reads, consensus, final_umis, umi_counts, low_gt, corrected, low_gt_corrected




######### Call ################


def run_dir_adj(bundle, threshold, stats, further_stats, mismatches, nprocs, gt_threshold):
    #threshold=1, stats=True, further_stats=True, mismatches=5
    umis = bundle.keys()
    # print(umis)
    # print(sorted(umis)[1:5])
    len_umis = [len(x) for x in umis]
    assert max(len_umis) == min(len_umis), (
        "not all umis are the same length(!):  %d - %d" % (
            min(len_umis), max(len_umis)))

    counts = {umi: bundle[umi]['count'] for umi in umis} #If two UMI's mergered, count will be taken only from the larger UMI group
    print('Getting directional adjacency list')
    adj_list = ham_adj_list_directional_adjacency(list(umis), counts, threshold)
    # print(len(adj_list), sorted(adj_list)[1:5])
    print('Getting connected components')
    clusters = get_connected_components_adjacency(umis, adj_list, counts)


    # print(len(clusters), sorted(clusters)[1:5])
    print('Resolving clusters')
    #gt_threshold implemented here too so will reduce the overall number of clusters
    rclusters = resolve_clusters(bundle, clusters, counts, mismatches, gt_threshold)
    print('rclutsers', len(rclusters))
    print('Reducing clusters')
    # if nprocs == 1: #if no additional cores available
    reads, consensus, final_umis, umi_counts, low_gt, corrected, low_gt_corrected =\
    reduce_clusters_single(bundle, rclusters, counts, stats, mismatches, gt_threshold)

    # else


    print('Unique:', len(set(reads)), 'All:', len(reads))
    #TODO: http://stackoverflow.com/questions/6974695/python-process-pool-non-daemonic/8963618#8963618
    # reads, consensus, final_umis, umi_counts, low_gt, corrected, low_gt_corrected =\
    # reduce_clusters_single_parallel(bundle, rclusters, counts, nprocs, stats, mismatches, gt_threshold)


    # (bundle, clusters, counts, stats, mismtch)

    if further_stats:
        topologies = collections.Counter()
        nodes = collections.Counter()

        if len(clusters) == len(umis):
            topologies["single node"] = len(umis)
            nodes[1] = len(umis)
        else:
            for cluster in clusters:
                if len(cluster) == 1:
                    topologies["single node"] += 1
                    nodes[1] += 1
                else:
                    most_con = max([len(adj_list[umi]) for umi in cluster])

                    if most_con == len(cluster):
                        topologies["single hub"] += 1
                        nodes[len(cluster)] += 1
                    else:
                        topologies["complex"] += 1
                        nodes[len(cluster)] += 1

    else:
        topologies = None
        nodes = None

    return reads, consensus, final_umis, umi_counts, low_gt, corrected, low_gt_corrected, topologies, nodes


# def count_change(labels, values):
#     '''Histogram low count cutoff auto <=
#     '''
#
#     if len(values) < 5:
#         return 1
#     current_value = 0
#     for position in range(len(values)):
#
#         if current_value == 5:
#             return 5
#         else:
#             if values[position] < values[position+1]:
#                 #if true check one more position
#                 if values[position+1] < values[position+2]:
#                     return labels[position]
#             else:
#                 continue
#         current_value += 1


def dir_adj_worker(bundle, threshold, stats, further_stats, mismatch, min_reads, nprocs, gt_threshold):
    ''' worker for dir_adj_bundle_parallel '''
    reads, consensus, final_umis, umi_counts, low_gt, corrected, low_gt_corrected, topologies, nodes =\
    run_dir_adj(bundle, threshold=threshold, stats=stats, further_stats=further_stats,
                mismatches=mismatch, nprocs=nprocs, gt_threshold=gt_threshold)
    # print(reads, umis, umi_counts, topologies, nodes)

    num_input = sum([bundle[umi]['count'] for umi in bundle])
    # collect pre-dudupe stats
    stats_pre_df_dict = {'UMI': [], 'counts': []}
    # pre_average_distance = ''
    if stats:
        stats_pre_df_dict['UMI'].extend(bundle) #umi + read
        stats_pre_df_dict['counts'].extend([bundle[UMI]['count'] for UMI in bundle]) #umi counts

        # pre_average_distance = get_average_umi_distance(bundle.keys()) #v_seq + umi

    return [reads, consensus, final_umis, umi_counts, low_gt, corrected, low_gt_corrected, topologies, nodes, num_input, stats_pre_df_dict] #, pre_average_distance]


def dir_adj_bundle_parallel(reads_dict, low_umi_out, out, threshold, min_reads, mismatch,
                   stats, further_stats, threads, pdf_out, gt_threshold, nprocs):
    '''
    :param reads_dict:
    :param min_reads: minimun number of reads required in a umi group [5]
    '''


    #list of result stats
    #REVIEW:Returns duplicates (splits file incorrectly?)
    dir_adj_results = Parallel(n_jobs=threads)(delayed(dir_adj_worker)(bundle, threshold, stats, further_stats,
    mismatch, min_reads, nprocs, gt_threshold) for bundle in reads_dict.values())



    stats_pre_df_dict_all = {'UMI': [], 'counts': []}
    stats_post_df_dict = {'UMI': [], 'counts': []}
    pre_cluster_stats = []
    post_cluster_stats = []

    topology_counts = collections.Counter()
    node_counts = collections.Counter()

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
        labels, values = zip(*collections.Counter(dir_adj_results[bundle][3]).items()) #*-operator to unpack the arguments out of a list or tuple
        non_dedup_values = tuple(l*v for l, v in zip(labels, values))

        if min_reads != None:
            # cut_point = count_change(labels, non_dedup_values)

            cut_point = min_reads

            plt.figure()
            plt.bar(labels, non_dedup_values)
            #extract title from read
            # j_nam = re.search('J\d+', dir_adj_results[bundle][0][0].split(' ')[0]).group(0)
            j_nam = dir_adj_results[bundle][0][0].split(' ')[0].split('_')[-3]
            barcode = dir_adj_results[bundle][0][0].split(' ')[0].split('_')[-2]
            plt.title(j_nam + ' ' + barcode + ' Cut point: ' + str(cut_point), ha='center') #need to put name to know which bundle J is being processed
            my_plot = plt.axvline(cut_point, linestyle='dashed', linewidth=2).get_figure()
            pdf_out.savefig(my_plot)


        num_input_all += dir_adj_results[bundle][9] #num_input

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
            stats_pre_df_dict_all['UMI'].extend(dir_adj_results[bundle][10]['UMI']) #stats_pre_df_dict umi + read
            stats_pre_df_dict_all['counts'].extend(dir_adj_results[bundle][10]['counts']) #stats_pre_df_dict umi counts

            # pre_cluster_stats.append(pre_average_distance)

            # collect post-dudupe stats
            #v_seq + umi
            post_cluster_umis = [qname.split(' ')[-1] for qname in dir_adj_results[bundle][0]] #reads are just qnames
            stats_post_df_dict['UMI'].extend(dir_adj_results[bundle][2]) #final_umis
            stats_post_df_dict['counts'].extend(dir_adj_results[bundle][3]) #umi_counts

            # post_average_distance = get_average_umi_distance(post_cluster_umis)
            # post_cluster_stats.append(post_average_distance)

        if further_stats:
            for c_type, count in dir_adj_results[bundle][7].most_common(): #from the most common to the least
                topology_counts[c_type] += count
            for c_type, count in dir_adj_results[bundle][8].most_common():
                node_counts[c_type] += count

    return [stats_pre_df_dict_all, stats_post_df_dict, pre_cluster_stats, post_cluster_stats,
    topology_counts, node_counts, num_input_all, num_output, low_gt_reads, corrected_reads,
    low_gt_corrected_reads, low_umi_count]





def dir_adj_bundle(reads_dict, low_umi_out, out, threshold, min_reads, mismatch,
                   stats, further_stats, threads, pdf_out, gt_threshold):
    '''
    :param reads_dict:
    :param min_reads: minimun number of reads required in a umi group [5]
    '''

    stats_pre_df_dict = {'UMI': [], 'counts': []}
    stats_post_df_dict = {'UMI': [], 'counts': []}
    pre_cluster_stats = []
    post_cluster_stats = []

    topology_counts = collections.Counter()
    node_counts = collections.Counter()

    num_input, num_output = 0, 0

    low_gt_reads = 0
    corrected_reads = 0
    low_gt_corrected_reads = 0
    low_umi_count = 0

    bundle_count = 1
    total_bundles = len(reads_dict.values())


    for bundle in reads_dict.values(): #bundle of v_seq + umi and read
        print('Processing bundle', bundle_count, 'out of', total_bundles)
        bundle_count += 1
        # print(list(bundle.keys())[-1])
        # print(bundle.keys()) #reads at location with different J / orientation
        reads, consensus, final_umis, umi_counts, low_gt, corrected, low_gt_corrected, topologies, nodes =\
        run_dir_adj(bundle, threshold=threshold, stats=stats, further_stats=further_stats,
                    mismatches=mismatch, nprocs=threads, gt_threshold=gt_threshold)
        # print(reads, umis, umi_counts, topologies, nodes)

        #plot umi_counts and determine cutoff
        labels, values = zip(*collections.Counter(umi_counts).items()) #*-operator to unpack the arguments out of a list or tuple
        non_dedup_values = tuple(l*v for l, v in zip(labels, values))

        if min_reads != None:
            # cut_point = count_change(labels, non_dedup_values)

            cut_point = min_reads

            plt.figure()
            plt.bar(labels, non_dedup_values)
            #extract title from read
            j_nam = re.search('J\d+', reads[0].split(' ')[0]).group(0)
            plt.title(j_nam + ' Cut point: ' + str(cut_point), ha='center') #need to put name to know which bundle J is being processed
            my_plot = plt.axvline(cut_point, linestyle='dashed', linewidth=2).get_figure()
            pdf_out.savefig(my_plot)


        num_input += sum([bundle[umi]["count"] for umi in bundle])

        assert len(reads)  == len(consensus) == len(umi_counts), 'Reads, consensus and counts differ in length'

        #remove low umi counts 1-5
        indx = 0
        for count in umi_counts:
            if count <= cut_point:
                #write out fasta
                low_umi_out.write(reads[indx].split(' ')[0].replace('@', '>') +'\n' + consensus[indx] + '\n')
                low_umi_count += 1
            else:
                #write out fasta
                out.write(reads[indx].split(' ')[0].replace('@', '>') + '\n' + consensus[indx] + '\n')
                num_output += 1
            indx += 1

        if stats:
            low_gt_reads += low_gt
            corrected_reads += corrected
            low_gt_corrected_reads += low_gt_corrected

            # collect pre-dudupe stats
            stats_pre_df_dict['UMI'].extend(bundle) #umi + read
            stats_pre_df_dict['counts'].extend([bundle[UMI]['count'] for UMI in bundle]) #umi counts

            pre_average_distance = get_average_umi_distance(bundle.keys()) #v_seq + umi
            pre_cluster_stats.append(pre_average_distance)

            # collect post-dudupe stats
            #v_seq + umi
            post_cluster_umis = [qname.split(' ')[-1] for qname in reads] #reads are just qnames
            stats_post_df_dict['UMI'].extend(final_umis)
            stats_post_df_dict['counts'].extend(umi_counts)

            post_average_distance = get_average_umi_distance(post_cluster_umis)
            post_cluster_stats.append(post_average_distance)

        if further_stats:
            for c_type, count in topologies.most_common(): #from the most common to the least
                topology_counts[c_type] += count
            for c_type, count in nodes.most_common():
                node_counts[c_type] += count

    return [stats_pre_df_dict, stats_post_df_dict, pre_cluster_stats, post_cluster_stats,
    topology_counts, node_counts, num_input, num_output, low_gt_reads, corrected_reads,
    low_gt_corrected_reads, low_umi_count]

######### Stats ################

# def get_average_umi_distance(umis):
#     '''Get average distance of dir_adj cluster UMI's
#     '''
#
#     if len(umis) == 1:
#         return -1
#
#     dists = [edit_distance(x.encode('utf-8'), y.encode('utf-8')) for
#              x, y in itertools.combinations(umis, 2)]
#
#     return float(sum(dists))/(len(dists))



def aggregateStatsDF(stats_df):
    ''' return a data from with aggregated counts per UMI'''

    agg_df_dict = {}

    agg_df_dict['total_counts'] = stats_df.pivot_table(
        columns="UMI", values="counts", aggfunc=np.sum)

    agg_df_dict['median_counts'] = stats_df.pivot_table(
        columns="UMI", values="counts", aggfunc=np.median)

    agg_df_dict['times_observed'] = stats_df.pivot_table(
        columns="UMI", values="counts", aggfunc=len)

    return pd.DataFrame(agg_df_dict)




def make_bundle(fastq, ignore_umi, spe, ignore_j, skip_unclear):
    '''bundle reads
    '''
    unclear_skip = 0
    #Deduplication without alignment
    reads_dict = collections.defaultdict(lambda: collections.defaultdict(dict))


    with general.file_open(fastq) as jv_in:
        lines = jv_in.read().splitlines()
        for item in general.fastq_parse(lines):
            qname = item[0]
            seq = item[1]
            thrd = item[2]
            qual = item[3]

            #Get UMI from qname
            if ignore_umi:
                umi = ''
            else:
                umi = qname.split(' ')[0].split('_')[-1]

            barcode = qname.split(' ')[0].split('_')[-2]

            v_seq = seq[-8:]

            if skip_unclear:
                if 'unclear' in qname:
                    unclear_skip += 1
                    continue

            j_idn = qname.split('_')[-3]


            if ignore_j:
                key = barcode
            else:
                key = j_idn + barcode

            # vs = v_seq + umi
            # my_set = {'TGAGGTCCCAGCGG', 'TGAGGTCCACCCGA', 'TGAGGTCCACGCGA', 'TGAGGTCCACGCGG', 'TGAGGTCCCCGCGG'}
            # if vs in my_set:
            #     print(key, vs, seq)

            try:
                reads_dict[key][v_seq + umi]['count'] += 1
                reads_dict[key][v_seq + umi]['seq'].update([seq]) #add all the seqs for consensus
                # read_groups_br1[reads_dict_br1[key][v_seq + umi]['read'].split(' ')[0]].add(seq)
            except KeyError:
                reads_dict[key][v_seq + umi]['count'] = 1
                reads_dict[key][v_seq + umi]['read'] = qname.split(' ')[0] + ' ' + v_seq #put v_seq into qname for stats
                reads_dict[key][v_seq + umi]['seq'] = collections.Counter([seq]) #add all the seqs for consensus
                # read_groups_br1[reads_dict_br1[key][v_seq + umi]['read'].split(' ')[0]].add(seq)

    return (reads_dict, unclear_skip)



class deduplicate:
    '''Deduplicate using J, V start and UMI
    barcode 1 GACTCGT  barcode 2 CTGCTCCT
    '''

    def __init__(self, file_directory, br1='GACTCGT', br2='CTGCTCCT', assembled=False):
        '''
        :param file_directory: where files are
        :param single: seperate V and J (not assembled)
        '''
        self.file_directory = file_directory
        self.br1 = br1
        self.br2 = br2
        self.out_dir = ''

        if assembled: #assembled reads (PEAR)
            self.jv_fastq_br1 = glob.glob(self.file_directory + '/*all_jv*' + br1)[0]
            self.jv_fastq_br2 = glob.glob(self.file_directory + '/*all_jv*' + br2)[0]
            self.jv_prefix = ''
            self.header_postfix_jv = ''
            # self.out_qnames = set()

        else:
            self.v_fastq_br1 = glob.glob(self.file_directory + '/*all_V*' + br1)[0] #unlist
            self.v_fastq_br2 = glob.glob(self.file_directory + '/*all_V*' + br2)[0]
            self.j_fastq_br1 = glob.glob(self.file_directory + '/*all_J*' + br1)[0]
            self.j_fastq_br2 = glob.glob(self.file_directory + '/*all_J*' + br2)[0]
            # self.tmp_dir = ''
            self.v_prefix = ''
            self.j_prefix = ''
            self.header_postfix_v = ''
            self.header_postfix_j = ''
            self.j_reads = collections.defaultdict(list)



    def create_dirs(self, out_dir=None):
        '''#Create directories
        '''

        dir_main = Path(os.path.abspath(self.v_fastq_br1)).parents[1] #1 dir up, create outside of preclean directory
        self.v_prefix = re.split('(_all)', os.path.basename(self.v_fastq_br1))[0]
        # self.v_prefix_br2 = re.split('(_all)', os.path.basename(self.v_fastq_br2))[0]
        self.j_prefix = re.split('(_all)', os.path.basename(self.j_fastq_br1))[0]


        #temp dir for bowtie output etc.
        # self.tmp_dir = str(dir_main) + '/' + str(self.v_prefix) + '_UMI_tmp'

        #final output dir
        if out_dir == None:
            self.out_dir = str(dir_main) + '/' + str(self.v_prefix) + '_Deduplicated'
        else:
            self.out_dir = os.path.abspath(out_dir)


        #Create directories
        # try:
        #     os.mkdir(self.tmp_dir)
        # except FileExistsError:
        #     print('Directory', self.tmp_dir, 'already exists')
        #     pass
        try:
            os.mkdir(self.out_dir)
        except FileExistsError:
            print('Directory', self.out_dir, 'already exists')
            pass


        #Get 1:N:0:GCCAAT from header
        with open(self.v_fastq_br1, 'r') as f:
                first_line = f.readline()
                self.header_postfix_v = first_line.split(' ')[1].rstrip()

        with open(self.j_fastq_br1, 'r') as f:
                first_line = f.readline()
                self.header_postfix_j = first_line.split(' ')[1].rstrip() #remove newline



    def create_dirs_assembled(self, out_dir=None):
        '''#Create directories
        '''

        dir_main = Path(os.path.abspath(self.jv_fastq_br1)).parents[1] #1 dir up, create outside of preclean directory
        self.jv_prefix = re.split('(_all)', os.path.basename(self.jv_fastq_br1))[0]


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
        with open(self.jv_fastq_br1, 'r') as f:
                first_line = f.readline()
                self.header_postfix_jv = first_line.split(' ')[1].rstrip()




    def load_j_reads(self):
        '''Put J reads into dictionary
        '''
        # print('Loading J reads')
        with general.file_open(self.j_fastq_br1) as j_br1, general.file_open(self.j_fastq_br2) as j_br2:
            lines = j_br1.read().splitlines()
            for item in general.fastq_parse(lines):
                title = item[0]
                seq = item[1]
                qual = item[3]

                self.j_reads[title.split(' ')[0][1:]] = [seq, qual]

            lines = j_br2.read().splitlines()
            for item in general.fastq_parse(lines):
                title = item[0]
                seq = item[1]
                qual = item[3]

                self.j_reads[title.split(' ')[0][1:]] = [seq, qual]



    def v_start_j_umi_dedup(self, cores, spe='mmu', ignore_umi=False, verbose=True, plot=False,
                            stats=True, dedup_all=False, ignore_j=False, mapq_filter=1, skip_unclear=True, keep=False):
        '''Determine start position of v reads
        some elemets inspired by umi_tools by tom smith cagt

        '''
        # print(self.v_fastq_br1, self.j_fastq_br1)
        # print(self.v_fastq_br2)
        # print('Starting deduplication')

        assert self.j_reads, 'Run load_j_reads first!'
        #Align v end

        igh = presets.prs(spe).igh() #exclude things mapping elsewhere in genome


        print('Aligning barcode 1 ', self.br1)
        br1_bowtie2 = general.bowtie2()
        #trim 7 barcode + 6 umi
        br1_bowtie2.align_single(fastq=self.v_fastq_br1, nthreads=cores, trim5='13', spe=spe, verbose=verbose, samtools_mapq=mapq_filter)

        if plot:
            br1_bowtie2.plot(plot_region=igh[0] + ':' + str(igh[1]) + '-' + str(igh[2]), spe=spe)

        if dedup_all:
            sam_v_br1 = br1_bowtie2.pysam_out(fetch=True)
        else:
            sam_v_br1 = br1_bowtie2.pysam_out(region=igh, fetch=True)

        sam_algn_v_br1 = br1_bowtie2.pysam_out(algn=True)

        br1_bowtie2.del_tmp()



        print('Aligning barcode 2 ', self.br2)
        br2_bowtie2 = general.bowtie2()
        #trim 8 barcode + 6 umi
        br2_bowtie2.align_single(fastq=self.v_fastq_br2, nthreads=cores, trim5='14', spe=spe, verbose=verbose, samtools_mapq=mapq_filter)

        if plot:
            br2_bowtie2.plot(plot_region=igh[0] + ':' + str(igh[1]) + '-' + str(igh[2]), spe=spe)

        if dedup_all:
            sam_v_br2 = br2_bowtie2.pysam_out(fetch=True)
        else:
            sam_v_br2 = br2_bowtie2.pysam_out(region=igh, fetch=True)

        # sam_algn_v_br2 = br2_bowtie2.pysam_out(algn=True)

        br2_bowtie2.del_tmp()


        #Don't have to worry about unmapped reads (bowtie2 --no-unal)/
        #paired reads (single end)/
        #low quality aligned reads (should set samtools MAPQ for that)


        #br1
        reads_dict_br1 = collections.defaultdict(lambda: collections.defaultdict(lambda: collections.defaultdict(dict)))

        for read in sam_v_br1:

            if read.is_unmapped:
                continue

            if ignore_umi:
                umi = ''
            else:
                umi = read.qname.split('_')[-1]

            j_idn = read.qname.split('_')[-3]

            if skip_unclear:
                if 'unclear' in j_idn:
                    continue

            #Cigar types
            # M 0 alignment match (can be a sequence match or mismatch)
            # I 1 insertion to the reference
            # D 2 deletion from the reference
            # N 3 skipped region from the reference
            # S 4 soft clipping (clipped sequences present in SEQ)
            # H 5 hard clipping (clipped sequences NOT present in SEQ)
            # P 6 padding (silent deletion from padded reference)
            # = 7 sequence match
            # X 8 sequence mismatch

            if read.is_reverse:
                pos = read.aend
                if read.cigar[-1][0] == 4: #cigar type = soft clipping
                    pos = pos + read.cigar[-1][1]

            else: #Skip forwards reads
                continue
                # pos = read.pos
                # if read.cigar[0][0] == 4:
                #     pos = pos - read.cigar[0][1]



            #Don't need to include reverse
            # key = (j_idn, read.is_reverse)
            if ignore_j:
                key = ''
            else:
                key = j_idn

            #Make bundles
            try:
                reads_dict_br1[pos][key][umi]['count'] += 1
            except KeyError:
                reads_dict_br1[pos][key][umi]['count'] = 1
                reads_dict_br1[pos][key][umi]['read'] = read

            # else:
            #
            #     if reads_dict[pos][key][umi]["read"].mapq > read.mapq:
            #         continue
            #
            #     if reads_dict[pos][key][umi]["read"].mapq < read.mapq:
            #         reads_dict[pos][key][umi]["read"] = read
            #         read_counts[pos][key][umi] = 0
            #         continue

        # if stats:
        # # set up arrays to hold stats data
        #     stats_pre_df_dict = {"UMI": [], "counts": []}
        #     stats_post_df_dict = {"UMI": [], "counts": []}
        #     pre_cluster_stats = []
        #     post_cluster_stats = []
        #     pre_cluster_stats_null = []
        #     post_cluster_stats_null = []
        #     topology_counts = collections.Counter()
        #     node_counts = collections.Counter()
        #     read_gn = random_read_generator(infile.filename, chrom=options.chrom)




        ##################

        #br2
        reads_dict_br2 = collections.defaultdict(lambda: collections.defaultdict(lambda: collections.defaultdict(dict)))

        for read in sam_v_br2:

            if read.is_unmapped:
                continue

            if ignore_umi:
                umi = ''
            else:
                umi = read.qname.split('_')[-1]

            j_idn = read.qname.split('_')[-3]

            if skip_unclear:
                if 'unclear' in j_idn:
                    continue

            if read.is_reverse:
                pos = read.aend
                if read.cigar[-1][0] == 4: #cigar type = soft clipping
                    pos = pos + read.cigar[-1][1]


            else: #Skip forward reads
                continue
                # pos = read.pos
                # if read.cigar[0][0] == 4:
                #     pos = pos - read.cigar[0][1]

            if ignore_j:
                key = ''
            else:
                key = j_idn
            #Make bundles
            try:
                reads_dict_br2[pos][key][umi]['count'] += 1
            except KeyError:
                reads_dict_br2[pos][key][umi]['count'] = 1
                reads_dict_br2[pos][key][umi]['read'] = read

        ########################

        #TODO: implement stats!!
        if stats:
            # collect pre-dudupe stats
            stats_pre_df_dict['UMI'].extend(bundle)
            stats_pre_df_dict['counts'].extend(
                [bundle[UMI]['count'] for UMI in bundle])

            # collect post-dudupe stats
            post_cluster_umis = [x.qname.split("_")[-1] for x in reads]
            stats_post_df_dict['UMI'].extend(umis)
            stats_post_df_dict['counts'].extend(umi_counts)

            average_distance = get_average_umi_distance(post_cluster_umis)
            post_cluster_stats.append(average_distance)

        #  if stats:
        #
        #         # collect pre-dudupe stats
        #         stats_pre_df_dict['UMI'].extend(bundle)
        #         stats_pre_df_dict['counts'].extend(
        #             [bundle[UMI]['count'] for UMI in bundle])
        #
        #         # collect post-dudupe stats
        #         post_cluster_umis = [x.qname.split("_")[-1] for x in reads]
        #         stats_post_df_dict['UMI'].extend(umis)
        #         stats_post_df_dict['counts'].extend(umi_counts)
        #
        #         average_distance = get_average_umi_distance(post_cluster_umis)
        #         post_cluster_stats.append(average_distance)
        #
        #         cluster_size = len(post_cluster_umis)
        #         random_umis = read_gn.getUmis(cluster_size)
        #         average_distance_null = get_average_umi_distance(random_umis)
        #         post_cluster_stats_null.append(average_distance_null)
        #
        #         if options.further_stats:
        #             for c_type, count in topologies.most_common():
        #                 topology_counts[c_type] += count
        #             for c_type, count in nodes.most_common():
        #                 node_counts[c_type] += count


        #Write bith barcodes into same file
        num_input_br1, num_output_br1 = 0, 0
        num_input_br2, num_output_br2 = 0, 0
        dj_end = presets.prs(spe).dj()[2]
        line_fmt = "@{0!s}\n{1!s}\n+\n{2!s}\n"
        dj_reads = 0
        v_reads = 0

        # print(reads_dict.keys())
        # with pysam.AlignmentFile(self.tmp_dir + '/' + self.v_prefix_br1 + '_' + self.br1 + '_dedup.bam', "wb", template=sam_algn_v_br1) as out_file:
        with open(self.out_dir + '/' + self.v_prefix + '_dedup_V.fastq', 'w') as v_out, \
        open(self.out_dir + '/' + self.v_prefix + '_dedup_DJ.fastq', 'w') as dj_out, \
        open(self.out_dir + '/' + self.j_prefix + '_dedup_V.fastq', 'w') as j_v_out, \
        open(self.out_dir + '/' + self.j_prefix + '_dedup_DJ.fastq', 'w') as j_dj_out, \
        pysam.AlignmentFile(self.out_dir + '/' + self.v_prefix + '_dedup_V.bam', 'wb', template=sam_algn_v_br1) as v_bam, \
        pysam.AlignmentFile(self.out_dir + '/' + self.v_prefix + '_dedup_DJ.bam', 'wb', template=sam_algn_v_br1) as dj_bam:
            for pos in reads_dict_br1:
                 # print(reads_dict[pos].keys())
                #br1
                for bundle in reads_dict_br1[pos].values(): #bundle of umi + read
                    # print(bundle.keys()) #reads at location with different J / orientation
                    reads, umis, umi_counts, topologies, nodes = run_dir_adj(bundle)
                    # print(reads, umis, umi_counts, topologies, nodes)


                    num_input_br1 += sum([bundle[umi]["count"] for umi in bundle])


                    #write out reads
                    for read in reads:
                        #don't need to worry about soft-clipping
                        if read.is_reverse: #double check?
                            pos = read.aend
                        else:
                            continue
                            # pos = read.pos

                        #merge barcodes into single file
                        if pos >= int(dj_end): #V read
                            v_out.write(line_fmt.format(read.qname + ' ' + self.header_postfix_v, general.reverse_complement(read.seq), read.qual[::-1]))
                            num_output_br1 += 1
                            v_reads += 1

                            if keep: #write out bam
                                v_bam.write(read)

                            try:
                                seq, qual = self.j_reads[read.qname]
                                j_v_out.write(line_fmt.format(read.qname + ' ' + self.header_postfix_j, seq, qual))
                            except KeyError:
                                pass

                        else: #DJ read
                            dj_out.write(line_fmt.format(read.qname + ' ' + self.header_postfix_v, general.reverse_complement(read.seq), read.qual[::-1]))
                            num_output_br1 += 1
                            dj_reads += 1

                            if keep: #write out bam
                                dj_bam.write(read)

                            try:
                                seq, qual = self.j_reads[read.qname]
                                j_dj_out.write(line_fmt.format(read.qname + ' ' + self.header_postfix_j, seq, qual))
                            except KeyError:
                                pass

                        # out_file.write(read)
            #br2
            for pos in reads_dict_br2:
                 # print(reads_dict[pos].keys())
                for bundle in reads_dict_br2[pos].values(): #bundle of umi + read
                    # print(bundle.keys()) #reads at location with different J / orientation
                    reads, umis, umi_counts, topologies, nodes = run_dir_adj(bundle)
                    # print(reads, umis, umi_counts, topologies, nodes)

                    num_input_br2 += sum([bundle[umi]['count'] for umi in bundle])

                    #write out reads
                    for read in reads:
                        #don't need to worry about soft-clipping
                        if read.is_reverse: #double check?
                            pos = read.aend
                        else:
                            continue
                            # pos = read.pos

                        #merge barcodes into single file
                        if pos >= int(dj_end): #V read
                            v_out.write(line_fmt.format(read.qname + ' ' + self.header_postfix_v, general.reverse_complement(read.seq), read.qual[::-1]))
                            num_output_br2 += 1
                            v_reads += 1

                            if keep: #write out bam
                                v_bam.write(read)

                            try:
                                seq, qual = self.j_reads[read.qname]
                                j_v_out.write(line_fmt.format(read.qname + ' ' + self.header_postfix_j, seq, qual))
                            except KeyError:
                                pass

                        else: #DJ read
                            dj_out.write(line_fmt.format(read.qname + ' ' + self.header_postfix_v, general.reverse_complement(read.seq), read.qual[::-1]))
                            num_output_br2 += 1
                            dj_reads += 1

                            if keep: #write out bam
                                dj_bam.write(read)

                            try:
                                seq, qual = self.j_reads[read.qname]
                                j_dj_out.write(line_fmt.format(read.qname + ' ' + self.header_postfix_j, seq, qual))
                            except KeyError:
                                pass

                        # out_file.write(read)


        print('Number of input reads barcode 1:', num_input_br1)
        print('Number of output reads barcode 1:', num_output_br1)

        print('Number of input reads barcode 2:', num_input_br2)
        print('Number of output reads barcode 2:', num_output_br2)

        print('Number of V reads:', v_reads)
        print('Number of DJ reads:', dj_reads)

        # print(reads_dict.keys())

        # with pysam.AlignmentFile(self.tmp_dir + '/' + self.v_prefix_br2 + '_' + self.br2 + '_dedup.bam', "wb", template=sam_algn_v_br2) as out_file:
        # with open(self.out_dir + '/' + self.v_prefix_br2 + '_' + self.br2 + '_dedup_V.fastq', 'w') as v_out, \
        # open(self.out_dir + '/' + self.v_prefix_br2 + '_' + self.br2 + '_dedup_DJ.fastq', 'w') as dj_out:
        #





    def v_start_j_umi_dedup_assembled(self, threshold, min_reads, threads, mismatch, gt_threshold, spe='mmu',
                                      ignore_umi=False, stats=False, ignore_j=False, skip_unclear=False,
                                      further_stats=False):
        '''Determine start position of v reads
        some elements inspired by umi_tools by tom smith cagt

        '''

        reads_dict_br1, unclear_skip_br1 = make_bundle(self.jv_fastq_br1, ignore_umi=ignore_umi, spe=spe, ignore_j=ignore_j, skip_unclear=skip_unclear)
        reads_dict_br2, unclear_skip_br2 = make_bundle(self.jv_fastq_br2, ignore_umi=ignore_umi, spe=spe, ignore_j=ignore_j, skip_unclear=skip_unclear)

        # merge the two dict_keys
        reads_dict = {**reads_dict_br1, **reads_dict_br2}

        # count = 0
        # print(self.jv_fastq_br2)
        # print(sorted(reads_dict_br2['J1'].keys())[1:5])
        ########################

        if stats:
            print('Unclear skiped:', unclear_skip_br1 + unclear_skip_br2)
            all_unclear =  unclear_skip_br1 + unclear_skip_br2
            logging.info('Unclear skiped:' +  str(all_unclear))
        # set up arrays to hold stats data
            stats_pre_df_dict_all = {'UMI': [], 'counts': []}
            stats_post_df_dict_all = {'UMI': [], 'counts': []}
            pre_cluster_stats_all = []
            post_cluster_stats_all = []
            # pre_cluster_stats_br2 = []
            # post_cluster_stats_br2 = []

            topology_counts_all = collections.Counter()
            node_counts_all = collections.Counter()
            # read_gn = random_read_generator(infile.filename, chrom=options.chrom)



            num_input_br1, num_output_br1 = 0, 0
            num_input_br2, num_output_br2 = 0, 0
            # line_fmt = "@{0!s}\n{1!s}\n+\n{2!s}\n"
            low_gt_reads_br1, low_gt_reads_br2 = 0, 0
            corrected_reads_br1, corrected_reads_br2 = 0, 0
            low_gt_corrected_reads_br1, low_gt_corrected_reads_br2 = 0, 0
            low_umi_count_br1, low_umi_count_br2 = 0, 0

        #Write both barcodes into same file
        #Can't split into DJ and V

        # print(reads_dict.keys())
        # with pysam.AlignmentFile(self.tmp_dir + '/' + self.v_prefix_br1 + '_' + self.br1 + '_dedup.bam', "wb", template=sam_algn_v_br1) as out_file:
        with open(self.out_dir + '/' + self.jv_prefix + '_dedup.fasta', 'w') as jv_out, \
        open(self.out_dir + '/' + self.jv_prefix + '_low_umi.fasta', 'w') as low_umi_out, \
        PdfPages(self.out_dir + '/' + self.jv_prefix + '_histogram.pdf') as pdf:

            #run br1 and br2 side by side
            if len(reads_dict) >= threads:
                nprocs = 1
            else:
                nprocs = int(threads/len(read_dict)) #how many unused cores are available?

            #br1
            stats_pre_df_dict, stats_post_df_dict, pre_cluster_stats, post_cluster_stats, \
            topology_counts, node_counts, num_input, num_output, low_gt_reads, corrected_reads, \
            low_gt_corrected_reads, low_umi_count=\
            dir_adj_bundle_parallel(reads_dict, low_umi_out, jv_out, threshold=threshold, min_reads=min_reads,
                           mismatch=mismatch, gt_threshold=gt_threshold,
                           stats=stats, further_stats=further_stats, threads=threads, pdf_out=pdf, nprocs=nprocs)
            # dir_adj_bundle(reads_dict_br1, low_umi_out, jv_out, threshold=threshold, min_reads=min_reads,
            #                mismatch=mismatch, gt_threshold=gt_threshold,
            #                stats=stats, further_stats=further_stats, threads=threads, pdf_out=pdf)

            #stats
            if stats:
                num_input_br1 += num_input
                num_output_br1 += num_output
                low_gt_reads_br1 += low_gt_reads
                corrected_reads_br1 += corrected_reads
                low_gt_corrected_reads_br1 += low_gt_corrected_reads
                low_umi_count_br1 += low_umi_count

                stats_pre_df_dict_all.update(stats_pre_df_dict)
                stats_post_df_dict_all.update(stats_post_df_dict)

                pre_cluster_stats_all.extend(pre_cluster_stats)
                post_cluster_stats_all.extend(post_cluster_stats)

                topology_counts_all.update(topology_counts)
                node_counts_all.update(node_counts)

            # #br2
            # stats_pre_df_dict, stats_post_df_dict, pre_cluster_stats, post_cluster_stats, \
            # topology_counts, node_counts, num_input, num_output, low_gt_reads, corrected_reads, \
            # low_gt_corrected_reads, low_umi_count=\
            # dir_adj_bundle_parallel(reads_dict_br2, low_umi_out, jv_out, threshold=threshold, min_reads=min_reads,
            #                mismatch=mismatch, gt_threshold=gt_threshold,
            #                stats=stats, further_stats=further_stats, threads=threads, pdf_out=pdf, nprocs=threads)
            # # dir_adj_bundle(reads_dict_br2, low_umi_out, jv_out, threshold=threshold, min_reads=min_reads,
            # #                mismatch=mismatch, gt_threshold=gt_threshold,
            # #                stats=stats, further_stats=further_stats, threads=threads, pdf_out=pdf)
            #
            #
            #
            # #stats
            # if stats:
            #     num_input_br2 += num_input
            #     num_output_br2 += num_output
            #     low_gt_reads_br2 += low_gt_reads
            #     corrected_reads_br2 += corrected_reads
            #     low_gt_corrected_reads_br2 += low_gt_corrected_reads
            #     low_umi_count_br2 += low_umi_count
            #
            #     stats_pre_df_dict_all.update(stats_pre_df_dict)
            #     stats_post_df_dict_all.update(stats_post_df_dict)
            #
            #     pre_cluster_stats_all.extend(pre_cluster_stats)
            #     post_cluster_stats_all.extend(post_cluster_stats)
            #
            #     topology_counts_all.update(topology_counts)
            #     node_counts_all.update(node_counts)


        if stats:
            print('Number of input reads:', num_input_br1)
            print('Number of output reads:', num_output_br1)
            logging.info('Number of input reads:' + str(num_input_br1))
            logging.info('Number of output reads:' + str(num_output_br1))
            # print('Number of input reads barcode 2:', num_input_br2)
            # print('Number of output reads barcode 2:', num_output_br2)

            print('Number of clusters with low ratio discarded:' + str(low_gt_reads_br1))
            # print('Number of clusters with low ratio discarded barcode 2:', low_gt_reads_br2)
            logging.info('Number of clusters with low ratio discarded:' + str(low_gt_reads_br1))
            print('Number of directional-adjacency corrected clusters:', corrected_reads_br1)
            # print('Number of directional-adjacency corrected clusters barcode 2:', corrected_reads_br2)
            logging.info('Number of directional-adjacency corrected clusters:' + str(corrected_reads_br1))
            print('Number of corrected clusters with low ratio discarded:', low_gt_corrected_reads_br1)
            # print('Number of corrected clusters with low ratio discarded barcode 2:', low_gt_corrected_reads_br2)
            logging.info('Number of corrected clusters with low ratio discarded:' +  str(low_gt_corrected_reads_br1))
            print('Number of low UMI count groups:', low_umi_count_br1)
            # print('Number of low UMI count groups barcode 2:', low_umi_count_br2)
            logging.info('Number of low UMI count groups:' + str(low_umi_count_br1))
            print('Topology:', topology_counts)
            print('Node count:', node_counts)
            logging.info('Topology:' + str(topology_counts))
            logging.info('Node count:' + str(node_counts))
            ##########################################

            #From UMI_tools
            stats_pre_df = pd.DataFrame(stats_pre_df_dict_all)
            stats_post_df = pd.DataFrame(stats_post_df_dict_all)

            # print(pd.DataFrame.from_dict(collections.Counter(stats_post_df_dict['counts']), orient='index').reset_index())

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
            # print(UMI_counts_df)
            # TS - if count value not observed either pre/post-dedup,
            # merge will leave an empty cell and the column will be cast as a float
            # see http://pandas.pydata.org/pandas-docs/dev/missing_data.html
            # --> Missing data casting rules and indexing
            # so, back fill with zeros and convert back to int
            UMI_counts_df = UMI_counts_df.fillna(0).astype(int)

            UMI_counts_df.to_csv(self.out_dir + '/' + self.jv_prefix + '_per_umi_per_position.tsv', sep='\t')

            ##########################

             # aggregate stats pre/post per UMI
            agg_pre_df = aggregateStatsDF(stats_pre_df)
            agg_post_df = aggregateStatsDF(stats_post_df)

            agg_df = pd.merge(agg_pre_df, agg_post_df, how='left',
                              left_index=True, right_index=True,
                              sort=True, suffixes=['_pre', '_post'])

            # TS - see comment above regarding missing values
            agg_df = agg_df.fillna(0).astype(int)
            agg_df.to_csv(self.out_dir + '/' + self.jv_prefix + '_per_umi.tsv', sep="\t")

            #########################

            # bin distances into integer bins (average_umi_distance)
            # max_ed = int(max(map(max, [pre_cluster_stats_all,
            #                            post_cluster_stats_all,
            #                            ])))
            #
            # cluster_bins = range(-1, int(max_ed) + 2)
            #
            # def bin_clusters(cluster_list, bins=cluster_bins):
            #     ''' take list of floats and return bins'''
            #     return np.digitize(cluster_list, bins, right=True)
            #
            # def tallyCounts(binned_cluster, max_edit_distance):
            #     ''' tally counts per bin '''
            #     return np.bincount(binned_cluster, minlength=max_edit_distance + 3)
            #
            # pre_cluster_binned = bin_clusters(pre_cluster_stats_all)
            # post_cluster_binned = bin_clusters(post_cluster_stats_all)
            #
            # edit_distance_df = pd.DataFrame({
            #     'unique': tallyCounts(pre_cluster_binned, max_ed),
            #     'dir_adj': tallyCounts(post_cluster_binned, max_ed),
            #     'edit_distance': cluster_bins})
            #
            # # TS - set lowest bin (-1) to "Single_UMI"
            # edit_distance_df['edit_distance'][0] = 'Single_UMI'
            #
            # edit_distance_df.to_csv(self.out_dir + '/' + self.jv_prefix + '_edit_distance.tsv', index=False, sep="\t")
        #End of from UMI_tools



    # def write_assembled(self):
    #     '''Write out deduplicated assembled reads
    #     '''
    #
    #     #Parse jv fastq and write out reads into same file
    #     with open(self.out_dir + '/' + self.jv_prefix + '_dedup.fastq', 'w') as jv_out:
    #         with general.file_open(self.jv_fastq_br1) as jv_br1:
    #             lines = jv_br1.read().splitlines()
    #             for item in general.fastq_parse(lines):
    #                 qname = item[0]
    #                 seq = item[1]
    #                 thrd = item[2]
    #                 qual = item[3]
    #
    #                 if qname.split(' ')[0] in self.out_qnames:
    #                     jv_out.write(qname + '\n' + seq + '\n' + thrd + '\n' + qual + '\n')
    #
    #         with general.file_open(self.jv_fastq_br2) as jv_br2:
    #             lines = jv_br2.read().splitlines()
    #             for item in general.fastq_parse(lines):
    #                 qname = item[0]
    #                 seq = item[1]
    #                 thrd = item[2]
    #                 qual = item[3]
    #
    #                 if qname.split(' ')[0] in self.out_qnames:
    #                     jv_out.write(qname + '\n' + seq + '\n' + thrd + '\n' + qual + '\n')



    def plot_after(self, spe='mmu'):
        '''Plot bam after deduplication
        '''

        plot_igh_br1 = ['Rscript', os.path.dirname(os.path.realpath(general.__file__)) +
        '/' + 'plot_igh.R','-o', self.out_dir + '/' + self.v_prefix + '_br1_dedup_coverage_V.pdf',
        '-n', 'br1', '--genome', presets.prs(spe).genome(), '-r', presets.prs(spe).v_region(),
        '-b', self.tmp_dir + '/' + self.v_prefix + '_' + self.br1 + '_sorted_dedup.bam']

        plot_igh_br2 = ['Rscript', os.path.dirname(os.path.realpath(general.__file__)) +
        '/' + 'plot_igh.R','-o', self.out_dir + '/' + self.v_prefix + '_br2_dedup_coverage_V.pdf',
        '-n', 'br2', '--genome', presets.prs(spe).genome(), '-r', presets.prs(spe).v_region(),
        '-b', self.tmp_dir + '/' + self.v_prefix + '_' + self.br2 + '_sorted_dedup.bam']

        subprocess.call(plot_igh_br1)
        subprocess.call(plot_igh_br2)


    # def delete_tmp(self, keep=False):
    #     '''Delete temporary files
    #     '''
    #     assert self.tmp_dir, 'Run create_dirs first!'
    #     if not keep:
    #         shutil.rmtree(self.tmp_dir)





def parse_args():
    parser = argparse.ArgumentParser(description='BabrahamLinkON Deduplicate')

    parser.add_argument('--input_dir', dest='in_dir', type=str, required=True, help='Input directory (created for/by preclean)')
    parser.add_argument('--cores', dest='nthreads', default=1, type=int, help='Number of cores to use (if aligning), default: 1')
    parser.add_argument('--species', dest='species', default='mmu', type=str, help='Which species (mmu, hsa), default: mmu')
    parser.add_argument('--br1', dest='br1', default='GACTCGT', type=str, help='Default: GACTCGT')
    parser.add_argument('--br2', dest='br2', default='CTGCTCCT', type=str, help='Default: CTGCTCCT')
    parser.add_argument('--verbose', action='store_true', help='Print detailed progress')
    parser.add_argument('--out', dest='out_dir', type=str, help='Output directory, default: creates Deduplicated in main directory')
    parser.add_argument('--ignore_umi', action='store_true', help='Deduplicate without using UMI')
    parser.add_argument('--ignore_j', action='store_true', help='Deduplicate without using J end identity')


    parser.add_argument('--keep', action='store_true', help='Keep temporary files (good for troubleshooting)')
    parser.add_argument('--mapq', dest='mapq', default=1, type=int, help='Skip alignments with MAPQ smaller than INT [1]')
    parser.add_argument('--plot', action='store_true', help='Plot V region before and after deduplication')
    parser.add_argument('--dedup_full', action='store_true', help='Do not exclude deduplication of reads outside of VDJ region')

    parser.add_argument('--stats', action='store_true', help='Output stats from UMI deduplication [False]')
    parser.add_argument('--further-stats', action='store_true', help='Output more stats from UMI deduplication [False]')
    parser.add_argument('--skip_unclear', action='store_true', help='Skip unclear J reads [False]')
    parser.add_argument('--assembled', action='store_true', help='Assembled reads are being provided as input (from PEAR) [False]')
    parser.add_argument('--mismatch', dest='mismatch', type=int, default=5, help='Number of mismatches allowed in consensus sequence comparison [5]')
    parser.add_argument('--threshold', dest='threshold', type=int, default=1, help='Number of mismatches allowed in UMI [1]')
    parser.add_argument('--min_reads', dest='minreads', type=int, default=5, help='Minimum number of reads in UMI group [5]')
    parser.add_argument('--gt_ratio', dest='gtratio', type=float, default=1, help='Ratio of good to total reads to mark UMI group as early PCR error 0-1 [1]')

    opts = parser.parse_args()

    return opts



def main():

    #argparse
    opts = parse_args()

    # old_stdout = sys.stdout
    # with open('message.log','w') as log_file:
    #     sys.stdout = log_file
    #     print "this will be written to message.log"
    #     sys.stdout = old_stdout



    dedup = deduplicate(file_directory=opts.in_dir, br1=opts.br1, br2=opts.br2, assembled=opts.assembled)
    if opts.assembled:
        dedup.create_dirs_assembled(out_dir=opts.out_dir)

        logging.basicConfig(level=logging.DEBUG, filename=dedup.out_dir +'/' + dedup.jv_prefix + '.log', filemode='a+',
                            format='%(asctime)-15s %(levelname)-8s %(message)s')

        print('Starting deduplication')
        dedup.v_start_j_umi_dedup_assembled(threshold=opts.threshold, min_reads=opts.minreads, threads=opts.nthreads,
                                            mismatch=opts.mismatch, gt_threshold=opts.gtratio, spe=opts.species, stats=opts.stats,
                                            further_stats=opts.further_stats, ignore_umi=opts.ignore_umi,
                                            ignore_j=opts.ignore_j, skip_unclear=opts.skip_unclear)
        # print('Writing out')
        # dedup.write_assembled()

    else:
        dedup.create_dirs(out_dir=opts.out_dir)

        logging.basicConfig(level=logging.DEBUG, filename=dedup.out_dir +'/' + dedup.v_prefix + '.log', filemode='a+',
                            format='%(asctime)-15s %(levelname)-8s %(message)s')

        print('Loading J reads')
        dedup.load_j_reads()
        print('Starting deduplication')
        dedup.v_start_j_umi_dedup(cores=opts.nthreads, spe=opts.species, verbose=opts.verbose,
                                  plot=opts.plot, stats=opts.stats, ignore_umi=opts.ignore_umi,
                                  dedup_all=opts.dedup_full, ignore_j=opts.ignore_j, mapq_filter=opts.mapq,
                                  skip_unclear=opts.skip_unclear, keep=opts.keep)
        # dedup.write_dj_region_file(spe=opts.species)
        # dedup.seperate_dj()

        # if opts.plot:
        #     dedup.plot_after(spe=opts.species)

        # dedup.delete_tmp(keep=opts.keep)


if __name__ == "__main__":
    main()
