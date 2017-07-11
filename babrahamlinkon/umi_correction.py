
from babrahamlinkon import deduplicate
import pyximport
from collections import defaultdict, Counter
import itertools
from copy import deepcopy

################################################################################
#Functions from UMI_tools to do UMI error corrections (some modifications)
################################################################################
# TODO: short option giving different results each run

'''
directional-adjacency
adapted from:
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

from babrahamlinkon._dedup_umi import edit_distance


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
    adj_list = defaultdict(list)
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

    found = set() #change to sets
    components = list()

    for node in sorted(graph, key=lambda x: counts[x], reverse=True):
        if node not in found:
            component = breadth_first_search(node, graph)
            found.update(component)
            components.append(component)

    return components





################################################################################


def consensus_difference(seq_counter_dict, no_msa=True, short=False):
    '''Resolve if umi in multiple networks, check head node to see where it fits best
    :param alignment: output from msa
    :return: number of differences between two umi group consensus sequences
    '''

    if no_msa:

        cntr_1, cntr_2 = list(seq_counter_dict.values())

        lst_lists_1 = [list(item) for item in cntr_1.elements()]
        lst_lists_2 = [list(item) for item in cntr_2.elements()]

        cons_seq_1 = deduplicate.consensus_unequal(lst_lists_1)
        cons_seq_2 = deduplicate.consensus_unequal(lst_lists_2)

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


    else:

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

        #without quality should be sufficient here
        cons_seq_1 = deduplicate.consensus(lst_lists_1)
        cons_seq_2 = deduplicate.consensus(lst_lists_2)


    num_diffs = edit_distance(cons_seq_1.encode('utf-8'), cons_seq_2.encode('utf-8'))

    return num_diffs



def merge_dict(in_dict):
    #Messes up in_dict, need to make a deepcopy
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

#
#
# import pickle
# import os
#
#
# def loop_writer(func): #for resolve_clusters
#     def func_wrapper(bundle, clusters, counts, mismatches, gt_threshold, msa, short):
#         count = len(os.listdir('/media/chovanec/My_Passport/test/'))
#         with open('/media/chovanec/My_Passport/test/save_file' + str(count) + '.pkl', 'wb') as out_pkl:
#             rcluster = func(bundle, clusters, counts, mismatches, gt_threshold, msa, short)
#             pickle.dump(rcluster, out_pkl)
#         return rcluster
#     return func_wrapper
#
#
#
# @loop_writer
def resolve_clusters(bundle, clusters, counts, differences, gt_threshold, no_msa=False, short=False):
    '''
    Which shared nodes belong to which head node (not required if not doing UMI error correction)
    '''
    single_components = []
    network_dict = defaultdict(set)
    cont_comp = 0

    for comp in clusters:
        if len(comp) == 1: #not a cluster
            single_components.append(comp)
        else:
            cont_comp += len(comp)
            ordered_network = sorted(comp, key=lambda x: counts[x], reverse=True)
            network_dict[ordered_network[0]] = set(ordered_network[1:])


    #TODO: change so don't need to make copy
    shared_list = merge_dict(deepcopy(network_dict))


    connect_comp = []
    added_key = set()
    duplicates_removed = 0
    head_node_removed = 0
    single_share = 0
    umbigous = 0


    remove_dict = defaultdict(set)
    #Process shared one at a time (create subset of network_dict)
    for shared in shared_list:
        #get all the networks that share components
        #skip next part if head node not sharing item with other head nodes
        if len(shared) > 1:

            remove_dict = defaultdict(set)

            for head_1, head_2 in itertools.combinations(list(shared),2):

                shared_values = network_dict[head_1].intersection(network_dict[head_2])

                if len(shared_values) > 0: #head nodes need to share values

                    to_1 = 0
                    to_2 = 0
                    for item in shared_values:

                        resolve_dict_1 = {head_1:bundle[head_1]['seq'], item:bundle[item]['seq']}
                        resolve_dict_2 = {head_2:bundle[head_2]['seq'], item:bundle[item]['seq']}

                        if no_msa:

                            diff_1 = consensus_difference(resolve_dict_1, short=short)
                            diff_2 = consensus_difference(resolve_dict_2, short=short)

                        else:

                            #Align head 1 and value
                            algn_1 = deduplicate.kalign_msa(resolve_dict_1)
                            diff_1 = consensus_difference(algn_1, no_msa=False)

                            #Align head 2 and value
                            algn_2 = deduplicate.kalign_msa(resolve_dict_2)
                            diff_2 = consensus_difference(algn_2, no_msa=False)


                        #which ever is lower asign to that
                        if diff_1 < diff_2:
                            remove_dict[head_2].update([item])
                        elif diff_1 > diff_2:
                            remove_dict[head_1].update([item])
                        elif diff_1 == diff_2: #remove from both
                            umbigous += 1
                            remove_dict[head_2].update([item])
                            remove_dict[head_1].update([item])

                        else:
                            print('Something else')

            for key, value in remove_dict.items():
                network_dict[key].difference_update(value)
        else:
            single_share += 1


    for key, value in network_dict.items():
        if key not in added_key:

            node = set([key]).union(value)
            connect_comp.append(node)
        else:
            duplicates_removed +=1

    out_components = single_components + connect_comp

    return out_components


################################################################################
# UMI error correction using directional adjaceny method
################################################################################



def run_dir_adj(bundle, threshold, stats, mismatches, nprocs, gt_threshold, qual_dict, no_msa, short, cons_no_qual):
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
    rclusters = resolve_clusters(bundle, clusters, counts, mismatches, gt_threshold, no_msa, short)


    print('Reducing clusters')
    # if nprocs == 1: #if no additional cores available
    reads, consensus_seqs, consensus_quals, final_umis, umi_counts, low_gt, corrected, low_gt_corrected, cons_dffs =\
    deduplicate.reduce_clusters_single(bundle, rclusters, counts, stats, mismatches, gt_threshold, qual_dict, no_msa, short, cons_no_qual)


    # print('Unique:', len(set(reads)), 'All:', len(reads))
    #TODO: http://stackoverflow.com/questions/6974695/python-process-pool-non-daemonic/8963618#8963618
    # reads, consensus, final_umis, umi_counts, low_gt, corrected, low_gt_corrected =\
    # reduce_clusters_single_parallel(bundle, rclusters, counts, nprocs, stats, mismatches, gt_threshold)

    # (bundle, clusters, counts, stats, mismtch)

    return reads, consensus_seqs, consensus_quals, final_umis, umi_counts, low_gt, corrected, low_gt_corrected, cons_dffs
