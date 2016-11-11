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
from babrahamlinkon import general
from skbio import DNA, TabularMSA
import multiprocessing
import math

import pyximport

### Modified from UMI tools: ###

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

#Use Cpython edit_distance
# def edit_dist(first, second):
#     ''' returns the edit distance/hamming distances between
#     its two arguements '''

#     dist = sum([not a == b for a, b in zip(first, second)])
#     return dist


def breadth_first_search(node, adj_list, removed):
    '''collapse key and value of dict
    '''
    searched = set()
    found = set()
    queue = set()
    queue.update((node,))
    found.update((node,))

    while len(queue)>0:
        node=(list(queue))[0]
        # Some nodes may have been removed in a previous search
        # through the network. If so, we don't want to include them in
        # multiple connected components. Quicker to check with list of
        # previously removed nodes than to re-generate the adj_list
        # for each step (important for consensus generation)
        if node not in removed:
            for connected_node in adj_list[node]:
                if connected_node not in removed:
                    found.update(adj_list[node])
                    queue.update(adj_list[node])

        searched.update((node,))
        queue.difference_update(searched)

    return found


######### read loss #########


def consensus_difference(seq_counter, differences=5):
    '''Read loss analysis
    :param seq_counter: Counter object with sequences
    :param differences: number of differences from consensus allowed
    '''

    #Convert Counter into fasta input for kalign
    seq_fasta = ''
    count = 0
    for item in seq_counter:
        clust_items = list(item.elements()) #Reverse Counter
        for seq in clust_items:
            seq_fasta = seq_fasta + '>seq' + str(count) + '\n' + seq + '\n'
            count += 1

    #Can't get consensus from single sequence, return single sequence
    if count == 1:
        return (1, list(seq_counter[0].elements())[0])

    #Multiple sequence alignment
    #http://msa.sbc.su.se/cgi-bin/msa.cgi
    kalign_cmd = ['kalign', '-f', 'fasta']


    p = subprocess.Popen(kalign_cmd, stdout=subprocess.PIPE, stdin=subprocess.PIPE, stderr=subprocess.STDOUT)
    kalign_stdout = p.communicate(input=seq_fasta.encode('utf-8'))[0]

    #parse kalign output and derive consensus seq
    head, sep, tail = kalign_stdout.partition(b'>') #remove head kalign intro

    alignment = sep.decode() + tail.decode()
    alignment = re.sub('\n', '', alignment)
    alignment = list(filter(None, re.split('(>seq\d+)', alignment)))


    #Convert into TabularMSA scikit bio object
    seq_dict = collections.defaultdict()

    for item in general.fasta_parse(alignment):
        qname = item[0]
        seq = item[1]
        # print(seq)
        seq_dict[qname] = DNA(seq)


    msa = TabularMSA.from_dict(seq_dict)

    cons_seq = str(msa.consensus()) #REVIEW: Experimental; keep an eye on development

    good = 0
    total = 0
    for seq in iter(msa):
        total += 1

        consensus_diff = edit_distance(str(seq).encode('utf-8'), cons_seq.encode('utf-8')) #how many mismatches present

        if consensus_diff <= int(differences):
            good += 1

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


def get_adj_list_directional_adjacency(umis, counts, threshold=1):
    ''' identify all umis within the hamming distance threshold (1 is best)
    and where the counts of the first umi is > (2 * second umi counts)-1'''

    return {umi: [umi2 for umi2 in umis if
                  edit_distance(umi.encode('utf-8'),
                                umi2.encode('utf-8')) == threshold and
                  counts[umi] >= (counts[umi2]*2)-1] for umi in umis}


######## "get_connected_components" methods ##########

def get_connected_components_adjacency(umis, graph, counts):
    ''' find the connected UMIs within an adjacency dictionary'''

    found = list()
    components = list()

    for node in sorted(graph, key=lambda x: counts[x], reverse=True):
        if node not in found:
            component = breadth_first_search(node, graph, found)
            found.extend(component)
            components.append(component)

    return components


######## "reduce_clusters" methods ##########


def reduce_clusters_single(bundle, clusters, counts, stats, mismtch):
    ''' collapse clusters down to the UMI which accounts for the cluster
    using the adjacency dictionary and return the list of final UMIs'''

    reads = []
    final_umis = []
    umi_counts = []
    low_gt = 0
    corrected = 0
    cons_not_match = 0
    for cluster in clusters:
        umi_in_cluster = 0
        #Consensus for read loss filter
        out = [] #This is what want a consensus of
        for umi in cluster: #if not corrected should have single UMI in cluster
            umi_in_cluster += 1
            try:
                out.append(bundle[umi]['seq'])
            except KeyError:
                print('UMI not in bundle')

        assert len(out) != 0, 'No sequence from umi'

        gt_ratio, consensus_seq = consensus_difference(out, differences=mismtch)


        #REVIEW: use consensus sequence or highest count sequence?
        if gt_ratio == 1: #100% agreement
            #Parent umi = highest count umi which account for the cluster
            parent_umi = get_best_higher_counts(cluster, counts)
            if consensus_seq != bundle[parent_umi]['seq'].most_common()[0][0]:
                cons_not_match += 1
                if cons_not_match < 100:
                    print(consensus_seq, bundle[parent_umi]['seq'].most_common()[0][0])
            reads.append(bundle[parent_umi]["read"])
        else:
            low_gt += 1
            if umi_in_cluster > 1: #contains umi with 1 error and sequence doesn't match
                corrected += 1
            continue

        if stats:
            final_umis.append(parent_umi)
            #Number of UMI's in the cluster (how many have been collapsed)
            umi_counts.append(sum([counts[x] for x in cluster]))
    print('Not match:', cons_not_match)
    #list of reads, final umi's used, list of umi counts within clusters
    return reads, final_umis, umi_counts, low_gt, corrected



def reduce_clusters_single_parallel(bundle, clusters, counts, nprocs, stats=False):
    ''' collapse clusters down to the UMI which accounts for the cluster
    using the adjacency dictionary and return the list of final UMIs'''

    def worker(bundle, clusters, counts, out_q):

        inter_results = {'reads':[], 'final_umis':[], 'umi_counts':[], 'low_gt':0, 'corrected':0}

        for cluster in clusters:
            umi_in_cluster = 0


            #Consensus for read loss filter
            out = [] #This is what want a consensus of
            for umi in cluster:
                umi_in_cluster += 1
                out.append(bundle[umi]['seq'])
            # print(out)
            gt_ratio, consensus_seq = consensus_difference(out)

            if umi_in_cluster > 1:
                    inter_results['corrected'] += 1
            #REVIEW: use consensus sequence or highest count sequence?
            if gt_ratio == 1:
                #Parent umi = highest count umi which account for the cluster
                parent_umi = get_best_higher_counts(cluster, counts)

                inter_results['reads'].append(bundle[parent_umi]["read"])
            else:
                inter_results['low_gt'] += 1

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
        p = multiprocessing.Process(target=worker, args=(bundle, clusters[cluster_chunk * i:cluster_chunk * (i+1)], counts, out_q))
        procs.append(p)
        p.start()

    # Collect all results into a single result dict
    reads = []
    final_umis = []
    umi_counts = []
    low_gt = 0
    corrected = 0
    for i in range(nprocs):
        out_dict = out_q.get()
        reads.extend(out_dict['reads'])
        final_umis.extend(out_dict['final_umis'])
        umi_counts.extend(out_dict['umi_counts'])
        low_gt += out_dict['low_gt']
        corrected += out_dict['corrected']

    # Wait for all worker processes to finish
    for p in procs:
        p.join()


    #list of reads, final umi's used, list of umi counts within clusters
    return reads, final_umis, umi_counts, low_gt, corrected



######### Call ################


def run_dir_adj(bundle, threshold, stats, further_stats, mismatches):
    #threshold=1, stats=True, further_stats=True, mismatches=5
    umis = bundle.keys()

    len_umis = [len(x) for x in umis]
    assert max(len_umis) == min(len_umis), (
        "not all umis are the same length(!):  %d - %d" % (
            min(len_umis), max(len_umis)))

    counts = {umi: bundle[umi]["count"] for umi in umis}

    adj_list = get_adj_list_directional_adjacency(umis, counts, threshold)

    clusters = get_connected_components_adjacency(umis, adj_list, counts)

    reads, final_umis, umi_counts, low_gt, corrected = reduce_clusters_single(
        bundle, clusters, counts, stats, mismatches)

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

    return reads, final_umis, umi_counts, low_gt, corrected, topologies, nodes



def get_average_umi_distance(umis):
    '''Get average distance of dir_adj cluster UMI's
    '''

    if len(umis) == 1:
        return -1

    dists = [edit_distance(x.encode('utf-8'), y.encode('utf-8')) for
             x, y in itertools.combinations(umis, 2)]

    return float(sum(dists))/(len(dists))



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


### End of UMI_tools ###


def bundle(fastq, ignore_umi, spe, ignore_j):
    '''bundle reads
    '''
    unclear_skip = 0
    #Deduplication without alignment
    reads_dict = collections.defaultdict(lambda: collections.defaultdict(dict))

    # read_groups_br1 = collections.defaultdict(set) #####################
    # read_groups_br2 = collections.defaultdict(set)

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

            #Use read length as proxy for v start
            # read_len = len(seq)

            v_seq = seq[-10:]

            # if skip_unclear:
            if 'unclear' in qname:
                unclear_skip += 1
                continue

            j_idn = qname.split('_')[-3]

            try: #REVIEW: This will also get rid of unclear, do we want junction sequence
                if spe == 'mmu':
                    j_len = len(general.species(spe).replace()[j_idn])

                    add = 5 #For J2 and J3
                    if j_idn == 'J1': #to get beyond primer sequence
                        add = 8
                    elif j_idn == 'J4':
                        add = 6

                    j_seq = seq[j_len+add:j_len+10+add]

                    trim_seq = seq[j_len+add:][:-10] #Used for consensus - trim J primer and 10bp v used with umi
                #TODO: throw exception here

            except KeyError:
                continue

            if ignore_j:
                key = ''
                # key = read_len
                # key = str(read_len)
            else:
                key = j_idn
                # key = str(read_len)

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
            self.out_qnames = set()

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
        print('Starting deduplication')

        assert self.j_reads, 'Run load_j_reads first!'
        #Align v end

        igh = general.species(spe).igh() #exclude things mapping elsewhere in genome


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
        dj_end = general.species(spe).dj()[2]
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





    def v_start_j_umi_dedup_assembled(self, threshold, spe='mmu', ignore_umi=False, verbose=False,
                            stats=False, ignore_j=False, skip_unclear=False, further_stats=False,
                            mismatch=5):
        '''Determine start position of v reads
        some elements inspired by umi_tools by tom smith cagt

        '''
        # print(self.v_fastq_br1, self.j_fastq_br1)
        # print(self.v_fastq_br2)
        # print('Starting deduplication')


        # igh = general.species(spe).igh() #exclude things mapping elsewhere in genome

        unclear_skip = 0

        #Deduplication without alignment
        #br1
        reads_dict_br1 = collections.defaultdict(lambda: collections.defaultdict(dict))

        read_groups_br1 = collections.defaultdict(set) #####################
        read_groups_br2 = collections.defaultdict(set)

        with general.file_open(self.jv_fastq_br1) as jv_in:
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

                #Use read length as proxy for v start
                read_len = len(seq)

                #Test with first 10bp

                # v_seq = seq[-10:]
                v_seq = seq[-10:]

                if skip_unclear:
                    if 'unclear' in qname:
                        unclear_skip += 1
                        continue

                j_idn = qname.split('_')[-3]

                try: #REVIEW: This will also get rid of unclear, do we want junction sequence
                    if spe == 'mmu':
                        j_len = len(general.species(spe).replace()[j_idn])

                        add = 5 #For J2 and J3
                        if j_idn == 'J1': #to get beyond primer sequence
                            add = 8
                        elif j_idn == 'J4':
                            add = 6

                        j_seq = seq[j_len+add:j_len+10+add]

                        trim_seq = seq[j_len+add:][:-10] #Used for consensus - trim J primer and 10bp v used with umi
                    #TODO: throw exception here

                except KeyError:
                    continue

                #Don't need to include reverse
                # key = (j_idn, read.is_reverse)
                if ignore_j:
                    key = ''
                    # key = read_len
                    # key = str(read_len)
                else:
                    key = j_idn
                    # key = str(read_len)

                #Make bundles
                # try:
                #     reads_dict_br1[key][j_seq + v_seq + umi]['count'] += 1
                # except KeyError:
                #     reads_dict_br1[key][j_seq + v_seq + umi]['count'] = 1
                #     reads_dict_br1[key][j_seq + v_seq + umi]['read'] = qname.split(' ')[0] + ' ' + j_seq + v_seq + umi #put v_seq into qname for stats

                try:
                    reads_dict_br1[key][v_seq + umi]['count'] += 1
                    reads_dict_br1[key][v_seq + umi]['seq'].update([trim_seq]) #add all the seqs for consensus
                    # read_groups_br1[reads_dict_br1[key][v_seq + umi]['read'].split(' ')[0]].add(seq)
                except KeyError:
                    reads_dict_br1[key][v_seq + umi]['count'] = 1
                    reads_dict_br1[key][v_seq + umi]['read'] = qname.split(' ')[0] + ' ' + v_seq #put v_seq into qname for stats
                    reads_dict_br1[key][v_seq + umi]['seq'] = collections.Counter([trim_seq]) #add all the seqs for consensus
                    # read_groups_br1[reads_dict_br1[key][v_seq + umi]['read'].split(' ')[0]].add(seq)
                # try:
                #     reads_dict_br1[key][umi]['count'] += 1
                # except KeyError:
                #     reads_dict_br1[key][umi]['count'] = 1
                #     reads_dict_br1[key][umi]['read'] = qname.split(' ')[0] + ' ' + v_seq #put v_seq into qname for stats



        ##################

        #br2
        reads_dict_br2 = collections.defaultdict(lambda: collections.defaultdict(dict))

        with general.file_open(self.jv_fastq_br2) as jv_in:
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

                #Use read length as proxy for v start
                read_len = len(seq)

                #Test with first 10bp
                # j_seq = seq[:10] #REVIEW: use 10bp beyond J?
                # v_seq = seq[-10:]
                v_seq = seq[-10:]

                if skip_unclear:
                    if 'unclear' in qname:
                        unclear_skip += 1
                        continue

                j_idn = qname.split('_')[-3]

                try:
                    if spe == 'mmu':
                        j_len = len(general.species(spe).replace()[j_idn])

                        #Goes last 5bp of J and 5bp junction/D
                        add = 5 #for J2 and J3
                        if j_idn == 'J1': #to get beyond primer sequence
                            add = 8
                        elif j_idn == 'J4':
                            add = 6

                        j_seq = seq[j_len+add:j_len+10+add]

                        trim_seq = seq[j_len+add:][:-10] #Used for consensus - trim J primer and 10bp v used with umi
                    #TODO: throw exception here if other than mmu

                except KeyError:
                    continue

                if ignore_j:
                    key = ''
                    # key = read_len
                    # key = str(read_len)
                else:
                    key = j_idn
                    # key = str(read_len)

                #Make bundles
                # try:
                #     reads_dict_br2[key][j_seq + v_seq + umi]['count'] += 1
                # except KeyError: #bundle
                #     reads_dict_br2[key][j_seq + v_seq + umi]['count'] = 1
                #     reads_dict_br2[key][j_seq + v_seq + umi]['read'] = qname.split(' ')[0] + ' ' + j_seq + v_seq + umi

                try:
                    reads_dict_br2[key][v_seq + umi]['count'] += 1
                    # read_groups_br2[reads_dict_br2[key][v_seq + umi]['read'].split(' ')[0]].add(seq)
                    reads_dict_br2[key][v_seq + umi]['seq'].update([trim_seq]) #add all the seqs for consensus
                except KeyError:
                    reads_dict_br2[key][v_seq + umi]['count'] = 1
                    reads_dict_br2[key][v_seq + umi]['read'] = qname.split(' ')[0] + ' ' + v_seq #put v_seq into qname for stats
                    reads_dict_br2[key][v_seq + umi]['seq'] = collections.Counter([trim_seq])
                    # read_groups_br2[reads_dict_br2[key][v_seq + umi]['read'].split(' ')[0]].add(seq)
                # try:
                #     reads_dict_br2[key][umi]['count'] += 1
                # except KeyError: #bundle
                #     reads_dict_br2[key][umi]['count'] = 1
                #     reads_dict_br2[key][umi]['read'] = qname.split(' ')[0] + ' ' + v_seq

        # print(reads_dict_br1)
        # print(reads_dict_br2)
        # for item in read_groups_br1.values():
        #     if len(item) > 1:
        #         print(item)
        # for item in read_groups_br2.values():
        #     if len(item) > 1:
        #         print(item)
        # print(read_groups_br1)
        # print(read_groups_br2)
        ########################

        if stats:
            print('Unclear skiped:', unclear_skip)
        # set up arrays to hold stats data
            stats_pre_df_dict = {'UMI': [], 'counts': []}
            stats_post_df_dict = {'UMI': [], 'counts': []}
            pre_cluster_stats = []
            post_cluster_stats = []

            topology_counts = collections.Counter()
            node_counts = collections.Counter()
            # read_gn = random_read_generator(infile.filename, chrom=options.chrom)


        #Write bith barcodes into same file
        #Can't split into DJ and V
        num_input_br1, num_output_br1 = 0, 0
        num_input_br2, num_output_br2 = 0, 0
        line_fmt = "@{0!s}\n{1!s}\n+\n{2!s}\n"
        low_gt_reads_br1, low_gt_reads_br2 = 0, 0
        corrected_reads_br1, corrected_reads_br2 = 0, 0

        # print(reads_dict.keys())
        # with pysam.AlignmentFile(self.tmp_dir + '/' + self.v_prefix_br1 + '_' + self.br1 + '_dedup.bam', "wb", template=sam_algn_v_br1) as out_file:
        with open(self.out_dir + '/' + self.jv_prefix + '_dedup.fastq', 'w') as jv_out, \
        open(self.out_dir + '/' + self.jv_prefix + '_low_umi.txt', 'w') as low_umi_out:
            #br1
            for bundle in reads_dict_br1.values(): #bundle of v_seq + umi and read

                # print(bundle.keys()) #reads at location with different J / orientation
                reads, umis, umi_counts, low_gt, corrected, topologies, nodes = run_dir_adj(bundle, threshold=threshold, stats=stats, further_stats=further_stats, mismatches=mismatch)
                # print(reads, umis, umi_counts, topologies, nodes)
                low_gt_reads_br1 += low_gt
                corrected_reads_br1 += corrected

                num_input_br1 += sum([bundle[umi]["count"] for umi in bundle])

                assert len(reads) == len(umi_counts), 'Reads and counts differ'
                #return reads with low umi counts 1-5 (What are these reads?) only print a list of read names to isolate in igblast
                indx = 0
                for count in umi_counts:
                    if count <= 5:
                        low_umi_out.write('Read:' + reads[indx].split(' ')[0] +'\n' + 'Count:' + str(count) + '\n' + 'UMI:' + umis[indx] + '\n')
                    indx += 1

                #qname into set
                for qname in reads:
                    self.out_qnames.add(qname.split(' ')[0])
                    num_output_br1 += 1

                if stats:
                    # collect pre-dudupe stats
                    stats_pre_df_dict['UMI'].extend(bundle) #umi + read
                    stats_pre_df_dict['counts'].extend([bundle[UMI]['count'] for UMI in bundle]) #umi counts

                    pre_average_distance = get_average_umi_distance(bundle.keys()) #v_seq + umi
                    pre_cluster_stats.append(pre_average_distance)

                    # collect post-dudupe stats
                    #v_seq + umi
                    post_cluster_umis = [qname.split(' ')[-1] for qname in reads] #reads are just qnames
                    stats_post_df_dict['UMI'].extend(umis)
                    stats_post_df_dict['counts'].extend(umi_counts)

                    post_average_distance = get_average_umi_distance(post_cluster_umis)
                    post_cluster_stats.append(post_average_distance)

                if further_stats:
                    for c_type, count in topologies.most_common(): #from the most common to the least
                        topology_counts[c_type] += count
                    for c_type, count in nodes.most_common():
                        node_counts[c_type] += count


            #br2
            for bundle in reads_dict_br2.values(): #bundle of umi + read
                # print(bundle.keys()) #reads at location with different J / orientation
                reads, umis, umi_counts, low_gt, corrected, topologies, nodes = run_dir_adj(bundle, threshold=threshold, stats=stats, further_stats=further_stats, mismatches=mismatch)
                # print(reads, umis, umi_counts, topologies, nodes)
                low_gt_reads_br2 += low_gt
                corrected_reads_br2 += corrected

                num_input_br2 += sum([bundle[umi]['count'] for umi in bundle])

                #return reads with low umi counts 1-5 (What are these reads?) only print a list of read names to isolate in igblast
                indx = 0
                for count in umi_counts:
                    if count <= 5:
                        low_umi_out.write('Read:' + reads[indx].split(' ')[0] +'\n' + 'Count:' + str(count) + '\n' + 'UMI:' + umis[indx] + '\n')
                    indx += 1


                #write out reads
                for qname in reads:
                    self.out_qnames.add(qname.split(' ')[0])
                    num_output_br2 += 1

                if stats:
                    # collect pre-dudupe stats
                    stats_pre_df_dict['UMI'].extend(bundle) #umi + read
                    stats_pre_df_dict['counts'].extend([bundle[UMI]['count'] for UMI in bundle]) #umi counts

                    pre_average_distance = get_average_umi_distance(bundle.keys()) #v_seq + umi
                    pre_cluster_stats.append(pre_average_distance)

                    # collect post-dudupe stats
                    #v_seq + umi
                    post_cluster_umis = [qname.split(' ')[-1] for qname in reads] #reads are just qnames
                    stats_post_df_dict['UMI'].extend(umis)
                    stats_post_df_dict['counts'].extend(umi_counts)

                    post_average_distance = get_average_umi_distance(post_cluster_umis)
                    post_cluster_stats.append(post_average_distance)

                if further_stats:
                    for c_type, count in topologies.most_common(): #from the most common to the least
                        topology_counts[c_type] += count
                    for c_type, count in nodes.most_common():
                        node_counts[c_type] += count

        if verbose:
            print('Number of input reads barcode 1:', num_input_br1)
            print('Number of output reads barcode 1:', num_output_br1)

            print('Number of input reads barcode 2:', num_input_br2)
            print('Number of output reads barcode 2:', num_output_br2)

            print('Number of clusters with low ratio discarded barcode 1:', low_gt_reads_br1)
            print('Number of clusters with low ratio discarded barcode 2:', low_gt_reads_br2)

            print('Number of corrected clusters with low ratio discarded barcode 1:', corrected_reads_br1)
            print('Number of corrected clusters with low ratio discarded barcode 2:', corrected_reads_br2)

        if stats:
            # print('Pre:', stats_pre_df_dict)
            # print('Post:', stats_post_df_dict)
            # print(pre_cluster_stats)
            # print(post_cluster_stats)
            # print(post_cluster_umis_all)

            print('Topology:', topology_counts)
            print('Node count:', node_counts)

            ##########################################

            #From UMI_tools
            stats_pre_df = pd.DataFrame(stats_pre_df_dict)
            stats_post_df = pd.DataFrame(stats_post_df_dict)

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
            max_ed = int(max(map(max, [pre_cluster_stats,
                                       post_cluster_stats,
                                       ])))

            cluster_bins = range(-1, int(max_ed) + 2)

            def bin_clusters(cluster_list, bins=cluster_bins):
                ''' take list of floats and return bins'''
                return np.digitize(cluster_list, bins, right=True)

            def tallyCounts(binned_cluster, max_edit_distance):
                ''' tally counts per bin '''
                return np.bincount(binned_cluster, minlength=max_edit_distance + 3)

            pre_cluster_binned = bin_clusters(pre_cluster_stats)
            post_cluster_binned = bin_clusters(post_cluster_stats)

            edit_distance_df = pd.DataFrame({
                'unique': tallyCounts(pre_cluster_binned, max_ed),
                'dir_adj': tallyCounts(post_cluster_binned, max_ed),
                'edit_distance': cluster_bins})

            # TS - set lowest bin (-1) to "Single_UMI"
            edit_distance_df['edit_distance'][0] = 'Single_UMI'

            edit_distance_df.to_csv(self.out_dir + '/' + self.jv_prefix + '_edit_distance.tsv', index=False, sep="\t")
        #End of from UMI_tools



    def write_assembled(self):
        '''Write out deduplicated assembled reads
        '''

        #Parse jv fastq and write out reads into same file
        with open(self.out_dir + '/' + self.jv_prefix + '_dedup.fastq', 'w') as jv_out:
            with general.file_open(self.jv_fastq_br1) as jv_br1:
                lines = jv_br1.read().splitlines()
                for item in general.fastq_parse(lines):
                    qname = item[0]
                    seq = item[1]
                    thrd = item[2]
                    qual = item[3]

                    if qname.split(' ')[0] in self.out_qnames:
                        jv_out.write(qname + '\n' + seq + '\n' + thrd + '\n' + qual + '\n')

            with general.file_open(self.jv_fastq_br2) as jv_br2:
                lines = jv_br2.read().splitlines()
                for item in general.fastq_parse(lines):
                    qname = item[0]
                    seq = item[1]
                    thrd = item[2]
                    qual = item[3]

                    if qname.split(' ')[0] in self.out_qnames:
                        jv_out.write(qname + '\n' + seq + '\n' + thrd + '\n' + qual + '\n')



    def plot_after(self, spe='mmu'):
        '''Plot bam after deduplication
        '''

        plot_igh_br1 = ['Rscript', os.path.dirname(os.path.realpath(general.__file__)) +
        '/' + 'plot_igh.R','-o', self.out_dir + '/' + self.v_prefix + '_br1_dedup_coverage_V.pdf',
        '-n', 'br1', '--genome', general.species(spe).genome(), '-r', general.species(spe).v_region(),
        '-b', self.tmp_dir + '/' + self.v_prefix + '_' + self.br1 + '_sorted_dedup.bam']

        plot_igh_br2 = ['Rscript', os.path.dirname(os.path.realpath(general.__file__)) +
        '/' + 'plot_igh.R','-o', self.out_dir + '/' + self.v_prefix + '_br2_dedup_coverage_V.pdf',
        '-n', 'br2', '--genome', general.species(spe).genome(), '-r', general.species(spe).v_region(),
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

    opts = parser.parse_args()

    return opts



def main():

    #argparse
    opts = parse_args()

    dedup = deduplicate(file_directory=opts.in_dir, br1=opts.br1, br2=opts.br2, assembled=opts.assembled)
    if opts.assembled:
        dedup.create_dirs_assembled(out_dir=opts.out_dir)
        print('Starting deduplication')
        dedup.v_start_j_umi_dedup_assembled(threshold=opts.threshold, spe=opts.species, verbose=opts.verbose,
                                  stats=opts.stats, further_stats=opts.further_stats, ignore_umi=opts.ignore_umi,
                                  ignore_j=opts.ignore_j, skip_unclear=opts.skip_unclear, mismatch=opts.mismatch)
        print('Writing out')
        dedup.write_assembled()

    else:
        dedup.create_dirs(out_dir=opts.out_dir)
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
