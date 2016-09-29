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
from babrahamlinkon import general

import pyximport

### From UMI tools: ###

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


def breadth_first_search(node, adj_list):
    '''collapse key and value of dict
    '''
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

#TODO: derive consensus sequence instead of highest count seq?
def get_best_higher_counts(cluster, counts):
    ''' return the UMI with the highest counts'''
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
            component = breadth_first_search(node, graph)
            found.extend(component)
            components.append(component)

    return components


######## "reduce_clusters" methods ##########


def reduce_clusters_single(bundle, clusters, adj_list, counts, stats=False):
    ''' collapse clusters down to the UMI which accounts for the cluster
    using the adjacency dictionary and return the list of final UMIs'''

    reads = []
    final_umis = []
    umi_counts = []

    for cluster in clusters:
        parent_umi = get_best_higher_counts(cluster, counts)
        reads.append(bundle[parent_umi]["read"])

        if stats:
            final_umis.append(parent_umi)
            umi_counts.append(sum([counts[x] for x in cluster]))

    return reads, final_umis, umi_counts


######### Call ################


def run_dir_adj(bundle, threshold=1, stats=True, further_stats=True):

    umis = bundle.keys()

    len_umis = [len(x) for x in umis]
    assert max(len_umis) == min(len_umis), (
        "not all umis are the same length(!):  %d - %d" % (
            min(len_umis), max(len_umis)))

    counts = {umi: bundle[umi]["count"] for umi in umis}

    adj_list = get_adj_list_directional_adjacency(umis, counts, threshold)

    clusters = get_connected_components_adjacency(umis, adj_list, counts)

    reads, final_umis, umi_counts = reduce_clusters_single(
        bundle, clusters, adj_list, counts, stats)

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

    return reads, final_umis, umi_counts, topologies, nodes


### End of UMI_tools ###


class deduplicate:
    '''Deduplicate using J, V start and UMI
    barcode 1 GACTCGT  barcode 2 CTGCTCCT
    '''

    def __init__(self, file_directory, br1='GACTCGT', br2='CTGCTCCT'):
        self.file_directory = file_directory
        self.br1 = br1
        self.br2 = br2
        self.v_fastq_br1 = glob.glob(self.file_directory + '/*all_V*' + br1)[0] #unlist
        self.v_fastq_br2 = glob.glob(self.file_directory + '/*all_V*' + br2)[0]
        self.j_fastq_br1 = glob.glob(self.file_directory + '/*all_J*' + br1)[0]
        self.j_fastq_br2 = glob.glob(self.file_directory + '/*all_J*' + br2)[0]
        self.tmp_dir = ''
        self.out_dir = ''
        self.v_prefix = ''
        self.j_prefix = ''
        # self.v_prefix_br2 = ''
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
        self.tmp_dir = str(dir_main) + '/' + str(self.v_prefix) + '_UMI_tmp'

        #final output dir
        if out_dir == None:
            self.out_dir = str(dir_main) + '/' + str(self.v_prefix) + '_Deduplicated'
        else:
            self.out_dir = os.path.abspath(out_dir)


        #Create directories
        try:
            os.mkdir(self.tmp_dir)
        except FileExistsError:
            print('Directory', self.tmp_dir, 'already exists')
            pass
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


    def load_j_reads(self):
        '''Put J reads into dictionary
        '''
        print('Loading J reads')
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


    def v_start_j_umi_dedup(self, cores, spe='mmu', ignore_umi=False, verbose=True, plot=False, stats=True, dedup_all=False):
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
        br1_bowtie2.align_single(fastq=self.v_fastq_br1, nthreads=cores, trim5='7', spe=spe, verbose=verbose)

        if plot:
            br1_bowtie2.plot(plot_region=igh[0] + ':' + str(igh[1]) + '-' + str(igh[2]), spe=spe)

        if dedup_all:
            sam_v_br1 = br1_bowtie2.pysam_out(fetch=True)
        else:
            sam_v_br1 = br1_bowtie2.pysam_out(region=igh, fetch=True)

        # sam_algn_v_br1 = br1_bowtie2.pysam_out(algn=True)

        br1_bowtie2.del_tmp()



        print('Aligning barcode 2 ', self.br2)
        br2_bowtie2 = general.bowtie2()
        br2_bowtie2.align_single(fastq=self.v_fastq_br2, nthreads=cores, trim5='8', spe=spe, verbose=verbose)

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

            if ignore_umi:
                umi = ''
            else:
                umi = read.qname.split('_')[-1]

            j_idn = read.qname.split('_')[-3]

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
                start = pos

            else:
                pos = read.pos
                if read.cigar[0][0] == 4:
                    pos = pos - read.cigar[0][1]
                start = pos


            # print(start)

            #Make bundles
            try:
                reads_dict_br1[start][(j_idn, read.is_reverse)][umi]['count'] += 1
            except KeyError:
                reads_dict_br1[start][(j_idn, read.is_reverse)][umi]['count'] = 1
                reads_dict_br1[start][(j_idn, read.is_reverse)][umi]['read'] = read

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

            if ignore_umi:
                umi = ''
            else:
                umi = read.qname.split('_')[-1]

            j_idn = read.qname.split('_')[-3]


            if read.is_reverse:
                pos = read.aend
                if read.cigar[-1][0] == 4: #cigar type = soft clipping
                    pos = pos + read.cigar[-1][1]
                start = pos

            else:
                pos = read.pos
                if read.cigar[0][0] == 4:
                    pos = pos - read.cigar[0][1]
                start = pos

            #Make bundles
            try:
                reads_dict_br2[start][(j_idn, read.is_reverse)][umi]['count'] += 1
            except KeyError:
                reads_dict_br2[start][(j_idn, read.is_reverse)][umi]['count'] = 1
                reads_dict_br2[start][(j_idn, read.is_reverse)][umi]['read'] = read

        ########################
        #Write bith barcodes into same file
        num_input_br1, num_output_br1 = 0, 0
        num_input_br2, num_output_br2 = 0, 0
        dj_end = general.species(spe).dj()[2]
        line_fmt = "@{0!s}\n{1!s}\n+\n{2!s}\n"

        # print(reads_dict.keys())
        # with pysam.AlignmentFile(self.tmp_dir + '/' + self.v_prefix_br1 + '_' + self.br1 + '_dedup.bam', "wb", template=sam_algn_v_br1) as out_file:
        with open(self.out_dir + '/' + self.v_prefix + '_dedup_V.fastq', 'w') as v_out, \
        open(self.out_dir + '/' + self.v_prefix + '_dedup_DJ.fastq', 'w') as dj_out, \
        open(self.out_dir + '/' + self.j_prefix + '_dedup_V.fastq', 'w') as j_v_out, \
        open(self.out_dir + '/' + self.j_prefix + '_dedup_DJ.fastq', 'w') as j_dj_out:
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
                        if read.is_reverse:
                            pos = read.aend
                        else:
                            pos = read.pos

                        #merge barcodes into single file
                        if pos >= int(dj_end): #V read
                            v_out.write(line_fmt.format(read.qname + ' ' + self.header_postfix_v, read.seq, read.qual))
                            num_output_br1 += 1

                            try:
                                seq, qual = self.j_reads[read.qname]
                                j_v_out.write(line_fmt.format(read.qname + ' ' + self.header_postfix_j, seq, qual))
                            except KeyError:
                                pass

                        else: #DJ read
                            dj_out.write(line_fmt.format(read.qname + ' ' + self.header_postfix_v,read.seq, read.qual))
                            num_output_br1 += 1

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
                        if read.is_reverse:
                            pos = read.aend
                        else:
                            pos = read.pos

                        #merge barcodes into single file
                        if pos >= int(dj_end): #V read
                            v_out.write(line_fmt.format(read.qname + ' ' + self.header_postfix_v, read.seq, read.qual))
                            num_output_br2 += 1

                            try:
                                seq, qual = self.j_reads[read.qname]
                                j_v_out.write(line_fmt.format(read.qname + ' ' + self.header_postfix_j, seq, qual))
                            except KeyError:
                                pass

                        else: #DJ read
                            dj_out.write(line_fmt.format(read.qname + ' ' + self.header_postfix_v, read.seq, read.qual))
                            num_output_br2 += 1

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

        # print(reads_dict.keys())

        # with pysam.AlignmentFile(self.tmp_dir + '/' + self.v_prefix_br2 + '_' + self.br2 + '_dedup.bam', "wb", template=sam_algn_v_br2) as out_file:
        # with open(self.out_dir + '/' + self.v_prefix_br2 + '_' + self.br2 + '_dedup_V.fastq', 'w') as v_out, \
        # open(self.out_dir + '/' + self.v_prefix_br2 + '_' + self.br2 + '_dedup_DJ.fastq', 'w') as dj_out:
        #




    def write_dj_region_file(self, spe='mmu'):
        '''Required by samtools view -L command
        '''
        #write v_region file, required for samtools -L
        with open(self.tmp_dir + '/' + 'DJ_region.txt', 'w') as v_region:
            dj_r = general.species(spe).dj()
            v_region.write(dj_r[0] + '\t' + dj_r[1] + '\t' + dj_r[2]) #last D + end of igh


    def seperate_dj(self):
        '''Seperate VDJ from DJ recombinations and convert bam to fastq
        '''

        #Sort and index bam files
        samtools_sort_br1 = ['samtools', 'sort', '-O', 'bam', '-o', self.tmp_dir + '/' + self.v_prefix + '_' + self.br1 + '_sorted_dedup.bam', '-T tmp',
                         self.tmp_dir + '/' + self.v_prefix + '_' + self.br1 + '_dedup.bam']
        samtools_sort_br2 = ['samtools', 'sort', '-O', 'bam', '-o', self.tmp_dir + '/' + self.v_prefix + '_' + self.br2 + '_sorted_dedup.bam', '-T tmp',
                         self.tmp_dir + '/' + self.v_prefix + '_' + self.br2 + '_dedup.bam']
        subprocess.call(samtools_sort_br1)
        subprocess.call(samtools_sort_br2)

        samtools_index_br1 = ['samtools', 'index', self.tmp_dir + '/' + self.v_prefix + '_' + self.br1 + '_sorted_dedup.bam']
        samtools_index_br2 = ['samtools', 'index', self.tmp_dir + '/' + self.v_prefix + '_' + self.br2 + '_sorted_dedup.bam']
        subprocess.call(samtools_index_br1)
        subprocess.call(samtools_index_br2)


        #Test presence of DJ_region.txt
        if os.path.isfile(self.tmp_dir + '/' + 'DJ_region.txt'):
            samtools_region_br1 = ['samtools', 'view', '-f 0x10',
                               self.tmp_dir + '/' + self.v_prefix + '_' + self.br1 + '_sorted_dedup.bam',
                              '-L', self.tmp_dir + '/' + 'DJ_region.txt',
                              '-o', self.tmp_dir + '/' + self.v_prefix + '_' + self.br1 + '_dedup_DJ.bam',
                              '-U', self.tmp_dir + '/' + self.v_prefix + '_' + self.br1 + '_dedup_V.bam'] #V with everything else, should be excluded by bowtie region

            subprocess.call(samtools_region_br1)

            samtools_region_br2 = ['samtools', 'view', '-f 0x10',
                               self.tmp_dir + '/' + self.v_prefix + '_' + self.br2 + '_sorted_dedup.bam',
                              '-L', self.tmp_dir + '/' + 'DJ_region.txt',
                              '-o', self.tmp_dir + '/' + self.v_prefix + '_' + self.br2 + '_dedup_DJ.bam',
                              '-U', self.tmp_dir + '/' + self.v_prefix + '_' + self.br2 + '_dedup_V.bam'] #V with everything else, should be excluded by bowtie region

            subprocess.call(samtools_region_br2)
        else:
            raise Exception('Need to run write_dj_region_file first')

        #bam to fastq
        # samtools_bam2fq_DJ = ['samtools', 'bam2fq', self.tmp_dir + '/' + self.v_prefix_br1 + '_' + self.br1 + '_dedup_DJ.bam']
        # samtools_bam2fq_V = ['samtools', 'bam2fq', self.tmp_dir + '/' + self.v_prefix_br1 + '_' + self.br1 + '_dedup_V.bam']
        #
        # with open(self.out_dir + '/' + self.v_prefix_br1 + '_' + self.br1 + '_dedup_V.fastq', 'wb') as out:
        #     subprocess.Popen(samtools_bam2fq_V, stdout=out)
        # with open(self.out_dir + '/' + self.v_prefix_br1 + '_' + self.br1 + '_dedup_DJ.fastq', 'wb') as out:
        #     subprocess.Popen(samtools_bam2fq_DJ, stdout=out)
        #
        #
        # samtools_bam2fq_DJ = ['samtools', 'bam2fq', self.tmp_dir + '/' + self.v_prefix_br2 + '_' + self.br2 + '_dedup_DJ.bam']
        # samtools_bam2fq_V = ['samtools', 'bam2fq', self.tmp_dir + '/' + self.v_prefix_br2 + '_' + self.br2 + '_dedup_V.bam']
        #
        # with open(self.out_dir + '/' + self.v_prefix_br2 + '_' + self.br2 + '_dedup_V.fastq', 'wb') as out:
        #     subprocess.Popen(samtools_bam2fq_V, stdout=out)
        # with open(self.out_dir + '/' + self.v_prefix_br2 + '_' + self.br2 + '_dedup_DJ.fastq', 'wb') as out:
        #     subprocess.Popen(samtools_bam2fq_DJ, stdout=out)




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


    def delete_tmp(self, keep=False):
        '''Delete temporary files
        '''
        assert self.tmp_dir, 'Run create_dirs first!'
        if not keep:
            shutil.rmtree(self.tmp_dir)





def parse_args():
    parser = argparse.ArgumentParser(description='BabrahamLinkON Deduplicate')

    parser.add_argument('--input_dir', dest='in_dir', type=str, required=True, help='Input directory (created for/by preclean)')
    parser.add_argument('--cores', dest='nthreads', default=1, type=int, help='Number of cores to use, default: 1')
    parser.add_argument('--species', dest='species', default='mmu', type=str, help='Which species (mmu, hsa), default: mmu')
    parser.add_argument('--keep', action='store_true', help='Keep temporary files (good for troubleshooting)')
    parser.add_argument('--br1', dest='br1', default='GACTCGT', type=str, help='Default: GACTCGT')
    parser.add_argument('--br2', dest='br2', default='CTGCTCCT', type=str, help='Default: CTGCTCCT')
    parser.add_argument('--verbose', action='store_true', help='Print detailed progress')
    parser.add_argument('--plot', action='store_true', help='Plot V region before and after deduplication')
    parser.add_argument('--out', dest='out_dir', type=str, help='Output directory, default: creates Deduplicated in main directory')
    parser.add_argument('--dedup_full', action='store_true', help='Do not exclude deduplication of reads outside of VDJ region')
    parser.add_argument('--ignore_umi', action='store_true', help='Deduplicate without using UMI')
    parser.add_argument('--stats', action='store_true', help='Output stats from UMI deduplication')


    opts = parser.parse_args()

    return opts



def main():

    #argparse
    opts = parse_args()

    dedup = deduplicate(file_directory=opts.in_dir, br1=opts.br1, br2=opts.br2)
    dedup.create_dirs(out_dir=opts.out_dir)
    dedup.load_j_reads()
    dedup.v_start_j_umi_dedup(cores=opts.nthreads, spe=opts.species, verbose=opts.verbose,
                              plot=opts.plot, stats=opts.stats, ignore_umi=opts.ignore_umi,
                              dedup_all=opts.dedup_full)
    # dedup.write_dj_region_file(spe=opts.species)
    # dedup.seperate_dj()

    # if opts.plot:
    #     dedup.plot_after(spe=opts.species)

    dedup.delete_tmp(keep=opts.keep)


if __name__ == "__main__":
    main()
