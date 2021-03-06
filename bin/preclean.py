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

import sys
import pysam
import subprocess
import os
import re
import argparse
import Levenshtein
import argparse
import glob
from babrahamlinkon import general, presets, bowtie2_wrapper, igblast_wrapper, mispriming_correction, mispriming_error_estimate
from collections import defaultdict
import logging
import shutil
import tempfile
from tqdm import tqdm
import json
from pathlib import Path
from contextlib import ExitStack
from babrahamlinkon.version import __version__


def assemble(V_region, J_region, out_dir, threads=1, prefix=None, short=False):
    '''Assemble reads using pear and join unassembled reads
    http://stackoverflow.com/questions/21953835/run-subprocess-and-print-output-to-logging
    '''

    if prefix==None:
        prefix = os.path.basename(V_region).split('.')[0]

    fp_V_region = os.path.abspath(V_region)
    fp_J_region = os.path.abspath(J_region)

    # dir_path =  os.path.dirname(fp_V_region)

    # prefix = ''.join(re.split('(L00[0-9])', os.path.basename(V_region))[:-1])

    pear_assemble = ['pear', '-f', fp_J_region, '-r', fp_V_region, '-o', out_dir + '/' + prefix, '-j', str(threads)]

    logging.info('Subprocess: "' + ' '.join(pear_assemble) + '"')


    try:
        command_line_process = subprocess.Popen(
            pear_assemble,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
        )

        process_output, _ =  command_line_process.communicate()

        # process_output is now a string, not a file,
        # you may want to do:
        # process_output = StringIO(process_output)
        logging.info(process_output.decode('utf-8'))

    except (OSError, CalledProcessError) as exception:
        logging.info('Exception occured: ' + str(exception))
        logging.info('Subprocess failed')
        return False
    else:
        # no exception was raised
        logging.info('Subprocess finished')



    #will output 4 files
    # subprocess.call(pear_assemble)

    # prefix.assembled.fastq
    # prefix.discarded.fastq
    # prefix.unassembled.forward.fastq
    # prefix.unassembled.reverse.fastq

    out_file_name = out_dir + '/' + prefix
    dict_out = {'assembled':[], 'unassembled':[]}

    if short:

        #Create merge file, basically do a fancy cat
        with open(out_dir + '/' + prefix + '.unassembled.forward.fastq', 'r') as forward, \
        open(out_dir + '/' + prefix + '.assembled.fastq', 'r') as assembled, \
        open(out_dir + '/' + prefix + '.all_J.fastq', 'w') as out_file, \
        open(out_dir + '/' + prefix + '.json', 'w') as json_out:

            #first write out assembled reads
            # out_file.write(assembled.read())

            #write out J unasembled reads
            # out_file.write(forward.read())

            #write out json dict with assembled and unassembled sets
            for qname, seq, thrd, qual in general.fastq_parse(assembled):
                out_file.write(qname + '\n' + seq + '\n' + thrd + '\n' + qual + '\n')
                dict_out['assembled'].append(qname.rstrip('\n').split(' ')[0])
            for qname, seq, thrd, qual in general.fastq_parse(forward):
                out_file.write(qname + '\n' + seq + '\n' + thrd + '\n' + qual + '\n')
                dict_out['unassembled'].append(qname.rstrip('\n').split(' ')[0])
            json.dump(dict_out, json_out) #output dict as json file so assembled and unassembed can be distinguished for merging

    return out_file_name, dict_out



#Don't need to do, will perform in deduplicate (useful for troubleshooting)
# def aligned_to_igh(V_region, thread_num, spe='mmu', write=False, plot=False, verbose=True):
    # '''Extract reads aligning to the IgH
    # :param V_region: R1 fastq file with V end sequences
    # :param spe: which organism
    # :return: set of read qnames that have aligned to the IgH
    # '''
    #
    # igh = general.species(spe).igh()
    # germ = general.species(spe).germline()
    #
    # fp_V_region = os.path.abspath(V_region)
    #
    # #Align with bowtie2
    # run_bowtie2 = bowtie2_wrapper.bowtie2()
    # run_bowtie2.align_single(fastq=fp_V_region, nthreads=thread_num, spe=spe, verbose=verbose)
    #
    # if plot:
    #     run_bowtie2.plot(plot_region=germ[0] + ':' + str(germ[1]) + '-' + str(germ[2]), spe=spe)
    #
    # sam_file = run_bowtie2.pysam_out(region=germ, fetch=True)
    #
    # if write:
    #     run_bowtie2.write(region=germ)
    #
    # run_bowtie2.del_tmp()


    # sam_file = bowtie2_wrapper.align_single(
    # fastq=fp_V_region,
    # nthreads=thread_num,
    # region=igh,
    # write=write,
    # spe=spe,
    # plot=plot,
    # plot_region=germ[0] + ':' + str(germ[1]) + '-' + str(germ[2]),
    # verbose=verbose)

    #Get read names into a set
    # igh_reads = set([read.qname for read in sam_file])
    #
    # return igh_reads


################################################################################
###### Data holder ######
################################################################################


class fastqHolder:
    '''
    Store fastq split by gene (J gene)
    '''
    def __init__(self):
        self.gene_split = defaultdict(set) #split into seperate J's and germline
        # self.preclean_v = defaultdict(set) #read has to be in both gene_split and preclean_v
        self.demultiplex = defaultdict(set)
        # self.germline = defaultdict(set)
        self.misprimed = defaultdict(lambda: dict())
        self.original_id = defaultdict(set)
        self.fastq_split = defaultdict(list)
        self.assembled = ''
        self.anchor_dict = defaultdict()
        self.umi_dict = defaultdict()
        self.j_umi_dict = defaultdict()
        self.rec_skip = set()

    # def add_to_gene(self, gene, qname):
    #     '''Add to gene_split dictionary
    #     :param gene: name of key
    #     :param qname: fastq qname to be stored under key
    #     '''
    #     self.gene_split[gene].add(qname)

    def add_to_misprimed(self, orig_gene, gene, cor_seq, qname):
        '''Add to gene_split dictionary
        :param orig_gene: initial identity of J before mispriming correction
        :param gene: name of key
        :param cor_seq: corrected misprimed sequence to be stored under key
        :param qname: fastq qname to be stored under key
        '''
        self.misprimed[gene][qname] = cor_seq
        self.original_id[qname] = orig_gene


    def write_demultiplex_unassigned(self, fastq_path, gene, out_path):
        '''Write to a fastq (for unasigned reads)
        :param fastq_path: fastq to subset
        :param gene: what to subset by (J gene in split_gene)
        :param out_path: write path
        :param j_end: are these J reads? if yes will correct misprimed sequences
        '''

        with general.file_open(fastq_path) as fq, \
        open(out_path, 'w') as out_file:
            for qname, seq, thrd, qual in general.fastq_parse(fq):
                if qname.split(' ')[0][1:] in self.demultiplex[gene]:
                    #Add J gene identity into qname
                    out_file.write(qname.split(' ')[0] + '_' + gene + ' ' +
                                   ''.join(qname.split(' ')[1:]) + '\n' + seq + '\n' + thrd + '\n' + qual + '\n')



    def write_preclean(self, fastq_path, gene, out_path):
        '''Write fastq file
        :param fastq: fastq to subset
        :param gene: what to subset by (J gene in split_gene)
        :param path: write path
        '''

        with general.file_open(fastq_path) as fq, \
        open(out_path, 'w') as out_file:
            for qname, seq, thrd, qual in general.fastq_parse(fq):
                if qname.split(' ')[0][1:] in self.gene_split[gene]:
                    out_file.write(qname.split(' ')[0] + '_' + gene + ' ' +
                                   ''.join(qname.split(' ')[1:]) + '\n' + seq + '\n' + thrd + '\n' + qual + '\n')



    def write_demultiplex_umi_extract_assembled(self, jv_fastq, gene_list, out_path, an1, an2,
                                                umi_len, beyond_anchor, q_score, verbose, j_len,
                                                in_len):
        '''Write everything into single fastq file and extract umi
        :param fastq: fastq to subset
        :param gene_list: what to subset by (J gene in split_gene)
        :param V_out_path: V out write path
        :param J_out_path: J out write path
        :param an_1: barcode 1 GACTCGT
        :param an_2: barcode 2 CTGCTCCT
        '''

        low_qual_UMI = 0
        an1_count = 0
        an2_count = 0
        gene_count = defaultdict(int)

        with general.file_open(jv_fastq) as jv_fq, \
        open(out_path + '_' + an1, 'w') as an1_out_file, \
        open(out_path + '_' + an2, 'w') as an2_out_file:

            for qname, seq, thrd, qual in general.fastq_parse(jv_fq):
                for gene in gene_list: #J1 J2 J3 J4...
                    if qname.split(' ')[0][1:] in self.demultiplex[gene]:

                        if an1 in gene:
                            anchor_len = 7
                        elif an2 in gene:
                            anchor_len = 8

                        umi = seq[-umi_len:] #last 6 bases NNNNNN (it is in reverse complement configuration compared to original V read)
                        umi += seq[-umi_len-anchor_len-beyond_anchor:-umi_len-anchor_len]
                        #If you don't have anchor you do not have a UMI, but if you have UMI you should not need J end (CDR3) sequences
                        #to further diversify the UMI so no need to add j_len parameter here (for long sequences without anchor UMI)
                        if j_len > 0:
                            umi += seq[:in_len][-j_len:]

                        # umi_dict[qname.split(' ')[0][1:]] = umi #add to dict so can add to J qname too
                        if self.misprimed: #empty dict evals to False
                            try:
                                base_gene = gene.split('_')[0] #get rid of barcode, mispriming dict doesn't have barcode
                                seq = self.misprimed[base_gene][qname.split(' ')[0][1:]] #overwrite seq
                                if len(qual)>=len(seq):
                                    qual = qual[len(qual)-len(seq):] #trim qual if it is longer than corrected seq
                                else:
                                    qual = 'I'*(len(seq)-len(qual)) + qual #if corrected seq longer, append 'I' at the beginning

                                #Don't write original identity, but identity after misprime correction
                                # gene = self.original_id[qname.split(' ')[0][1:]]

                                assert len(seq) == len(qual), 'Sequence and quality length do not match!'
                            except KeyError:
                                pass

                        #if taking sequences beyond the umi, recheck the quality
                        if beyond_anchor > 0:
                            umi_qual = qual[-umi_len:]
                            umi_qual += qual[-umi_len-anchor_len-beyond_anchor:-umi_len-anchor_len]

                            if general.check_qual(umi_qual, q_score):
                                #poor quality skip record
                                low_qual_UMI += 1
                                continue #skip low quality record

                        #Remove UMI and barcode from seq (keeping it in read qname only)
                        #Do not remove sequnce beyond the anchor as it holds value in annotation stage
                        if an1 in gene:
                            an1_out_file.write(qname.split(' ')[0] + '_' + gene + '_' + umi + ' ' +
                            ''.join(qname.split(' ')[1:]) + '\n' + seq[:-(anchor_len+umi_len)] +
                            '\n' + thrd + '\n' + qual[:-(anchor_len+umi_len)] + '\n')
                            an1_count += 1
                            gene_count[gene] += 1
                        else:
                            an2_out_file.write(qname.split(' ')[0] + '_' + gene + '_' + umi + ' ' +
                            ''.join(qname.split(' ')[1:]) + '\n' + seq[:-(anchor_len+umi_len)] +
                            '\n' + thrd + '\n' + qual[:-(anchor_len+umi_len)] + '\n')
                            an2_count += 1
                            gene_count[gene] += 1

        if verbose:
            # print('Low quality UMIs:', low_qual_UMI)
            print(an1, 'reads written out:', an1_count)
            print(an2, 'reads written out:', an2_count)

            for j in gene_count.keys():
                print('Number of', j, 'reads written out:', gene_count[j])

        return [low_qual_UMI, an1_count, an2_count]


    #for short sequences only --------------------------------------------------
    def j_end_umi(self, fastq_path_jv, j_len, in_len):

        with general.file_open(fastq_path_jv) as j_fq:
            for qname, seq, thrd, qual in general.fastq_parse(j_fq):
                self.j_umi_dict[qname.split(' ')[0]] = seq[:in_len][-j_len:]


    def process_v_end(self, fastq_path_v, umi_len, out_path, v_iden_dict, genes, q_score, anchor, beyond_anchor):

        low_qual_UMI = 0
        no_anchor = 0
        v_short = 0
        assembled = 0

        with general.file_open(fastq_path_v) as v_fq:
            for qname, seq, thrd, qual in tqdm(general.fastq_parse(v_fq)):
                if anchor: #if anochor present verify it is correct
                    if seq[umi_len:umi_len+7] == 'GACTCGT':
                        self.anchor_dict[qname.split(' ')[0]] = 'GACTCGT'
                        anchor_len = 7
                    elif seq[umi_len:umi_len+8] == 'CTGCTCCT':
                        self.anchor_dict[qname.split(' ')[0]] = 'CTGCTCCT'
                        anchor_len = 8
                    else:
                        self.rec_skip.add(qname.split(' ')[0][1:])
                        no_anchor += 1
                        continue


                #take sequence beyond the anchor (V start to diversify short UMI)
                if beyond_anchor:
                    umi_seq = seq[:umi_len]
                    umi_seq += seq[umi_len+anchor_len:umi_len+anchor_len+beyond_anchor]
                    umi_qual = qual[:umi_len]
                    #check qual of umi_len (this is for cases when you want to take part of V read beyond anchor)
                    umi_qual += qual[umi_len+anchor_len:umi_len+anchor_len+beyond_anchor]

                    if general.check_qual(umi_qual, q_score):
                        #poor quality skip record
                        low_qual_UMI += 1
                        self.rec_skip.add(qname.split(' ')[0][1:])
                    else:
                        if len(umi_seq) == umi_len+beyond_anchor: #if V read is too short skip it
                            self.umi_dict[qname.split(' ')[0]] = umi_seq
                        else:
                            v_short += 1
                            self.rec_skip.add(qname.split(' ')[0][1:])
                else:
                    umi_seq = seq[:umi_len]
                    umi_qual = qual[:umi_len]
                    if general.check_qual(umi_qual, q_score):
                        #poor quality skip record
                        low_qual_UMI += 1
                        self.rec_skip.add(qname.split(' ')[0][1:])
                    else:
                        if len(umi_seq) == umi_len: #if V read is too short skip it
                            self.umi_dict[qname.split(' ')[0]] = umi_seq
                        else:
                            v_short += 1
                            self.rec_skip.add(qname.split(' ')[0][1:])

        #write out v end for non-assembled reads
        #loop through file again, TODO: figure out a better way

        #remove @ from self.assembled
        assembled_seqs = set([i[1:] for i in self.assembled['assembled']])
        #make single set with all seqs to be written out
        skip_seqs = set.union(self.rec_skip, assembled_seqs)

        with ExitStack() as stack:
            v_fq = stack.enter_context(general.file_open(fastq_path_v))
            if anchor:
                out_GACTCGT = stack.enter_context(open(out_path + '_GACTCGT', 'w'))
                out_CTGCTCCT = stack.enter_context(open(out_path + '_CTGCTCCT', 'w'))
            else:
                v_fq_filt = stack.enter_context(open(out_path, 'w'))

            for qname, seq, thrd, qual in tqdm(general.fastq_parse(v_fq)):

                if len(self.j_umi_dict) > 0:
                    umi = self.j_umi_dict.get(qname.split(' ')[0], 'not_present') + \
                          self.umi_dict.get(qname.split(' ')[0], 'not_present')
                else:
                    umi = self.umi_dict.get(qname.split(' ')[0], 'not_present')

                #skip sequences without a umi
                if 'not_present' in umi:
                    continue
                #UMI is not removed from sequence as it hold value at the annotation stage
                for gene in genes:
                    if 'germline' not in gene and 'other' not in gene:
                        if qname.split(' ')[0][1:] in self.gene_split[gene] or \
                           qname.split(' ')[0][1:] in self.misprimed[gene]:

                            v_iden_out = v_iden_dict.get(qname.split(' ')[0][1:], '')

                            # if qname.split(' ')[0][1:] not in self.rec_skip and qname.split(' ')[0] not in self.assembled['assembled']:
                            if qname.split(' ')[0][1:] not in skip_seqs:
                                if anchor:
                                    if self.anchor_dict.get(qname.split(' ')[0], '') == 'GACTCGT':
                                        out_GACTCGT.write(qname.split(' ')[0] + '_' + gene + '_' + + v_iden_out.upper() + '_GACTCGT_' +
                                                          umi + ' ' + ''.join(qname.split(' ')[1:]) +
                                                          '\n' + seq + '\n' + thrd + '\n' + qual + '\n')
                                    elif self.anchor_dict.get(qname.split(' ')[0],'') == 'CTGCTCCT':
                                        out_CTGCTCCT.write(qname.split(' ')[0] + '_' + gene + '_' + + v_iden_out.upper() + '_CTGCTCCT_' +
                                                           umi + ' ' + ''.join(qname.split(' ')[1:]) +
                                                           '\n' + seq + '\n' + thrd + '\n' + qual + '\n')
                                else:
                                    v_fq_filt.write(qname.split(' ')[0] + '_' + gene + '_' + v_iden_out.upper() + '_' + umi +
                                                    ' ' + ''.join(qname.split(' ')[1:]) +
                                                    '\n' + seq + '\n' + thrd + '\n' + qual + '\n')
                            elif qname.split(' ')[0][1:] in assembled_seqs:
                                assembled += 1


        return low_qual_UMI, no_anchor, v_short, assembled

    #make a decorator to gzip files if specified (TODO: when you have more time)
    # def _write_gzip_or_not(func):
    #     def my_wrapper(self, fastq_path_jv, genes, out_path, v_iden_dict,
    #                              merge, anchor):
    #https://stackoverflow.com/questions/1263451/python-decorators-in-classes

    #TODO: Don't hardcode the anchor sequences
    def write_preclean_short(self, fastq_path_jv, genes, out_path, v_iden_dict,
                             merge, anchor):
        '''Write fastq file and extract UMI (V start) from V fastq
        assumes UMI before anchor is 6N
        :param fastq: fastq to subset (V for UMI extraction)
        :param gene: what to subset by (J gene in split_gene)
        :param path: write path
        '''
        out_reads = 0

        if merge:

            with ExitStack() as stack:
                fq = stack.enter_context(general.file_open(fastq_path_jv))
                if anchor:
                    out_an1 = stack.enter_context(open(out_path + '_GACTCGT', 'w'))
                    out_an2 = stack.enter_context(open(out_path + '_CTGCTCCT', 'w'))
                else:
                    out_file = stack.enter_context(open(out_path, 'w'))

                for qname, seq, thrd, qual in tqdm(general.fastq_parse(fq)):
                    if qname.split(' ')[0][1:] in self.rec_skip:
                        continue

                    v_iden_out = v_iden_dict.get(qname.split(' ')[0][1:], '')

                    if len(self.j_umi_dict) > 0:
                        umi = self.j_umi_dict.get(qname.split(' ')[0], 'not_present') + \
                              self.umi_dict.get(qname.split(' ')[0], 'not_present')
                    else:
                        umi = self.umi_dict.get(qname.split(' ')[0], 'not_present')

                    #skip sequences without a umi
                    if 'not_present' in umi:
                        continue

                    for gene in genes:
                        if 'germline' not in gene and 'other' not in gene:
                            if qname.split(' ')[0][1:] in self.gene_split[gene] or qname.split(' ')[0][1:] in self.misprimed[gene]:
                                if anchor:
                                    if self.anchor_dict.get(qname.split(' ')[0], '') == 'GACTCGT':
                                        out_an1.write(qname.split(' ')[0] + '_' + gene + '_' + v_iden_out.upper() + '_' +
                                                       self.anchor_dict.get(qname.split(' ')[0]) + '_' + umi + ' ' +
                                                       ''.join(qname.split(' ')[1:]) + '\n' + seq + '\n' + thrd + '\n' + qual + '\n')

                                    elif self.anchor_dict.get(qname.split(' ')[0],'') == 'CTGCTCCT':
                                        out_an2.write(qname.split(' ')[0] + '_' + gene + '_' + v_iden_out.upper() + '_' +
                                                       self.anchor_dict.get(qname.split(' ')[0]) + '_' + umi + ' ' +
                                                       ''.join(qname.split(' ')[1:]) + '\n' + seq + '\n' + thrd + '\n' + qual + '\n')
                                else:
                                    out_file.write(qname.split(' ')[0] + '_' + gene + '_' + v_iden_out.upper() + '_' + umi + ' ' +
                                                   ''.join(qname.split(' ')[1:]) + '\n' + seq + '\n' + thrd + '\n' + qual + '\n')
                                out_reads += 1


        else:
            if isinstance(genes, list):
                raise ValueError('Multiple genes supplied without wanting to merge')

            with general.file_open(fastq_path_jv) as fq, open(out_path, 'w') as out_file:
                for qname, seq, thrd, qual in general.fastq_parse(fq):
                    # try:
                    if len(self.j_umi_dict) > 0:
                        umi = self.j_umi_dict.get(qname.split(' ')[0], 'not_present') + \
                              self.umi_dict.get(qname.split(' ')[0], 'not_present')
                    else:
                        umi = self.umi_dict.get(qname.split(' ')[0], 'not_present')

                    #skip sequences without a umi
                    if 'not_present' in umi:
                        continue

                    if qname.split(' ')[0][1:] in self.rec_skip:
                        continue
                    if qname.split(' ')[0][1:] in self.gene_split[genes]:
                        out_file.write(qname.split(' ')[0] + '_' + genes + '_' + umi + ' ' +
                                       ''.join(qname.split(' ')[1:]) + '\n' + seq + '\n' + thrd + '\n' + qual + '\n')
                        out_reads += 1

        return out_reads




################################################################################


#Don't need to remove germline other than to reduce file size (speed up)
#Almost half the data is germline, required for germline mispriming.
def germline(J_region, thread_num, spe='mmu', plot=False, write=False, verbose=True, if_assembled=None):
    '''Extract reads aligning to J genes i.e. germline (used to correct mispriming)
    :param J_region: R2 fastq file with J end sequences
    :param spe: which species
    :param plot: generate pileup plot; default False
    :return: set of read qnames that have aligned to J genes
    '''

    germ = presets.prs(spe).germline() #['chr12', 113428237, 113430474]
    ref_index = presets.prs(spe).bowtie_index()

    fq_data = fastqHolder()

    #add if read has been assembled or not (only for short)
    fq_data.assembled = if_assembled

    #for mm take only first 21 nts (smallest J)
    j_size_ord = sorted(presets.prs(spe).J_seq().values(), key=len)
    shortest_J = len(j_size_ord[0])
    # location_j = presets.prs(spe).J_location()

    fp_J_region = os.path.abspath(J_region)

    #Don't need to align V reads, if germline present will get discarded in deduplicaiton DJ V seperation
    run_bowtie2 = bowtie2_wrapper.bowtie2()
    run_bowtie2.align_single(fastq=fp_J_region, nthreads=thread_num, trim5=str(shortest_J), spe=spe, verbose=verbose)

    if plot:
        plot_file = run_bowtie2.plot(plot_region=germ[0] + ':' + str(germ[1]) + '-' + str(germ[2]), spe=spe)
        if not Path(plot_file).is_file():
            print('No plot produced')

    sam_file_j = run_bowtie2.pysam_out(region=germ, fetch=True)

    if write:
        run_bowtie2.write(region=germ)

    run_bowtie2.del_tmp()

    #get rid of everything aligning to J genes
    count = 0
    for read in sam_file_j:
        count += 1
        if not read.is_reverse: #should get excluded by region constrainst
            fq_data.gene_split['germline_J'].add(read.qname) #add to same dict as normal J's
        elif spe == 'mmuk': #kappa are in different orientation
            fq_data.gene_split['germline_J'].add(read.qname) #add to same dict as normal J's

    if verbose:
        print('Germline reads', len(fq_data.gene_split['germline_J']))
        print('Germline reverse and forward reads', count)
    #Get read names into a set
    # germline_reads = set([read.qname for read in region_reads])

    return fq_data




def preclean_assembled(jv_region, fq_dict_germ, q_score, umi_len, spe='mmu', verbose = True,
                       no_misprime_correct=False, discard_germline=True, fast=False, short=False):
    '''Seperate out J reads and clean up low quality/unknown reads
    :param jv_region: assembled J and V sequences (from PEAR)
    :param fq_dict_germ: fastqHolder object from germline
    :param q_score: Minimum Q score within UMI
    :param umi_len: length of the UMI
    :param spe: which species
    :param verbose: output run startswith
    :param no_misprime_correct: do not run mispriming correction
    :param fast: use only quick initial J identification (without SSW and mispriming correction)
    :return: fastqHolder object
    '''

    #Process jv reads now

    other_J = 0
    unclear_J = 0
    count = 0
    low_qual_seq = 0
    # mapped_to_other = 0
    germ_count = 0

    j_align = mispriming_correction.SSW_align()

    Js = presets.prs(spe).J_seq()
    ref = j_align.reference(spe, quick_align=fast)

    #Germline set
    germline = fq_dict_germ.gene_split['germline_J']


    # j_size_ord = sorted(presets.prs(spe).J_seq().values(), key=len)
    # shortest_J = len(j_size_ord[0])
    # print(germline)
    fp_jv_region = os.path.abspath(jv_region)

    misprimed_corrected = 0
    primer_corrected = 0
    first_filter = 0
    low_qual_UMI = 0
    with general.file_open(fp_jv_region) as jvr:
        for item in tqdm(general.fastq_parse(jvr)):
            title = item[0]
            seq = item[1]
            thrd = item[2]
            qual = item[3]
            count +=1
            # if verbose:
            #     if count % 1000000 == 0:
            #         print('Processed', count, 'J sequences')

            #Remove low quality and short reads and reads that were low quality
            if short:
                min_seq_len = 60
            else:
                longest_j = max([len(v) for v in presets.prs(spe).J_seq().values()])
                min_seq_len = 60+8+longest_j+umi_len  #+ size of anchor and UMI 13 or 14 + longest J


            if '#' in qual or len(seq) < min_seq_len:
                low_qual_seq +=1
                fq_dict_germ.gene_split['other_J'].add(title.split(' ')[0][1:]) #removes @ as well
                continue
            #Skip germline reads
            elif discard_germline:
                if title.split(' ')[0][1:] in germline:
                    germ_count += 1
                    # fq_dict_germ.split_gene['germline_J'].add(title.split(' ')[0][1:])
                    continue

            #Check quality of UMI Q30, skip low qual UMI reads
            #Demultiplex will make sure the last 6 bases are actually the UMI!
            if not short:
                umi_qual = qual[-umi_len:]
                if general.check_qual(umi_qual, q_score=q_score):
                    low_qual_UMI += 1
                    fq_dict_germ.gene_split['other_J'].add(title.split(' ')[0][1:])
                    continue

            first_filter += 1
            #allows 2 mismatches for 21nt seq
            #Will return touple is misprimed
            J = j_align.align(ref, seq, no_misprim_cor=no_misprime_correct, quick_align=fast, spe=spe) #might need to adjust
            #@HWI-1KL136:214:D1MR5ACXX:5:1101:1821:2154 1:N:0:TGACCA remove 1:N:0:TGACCA

            #Output from align
            # ['J1', 'J1', 'CCCTGTGCCCCAGACATCGAAGTACCAGTA']
            # ['J2', 'J2', 'AGTGGTGCCTTGGCCCCAGTAGTCAAA']
            # [initial, misprime corrected, corrected seq]

            if isinstance(J, list): #If a list, then it is a useful J, else a string returned which goes into discard pile
                if len(J) == 1: #not misprime correcting/fast align
                    fq_dict_germ.gene_split[J[0]].add(title.split(' ')[0][1:])
                else:
                    # if J[0] == J[1]:#don't need to do anything
                    if seq == J[2]: #seq == corrected_seq
                        fq_dict_germ.gene_split[J[0]].add(title.split(' ')[0][1:]) #without seq
                    else: #new J, old J, qname, corrected_seq, 5bp beyond primer
                        if J[0] == J[1]: #correct sequences with error in the primer sequence here
                            fq_dict_germ.add_to_misprimed(J[0], J[1], J[2], title.split(' ')[0][1:]) #with seq also removes @
                            primer_corrected += 1
                        else:
                            misprimed_corrected += 1
                            fq_dict_germ.add_to_misprimed(J[0], J[1], J[2], title.split(' ')[0][1:]) #with seq also removes @


            else:
                if J is 'other':
                    other_J += 1
                    fq_dict_germ.gene_split['other_J'].add(title.split(' ')[0][1:])
                elif 'unclear' in J:
                    unclear_J += 1
                    fq_dict_germ.gene_split[J].add(title.split(' ')[0][1:])



    #TODO: clean up whats printed
    if verbose:
        for key in sorted(fq_dict_germ.gene_split.keys()): #order should be always the same
            print('Number of reads in', key, ':', len(fq_dict_germ.gene_split[key]))

        print('Number of mispriming corrected reads:', misprimed_corrected)
        print('Number of primer corrected reads:', primer_corrected)
        print('Number of reads processed:', count)
        print('Number of total unclear JV reads:', unclear_J)
        print('Number of total mispriming other JV reads:', other_J)
        print('Number of reads passing first filter:', first_filter)
        print('Number of low quality reads:', low_qual_seq)
        print('Number of low quality UMIs:', low_qual_UMI)
        print('Number of germline reads removed:', germ_count)

    logging.info('Number of mispriming corrected reads:' + str(misprimed_corrected))
    logging.info('Number of primer corrected reads:' + str(primer_corrected))
    logging.info('Number of reads processed:' + str(count))
    logging.info('Number of total unclear JV reads:' + str(unclear_J))
    logging.info('Number of total mispriming other JV reads:' + str(other_J))
    logging.info('Number of reads passing first filter:' + str(first_filter))
    logging.info('Number of low quality reads:' + str(low_qual_seq))
    if not short:
        logging.info('Number of low quality UMIs:' + str(low_qual_UMI))
    logging.info('Number of germline reads removed:' + str(germ_count))

    return fq_dict_germ



def demultiplex_assembled(jv_region, fq_dict_pcln, umi_len, anchor_1='GACTCGT', anchor_2='CTGCTCCT', verbose=True):
    '''Demultiplex V reads and extract the umi from assembled reads
    :param V_region: R1 fastq file with V end sequences
    :param fq_dict_pcln: fastqHolder object from preclean
    :param umi_len: length of the umi
    :param anchor_1: anchor 1 GACTCGT
    :param anchor_2: anchor 2 CTGCTCCT
    :param verbose: print stats
    '''

    #Test if preclean run on data

    #anchor 1 and 2 need reverse complement
    rv_anchor_1 = general.reverse_complement(anchor_1)
    rv_anchor_2 = general.reverse_complement(anchor_2)

    fp_jv_region = os.path.abspath(jv_region)

    req_files = fq_dict_pcln.gene_split.keys()
    mis_files = fq_dict_pcln.misprimed.keys()


    with general.file_open(fp_jv_region) as jvr:
        for item in general.fastq_parse(jvr):
            title = item[0]
            seq = item[1]
            thrd = item[2]
            qual = item[3]
            for key in req_files:
                if 'germline' not in key and 'other' not in key: #Don't need to demultiplex germline reads!
                    if title.split(' ')[0][1:] in fq_dict_pcln.gene_split[key]:

                        #perfect match wanted just so indels don't mess up UMI
                        if seq[-umi_len-len(anchor_1):-umi_len] == rv_anchor_1:
                            fq_dict_pcln.demultiplex[key + '_' + anchor_1].add(title.split(' ')[0][1:])
                        elif seq[-umi_len-len(anchor_2):-umi_len] == rv_anchor_2:
                            fq_dict_pcln.demultiplex[key + '_' + anchor_2].add(title.split(' ')[0][1:])
                        else:
                            fq_dict_pcln.demultiplex['unassigned_J'].add(title.split(' ')[0][1:])

            #Demultiplex misprimed reads
            for key in mis_files:
                if title.split(' ')[0][1:] in fq_dict_pcln.misprimed[key]: #or fq_dict_pcln.misprimed[key].keys() does same thing

                    #perfect match wanted just so indels don't mess up UMI
                    if seq[-umi_len-len(anchor_1):-umi_len] == rv_anchor_1:
                        fq_dict_pcln.demultiplex[key + '_' + anchor_1].add(title.split(' ')[0][1:])
                    elif seq[-umi_len-len(anchor_2):-umi_len] == rv_anchor_2:
                        fq_dict_pcln.demultiplex[key + '_' + anchor_2].add(title.split(' ')[0][1:])
                    else:
                        fq_dict_pcln.demultiplex['unassigned_J'].add(title.split(' ')[0][1:])

    if verbose:
        for key in fq_dict_pcln.demultiplex.keys():
            print('Number of reads in', key, ':', len(fq_dict_pcln.demultiplex[key]))

    return fq_dict_pcln




def write_assembled(jv_region, fq_dict_demult, umi_len, prefix=None, out_dir=None,
                    anchor_1='GACTCGT', anchor_2='CTGCTCCT', beyond_anchor=0, q_score=30,
                    verbose=False, j_len=0, in_len=0):
    '''Write out precleaned fastq files
    :param jv_region: (PEAR) assembled fastq file with J and V end sequences
    :param fq_dict: a fastqHolder object from demultiplex
    :param umi_len: length of UMI
    :param q_score: min phred quality for bases in the UMI
    :param prefix: prefix for out files, default is basename of input_J_R2
    :param out_dir: output directory, default is preclean folder created in input directory
    :param anchor_1: anchor 1 GACTCGT
    :param anchor_2: anchor 2 CTGCTCCT
    '''
    # J_to_V = {'unclear_J':'unclear_V', 'germline_J':'germline_V', 'other_J':'other_V', 'J1':'V1', 'J2':'V2', 'J3':'V3', 'J4':'V4'} #generate automatically

    req_files = fq_dict_demult.demultiplex.keys()
    germ_files = fq_dict_demult.gene_split.keys()

    fp_jv_region = os.path.abspath(jv_region)
    dir_nam = os.path.dirname(fp_jv_region)

    if prefix == None:
        prefix_jv = os.path.basename(jv_region).split('.')[0]
    else:
        prefix_jv = prefix

    print('Writing demultiplexed files to:', out_dir)

    #write everything into a single file
    key_list = []
    for key in req_files:
        if 'unassigned' not in key:
            key_list.append(key)
        else: #write unassigned into a seperate file
            fq_dict_demult.write_demultiplex_unassigned(fp_jv_region, key, out_dir + '/' + prefix_jv + '_' + key)

    if j_len > 0:
        fq_dict_pcln.j_end_umi(jv_region, j_len, in_len)
    #Write out everything else
    low_qual, an1_count, an2_count = fq_dict_demult.write_demultiplex_umi_extract_assembled(fp_jv_region, key_list,
                                                           out_path=out_dir + '/' + prefix_jv + '_' + 'all_jv',
                                                           an1=anchor_1, an2=anchor_2, umi_len=umi_len, beyond_anchor=beyond_anchor,
                                                           q_score=q_score, verbose=verbose, j_len=j_len, in_len=in_len)

    logging.info('Number of ' + anchor_1 + ' reads written out:' + str(an1_count))
    logging.info('Number of ' + anchor_2 + ' reads written out:' + str(an2_count))
    if beyond_anchor > 0:
        logging.info('Number of additional low umi quality reads filtered:' + str(low_qual))

    #write out germline and other
    for key in germ_files:
        if 'germline' in key or 'other' in key:
            fq_dict_demult.write_preclean(fp_jv_region, key, out_dir + '/' + prefix_jv + '_' + key)




def write_short(V_region, jv_region, fq_dict_pcln, v_iden_out, umi_len, j_len, in_len,
                prefix=None, out_dir=None, q_score=30, anchor=False, beyond_anchor=0):
    '''Write out short reads
    :param V_region: path to v end fastq
    :param JV_region: path to jv assembled fastq
    :param fq_dict_pcln: fastqHolder object from preclean function
    :param umi_len: how many bases to use from V read as the UMI
    :param j_len: how many bases to use from J read as the UMI
    :param in_len: how many bases into the J read to take bases from
    '''

    genes = set(list(fq_dict_pcln.gene_split.keys()) + list(fq_dict_pcln.misprimed.keys()))

    fp_v_region = os.path.abspath(V_region)
    fp_jv_region = os.path.abspath(jv_region)

    if prefix == None:
        prefix_jv = os.path.basename(jv_region).split('.')[0]
    else:
        prefix_jv = prefix

    #get umi from j end
    if j_len > 0:
        fq_dict_pcln.j_end_umi(jv_region, j_len, in_len)

    #process v end first; V end file will be smaller than J end as it doesn't contain assembled reads
    #need v_end as deduplication removes all_v
    low_qual_UMI, no_anchor, v_short, assembled = fq_dict_pcln.process_v_end(fp_v_region, umi_len,
                                                                  out_dir + '/' + prefix_jv + '_' + 'v_end_all_v',
                                                                  v_iden_out, list(genes),
                                                                  q_score=q_score, anchor=anchor,
                                                                  beyond_anchor=beyond_anchor)

    #write out, get umi from v end
    key_out_reads = defaultdict()
    for key in genes:
        if 'germline' in key or 'other' in key:
            out_reads_n = fq_dict_pcln.write_preclean_short(fp_jv_region, key,
                                                            out_dir + '/' + prefix_jv + '_' + key,
                                                            v_iden_out, merge=False,
                                                            anchor=anchor)
            key_out_reads[key] = out_reads_n
    #else write everything else in the same file
    out_reads = fq_dict_pcln.write_preclean_short(fp_jv_region, list(genes),
                                                  out_dir + '/' + prefix_jv + '_' + 'all_j',
                                                  v_iden_out, merge=True, anchor=anchor)

    other_stats = ', '.join("{!s}: {!r}".format(key,val) for (key,val) in key_out_reads.items())

    print('Number of low quality UMIs: ', low_qual_UMI)
    print('Number of missing anchors: ', no_anchor)
    print('Number of short V reads: ', v_short)
    print('Number of J reads written out', out_reads)
    print('Number of germline and other reads written out - ', other_stats)
    print('Number of assembled read not included in V end ', assembled)
    logging.info('Number of low quality UMIs:' + str(low_qual_UMI))
    logging.info('Number of missing anchors:' + str(no_anchor))
    logging.info('Number of short V reads:' + str(v_short))
    logging.info('Number of J reads written out: ' + str(out_reads))
    logging.info('Number of germline and other reads written out -' + other_stats)
    logging.info('Number of assembled read not included in V end' + str(assembled))


def gemline_removed_qc(V_region, out_dir, spe='mmu', prefix=None, thread_num=8, verbose=False):
    '''Realign and plot J region to determine success of germline removal

    '''
    germ = presets.prs(spe).germline() #['chr12', 113428237, 113430474]

    fp_V_region = os.path.abspath(V_region)
    dir_nam = os.path.dirname(fp_V_region)

    if prefix == None:
        prefix_V = os.path.basename(V_region).split('.')[0]
    else:
        prefix_V = prefix

    if out_dir == None:
        out_dir = dir_nam + '/' + prefix_V + '_preclean'
    print(out_dir)
    #for mm take only first 21 nts (smallest J)
    j_size_ord = sorted(presets.prs(spe).J_seq().values(), key=len)
    shortest_J = len(j_size_ord[0])
    location_j = presets.prs(spe).J_location()

    for name in glob.glob(out_dir + '/*J[0-9]_*'):
        print('Processing', name)

        run_bowtie2 = bowtie2_wrapper.bowtie2()
        run_bowtie2.align_single(fastq=name, nthreads=thread_num, trim5=str(shortest_J), spe=spe, verbose=verbose)

        run_bowtie2.plot(plot_region=germ[0] + ':' + str(germ[1]) + '-' + str(germ[2]), spe=spe)

        run_bowtie2.write(region=germ)

        run_bowtie2.del_tmp()




def v_end_identity(igh_ref, V_region, thread_num, spe='mmu'):
    #Get position of all V genes

    v_genes = defaultdict()
    with open(igh_ref, 'r') as annotation:
        for line in annotation:
            split_line = line.split('\t')
            if 'IGHV' in split_line[3].upper():
                # v_genes[split_line[3].replace('"', '')] = [split_line[4], split_line[5], split_line[6]]
                v_genes[split_line[5]] = split_line[3]


    start_pos = sorted(v_genes.keys())

    #v bins from rss until the next rss (last gene +2000)
    v_bins = defaultdict()
    for i in range(len(start_pos)):
        try:
            v_bins[start_pos[i] + '_' + start_pos[i+1]] = v_genes[start_pos[i]]
        except IndexError:
            v_bins[start_pos[i] + '_' + str(int(start_pos[i])+2000)] = v_genes[start_pos[i]]

    #Align V end with bowtie2 and put idenity of V into qname
    run_bowtie2 = bowtie2_wrapper.bowtie2()
    run_bowtie2.align_single(fastq=V_region, nthreads=thread_num, spe=spe, verbose=True)

    #fetch only within the IG region
    ig = presets.prs(spe).ig()

    sam_file_v = run_bowtie2.pysam_out(region=ig, fetch=True)

    # if write:
    #     run_bowtie2.write(region=igh) #Doesn't work with Hydrogen for some reason!

    run_bowtie2.del_tmp()

    count = 0
    v_iden_dict = defaultdict()
    for read in sam_file_v:
        count += 1
        if read.is_reverse: #should get excluded by region constrainst

            for key in v_bins:
                #start-100<= pos < start+1000
                if int(key.split('_')[0])-100 <= read.pos and int(key.split('_')[0])+1000 > read.pos:
                    v_iden_dict[read.qname] = v_bins[key]

    return v_iden_dict


#
# def v_identity_igblast(V_region, thread_num, spe='mmu'):
#
#     #fastq to fasta
#     fasta = ''
#     for record in general.fastq_to_fasta_iter(V_region):
#         fasta += record
#
#     #make tmp directory with igblast run files
#     tmp_dir = tempfile.mkdtemp()
#     tmp_fmt = os.path.join(tmp_dir, "igblast.fmt7")
#
#     #need to write out fasta for changeo (find alternative which accepts stdin?)
#     with open(tmp_dir + '/igblast.fasta', 'w') as fa_out:
#         fa_out.write(fasta)
#
#     igblast_wrapper.run_igblast(tmp_dir + '/igblast.fasta', tmp_fmt, 10000, spe, thread_num, aux_file=None, additional_flags=['-num_alignments_V', '1'])
#     igblast_wrapper.parse_igblast(tmp_fmt, tmp_dir + '/igblast.fasta', spe)
#     #make v_identity dict key=qname value=idenity
#
#     #need to find the output of changeo
#     tmp_tab = glob.glob(tmp_dir + '/*.tab')
#
#     v_iden_dict = defaultdict()
#     with open(tmp_tab[0], 'r') as tab_file:
#         header = tab_file.readline()
#
#         for line in tab_file:
#             sp_line = line.rstrip('\n').split('\t')
#             SEQUENCE_ID = sp_line[0]
#             V_CALL = sp_line[7]
#             V_SCORE = sp_line[37]
#
#             v_iden_dict[SEQUENCE_ID] = V_CALL + '.' + V_SCORE
#
#     #Delete temporary files
#     shutil.rmtree(tmp_dir)
#
#     return v_iden_dict



def parse_args():
    '''Set up parser
    '''
    parser = argparse.ArgumentParser(description='BabrahamLinkON Preclean')

    parser.add_argument('--version', action='version', version='%(prog)s ' + __version__)

    sub = parser.add_subparsers(dest='action', description='Choose pipeline')

    sp1 = sub.add_parser('umi')
    sp2 = sub.add_parser('short')
    sp3 = sub.add_parser('short_anchor')
    sp4 = sub.add_parser('mispriming_error')
    # sp4 = sub.add_parser('no_anchor')

    # group = parser.add_mutually_exclusive_group(required=True)
    #TODO:remove files not required for deduplication like output from PEAR

    for sp in [sp1, sp2, sp3, sp4]:

        sp.add_argument('--species', dest='species', default='mmu', type=str, help='Which species/locus (mmu, hsa, mmuk), [mmu]')
        sp.add_argument('-t', '--threads', dest='nthreads', default=1, type=int, help='Number of threads to use, [1]')
        # parser.add_argument('-v', '--V_r1', dest='input_V', type=str, metavar='v.fastq', nargs='+', help='Input fastq file(s) with V end sequences')
        # parser.add_argument('-j', '--J_r2', dest='input_J', type=str, metavar='j.fastq', nargs='+', help='Input fastq file(s) with J end sequences')
    for sp in [sp1, sp2, sp3]:

        sp.add_argument('-v', '--V_r1', dest='input_V', type=str, help='Input fastq file with V end sequences')
        sp.add_argument('-j', '--J_r2', dest='input_J', type=str, help='Input fastq file with J end sequences') #required=True

        # parser.add_argument('-jv', '--jv', dest='input_jv', type=str, metavar='jv.fastq', nargs='+', help='Input fastq file(s) from PEAR with J (forward) end and V (reverse) end sequences')

        sp.add_argument('--fast', action='store_true', help='Perform fast inaccurate J identification (additionally need to use --no_mispriming)')
        sp.add_argument('--no_mispriming', action='store_true', help='Don\'t perform mispriming correction ')
        # parser.add_argument('--prefix', dest='prefix', type=str, metavar='N', nargs='+', help='Prefix of the output file (need to provide one for each input)')
        sp.add_argument('--prefix', dest='prefix', type=str, help='Prefix of the output file (need to provide one for each input)')

        sp.add_argument('--out', dest='out_dir', type=str, help='Output direcotry')

        sp.add_argument('--verbose', action='store_true', help='Print detailed progress')
        sp.add_argument('--plot', action='store_true', help='Plot alignments')
        sp.add_argument('-q', '--q_score', dest='q_score', type=int, default=30, help='Minimum Phred quality score for bases in UMI')
        sp.add_argument('-ul', '--umi_len', dest='umi_len', type=int, default=12, help='Length of the UMI')
        sp.add_argument('--j_len', dest='j_len', type=int, default=0, help='Length of J end sequence, 50bp into read (-), \
                                                                            to add to the UMI to further diversify the UMI [0]')
        sp.add_argument('--in_len', dest='in_len', type=int, default=50, help='Length into the J end sequence to go, xbp into read (+umi) [50]')

        sp.add_argument('--an1', dest='an1', default='GACTCGT', type=str, help='Adaptor 1 anchor sequence [GACTCGT]')
        sp.add_argument('--an2', dest='an2', default='CTGCTCCT', type=str, help='Adaptor 2 anchor sequence [CTGCTCCT]')

        sp.add_argument('--keep_germline', action='store_false', help='Skip germline removal step')
        sp.add_argument('--keep_pear', action='store_false', help='Don\'t delete output files from pear')

    for sp in [sp2, sp3]:

        #If analysing short sequences
        # sp.add_argument('--short', action='store_true', help='If using short reads <250bp with no anchor+umi')
        # sp.add_argument('--anchor', action='store_true', help='If using short reads <250bp with anchor+umi')
        sp.add_argument('--ref', dest='ref_path', default=None, type=str, help='Igh reference files path')

    for sp in [sp1, sp3]:
        sp.add_argument('-ba', '--beyond_anchor', dest='beyond_anchor', type=int, default=0, help='Length of V end to take beyond the anchor sequence')

    sp4.add_argument('--input_dir', dest='in_dir', type=str, required=True, help='Input directory (created for/by preclean)')
    # parser.add_argument('--plot_QC', action='store_true', help='QC plot showing if all germline reads were removed (few will be present J-J rearrangements)')

    sp1.set_defaults(short=False, mispriming_error=False)
    sp2.set_defaults(short=True, anchor=False, mispriming_error=False, beyond_anchor=0)
    sp3.set_defaults(short=True, anchor=True, mispriming_error=False)
    sp4.set_defaults(short=False, mispriming_error=True)

    opts = parser.parse_args()

    return opts


#TODO: compress output?
def main():


    #argparse
    opts = parse_args()

    if opts.mispriming_error:

        #variables
        germline_fq = glob.glob(opts.in_dir + '/*germline*')[0]

        verbose=True
        plot=False
        write=False


        all_seqs = mispriming_error_estimate.align_germline(germline_fq, opts.species, plot, write, verbose, opts.nthreads)
        seq_dict = mispriming_error_estimate.extract_j_genes(all_seqs, opts.species)

        before, after = mispriming_error_estimate.mispriming_matrix(seq_dict, opts.species)
        mispriming_error_estimate.make_table(before, after, opts.species, germline_fq + '_mispriming_error.csv')

    else:

        if opts.input_V == None and opts.input_J == None and len(opts.input_V) != len(opts.input_J):
            raise FileNotFoundError('No input files supplied!')

        #Repeat, maybe condense
        if opts.prefix==None:
            prefix = os.path.basename(opts.input_V).split('.')[0]

        fp_v_region = os.path.abspath(opts.input_V)
        dir_nam = os.path.dirname(fp_v_region)

        if prefix == None:
            prefix_jv = os.path.basename(fp_v_region).split('.')[0]
        else:
            prefix_jv = prefix


        #add prefix to out directory
        if opts.out_dir == None:
            out_dir = dir_nam + '/' + prefix_jv + '_preclean'
            print('Output directory:', out_dir)
            try:
                os.mkdir(out_dir)
            except FileExistsError:
                print('Default directory', out_dir, 'already exists. Might overwrite files!')
        else:
            out_dir = opts.out_dir
            print('Output directory:', out_dir)
            try:
                os.mkdir(out_dir)
            except FileExistsError:
                print('Default directory', out_dir, 'already exists. Might overwrite files!')


        logging.basicConfig(level=logging.DEBUG, filename=out_dir + '/' + prefix + '_preclean.log', filemode='a+',
                            format='%(asctime)-15s %(levelname)-8s %(message)s')


        assembled_file, if_assembled = assemble(opts.input_V, opts.input_J, out_dir, threads=opts.nthreads,
                                                prefix=opts.prefix, short=opts.short)

        # os.chdir('/media/chovanec/My_Passport/Old_vdj_seq_data/')
        # assembled_file = assemble('/media/chovanec/My_Passport/Old_vdj_seq_data/lane5_TGACCA_WT_BC_L005_R1_val_1.fq.gz',
        # '/media/chovanec/My_Passport/Old_vdj_seq_data/lane5_TGACCA_WT_BC_L005_R3_val_2.fq.gz', threads=6, short=True)
        #
        # os.getcwd()

        # prefix.assembled.fastq
        # prefix.discarded.fastq
        # prefix.unassembled.forward.fastq
        # prefix.unassembled.reverse.fastq
        # prefix.all_assembled.fastq

        if opts.short:
            germ_assembled = germline(assembled_file + '.all_J.fastq', spe=opts.species, thread_num=opts.nthreads,
                                      plot=opts.plot, verbose=opts.verbose, if_assembled=if_assembled)

            #Merge the two files into one (pairing unassembled reads)
            fq_clean = preclean_assembled(assembled_file + '.all_J.fastq', germ_assembled, q_score=opts.q_score,
                                          umi_len=opts.umi_len, spe=opts.species, verbose=opts.verbose,
                                          no_misprime_correct=opts.no_mispriming, discard_germline=opts.keep_germline,
                                          fast=opts.fast, short=opts.short)

            #get identity of V end using bowtie2 alignment, if ref not specified skip this step.
            if opts.ref_path == None:
                v_iden_dict = dict()
            else:
                v_iden_dict = v_end_identity(opts.ref_path, opts.input_V, thread_num=opts.nthreads, spe=opts.species)
            #get identity of V end using igblast
            # v_iden_dict = v_identity_igblast(opts.input_V[0], thread_num=opts.nthreads, spe=opts.species)

            #Old short reads don't have any anchor (short reads with anchor ignore for now)
            write_short(opts.input_V, assembled_file + '.all_J.fastq', fq_clean,
                        v_iden_dict, umi_len=opts.umi_len,
                        j_len=opts.j_len, in_len=opts.in_len,
                        prefix=opts.prefix, out_dir=out_dir, q_score=opts.q_score,
                        anchor=opts.anchor, beyond_anchor=opts.beyond_anchor)

            if opts.keep_pear:
                print('Removing files')
                os.remove(assembled_file + '.all_J.fastq')
                os.remove(assembled_file + '.assembled.fastq')
                os.remove(assembled_file + '.discarded.fastq')
                os.remove(assembled_file + '.unassembled.forward.fastq')
                os.remove(assembled_file + '.unassembled.reverse.fastq')


        else:
            germ_assembled = germline(assembled_file + '.assembled.fastq', spe=opts.species, thread_num=opts.nthreads,
                                      plot=opts.plot, verbose=opts.verbose, if_assembled=if_assembled)

            fq_clean = preclean_assembled(assembled_file + '.assembled.fastq', germ_assembled, q_score=opts.q_score,
                                          umi_len=opts.umi_len, spe=opts.species, verbose=opts.verbose,
                                          no_misprime_correct=opts.no_mispriming, discard_germline=opts.keep_germline,
                                          fast=opts.fast)

            fq_demultiplex = demultiplex_assembled(assembled_file + '.assembled.fastq', fq_clean, umi_len=opts.umi_len,
                                                   anchor_1=opts.an1, anchor_2=opts.an2, verbose=opts.verbose)

            write_assembled(assembled_file + '.assembled.fastq', fq_demultiplex, umi_len=opts.umi_len, prefix=opts.prefix,
                            out_dir=out_dir, beyond_anchor=opts.beyond_anchor, q_score=opts.q_score, verbose=opts.verbose,
                            j_len=opts.j_len, in_len=opts.in_len)

            if opts.keep_pear:
                print('Removing files')
                os.remove(assembled_file + '.assembled.fastq')
                os.remove(assembled_file + '.discarded.fastq')
                os.remove(assembled_file + '.unassembled.forward.fastq')
                os.remove(assembled_file + '.unassembled.reverse.fastq')



if __name__ == "__main__":
    main()
