from collections import defaultdict, OrderedDict
import gzip
import re
import numpy as np
import pandas as pd
import os
import Levenshtein

import subprocess
import pysam
import tempfile
import shutil
import shlex
from babrahamlinkon import presets
from itertools import groupby
import tempfile


try:
    from skbio.alignment import StripedSmithWaterman #put into setup.py install
except:
    print('ERROR: This script requires scikit-bio: pip install scikit-bio')
    sys.exit()



def file_open(filename):
    """
    Open as normal or as gzip
    Faster using zcat?
    """
    #does file exist?
    f = open(filename,'rb')
    if (f.read(2) == b'\x1f\x8b'): #compressed alsways start with these two bytes
        f.seek(0) #return to start of file
        return gzip.GzipFile(fileobj=f, mode='rb')
    else:
        f.seek(0)
        return f



def fastq_parse(fp):
    """
    Parse fastq file.
    """
    linecount = 0
    name, seq, thrd, qual = [None] * 4
    for line in fp:

        linecount += 1
        if linecount % 4 == 1:
            try:
                name = line.decode('UTF-8').rstrip()
            except AttributeError:
                name = line.rstrip()
            assert name.startswith('@'),\
                   "ERROR: The 1st line in fastq element does not start with '@'.\n\
                   Please check FastQ file near line number %s" % (linecount)
        elif linecount % 4 == 2:
            try:
                seq = line.decode('UTF-8').rstrip()
            except AttributeError:
                seq = line.rstrip()
        elif linecount % 4 == 3:
            try:
                thrd = line.decode('UTF-8').rstrip()
            except AttributeError:
                thrd = line.rstrip()
            assert thrd.startswith('+'),\
                   "ERROR: The 3st line in fastq element does not start with '+'.\n\
                   Please check FastQ file near line number %s" % (linecount)
        elif linecount % 4 == 0:
            try:
                qual = line.decode('UTF-8').rstrip()
            except AttributeError:
                qual = line.rstrip()
            assert len(seq) == len(qual),\
                    "ERROR: The length of Sequence and Quality aren't equal.\n\
                    Please check FastQ file near line number %s" % (linecount)

            yield name, seq, thrd, qual,
            name, seq, thrd, qual = [None] * 4




def fasta_parse(fp):
    '''
    Parse fasta file
    '''
    linecount = 0
    name, seq = [None] * 2
    for line in fp:

        linecount += 1
        if linecount % 2 == 1:

            try:
                name = line.decode('UTF-8').rstrip()
            except AttributeError:
                name = line.rstrip()

            assert name.startswith('>'),\
                   "ERROR: The 1st line in fasta element does not start with '>'.\n\
                   Please check FastQ file near line number %s" % (linecount)
        elif linecount % 2 == 0:
            try:
                seq = line.decode('UTF-8').rstrip()
            except AttributeError:
                seq = line.rstrip()

            yield name, seq
            name, seq = [None] * 2

#TODO: better fasta parser in case sequence is not all on the same line
def fasta_iter(fasta_name):
    '''
    given a fasta file. yield tuples of header, sequence
    https://drj11.wordpress.com/2010/02/22/python-getting-fasta-with-itertools-groupby/
    https://www.biostars.org/p/710/
    '''
    with open(fasta_name) as fasta_in:
        # ditch the boolean (x[0]) and just keep the header or sequence since we know they alternate.
        faiter = (x[1] for x in groupby(fasta_in, lambda line: line[0] == '>'))
        for header in faiter:
            # drop the ">"
            header = header.__next__()[1:].strip() #first group
            # join all sequence lines to one.
            seq = "".join(s.strip() for s in faiter.__next__()) #second group
            yield header, seq



def trim_fastq(J_region, out_file, trim_size=100):
    '''Trim all reads to a set length
    :param J_region: J fastq
    :param out_file: path of output file
    :param trim_size: what size to trim all the reads down to (default: 100)
    :return: fastq with trimmed sequences
    '''
    with open(out_file, 'w') as out, file_open(J_region) as Jr:
            lines = Jr.read().splitlines()
            for item in fastq_parse(lines):
                title = item[0]
                seq = item[1]
                thrd = item[2]
                qual = item[3]

                if len(seq) <= trim_size:
                    out.write(title + '\n' + seq + '\n' + thrd + '\n' + qual + '\n')
                else:
                    seq = seq[:trim_size]
                    qual = qual[:trim_size]
                    out.write(title + '\n' + seq + '\n' + thrd + '\n' + qual + '\n')



def reverse_complement(seq):
    '''Return reverse complement of DNA sequence
    '''
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N':'N'}
    return ''.join([complement[base] for base in seq.upper()[::-1]])


def reverse_string(my_str):
    '''
    Reverse the qual of a sequence
    '''
    return ''.join(let for let in my_str.upper()[::-1])


def check_qual(umi_qual, q_score=30):
    for val in umi_qual:
        phred = ord(val)-33
        assert phred <= 41 and phred >= 0, 'Phred score out side of range 0-41'
        if phred < q_score:
            return True
        else:
            return False


################################################################################
###### Data holder ######
################################################################################


class fastqHolder:
    '''
    Store fastq split by gene (J gene)
    '''
    def __init__(self):
        self.gene_split = defaultdict(set) #split into seperate J's and germline
        self.preclean_v = defaultdict(set) #read has to be in both gene_split and preclean_v
        self.demultiplex = defaultdict(set)
        # self.germline = defaultdict(set)
        self.misprimed = defaultdict(lambda: dict())
        self.original_id = defaultdict(set)
        self.fastq_split = defaultdict(list)

    # def add_to_gene(self, gene, qname):
    #     '''Add to gene_split dictionary
    #     :param gene: name of key
    #     :param qname: fastq qname to be stored under key
    #     '''
    #     self.gene_split[gene].add(qname)

    def add_to_misprimed(self, orig_gene, gene, cor_seq, qname):
        '''Add to gene_split dictionary
        :param gene: name of key
        :param orig_gene: initial identity of J before mispriming correction
        :param qname: fastq qname to be stored under key
        :param cor_seq: corrected misprimed sequence to be stored under key
        '''
        self.misprimed[gene][qname] = cor_seq
        self.original_id[qname] = orig_gene

    def add_to_fastq(self, gene, qname, seq, third, qual):
        self.fastq_split[gene].append([qname, seq, third, qual])

    # def write_to_file(self, gene, path):
    #     '''Write fastq_split to file
    #     '''
    #     with open(path, 'w') as out_file:
    #         for qname, seq, thrd, qual in self.fastq_split[gene]:
    #             out_file.write(qname + '\n' + seq + '\n' + thrd + '\n' + qual + '\n')

    def write_demultiplex_unassigned(self, fastq_path, gene, out_path):
        '''Write to a fastq (for unasigned reads)
        :param fastq_path: fastq to subset
        :param gene: what to subset by (J gene in split_gene)
        :param out_path: write path
        :param j_end: are these J reads? if yes will correct misprimed sequences
        '''

        with file_open(fastq_path) as fq, \
        open(out_path, 'w') as out_file:
            for qname, seq, thrd, qual in fastq_parse(fq):
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

        with file_open(fastq_path) as fq, \
        open(out_path, 'w') as out_file:
            for qname, seq, thrd, qual in fastq_parse(fq):
                if qname.split(' ')[0][1:] in self.gene_split[gene]:
                    out_file.write(qname.split(' ')[0] + '_' + gene + ' ' +
                                   ''.join(qname.split(' ')[1:]) + '\n' + seq + '\n' + thrd + '\n' + qual + '\n')


    def write_preclean_short(self, fastq_path_v, fastq_path_jv, genes, out_path, umi_len, v_iden_dict, merge, q_score):
        '''Write fastq file
        :param fastq: fastq to subset
        :param gene: what to subset by (J gene in split_gene)
        :param path: write path
        '''

        umi_dict = defaultdict()
        rec_skip = set()
        low_qual_UMI = 0
        if merge:
            with file_open(fastq_path_jv) as fq, open(out_path, 'w') as out_file, file_open(fastq_path_v) as v_fq:
                for qname, seq, thrd, qual in fastq_parse(v_fq):
                    #check qual of umi_len
                    umi_qual = qual[-umi_len:]
                    if check_qual(umi_qual, q_score):
                        #poor quality skip record
                        low_qual_UMI += 1
                        rec_skip.add(qname.split(' ')[0][1:])
                    else:
                        umi_dict[qname.split(' ')[0]] = seq[-umi_len:]
                for qname, seq, thrd, qual in fastq_parse(fq):
                    if qname.split(' ')[0][1:] in rec_skip:
                        continue
                    try:
                        v_iden_out = v_iden_dict[qname.split(' ')[0][1:]]
                    except KeyError:
                        v_iden_out = ''
                    for gene in genes:
                        if 'germline' not in gene and 'other' not in gene:
                            if qname.split(' ')[0][1:] in self.gene_split[gene]:
                                out_file.write(qname.split(' ')[0] + '_' + gene + '_' + v_iden_out.upper() + '_' + umi_dict[qname.split(' ')[0]] + ' ' +
                                               ''.join(qname.split(' ')[1:]) + '\n' + seq + '\n' + thrd + '\n' + qual + '\n')

        else:
            if isinstance(genes, list):
                raise ValueError('Multiple genes supplied without wanting to merge')

            with file_open(fastq_path_jv) as fq, \
            open(out_path, 'w') as out_file, \
            file_open(fastq_path_v) as v_fq:
                for qname, seq, thrd, qual in fastq_parse(v_fq):
                    #check qual of umi_len
                    umi_qual = qual[-umi_len:]
                    if check_qual(umi_qual, q_score):
                        #poor quality skip record
                        low_qual_UMI += 1
                        rec_skip.add(qname.split(' ')[0][1:])
                    else:
                        umi_dict[qname.split(' ')[0]] = seq[-umi_len:]
                for qname, seq, thrd, qual in fastq_parse(fq):
                    if qname.split(' ')[0][1:] in rec_skip:
                        continue
                    if qname.split(' ')[0][1:] in self.gene_split[genes]:
                        out_file.write(qname.split(' ')[0] + '_' + genes + '_' + umi_dict[qname.split(' ')[0]] + ' ' +
                                       ''.join(qname.split(' ')[1:]) + '\n' + seq + '\n' + thrd + '\n' + qual + '\n')
        return low_qual_UMI


    def write_demultiplex_umi_extract_assembled(self, jv_fastq, gene_list, out_path, an1, an2, umi_len):
        '''Write everything into single fastq file and extract umi
        :param fastq: fastq to subset
        :param gene_list: what to subset by (J gene in split_gene)
        :param V_out_path: V out write path
        :param J_out_path: J out write path
        :param an_1: barcode 1 GACTCGT
        :param an_2: barcode 2 CTGCTCCT
        '''
        # umi_dict = defaultdict()



        low_qual_UMI = 0
        an1_count = 0
        an2_count = 0
        gene_count = defaultdict(int)

        with file_open(jv_fastq) as jv_fq, open(out_path + '_' + an1, 'w') as an1_out_file, open(out_path + '_' + an2, 'w') as an2_out_file:
            for qname, seq, thrd, qual in fastq_parse(jv_fq):
                for gene in gene_list: #J1 J2 J3 J4...
                    if qname.split(' ')[0][1:] in self.demultiplex[gene]:

                        umi = seq[-umi_len:] #last 6 bases NNNNNN (it is in reverse complement configuration compared to original V read)

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
                        #Remove UMI and barcode from seq (keeping it in read qname only)
                        if an1 in gene:
                            an1_out_file.write(qname.split(' ')[0] + '_' + gene + '_' + umi + ' ' +
                            ''.join(qname.split(' ')[1:]) + '\n' + seq[:-(7+umi_len)] + '\n' + thrd + '\n' + qual[:-(7+umi_len)] + '\n')
                            an1_count += 1
                            gene_count[gene] += 1
                        else:
                            an2_out_file.write(qname.split(' ')[0] + '_' + gene + '_' + umi + ' ' +
                            ''.join(qname.split(' ')[1:]) + '\n' + seq[:-(8+umi_len)] + '\n' + thrd + '\n' + qual[:-(8+umi_len)] + '\n')
                            an2_count += 1
                            gene_count[gene] += 1

        # print('Low quality UMIs:', low_qual_UMI)
        print(an1, 'reads written out:', an1_count)
        print(an2, 'reads written out:', an2_count)

        for j in gene_count.keys():
            print('Number of', j, 'reads written out:', gene_count[j])











################################################################################
###### J alignment for mispriming ######
################################################################################

class SSW_align:
    '''Align reads using StripedSmithWaterman from scikit-bio
    http://scikit-bio.org/docs/latest/generated/skbio.alignment.StripedSmithWaterman.html#skbio.alignment.StripedSmithWaterman
    '''
    #align_coord = {'J1':0,'J2':27, 'J3':55, 'J4':85}
    def __init__(self):
        self.alignments = []
        # self.print_full = []
        self.ref_overlap = defaultdict(list) #created by reference function
        self.align_coord = defaultdict(int) #created by reference function

    def align(self, target, qry, score=0, misprim_cor=True, mismatch_allow=1, spe='mmu', quick_align=False):
        '''Align against reference and return which J sequence aligns to
        :param target: reference
        :param qry: sequence to align against reference
        :param score: score filter, min score to pass as an alignment
        :param misprim_cor: perform mispriming correction
        :param mismatch_allow: number of mismatches to allow in sequence beyond J primer
        :param spe: which organism
        :param quick_align: skip SSW and mispriming and return initial identification using only levenshtein distance on first 5bps
        :return: read identity (store alignment in alignments)
        '''
        #Query needs to be a single string
        assert isinstance(qry, str), 'Align query should not be a list'

        #use 5 bp to do initial identification
        ini_len = 5
        restart = True
        while restart:
            restart = False
            J_identities = defaultdict(int)
            for J in presets.prs(spe).J_seq().keys():
                dist = Levenshtein.distance(qry[:ini_len], presets.prs(spe).J_seq()[J][:ini_len])
                J_identities[J] = dist #if incorrect will catch later
            #chose J with lowest score
            if J_identities and ini_len <= 10: #if dict not empty and iteration is less then 10bp
                min_dist = min(J_identities.values())
                min_dist_key = [k for k, v in J_identities.items() if v == min_dist]
                if len(min_dist_key) > 1: #if two keys equal score increase size and repeat until only one hit obtained
                    restart = True
                    ini_len += 1
                else:
                    initial_identity = min_dist_key[0]
            else:
                #if all the same keys are haplotypes, take one of them
                min_dist = min(J_identities.values())
                min_dist_key = [k for k, v in J_identities.items() if v == min_dist]
                keys = set()
                for key in min_dist_key:
                    keys.add(key.split('.')[0])
                if len(keys) > 1:
                    return 'other' #no hits
                else:
                    initial_identity = min_dist_key[0]

        if quick_align: #in case only interested in V genes skipp SSW and mispriming
            return [initial_identity]

        #while loop in case initial identity is wrong
        redo = True
        add = 0
        not_added = True #prevent infinit loop
        not_switched = True
        while redo:
            redo = False

            #Use initial identity to align x bp
            align_size = len(presets.prs(spe).J_for_mispriming()[initial_identity]) + add

            #allows 1 error per 5 bases
            qry_err = align_size-round(align_size/5)

            query_out = StripedSmithWaterman(qry[:align_size], match_score=1, mismatch_score=0,
                                                 gap_open_penalty=1, gap_extend_penalty=1,
                                                 score_filter=qry_err)
            #score_filter filter out low score alignments (lowers computational time)

            alignment = query_out(target)

            if alignment['target_begin'] != -1: #target_begin and cigar not returned if score low


                #Which J does alignment correspond to
                #Overlap
                try:
                    read_identity = self.ref_overlap[alignment['target_begin']]
                except KeyError:
                    print('Need to run reference function first!')


                #if read doesn't match initial identification
                if read_identity != initial_identity:
                    if not_switched:
                        redo = True
                        initial_identity = read_identity
                        not_switched = False
                        continue
                    else: #initial and after alignment identity didn't match for second time
                        return 'other'
                else:

                    #If mispriming include original read identity in header; output: (before, after, correct_seq)
                    if misprim_cor:
                        self.alignments.append(alignment)


                        s = alignment['cigar']
                        matches = re.findall(r'(\d+)([A-Z]{1})', s)
                        cigar_dict = ([{'type':m[1], 'length':int(m[0])} for m in matches])

                        start = alignment['target_begin']
                        end = alignment['target_end_optimal']+1
                        q_start = alignment['query_begin']
                        q_end = alignment['query_end']+1

                        print_seq = '-'*start
                        poscounter = q_start

                        for cig in cigar_dict:
                            if cig['type'] == 'M':
                                #match just print seq
                                c_start = poscounter
                                poscounter += cig['length']
                                print_seq = print_seq + alignment['query_sequence'][c_start:poscounter]
                            elif cig['type'] == 'D':
                                print_seq = print_seq + '-'*cig['length']
                            elif cig['type'] == 'I' and q_end-poscounter < 5: #insertions may cause aligning sequence to be too short!
                                #don't skip insertions near the end
                                c_start = poscounter
                                poscounter += cig['length']  #skip over insertion
                                print_seq = print_seq + alignment['query_sequence'][c_start:poscounter]
                            elif cig['type'] == 'I':
                                poscounter += cig['length']  #skip over insertion

                        #get sequence beyond reference/alignment
                        #know where reference ends so can go 5bp back even if the end did not align / add the full query sequcence back
                        print_seq_full = print_seq + qry[poscounter:] #+ '-'*(len(alignment['target_sequence'])-end-len(qry[poscounter:]))
                        # CCCTGTGCCCCAGACATCNNNNNNAGTGGTGCCTTGGCCCCAGTAGTCAAANNNNNNACCAGAGTCCCTTGGCCCCAGTAAGCAAANNNNNNTGAGGTTCCTTGACCCCAGTAGTCCAT
                        # CCCCTGTGCCCCAGACATCGAAGTACCAGTA raw read
                        # CCTTGTGCCCCAGACATCGAAGTACCAGTA-----------------------------------------------------------------------------------------
                        # print_seq = print_seq + '-'*(len(alignment['target_sequence'])-end)
                        # CCTTGTGCCCCAGACAT------------------------------------------------------------------------------------------------------
                        # self.print_full = print_seq_full

                        ref_end = self.align_coord[read_identity]['end'] #End of J in the reference


                        differences = defaultdict(int)
                        primer_end = defaultdict(int)
                        for key in presets.prs(spe).mispriming().keys():
                            seq_chk_lst = []
                            try: #if in offset dict do ... else do ...
                                #Special case for J1 (rare)

                                offst = presets.prs(spe).offset()[read_identity][key]
                                for off in offst:
                                    seq_chk_lst.append(print_seq_full[ref_end-off:ref_end-off+5])
                                    primer_end[key] = ref_end-off
                                # read_offset = True
                            except KeyError: #if don't need to offset using offset dict then just use last 5bp
                                #reference end - 5 : reference end - distance from actual end of reference + last letters that have been cut off
                                seq_chk_lst.append(print_seq_full[ref_end-5:ref_end])
                                primer_end[key] = ref_end-5
                                # read_offset = False
                                # pass

                            diffs = []
                            for seq_chk in seq_chk_lst:
                                diffs.append(Levenshtein.distance(seq_chk, presets.prs(spe).mispriming()[key][:len(seq_chk)]))
                            #Pick smallest diff
                            differences[key] = min(diffs)

                            # if len(seq_chk) >= 4:
                            #     # print('seq_chk', seq_chk)
                            #     differences[key] = (Levenshtein.distance(seq_chk, presets.prs(spe).mispriming()[key][:len(seq_chk)]))
                            # #REVIEW: Can skip this now?
                            # elif not_added: #if sequence beyond primer is shorter than 4bp redo adding x bp to initial alignment
                            #     add = 5-len(seq_chk)
                            #     redo=True
                            #     not_added = False
                            # else: #shouldn't return anything here
                            #     return 'unclear-' + read_identity + '-bpa'

                        if redo: #redo with a longer read if too many insertions present
                            continue

                        min_val = min(differences.values())

                        if min_val <= 0: #allow only 0 mismatches by default in beyond J seq
                            #retrieve key of min val (J)
                            min_val_key = [k for k, v in differences.items() if v == min_val]

                            if len(min_val_key) > 1: #if there are multiple values with same score (shouldn't happen)
                                return 'unclear-' + read_identity + '-mh'
                            else:
                                replace_seq = presets.prs(spe).replace()[min_val_key[0]]

                                corrected_seq = replace_seq + print_seq_full[primer_end[min_val_key[0]]:]
                                # if '-' in corrected_seq:
                                #     print(min_val_key, differences)
                                #     print(corrected_seq)
                                #     print(print_seq_full)
                                #     corrected_seq.replace('-', '')
                                assert '-' not in corrected_seq, 'Unknown space in corrected sequence'

                                #Will correct qual when writing out fastq (fastqHolder), just trim start or add 'I'
                                # print(read_identity, min_val_key[0], corrected_seq)
                                return [read_identity, min_val_key[0], corrected_seq] #identity of J and corrected sequence

                        else: #if more than allowed mismatches (chewing back?)
                            return 'unclear-' + read_identity + '-chw'
                            # #return 15bp beyond primer for lowest score
                            # min_val_key = [k for k, v in differences.items() if v == min_val]
                            # if len(min_val_key) > 1: #if there are multiple values with same score
                            #     return 'unclear-' + read_identity + '-mh'
                            # else:
                            #     return 'unclear-' + read_identity + '-' + print_seq_full[primer_end[min_val_key[0]]:primer_end[min_val_key[0]]+15]

                    else: #if not misprime correcting return identity of primer seq
                        return [read_identity]
            else: #doesn't match primer seq within acceptable paramenters
                return 'other'


    def print_align(self, num=100, print_out=True):
        '''Print pretty alignemnt
        :param num: how many alignments to print
        '''
        #Target sequence should be the same for all (reference)
        seqs = []
        if print_out:
            print(self.alignments[0]['target_sequence'])

        seqs.append(self.alignments[0]['target_sequence'])
        count = 0

        for item in self.alignments:
            if count <= num: #print only num alignments
                s = item['cigar']
                matches = re.findall(r'(\d+)([A-Z]{1})', s)
                cigar_dict = ([{'type':m[1], 'length':int(m[0])} for m in matches])

                start = item['target_begin']
                end = item['target_end_optimal']+1
                q_start = item['query_begin']
                q_end = item['query_end']+1

                print_seq = '-'*start
                poscounter = q_start

                for cig in cigar_dict:
                    if cig['type'] == 'M':
                        #match just print seq
                        c_start = poscounter
                        poscounter += cig['length']
                        print_seq = print_seq + item['query_sequence'][c_start:poscounter]
                    elif cig['type'] == 'D':
                        #deletion print '-'
                        print_seq = print_seq + '-'*cig['length']
                    elif cig['type'] == 'I':
#                         print_seq = print_seq + '^'*cig['length']
                        poscounter += cig['length']  #skip over insertion

                print_seq = print_seq + '-'*(len(item['target_sequence'])-end)
                seqs.append(print_seq)
                if print_out:
                    print(print_seq)
                count += 1
            else:
                break

        return(seqs)



    def reference(self, spe='mmu', quick_align=False):
        '''Create reference from J sequences supplied (will seperate each J by 6 Ns)
        :param spe: which organism
        '''

        ref = ''

        if quick_align:
            return ref

        for key in sorted(presets.prs(spe).J_for_mispriming().keys()):
            start = len(ref)
            ref = ref + presets.prs(spe).J_for_mispriming()[key] + 'N'*6
            end = len(ref)-6
            if start == 0:
                self.align_coord[key] = {'start':start, 'end':end}
            else:
                self.align_coord[key] = {'start':start, 'end':end} #0 based coordinates

            for val in range(start, end):
                self.ref_overlap[val] = key
        ref = ref[:-6]

        return ref




################################################################################
###### Bowtie2 ######
################################################################################

class bowtie2:
    """Call bowtie2-align on single read data

    """

    def __init__(self):
        # # check bowtie2 and samtools is accessible
        # try:
        #     subprocess.check_output(['bowtie2', '-h'])
        #     subprocess.check_output(['samtools', '--help'])
        # except OSError:
        #     raise RuntimeError('bowtie2/samtools not found; put directory in $PATH\n')

        # Create tmp files
        self.tmp_dir = tempfile.mkdtemp()
        self.tmp_prefix = os.path.join(self.tmp_dir, "bowtiepipe")

        self.dir_nam = ''
        self.prefix = ''
        self.sam_algn = ''


    def align_single(self, fastq, nthreads, trim5='0', score='L,0,-1', flags=('--very-sensitive', '--no-unal', '--quiet'),
                     spe='mmu', samtools_mapq=0, verbose=True, out_dir=None):
        '''Use bowtie2 to perform single end alignment
        :param fastq: full path to fastq to aligned
        :param nthreads: number of cores to Use
        :param trim: trim the barcode from 5' end
        :param score: score min to be used  by bowtie2
        :param flags: additional bowtie2 flags, remove reads that failed to align, prints only serious errors (won't get alignment rate :( )
        :param spe: species
        :param samtools_mapq: set MAPQ score on which to filter
        :param verbose: Output some stats
        :param out_dir: custom output directory
        '''

        ref_index = presets.prs(spe).bowtie_index()

        #Set up variables for write functions
        if out_dir == None:
            self.dir_nam = os.path.dirname(fastq)
        else:
            self.dir_nam = out_dir

        self.prefix = os.path.basename(fastq).split('.')[0]


        # stream output from bowtie2
        bowtie2_args = ['bowtie2', \
        '-x', ref_index, \
        '-U', fastq, \
        '-p', str(nthreads), \
        '--trim5', str(trim5), \
        '--score-min', str(score)] + \
        list(flags)


        samtools_bam = shlex.split('samtools view -bS -q ' + str(samtools_mapq) + ' -') #default = 0

        bowtie2_out = subprocess.Popen(bowtie2_args, stdout=subprocess.PIPE)
        samtools_out = subprocess.Popen(samtools_bam, stdin=bowtie2_out.stdout,
                                        stdout=open(self.tmp_prefix + '.bam', 'wb')).communicate() #avoid communicate

        #If bowtie returns error, throw exception
        # out, err = bowtie2_out.communicate()
        # if b'(ERR)' in err:
        #     raise Exception(err.decode("utf-8"))

        bowtie2_out.stdout.close()

        #unique temp name
        tf = tempfile.NamedTemporaryFile()

        #Sort and index bam file
        samtools_sort = ['samtools', 'sort', '-O', 'bam', '-o', self.tmp_prefix + '_sorted.bam',
                         '-T', os.path.basename(tf.name), self.tmp_prefix + '.bam']
        samtools_index = ['samtools', 'index', self.tmp_prefix + '_sorted.bam']

        subprocess.call(samtools_sort)
        subprocess.call(samtools_index)


        if verbose:
            # print(err.decode('utf-8'))
            samtools_count = shlex.split('samtools view -F 0x904 -c ' + self.tmp_prefix + '_sorted.bam')
            read_count = subprocess.check_output(samtools_count)
            print('Number of reads that aligned and passed filter:', read_count.decode("utf-8"))



    def plot(self, plot_region='', spe='mmu'):
        '''Plot bam file using rscript
        '''
        assert self.dir_nam, 'Run align_single first!'

        print('Saving coverage plot to', self.dir_nam)

        #location of current scripts
        plot_igh_bef = ['Rscript', os.path.dirname(os.path.realpath(__file__)) +
        '/plot_igh.R','-o', self.dir_nam + '/' + self.prefix + '.pdf',
        '-n', 'Coverage','--genome', presets.prs(spe).genome(), '-r', plot_region,
        '-b', self.tmp_prefix + '_sorted.bam']

        #TODO: capture error here
        subprocess.call(plot_igh_bef)


    def pysam_out(self, region=None, algn=False, fetch=False):
        '''Pysam objects out
        :param algn: output AlignmentFile
        :param fetch: output fetch object
        :return: fetch pysam object or AlignmentFile pysam object
        '''
        # Read from tmp file
        self.sam_algn = pysam.AlignmentFile(self.tmp_prefix + '_sorted.bam', "rb")


        # sam_pileup = sam_algn.pileup(region[0], region[1], region[2])

        if algn:
            return self.sam_algn
        elif fetch:
            if region == None:
                sam_fetch = self.sam_algn.fetch()
            else:
                try: #chr12 or 12
                    sam_fetch = self.sam_algn.fetch(region[0], region[1], region[2])
                except ValueError:
                    sam_fetch = self.sam_algn.fetch(region[0].split('chr')[1], region[1], region[2])
            return sam_fetch


    def write(self, region=None):
        '''Write out tmp bam files
        :param region: bam region to output
        '''
        assert self.dir_nam, 'Run align_single first!'
        assert self.sam_algn, 'Run pysam_out first!'

        name_out = self.dir_nam + '/' + self.prefix + '.bam' #output same location as input
        print('Saving bam', self.prefix + '.bam', 'to', self.dir_nam)

        with pysam.AlignmentFile(name_out, 'wb', template=self.sam_algn) as write_reads:

            if region == None:
                for read in self.sam_algn.fetch():
                    write_reads.write(read)
            else:
                try:
                    sam_fetch = self.sam_algn.fetch(region[0], region[1], region[2])
                except ValueError:
                    sam_fetch = self.sam_algn.fetch(region[0].split('chr')[1], region[1], region[2])

                for read in sam_fetch:
                    write_reads.write(read)

    def del_tmp(self):
        '''Clean up the named pipe and associated temp directory
        '''
        assert self.tmp_dir, 'Run align_single first!'

        shutil.rmtree(self.tmp_dir)
