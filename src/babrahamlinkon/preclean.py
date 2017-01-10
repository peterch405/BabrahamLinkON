#!/usr/bin/env python3

import sys
import pysam
import subprocess
import os
import re
import argparse
import Levenshtein
import argparse
import glob
from babrahamlinkon import general, presets
from collections import defaultdict
import logging



def assemble(V_region, J_region, threads=1, prefix=None, short=False):
    '''Assemble reads using pear and join unassembled reads
    http://stackoverflow.com/questions/21953835/run-subprocess-and-print-output-to-logging
    '''

    if prefix==None:
        prefix = os.path.basename(V_region).split('.')[0]

    fp_V_region = os.path.abspath(V_region)
    fp_J_region = os.path.abspath(J_region)

    dir_path =  os.path.dirname(fp_V_region)

    # prefix = ''.join(re.split('(L00[0-9])', os.path.basename(V_region))[:-1])

    pear_assemble = ['pear', '-f', fp_J_region, '-r', fp_V_region, '-o', prefix, '-j', str(threads)]

    logging.info('Subprocess: "' + ' '.join(pear_assemble) + '"')

    #TODO: do similar logging for bowtie?
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

    out_file_name = dir_path + '/' + prefix #+ '.all_assembled.fastq'

    if short:
        #Create merge file (won't be used for germline)
        reverse_dict = defaultdict()

        with open(dir_path + '/' + prefix + '.unassembled.forward.fastq', 'r') as forward, \
        open(dir_path + '/' + prefix + '.unassembled.reverse.fastq', 'r') as reverse, \
        open(dir_path + '/' + prefix + '.assembled.fastq', 'r') as assembled, \
        open(dir_path + '/' + prefix + '.all_assembled.fastq', 'w') as out_file:

            #first write out assembled reads
            out_file.write(assembled.read())

            #write out fused unasembled reads
            for item in general.fastq_parse(reverse):
                reverse_dict[item[0].split(' ')[0]] = item

            for item in general.fastq_parse(forward):
                qname = item[0]
                seq = item[1]
                thrd = item[2]
                qual = item[3]

                #PEAR already reverse complements reads unless --keep-originl
                # try:
                #     rv_seq = general.reverse_complement(reverse_dict[qname.split(' ')[0]][1])
                #     rv_qual = general.reverse_string(reverse_dict[qname.split(' ')[0]][3])
                # except KeyError:
                #     print("Fastq pairs don't match:", qname)

                # out_file.write(qname + '\n' + seq + rev_seq + '\n' + thrd + '\n' + qual + rv_qual + '\n')
                v_len = len(reverse_dict[qname.split(' ')[0]][1])
                j_len = len(seq)
                qname_sp = qname.split(' ')
                out_file.write(qname_sp[0] + '_V'+ str(v_len) + '-J' + str(j_len) + ' ' + qname_sp[1] + \
                               '\n' + seq + reverse_dict[qname.split(' ')[0]][1] + \
                               '\n' + thrd + '\n' + qual + reverse_dict[qname.split(' ')[0]][3] + '\n')


    return out_file_name





#Don't need to do, will perform in deduplicate (useful for troubleshooting)
# def aligned_to_igh(V_region, cores_num, spe='mmu', write=False, plot=False, verbose=True):
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
    # run_bowtie2 = general.bowtie2()
    # run_bowtie2.align_single(fastq=fp_V_region, nthreads=cores_num, spe=spe, verbose=verbose)
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
    # nthreads=cores_num,
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

#Don't need to remove germline other than to reduce file size (speed up) (MiXCR removes germline as it isn't going to have a V)
#Almost half the data is germline, required for germline mispriming.
def germline(J_region, cores_num, spe='mmu', plot=False, write=False, verbose=True):
    '''Extract reads aligning to J genes i.e. germline (used to correct mispriming)
    :param J_region: R2 fastq file with J end sequences
    :param spe: which species
    :param plot: generate pileup plot; default False
    :return: set of read qnames that have aligned to J genes
    '''

    germ = presets.prs(spe).germline() #['chr12', 113428237, 113430474]
    ref_index = presets.prs(spe).bowtie_index()

    fq_data = general.fastqHolder()

    #for mm take only first 21 nts (smallest J)
    j_size_ord = sorted(presets.prs(spe).J_seq().values(), key=len)
    shortest_J = len(j_size_ord[0])
    # location_j = presets.prs(spe).J_location()

    fp_J_region = os.path.abspath(J_region)

    #Don't need to align V reads, if germline present will get discarded in deduplicaiton DJ V seperation
    run_bowtie2 = general.bowtie2()
    run_bowtie2.align_single(fastq=fp_J_region, nthreads=cores_num, trim5=str(shortest_J), spe=spe, verbose=verbose)

    if plot:
        run_bowtie2.plot(plot_region=germ[0] + ':' + str(germ[1]) + '-' + str(germ[2]), spe=spe)

    sam_file_j = run_bowtie2.pysam_out(region=germ, fetch=True)

    if write:
        run_bowtie2.write(region=germ)

    run_bowtie2.del_tmp()

    # sam_file_j = bowtie2_wrapper.align_single(
    # fastq=fp_J_region,
    # nthreads=cores_num,
    # region=germ,
    # trim5=str(shortest_J), #REVIEW: should use longest J?
    # write=write,
    # spe=spe,
    # plot=plot,
    # plot_region=germ[0] + ':' + str(germ[1]) + '-' + str(germ[2]),
    # verbose=verbose) #REVIEW: Filter on MAPQ?

    #get rid of everything aligning to J genes
    count = 0
    for read in sam_file_j:
        count += 1
        if not read.is_reverse: #should get excluded by region constrainst
            fq_data.gene_split['germline_J'].add(read.qname) #add to same dict as normal J's

    if verbose:
        print('Germline reverse reads', len(fq_data.gene_split['germline_J']))
        print('Germline reverse and forward reads', count)
    #Get read names into a set
    # germline_reads = set([read.qname for read in region_reads])

    return fq_data


def preclean(V_region, J_region, fq_dict_germ, spe='mmu', verbose = True, misprime_correct=True, discard_germline=True, fast=False):
    '''Seperate out J reads and clean up low quality/unknown reads
    :param V_region: R1 fastq file with V end sequences
    :param J_region: R2 fastq file with J end sequences
    :param fq_dict_germ: fastqHolder object from germline
    :param spe: which species
    :param verbose: output run startswith
    :param misprime_correct: run mispriming correction
    :param fast: use only quick initial J identification (without SSW and mispriming correction)
    :return: fastqHolder object
    '''

    fp_V_region = os.path.abspath(V_region)
    low_qual_V = set()
    low_qual_count_V = 0
    count_v = 0
    #Remove low quality reads from V file
    with general.file_open(fp_V_region) as Vr:
        lines = Vr.read().splitlines()
        for item in general.fastq_parse(lines):
            title = item[0]
            seq = item[1]
            thrd = item[2]
            qual = item[3]
            count_v +=1
            if verbose:
                if count_v % 1000000 == 0:
                    print('Processed', count_v, 'V sequences')

            #Remove low quality and short reads
            if '#' in qual:
                low_qual_count_V +=1
                low_qual_V.add(title.split(' ')[0][1:]) #removes @ as well
                continue

    #Process J reads now

    other_J = 0
    unclear_J = 0
    count = 0
    low_qual_seq = 0
    # mapped_to_other = 0
    germ_count = 0

    j_align = general.SSW_align()

    Js = presets.prs(spe).J_seq()
    ref = j_align.reference(spe, quick_align=fast)

    #Germline set
    germline = fq_dict_germ.gene_split['germline_J']


    # j_size_ord = sorted(presets.prs(spe).J_seq().values(), key=len)
    # shortest_J = len(j_size_ord[0])
    # print(germline)
    fp_J_region = os.path.abspath(J_region)

    misprimed_corrected = 0
    primer_corrected = 0
    first_filter = 0
    with general.file_open(fp_J_region) as Jr:
        lines = Jr.read().splitlines()
        for item in general.fastq_parse(lines):
            title = item[0]
            seq = item[1]
            thrd = item[2]
            qual = item[3]
            count +=1
            if verbose:
                if count % 1000000 == 0:
                    print('Processed', count, 'J sequences')
            #Remove reads that don't map to IgH
            # if title.split(' ')[0][1:] not in igh_reads:
            #     mapped_to_other += 1
            #     continue
            #Remove low quality and short reads and reads that were low quality in V file
            if '#' in qual or len(seq) < 60 or title.split(' ')[0][1:] in low_qual_V:
                low_qual_seq +=1
                fq_dict_germ.gene_split['other_J'].add(title.split(' ')[0][1:]) #removes @ as well
                continue
            #Skip germline reads
            elif discard_germline:
                if title.split(' ')[0][1:] in germline:
                    germ_count += 1
                    # fq_dict_germ.split_gene['germline_J'].add(title.split(' ')[0][1:])
                    continue
            first_filter += 1
            #allows 2 mismatches for 21nt seq
            #Will return touple is misprimed
            J = j_align.align(ref, seq, misprim_cor=misprime_correct, quick_align=fast, spe=spe) #might need to adjust
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
                        if J[0] == J[1]:
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

        print('Number of low quality V reads:', low_qual_count_V)
        print('Number of mispriming corrected reads:', misprimed_corrected)
        print('Number of primer corrected reads:', primer_corrected)
        print('Number of reads processed:', count)
        print('Number of total unclear J reads:', unclear_J)
        print('Number of total other J reads:', other_J)
        print('Number of reads passing first filter', first_filter)
        print('Number of low quality reads:', low_qual_seq)
        print('Number of germline reads removed:', germ_count)

        # print('Total read count: {cnt}'.format(cnt=count))
        # print('Number of reads not assigned to a J: {num}'.format(num=unassigned_J))
        # print('Number of reads with unclear J: {unc}'.format(unc=unclear_J))
        # print('Number of low quality reads: {qual}'.format(qual=low_qual_seq))
        # # print('Number of reads mapped elsewhere: {els}'.format(els=mapped_to_other))
        # print('Number of germline reads: {germ}'.format(germ=germ_count))

    return fq_dict_germ


def check_qual(umi_qual, q_score=30):
    for val in umi_qual:
        phred = ord(val)-33
        assert phred <= 41 and phred >= 0, 'Phred score out side of range 0-41'
        if phred < q_score:
            return True
        else:
            return False



def preclean_assembled(jv_region, fq_dict_germ, q_score, spe='mmu', verbose = True,
                       misprime_correct=True, discard_germline=True, fast=False):
    '''Seperate out J reads and clean up low quality/unknown reads
    :param jv_region: assembled J and V sequences (from PEAR)
    :param fq_dict_germ: fastqHolder object from germline
    :param spe: which species
    :param verbose: output run startswith
    :param misprime_correct: run mispriming correction
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

    j_align = general.SSW_align()

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
        lines = jvr.read().splitlines()
        for item in general.fastq_parse(lines):
            title = item[0]
            seq = item[1]
            thrd = item[2]
            qual = item[3]
            count +=1
            if verbose:
                if count % 1000000 == 0:
                    print('Processed', count, 'J sequences')

            #Remove low quality and short reads and reads that were low quality
            if '#' in qual or len(seq) < 60+14+28: #+ size of barcode and UMI 13 or 14 + longest J
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
            umi_qual = qual[-6:]
            if check_qual(umi_qual, q_score=q_score):
                low_qual_UMI += 1
                continue

            first_filter += 1
            #allows 2 mismatches for 21nt seq
            #Will return touple is misprimed
            J = j_align.align(ref, seq, misprim_cor=misprime_correct, quick_align=fast, spe=spe) #might need to adjust
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
        print('Number of reads passing first filter', first_filter)
        print('Number of low quality reads:', low_qual_seq)
        print('Number of low quality UMIs:', low_qual_UMI)
        print('Number of germline reads removed:', germ_count)

    logging.info('Number of mispriming corrected reads:' + str(misprimed_corrected))
    logging.info('Number of primer corrected reads:' + str(primer_corrected))
    logging.info('Number of reads processed:' + str(count))
    logging.info('Number of total unclear JV reads:' + str(unclear_J))
    logging.info('Number of total mispriming other JV reads:' + str(other_J))
    logging.info('Number of reads passing first filter' + str(first_filter))
    logging.info('Number of low quality reads:' + str(low_qual_seq))
    logging.info('Number of low quality UMIs:' + str(low_qual_UMI))
    logging.info('Number of germline reads removed:' + str(germ_count))

    return fq_dict_germ



#
# def demultiplex(V_region, fq_dict_pcln, barcode_1='GACTCGT', barcode_2='CTGCTCCT', verbose=True):
#     '''Demultiplex V reads and extract the umi
#     :param V_region: R1 fastq file with V end sequences
#     :param fq_dict_pcln: fastqHolder object from preclean
#     :param barcode_1: barcode 1 GACTCGT
#     :param barcode_2: barcode 2 CTGCTCCT
#     :param verbose: print stats
#     '''
#
#     #Test if preclean run on data
#
#     fp_V_region = os.path.abspath(V_region)
#
#     req_files = fq_dict_pcln.gene_split.keys()
#     mis_files = fq_dict_pcln.misprimed.keys()
#
#     with general.file_open(fp_V_region) as Vr:
#         lines = Vr.read().splitlines()
#         for item in general.fastq_parse(lines):
#             title = item[0]
#             seq = item[1]
#             thrd = item[2]
#             qual = item[3]
#             for key in req_files:
#                 if 'germline' not in key and 'other' not in key: #Don't need to demultiplex germline reads!
#                     if title.split(' ')[0][1:] in fq_dict_pcln.gene_split[key]: #missing missprimed!!
#
#                         # if Levenshtein.distance(barcode_1, barcode_2) <= 2:
#                         #     raise Exception('Barcodes not unique enough!')
#
#                         #Allowing 2 mismatches, first 6 bases UMI random N bases
#                         # if Levenshtein.distance(seq[6:6+len(barcode_1)], barcode_1) <= 2:
#                         #     fq_dict_pcln.demultiplex[key + '_' + barcode_1].add(title.split(' ')[0][1:])
#                         # elif Levenshtein.distance(seq[6:6+len(barcode_2)], barcode_2) <= 2:
#                         #     fq_dict_pcln.demultiplex[key + '_' + barcode_2].add(title.split(' ')[0][1:])
#                         # else:
#                         #     fq_dict_pcln.demultiplex['unassigned_J'].add(title.split(' ')[0][1:])
#
#                         #perfect match wanted just so indels don't mess up UMI
#                         if seq[6:6+len(barcode_1)] == barcode_1:
#                             fq_dict_pcln.demultiplex[key + '_' + barcode_1].add(title.split(' ')[0][1:])
#                         elif seq[6:6+len(barcode_2)] == barcode_2:
#                             fq_dict_pcln.demultiplex[key + '_' + barcode_2].add(title.split(' ')[0][1:])
#                         else:
#                             fq_dict_pcln.demultiplex['unassigned_J'].add(title.split(' ')[0][1:])
#
#             #Demultiplex misprimed reads
#             for key in mis_files:
#                 if title.split(' ')[0][1:] in fq_dict_pcln.misprimed[key]: #or fq_dict_pcln.misprimed[key].keys() does same thing
#                     #perfect match wanted just so indels don't mess up UMI
#                     if seq[6:6+len(barcode_1)] == barcode_1:
#                         fq_dict_pcln.demultiplex[key + '_' + barcode_1].add(title.split(' ')[0][1:])
#                     elif seq[6:6+len(barcode_2)] == barcode_2:
#                         fq_dict_pcln.demultiplex[key + '_' + barcode_2].add(title.split(' ')[0][1:])
#                     else:
#                         fq_dict_pcln.demultiplex['unassigned_J'].add(title.split(' ')[0][1:])
#
#
#     if verbose:
#         for key in fq_dict_pcln.demultiplex.keys():
#             print('Number of reads in', key, ':', len(fq_dict_pcln.demultiplex[key]))
#
#     return fq_dict_pcln


def demultiplex_assembled(jv_region, fq_dict_pcln, barcode_1='GACTCGT', barcode_2='CTGCTCCT', verbose=True):
    '''Demultiplex V reads and extract the umi from assembled reads
    :param V_region: R1 fastq file with V end sequences
    :param fq_dict_pcln: fastqHolder object from preclean
    :param barcode_1: barcode 1 GACTCGT
    :param barcode_2: barcode 2 CTGCTCCT
    :param verbose: print stats
    '''

    #Test if preclean run on data

    #Barcode 1 and 2 need reverse complement
    rv_barcode_1 = general.reverse_complement(barcode_1)
    rv_barcode_2 = general.reverse_complement(barcode_2)

    fp_jv_region = os.path.abspath(jv_region)

    req_files = fq_dict_pcln.gene_split.keys()
    mis_files = fq_dict_pcln.misprimed.keys()


    with general.file_open(fp_jv_region) as jvr:
        lines = jvr.read().splitlines() #REVIEW: don't need to do this
        for item in general.fastq_parse(lines):
            title = item[0]
            seq = item[1]
            thrd = item[2]
            qual = item[3]
            for key in req_files:
                if 'germline' not in key and 'other' not in key: #Don't need to demultiplex germline reads!
                    if title.split(' ')[0][1:] in fq_dict_pcln.gene_split[key]:

                        #perfect match wanted just so indels don't mess up UMI
                        if seq[-6-len(barcode_1):-6] == rv_barcode_1:
                            fq_dict_pcln.demultiplex[key + '_' + barcode_1].add(title.split(' ')[0][1:])
                        elif seq[-6-len(barcode_2):-6] == rv_barcode_2:
                            fq_dict_pcln.demultiplex[key + '_' + barcode_2].add(title.split(' ')[0][1:])
                        else:
                            fq_dict_pcln.demultiplex['unassigned_J'].add(title.split(' ')[0][1:])

            #Demultiplex misprimed reads
            for key in mis_files:
                if title.split(' ')[0][1:] in fq_dict_pcln.misprimed[key]: #or fq_dict_pcln.misprimed[key].keys() does same thing

                    #perfect match wanted just so indels don't mess up UMI
                    if seq[-6-len(barcode_1):-6] == rv_barcode_1:
                        fq_dict_pcln.demultiplex[key + '_' + barcode_1].add(title.split(' ')[0][1:])
                    elif seq[-6-len(barcode_2):-6] == rv_barcode_2:
                        fq_dict_pcln.demultiplex[key + '_' + barcode_2].add(title.split(' ')[0][1:])
                    else:
                        fq_dict_pcln.demultiplex['unassigned_J'].add(title.split(' ')[0][1:])

    if verbose:
        for key in fq_dict_pcln.demultiplex.keys():
            print('Number of reads in', key, ':', len(fq_dict_pcln.demultiplex[key]))

    return fq_dict_pcln



#REVIEW: that misprimed reads are being written out!
# def write_fastq(V_region, J_region, fq_dict_demult, q_score, prefix=None, out_dir=None,
#                 barcode_1='GACTCGT', barcode_2='CTGCTCCT'):
#     '''Write out precleaned fastq files
#     :param V_region: R1 fastq file with V end sequences
#     :param J_region: R2 fastq file with J end sequences
#     :param fq_dict_demult: a fastqHolder object from demultiplex
#     :param q_score: min phred quality for bases in the UMI
#     :param prefix: prefix for out files, default is basename of input_J_R2
#     :param out_dir: output directory, default is preclean folder created in input directory
#
#     :param barcode_1: barcode 1 GACTCGT
#     :param barcode_2: barcode 2 CTGCTCCT
#     '''
#     # J_to_V = {'unclear_J':'unclear_V', 'germline_J':'germline_V', 'other_J':'other_V', 'J1':'V1', 'J2':'V2', 'J3':'V3', 'J4':'V4'} #generate automatically
#
#     req_files = fq_dict_demult.demultiplex.keys()
#     germ_files = fq_dict_demult.gene_split.keys()
#
#     fp_V_region = os.path.abspath(V_region)
#     fp_J_region = os.path.abspath(J_region)
#     dir_nam = os.path.dirname(fp_J_region)
#
#     if prefix == None:
#         prefix_J = os.path.basename(J_region).split('.')[0]
#         prefix_V = os.path.basename(V_region).split('.')[0]
#     else:
#         prefix_J = prefix
#         prefix_V = prefix
#
#     #add prefix to out directory
#     if out_dir == None:
#         out_dir = dir_nam + '/' + prefix_V + '_preclean'
#         print(out_dir)
#         try:
#             os.mkdir(out_dir)
#         except FileExistsError:
#             print('Default directory', out_dir, 'already exists. Might overwrite files!')
#
#
#     print('Writing demultiplexed files to:', out_dir)
#
#     # if seperate_files:
#     #     #unassigned, J_barcode
#     #     for key in req_files:
#     #         fq_dict_demult.write_demultiplex(fp_J_region, key, out_dir + '/' + prefix_J + '_' + key, j_end=True)
#     #         fq_dict_demult.write_demultiplex(fp_V_region, key, out_dir + '/' + prefix_V + '_' + key.replace('J', 'V'))
#
#     # else:
#
#     #write everything into a single file
#     key_list = []
#     for key in req_files:
#         if 'unassigned' not in key:
#             key_list.append(key)
#         else: #write unassigned into a seperate file
#             fq_dict_demult.write_demultiplex_unassigned(fp_J_region, key, out_dir + '/' + prefix_J + '_' + key, j_end=True)
#             fq_dict_demult.write_demultiplex_unassigned(fp_V_region, key, out_dir + '/' + prefix_V + '_' + key.replace('J', 'V'))
#
#     #Write out everything else
#     fq_dict_demult.write_demultiplex_umi_extract(fp_V_region, fp_J_region, key_list,
#                                                  q_score=q_score, V_out_path=out_dir + '/' + prefix_V + '_' + 'all_V',
#                                                  J_out_path=out_dir + '/' + prefix_J + '_' + 'all_J', br1=barcode_1, br2=barcode_2)
#
#     #write out germline and other
#     for key in germ_files:
#         if 'germline' in key or 'other' in key:
#             fq_dict_demult.write_preclean(fp_J_region, key, out_dir + '/' + prefix_J + '_' + key)
#             fq_dict_demult.write_preclean(fp_V_region, key, out_dir + '/' + prefix_V + '_' + key.replace('J', 'V'))
#

def write_assembled(jv_region, fq_dict_demult, prefix=None, out_dir=None, barcode_1='GACTCGT', barcode_2='CTGCTCCT'):
    '''Write out precleaned fastq files
    :param jv_region: (PEAR) assembled fastq file with J and V end sequences
    :param fq_dict: a fastqHolder object from demultiplex
    :param q_score: min phred quality for bases in the UMI
    :param prefix: prefix for out files, default is basename of input_J_R2
    :param out_dir: output directory, default is preclean folder created in input directory
    :param barcode_1: barcode 1 GACTCGT
    :param barcode_2: barcode 2 CTGCTCCT
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


    #add prefix to out directory
    if out_dir == None:
        out_dir = dir_nam + '/' + prefix_jv + '_preclean'
        print(out_dir)
        try:
            os.mkdir(out_dir)
        except FileExistsError:
            print('Default directory', out_dir, 'already exists. Might overwrite files!')


    print('Writing demultiplexed files to:', out_dir)

    #write everything into a single file
    key_list = []
    for key in req_files:
        if 'unassigned' not in key:
            key_list.append(key)
        else: #write unassigned into a seperate file
            fq_dict_demult.write_demultiplex_unassigned(fp_jv_region, key, out_dir + '/' + prefix_jv + '_' + key, j_end=True)

    #Write out everything else
    fq_dict_demult.write_demultiplex_umi_extract_assembled(fp_jv_region, key_list,
                                                           out_path=out_dir + '/' + prefix_jv + '_' + 'all_jv',
                                                           br1=barcode_1, br2=barcode_2)

    #write out germline and other
    for key in germ_files:
        if 'germline' in key or 'other' in key:
            fq_dict_demult.write_preclean(fp_jv_region, key, out_dir + '/' + prefix_jv + '_' + key)





def gemline_removed_qc(V_region, out_dir, spe='mmu', prefix=None, cores_num=8, verbose=False):
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

        run_bowtie2 = general.bowtie2()
        run_bowtie2.align_single(fastq=name, nthreads=cores_num, trim5=str(shortest_J), spe=spe, verbose=verbose)

        run_bowtie2.plot(plot_region=germ[0] + ':' + str(germ[1]) + '-' + str(germ[2]), spe=spe)

        run_bowtie2.write(region=germ)

        run_bowtie2.del_tmp()


        # sam_file = bowtie2_wrapper.align_single(
        # fastq=name,
        # nthreads=cores_num,
        # region=germ,
        # trim5=str(shortest_J),
        # write=True,
        # spe=spe,
        # plot=True,
        # samtools_mapq=1, #no multimapping reads
        # plot_region=germ[0] + ':' + str(germ[1]) + '-' + str(germ[2]),
        # verbose=verbose,
        # flags=['--very-sensitive']) #Keep unaligned reads (the V end might align!)




def parse_args():
    '''Set up parser
    '''
    parser = argparse.ArgumentParser(description='BabrahamLinkON Preclean')

    # group = parser.add_mutually_exclusive_group(required=True)
    parser.add_argument('-v', '--V_r1', dest='input_V', type=str, metavar='v.fastq', nargs='+', help='Input fastq file(s) with V end sequences')
    parser.add_argument('-j', '--J_r2', dest='input_J', type=str, metavar='j.fastq', nargs='+', help='Input fastq file(s) with J end sequences')
    parser.add_argument('-jv', '--jv', dest='input_jv', type=str, metavar='jv.fastq', nargs='+', help='Input fastq file(s) from PEAR with J (forward) end and V (reverse) end sequences')
    parser.add_argument('--short', action='store_true', help='Pair if using short reads <250bp')
    parser.add_argument('--species', dest='species', default='mmu', type=str, help='Which species (mmu hsa), default: mmu')
    parser.add_argument('--fast', action='store_true', help='Perform fast inaccurate J identification (not recommended if deduplicating using J)')
    parser.add_argument('--mispriming', action='store_true', help='Perform mispriming correction (not compatible with --fast)')
    parser.add_argument('--prefix', dest='prefix', type=str, metavar='N', nargs='+', help='Prefix of the output file (need to provide one for each input)')
    parser.add_argument('--out', dest='out_dir', type=str, help='Output direcotry')
    parser.add_argument('--cores', dest='nthreads', default=1, type=int, help='Number of cores to use, default: 1')
    parser.add_argument('--br1', dest='br1', default='GACTCGT', type=str, help='Default: GACTCGT')
    parser.add_argument('--br2', dest='br2', default='CTGCTCCT', type=str, help='Default: CTGCTCCT')
    parser.add_argument('--verbose', action='store_true', help='Print detailed progress')
    parser.add_argument('--plot', action='store_true', help='Plot alignments')

    parser.add_argument('--keep_germline', action='store_false', help='Skip germline removal step')
    parser.add_argument('-q', '--q_score', dest='q_score', type=int, default=30, help='Minimum Phred quality score for bases in UMI')

    parser.add_argument('--plot_QC', action='store_true', help='QC plot showing if all germline reads were removed (few will be present J-J rearrangements)')

    opts = parser.parse_args()

    return opts



def main():


    #argparse
    opts = parse_args()

    if opts.input_V != None and opts.input_J != None and len(opts.input_V) == len(opts.input_J):

        #Repeat, maybe codense
        if opts.prefix==None:
            prefix = os.path.basename(opts.input_V[0]).split('.')[0]

        fp_V_region = os.path.abspath(opts.input_V[0])
        dir_path =  os.path.dirname(fp_V_region)

        logging.basicConfig(level=logging.DEBUG, filename=dir_path + '/' + prefix + '_preclean/' + prefix + '_preclean.log', filemode='a+',
                            format='%(asctime)-15s %(levelname)-8s %(message)s')


        assembled_file = assemble(opts.input_V[0], opts.input_J[0], threads=opts.nthreads, prefix=opts.prefix, short=opts.short)

        # prefix.assembled.fastq
        # prefix.discarded.fastq
        # prefix.unassembled.forward.fastq
        # prefix.unassembled.reverse.fastq
        # prefix.all_assembled.fastq

        germ_assembled = germline(assembled_file + '.assembled.fastq', spe=opts.species, cores_num=opts.nthreads, plot=opts.plot, verbose=opts.verbose)

        if opts.short:
            #TODO: running twice will overwrite the output plot!!
            germ_unassembled = germline(assembled_file + '.unassembled.forward.fastq', spe=opts.species, cores_num=opts.nthreads, plot=opts.plot, verbose=opts.verbose)

            #will results in fq_data.gene_split['germline_J']) merge?
            germ_assembled.gene_split['germline_J'].update(germ_unassembled.gene_split['germline_J']) #set

            #Merge the two files into one (pairing unassembled reads)
            fq_clean = preclean_assembled(assembled_file + '.all_assembled.fastq', germ_assembled, q_score=opts.q_score, spe=opts.species, verbose=opts.verbose,
                                          misprime_correct=opts.mispriming, discard_germline=opts.keep_germline, fast=opts.fast)

            fq_demultiplex = demultiplex_assembled(assembled_file + '.all_assembled.fastq', fq_clean, barcode_1=opts.br1, barcode_2=opts.br2, verbose=opts.verbose)

            write_assembled(assembled_file + '.all_assembled.fastq', fq_demultiplex, prefix=opts.prefix, out_dir=opts.out_dir)
        else:
            fq_clean = preclean_assembled(assembled_file + '.assembled.fastq', germ_assembled, q_score=opts.q_score, spe=opts.species, verbose=opts.verbose,
                                          misprime_correct=opts.mispriming, discard_germline=opts.keep_germline, fast=opts.fast)

            fq_demultiplex = demultiplex_assembled(assembled_file + '.assembled.fastq', fq_clean, barcode_1=opts.br1, barcode_2=opts.br2, verbose=opts.verbose)

            write_assembled(assembled_file + '.assembled.fastq', fq_demultiplex, prefix=opts.prefix, out_dir=opts.out_dir)


        # for item in range(len(opts.input_V)):
        #
        #     if not os.path.isfile(opts.input_V[item]) or not os.path.isfile(opts.input_J[item]):
        #         raise FileNotFoundError('File does not exit!')
        #
        #     germ = germline(opts.input_J[item], spe=opts.species, cores_num=opts.nthreads, plot=opts.plot, verbose=opts.verbose)
        #
        #     fq_clean = preclean(opts.input_V[item], opts.input_J[item], germ, spe=opts.species, verbose=opts.verbose,
        #                         misprime_correct=opts.mispriming, discard_germline=opts.keep_germline, fast=opts.fast)
        #
        #     fq_demultiplex = demultiplex(opts.input_V[item], fq_clean, barcode_1=opts.br1, barcode_2=opts.br2, verbose=opts.verbose)
        #
        #     write_fastq(opts.input_V[item], opts.input_J[item], fq_demultiplex, prefix=opts.prefix, out_dir=opts.out_dir)
        #
        #     if opts.plot_QC:
        #         gemline_removed_qc(opts.input_V[item], out_dir=opts.out_dir, spe=opts.species, prefix=opts.prefix, cores_num=opts.nthreads, verbose=opts.verbose)

    elif opts.input_jv != None: #using PEAR assembled
        for item in range(len(opts.input_jv)):

            if not os.path.isfile(opts.input_jv[item]):
                raise FileNotFoundError('File does not exit!')

            germ = germline(opts.input_jv[item], spe=opts.species, cores_num=opts.nthreads, plot=opts.plot, verbose=opts.verbose)

            fq_clean = preclean_assembled(opts.input_jv[item], germ, q_score=opts.q_score, spe=opts.species, verbose=opts.verbose,
                                          misprime_correct=opts.mispriming, discard_germline=opts.keep_germline, fast=opts.fast)

            fq_demultiplex = demultiplex_assembled(opts.input_jv[item], fq_clean, barcode_1=opts.br1, barcode_2=opts.br2, verbose=opts.verbose)

            write_assembled(opts.input_jv[item], fq_demultiplex, prefix=opts.prefix, out_dir=opts.out_dir)

    else: #No input supplied
        raise FileNotFoundError('No input files supplied!')



    # V_igh_reads = aligned_to_igh('/media/chovanec/My_Passport/Sync/BabrahamLinkON/tests/test1.fq.gz', plot=True)
    #
    # fq_clean = preclean('/media/chovanec/My_Passport/Sync/BabrahamLinkON/tests/test2.fq.gz', V_igh_reads)
    #
    # write_preclean('/media/chovanec/My_Passport/Sync/BabrahamLinkON/tests/test1.fq.gz',
    # '/media/chovanec/My_Passport/Sync/BabrahamLinkON/tests/test2.fq.gz', fq_clean)

    # germ = germline('/media/chovanec/My_Passport/Sync/BabrahamLinkON/tests/lane3_WT_FrBC_1_GCCAAT_L003_R2_val_2_40k.fq', cores_num=8, plot=True, write=True)

    # V_igh_reads = aligned_to_igh('/media/chovanec/My_Passport/Sync/BabrahamLinkON/tests/lane3_WT_FrBC_1_GCCAAT_L003_R1_val_1_40k.fq', cores_num=8, plot=True)

    # fq_clean = preclean('/media/chovanec/My_Passport/Sync/BabrahamLinkON/tests/lane3_WT_FrBC_1_GCCAAT_L003_R2_val_2_40k.fq', germ)
    #
    # fq_demultiplex = demultiplex('/media/chovanec/My_Passport/Sync/BabrahamLinkON/tests/lane3_WT_FrBC_1_GCCAAT_L003_R1_val_1_40k.fq', fq_clean)
    #
    # write_fastq('/media/chovanec/My_Passport/Sync/BabrahamLinkON/tests/lane3_WT_FrBC_1_GCCAAT_L003_R1_val_1_40k.fq',
    # '/media/chovanec/My_Passport/Sync/BabrahamLinkON/tests/lane3_WT_FrBC_1_GCCAAT_L003_R2_val_2_40k.fq', fq_demultiplex)



    # V_igh_reads = aligned_to_igh('/media/chovanec/My_Passport/CS_VDJ_seq/lane4850_CGATGT_1_BC_L001_R1_val_1_20k.fq', cores_num=8, spe='mm', write=True, plot=True)
    # V_igh_reads = aligned_to_igh('/media/chovanec/My_Passport/CS_VDJ_seq/lane4850_CGATGT_1_BC_L001_R1_val_1_20k.fq', cores_num=8, spe='hs', write=True)
    # V_igh_reads = aligned_to_igh('/media/chovanec/My_Passport/CS_VDJ_seq/lane4850_ACAGTG_5_BC_L001_R1_val_1_20k.fq', cores_num=8, spe='hs', write=True)

    #45.72

if __name__ == "__main__":
    main()
