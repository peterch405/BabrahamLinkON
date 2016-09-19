#!/usr/bin/env python3

from collections import defaultdict, Counter
import pandas as pd
import numpy as np
import argparse
from babrahamlinkon import general, bowtie2_wrapper

# from general import SSW_align, species
# from bowtie2_wrapper import align_single

#TODO: Test, produce mispriming meaningful output

def mispriming_germline(fastq_germline_J, V_fastq, aligned_j, spe='mmu', out_file=None):
    '''
    '''

    germ = general.species(spe).germline()

    fp_V_region = os.path.abspath(V_region)

    sam_file_v = bowtie2_wrapper.align_single(
    fastq=fp_V_region,
    nthreads=cores_num,
    region=germ,
    write=write,
    spe=spe,
    plot=plot,
    plot_region=germ[0] + ':' + str(germ[1]) + '-' + str(germ[2]),
    verbose=verbose)

    v_germline_reads = set([read.qname for read in sam_file_v])

    sam_file = bowtie2_wrapper.align_single(
        spe=spe,
        fastq=fastq_germline_J,
        nthreads=2,
        region=germ,
        trim5=0)


    region_reads = [] #can be read multiple times
    for read in sam_file:
        if read in v_germline_reads:
            region_reads.append(read)


    j_size_ord = sorted(general.species(spe).J_seq().values(), key=len)
    shortest_J = len(j_size_ord[0])


    location_j = general.species(spe).J_location()
    seq_dict = defaultdict(int)
    #Accumulate germline reads in dictionary
    for read in region_reads:
            if read.is_reverse == False: #should get excluded by region constrainst
                if read.qname in v_germline_reads: #true germline
    #             seq_dict[read.seq] += 1
                    if read.pos >=(location_j[aligned_j]['start']-50) and read.pos <=(location_j[aligned_j]['end']+100):
                        seq_dict[read.seq[:31]] += 1 #21+10 = full J



    #Find the most frequent reads
    print('Major reads present:')
    #10 most common reads
    most_common = dict(Counter(seq_dict).most_common(10))
    major_reads = []
    for val in range(0,len(most_common)):
        # print(list(most_common.keys())[val], list(most_common.values())[val])
        major_reads.append(list(most_common.keys())[val])

    df = pd.DataFrame.from_dict(most_common, orient='index')
    df.columns = ['Count']


    iden = []

    print('Indentify priming sequence:')
    align_J = general.SSW_align()
    ref = align_J.reference(spe)
    J_len = len(general.species(spe).J_seq()[aligned_j])

    for read in major_reads:
        identity_j = align_J.align(ref, read)
        if isinstance(identity_j, tuple):
            if identity_j[0] != aligned_j:
                iden.append((identity_j[0],read))
            else:
                iden.append(identity_j[0])
        else:
            if identity_j != aligned_j:
                iden.append((identity_j,read))
            else:
                iden.append(identity_j)
    align_J.print_align(print_out=True)

    df['Identity'] = iden
    # df[algn[0]] = algn[1:]

    print(df)
    # with open(out_file, 'w') as out:
    #     df.to_html(out)



def parse_args():
    parser = argparse.ArgumentParser(description='BabrahamLinkON Germline mispriming')

    parser.add_argument('-j', '--J_germline', dest='input_J', type=str, required=True, help='Input germline file with J end sequences')
    parser.add_argument('-v', '--V_fastq', dest='input_V', type=str, required=True, help='Input R1 file with V end sequences')
    parser.add_argument('--which_J', dest='which_J', type=str, required=True, help='Identity of germline sequences')
    parser.add_argument('--species', dest='species', default='mmu', type=str, help='Which species (mmu, hsa), default: mmu')

    opts = parser.parse_args()

    return opts



def main():

    # if sys.version_info.major != 3:
    #     raise RuntimeError('run_mixcr requires Python 3')

    #argparse
    opts = parse_args()

    mispriming_germline(opts.input_J, opts.input_V, opts.which_J, opts.species)

    # mispriming_germline('/media/chovanec/My_Passport/Sync/BabrahamLinkON/tests/preclean/lane3_WT_FrBC_1_GCCAAT_L003_R2_val_2_40k_germline_J1', 'J1')
    # mispriming_germline('/media/chovanec/My_Passport/Sync/BabrahamLinkON/tests/preclean/lane3_WT_FrBC_1_GCCAAT_L003_R2_val_2_40k_germline_J2', 'J2')
    # mispriming_germline('/media/chovanec/My_Passport/Sync/BabrahamLinkON/tests/preclean/lane3_WT_FrBC_1_GCCAAT_L003_R2_val_2_40k_germline_J3', 'J3')
    # mispriming_germline('/media/chovanec/My_Passport/Sync/BabrahamLinkON/tests/preclean/lane3_WT_FrBC_1_GCCAAT_L003_R2_val_2_40k_germline_J4', 'J4')


if __name__ == "__main__":
    main()
