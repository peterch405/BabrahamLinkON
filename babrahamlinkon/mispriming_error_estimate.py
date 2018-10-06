
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

from babrahamlinkon import general, presets, bowtie2_wrapper, mispriming_correction
from collections import defaultdict
import os
from tqdm import tqdm
import pandas as pd
import numpy as np

#read in first 10 lines of germline fastq to determine if it is long reads or short, do they need trimming?


def align_germline(jv_region, spe, plot, write, verbose, threads_num):
    j_size_ord = sorted(presets.prs(spe).J_seq().values(), key=len)
    shortest_J = len(j_size_ord[0])

    germ = presets.prs(spe).germline()

    fp_jv_region = os.path.abspath(jv_region)


    run_bowtie2 = bowtie2_wrapper.bowtie2()
    run_bowtie2.align_single(fastq=fp_jv_region, nthreads=threads_num, spe=spe, verbose=verbose)

    if plot:
        run_bowtie2.plot(plot_region=germ[0] + ':' + str(germ[1]) + '-' + str(germ[2]), spe=spe)

    sam_file_jv = run_bowtie2.pysam_out(region=germ, fetch=True)

    if write:
        run_bowtie2.write(region=germ)

    run_bowtie2.del_tmp()

    region_reads = [] #can be read multiple times
    for read in sam_file_jv:
            region_reads.append(read)

    return region_reads



def extract_j_genes(region_reads, spe):
    seq_dict = defaultdict(lambda: defaultdict(int))
    location_j = presets.prs(spe).J_location()
    for read in region_reads:
        if spe == 'mmuk': #is opposite strand than Igh
            if read.is_reverse == True:
                for k,v in presets.prs(spe).J_location().items():
                    if read.pos >= (location_j[k]['start']-100) and read.pos <=(location_j[k]['end']+50):
                        seq_dict[k][general.reverse_complement(read.seq)] += 1 #are the first 4 bps removed?
        else:
            if read.is_reverse == False:
                for k,v in presets.prs(spe).J_location().items():
                    if read.pos >= (location_j[k]['start']-50) and read.pos <=(location_j[k]['end']+100):
                        seq_dict[k][read.seq] += 1
    return seq_dict


def mispriming_matrix(seq_dict, spe):
    before = defaultdict(lambda: defaultdict(int))
    after = defaultdict(lambda: defaultdict(int))
    location_j = presets.prs(spe).J_location()

    for j in location_j.keys():
        major_reads = []
        major_reads_count = []
        for key, val in tqdm(seq_dict[j].items()):
        # for val in tqdm(range(0,len(seq_dict))):
            # if list(seq_dict.values())[val] > 1:
            major_reads.append(key)
            major_reads_count.append(val)



        print('Indentify priming sequence:')
        align_J = mispriming_correction.SSW_align()
        ref = align_J.reference(spe)
        # J_len = len(presets.prs(spe).J_seq()[j])

        # after_j_counts = defaultdict(int)
        # before_j_counts = defaultdict(int)

        # [initial, misprime corrected, corrected seq]

        for read in major_reads:

            identity_j = align_J.align(ref, read, no_misprim_cor=False, quick_align=False, spe=spe)

            if isinstance(identity_j, list):
                after[j][identity_j[1]] += 1
                before[j][identity_j[0]] += 1

            else:
                after[j][identity_j] += 1
                before[j][identity_j] += 1


    return(before, after)

# test_df = pd.read_csv('/media/chovanec/My_Passport/CS_VDJ_seq/lane4850_ACAGTG_5_BC_L001_R1_val_1_preclean/tmp.csv', index_col=0)
# df_after = pd.read_csv('/media/chovanec/My_Passport/CS_VDJ_seq/df_after_sub.csv', index_col=0)
# df_before = pd.read_csv('/media/chovanec/My_Passport/CS_VDJ_seq/df_before_sub.csv', index_col=0)
#
# test_df_sub = test_df.loc[:,test_df.columns[~test_df.columns.str.contains('other|unclear')]]
# test_df_sub.fillna(float(0), inplace=True)
# df_after_sub = variant_merge(df_after, 'hsa')
# df_before_sub = variant_merge(df_before, 'hsa')
# spe = 'hsa'
#
# sum_df = defaultdict()
# for k in df_after.index:
#     to_merge = df_after.columns.str.contains('^' + '$|^'.join(j_merge.get(k)) + '$')
#     sum_df[k] = df_after.loc[:,to_merge].sum(axis=1)
#
# #combine all dfs
# simple_df = pd.DataFrame.from_dict(sum_df)


def variant_merge(m_df, spe):
    #for human, the mispriming is more complicated with multiple variants
    #in order to get a value merge variants into groups to reflect the index
    all_js = sorted(list(presets.prs(spe).J_seq().keys()))
    #make merge dict
    j_merge = defaultdict(set)
    for i in all_js:
        if '.' in i:
            j_merge[i.split('.')[0]].add(i)
        else:
            j_merge[i].add(i)
    #also add J2P into dict
    j_merge['J2P'] = {'J2P'}

    #sum rows of columns to merge
    sum_df = defaultdict()
    for k in m_df.index:
        if j_merge.get(k) is None:
            sum_df[k] = 0
        else:
            to_merge = m_df.columns.str.contains('^' + '$|^'.join(j_merge.get(k)) + '$')
            sum_df[k] = m_df.loc[:,to_merge].sum(axis=1)

    #combine all dfs
    simple_df = pd.DataFrame.from_dict(sum_df)

    return(simple_df)


#make a table of clear and unclear, before and after
def make_table(before, after, spe, out_file):


    all_js = sorted(list(presets.prs(spe).J_location().keys()))

    df_before = pd.DataFrame.from_dict(before, orient='index')
    df_before_cols = sorted(df_before.columns.tolist())
    df_before = df_before[df_before_cols]

    df_after = pd.DataFrame.from_dict(after, orient='index')
    df_after_cols = sorted(df_after.columns.tolist())
    df_after = df_after[df_after_cols]

    #for small datasets if NaN present change to 0
    df_before.fillna(float(0), inplace=True)
    df_after.fillna(float(0), inplace=True)

    if spe == 'hsa':
        df_after_sub = variant_merge(df_after, spe)
        df_before_sub = variant_merge(df_before, spe)
    else:
        df_after_sub = df_after[all_js]
        df_before_sub = df_before[all_js]

    #subset table to only properly identified J's
    #after

    df_after_sum = pd.DataFrame(df_after_sub.sum())

    diag = pd.DataFrame(np.diag(df_after_sub), index=df_after_sub.index)

    out_df_after = pd.concat([df_after_sum, abs(df_after_sum.subtract(diag))], axis=1)
    out_df_after.columns = ['Total', 'Misprimed']
    out_df_after['Percent error'] = (out_df_after['Misprimed']/out_df_after['Total'])*100

    #before

    df_before_sum = pd.DataFrame(df_before_sub.sum())

    diag_before = pd.DataFrame(np.diag(df_before_sub), index=df_before_sub.index)

    out_df_before = pd.concat([df_before_sum, abs(df_before_sum.subtract(diag_before))], axis=1)
    out_df_before.columns = ['Total', 'Misprimed']
    out_df_before['Percent error'] = (out_df_before['Misprimed']/out_df_before['Total'])*100


    #concat all tables into one data.frame
    with open(out_file, 'a') as out_f:
        out_f.write('Before mispriming correction \n')
        df_before.to_csv(out_f, header=True)
        out_f.write('\n')
        out_f.write('After mispriming correction \n')
        df_after.to_csv(out_f, mode='a', header=True)
        out_f.write('\n')
        out_f.write('Before mispriming error \n')
        out_df_before.to_csv(out_f, mode='a', header=True)
        out_f.write('\n')
        out_f.write('After mispriming error \n')
        out_df_after.to_csv(out_f, mode='a', header=True)
