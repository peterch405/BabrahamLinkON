#!/usr/bin/env python3

from collections import defaultdict, Counter
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import pandas as pd
import numpy as np
import Levenshtein
from babrahamlinkon import deduplicate
import os
import argparse
import glob

# %matplotlib inline
# %config InlineBackend.figure_format = 'svg'
# '/media/chovanec/My_Passport/Dan_VDJ-seq_cycles_new/J_merged_1c_Deduplicated_test/J_merged_1c_dedup.0.fasta_db-pass.tab'


def ambigious_calls(item, full_name=False):
    '''remove v calls that are ambigous'''
    v_call = item.split(',')
    if len(v_call) >1:
        first_call = v_call[0].split('*')[0]
        for calls in range(1,len(v_call)):
            if v_call[calls].split('*')[0] == first_call:
                continue
            else:
                return None
        if full_name:
            return v_call[0]
        else:
            return first_call
    else:
        if full_name:
            return v_call[0]
        else:
            return v_call[0].split('*')[0]


# def lev_adj_list_directional_adjacency(umis, counts, threshold=1):
#     ''' identify all umis within the levenshtein distance threshold
#     and where the counts of the first umi is > (2 * second umi counts)-1
#     will have duplicates'''
#
#     #should be 50% faster
#     adj_list = defaultdict(list)
#     for i,umi in enumerate(umis):
#         a1 = adj_list[umi]
#         c1 = counts[umi]
#         for j in range(i+1,len(umis)):
#             umi2 = umis[j] #dict_keys object doesn't support indexing
#
#             if Levenshtein.distance(umi, umi2) == threshold:
#                 c2 = counts[umi2]
#                 if c1 >= (c2*2)-1:
#                     adj_list[umi].append(umi2)
#                 if c2 >= (c1*2)-1:
#                     adj_list[umi2].append(umi)
#
#     return adj_list


def lev_adj_list_adjacency(umis, counts, threshold=1):
    ''' identify all umis within the levenshtein distance threshold'''

    #should be 50% faster
    adj_list = defaultdict(list)
    for i,umi in enumerate(umis):
        a1 = adj_list[umi]
        for j in range(i+1,len(umis)):
            umi2 = umis[j] #dict_keys object doesn't support indexing
            if Levenshtein.distance(umi, umi2) == threshold:
                adj_list[umi].append(umi2)

    return adj_list


def change_v_call(row):
    count = 0
    if row['V_SCORE'] <= 50:
        count += 1
        v_gene = row['SEQUENCE_ID'].split('_')[-3].upper() #would normally be the barcode!
        #replace v gene call with bowtie alignment call

        if len(v_gene) < 1:
            return np.NaN
    else:
        v_gene = row['V_CALL']
    return v_gene


def read_changeo_out(tab_file, out, prefix, plot=False, retain_nam=False, minimal=False, short=False):

    igblast_out = pd.DataFrame()
    df_list = []
    for f in tab_file:
        df = pd.read_table(f, header=0)
        #add column with file identity
        id_label = os.path.basename(f).split('.')[0]
        df['file_ID'] = pd.Series(np.repeat(id_label, len(df.index)))
        df_list.append(df)
    igblast_out = pd.concat(df_list)
    igblast_out.reset_index(drop=True, inplace=True)
    # len(igblast_out)

    if plot:
        with PdfPages(out + '/' + prefix + '_score_plots.pdf') as pdf_out:
        #Plot V and J scores
            labels, values = zip(*Counter(igblast_out['V_SCORE']).items())
            non_dedup_values = tuple(l*v for l, v in zip(labels, values))

            plt.figure()
            plt.bar(labels, non_dedup_values)
            my_plot = plt.axvline(50, linestyle='dashed', linewidth=2).get_figure()
            pdf_out.savefig(my_plot)

            labels, values = zip(*Counter(igblast_out['J_SCORE']).items())
            non_dedup_values = tuple(l*v for l, v in zip(labels, values))

            plt.figure()
            plt.bar(labels, non_dedup_values)
            my_plot = plt.axvline(35, linestyle='dashed', linewidth=2).get_figure()
            pdf_out.savefig(my_plot)

    #if low quality replace by bowtie call in qname
    if short:
        # HWI-1KL136:214:D1MR5ACXX:5:1103:17395:138047_J4_Ighv13-2_GTGTCTAC_11

        igblast_out_v = igblast_out.apply(change_v_call, axis=1)
        igblast_out['V_CALL'] = igblast_out_v
        #drop low scoring V genes that don't have bowtie idenitity
        igblast_out_hsv = igblast_out.dropna(subset = ['V_CALL'])
        #drop low quality J calls
        igblast_out_hs = igblast_out_hsv[(igblast_out_hsv['J_SCORE'] > 35)]
    else:
        #Filter out low quality scores
        igblast_out_hs = igblast_out[(igblast_out['V_SCORE'] > 50) & (igblast_out['J_SCORE'] > 35)]


    #Filter mutiple different V calls
    # igblast_out_hs['V_CALL'].str.split('(,|\*)').str[0]
    # igblast_out_hs['V_CALL'].str.split('(,|\*)').str[4]

    pd.options.mode.chained_assignment = None #Turn of warning http://stackoverflow.com/questions/20625582/how-to-deal-with-settingwithcopywarning-in-pandas
    igblast_out_hs.loc[:,'V_CALL'] = igblast_out_hs['V_CALL'].apply(ambigious_calls, args=(retain_nam,))
    igblast_out_hs.loc[:,'J_CALL'] = igblast_out_hs['J_CALL'].apply(ambigious_calls, args=(retain_nam,))

    # Counter(igblast_out_hs['V_CALL'])
    # Counter(igblast_out_hs['J_CALL'])

    #Remove rows with None call
    igblast_out_hs = igblast_out_hs[igblast_out_hs['V_CALL'].notnull()]
    igblast_out_hs = igblast_out_hs[igblast_out_hs['J_CALL'].notnull()]

    # len(igblast_out_hs)

    #Seperate into functional non-functional
    # functional = igblast_out_hs[igblast_out_hs['FUNCTIONAL'] == 'T']
    # non_functional = igblast_out_hs[igblast_out_hs['FUNCTIONAL'] == 'F']

    # len(functional)
    # len(non_functional)

    #If not cdr3 present drop record
    igblast_out_hs_cln = igblast_out_hs.dropna(subset = ['CDR3_IMGT'])
    # functional_cln = functional.dropna(subset = ['CDR3_IMGT'])
    # non_functional_cln = non_functional.dropna(subset = ['CDR3_IMGT'])

    #only output a minimal table
    if minimal:
        igblast_out_hs_cln = igblast_out_hs_cln[['SEQUENCE_ID', 'SEQUENCE_INPUT', 'V_CALL', 'D_CALL', 'J_CALL', 'CDR3_IMGT', 'file_ID']]
    # return (functional_cln, non_functional_cln)
    return igblast_out_hs_cln

# functional_cln['SEQUENCE_ID']['HWI-M02293:218:000000000-AKGG1:1:1101:8139:9569_J3_CTGCTCCT_AGCGGA_5']
# functional_cln.SEQUENCE_ID[functional_cln.SEQUENCE_ID == 'HWI-M02293:218:000000000-AKGG1:1:1101:8139:9569_J3_CTGCTCCT_AGCGGA_5'].index.tolist()[0]

def make_bundle(pd_data_frame, only_v=False):
    '''Make dictionary of V-CRD3-J (bundle)'''
    clonotype_dict = defaultdict(lambda: defaultdict(dict))

    assert 'V_CALL' in pd_data_frame.columns and 'J_CALL' in pd_data_frame.columns \
    and 'CDR3_IMGT' in pd_data_frame.columns and 'SEQUENCE_ID' in pd_data_frame.columns \
    and 'SEQUENCE_INPUT' in pd_data_frame.columns, 'Requried columns not in data frame'
    # print(pd_data_frame)
    for line in pd_data_frame.index:
        if only_v:
            v_j = pd_data_frame['V_CALL'][line]
        else:
            v_j = pd_data_frame['V_CALL'][line] + '_' + pd_data_frame['J_CALL'][line]

        cdr3 = pd_data_frame['CDR3_IMGT'][line]
        try:
            clonotype_dict[v_j][cdr3]['qname'].append(pd_data_frame['SEQUENCE_ID'][line])
            clonotype_dict[v_j][cdr3]['read'].update([pd_data_frame['SEQUENCE_INPUT'][line]])
            clonotype_dict[v_j][cdr3]['count'] += 1
        except KeyError:
            clonotype_dict[v_j][cdr3]['qname'] = [pd_data_frame['SEQUENCE_ID'][line]]
            clonotype_dict[v_j][cdr3]['read'] = Counter([pd_data_frame['SEQUENCE_INPUT'][line]])
            clonotype_dict[v_j][cdr3]['count'] = 1

    return clonotype_dict

# make_bundle(functional_cln)

# Assemble clones allowing 1 mismatch in CDR3 (how many times does the same recombination take place?)
def assemble_colonotype(pd_data_frame, bundles, threshold):
    pd_data_frame['clonotype'] = ''

    bundle_count = 0
    for bundle in bundles.values():
        umis = bundle.keys()

        # len_umis = [len(x) for x in umis]
        # assert max(len_umis) == min(len_umis), ('not all umis are the same length(!):  %d - %d' % (min(len_umis), max(len_umis)))

        counts = {umi: bundle[umi]['count'] for umi in umis} #If two UMI's mergered, count will be taken only from the larger UMI group
        # print('Getting directional adjacency list')
        adj_list = lev_adj_list_adjacency(list(umis), counts, threshold)

        # print(len(adj_list), sorted(adj_list)[1:5])
        # print('Getting connected components')
        clusters = deduplicate.get_connected_components_adjacency(umis, adj_list, counts)
        # print(clusters)
        #Assign clonotypes
        cluster_count = 0
        for cluster in clusters:
            for cdr3 in cluster:
                #write into original pandas table
                for qname in bundle[cdr3]['qname']:
                    row_loc = pd_data_frame.SEQUENCE_ID[pd_data_frame.SEQUENCE_ID == qname].index.tolist()[0]
                    pd_data_frame['clonotype'][row_loc] = str(bundle_count) + '_' + str(cluster_count)
            cluster_count += 1
        bundle_count += 1

    return pd_data_frame

# assemble_colonotype(functional_cln, clonotype_dict)

#Write out pandas table
def write_out(pd_data_frame, out):
    pd_data_frame.to_csv(out, sep='\t')

# pd_data_frame.to_csv('/media/chovanec/My_Passport/Dan_VDJ-seq_cycles_new/J_merged_1c_Deduplicated_test/J_merged_1c_dedup.0.fasta_db-pass_assembled.tab', sep='\t')


#Assemble clones comparing entire sequence (allowing x mismatches?) msa? different V starts makes this complicated



#parser
def parse_args():
    parser = argparse.ArgumentParser(description='BabrahamLinkON Assemble Clones')

    parser.add_argument('--tab_file', dest='in_file', type=str, required=True, help='Input tab file from changeo IgBlast MakeDb (or file wildcard)')
    parser.add_argument('--plot', action='store_true', help='Plot V and J scores with cutoff')
    parser.add_argument('--out', dest='out_dir', type=str, help='Output directory, default: creates Deduplicated in main directory')
    parser.add_argument('--threshold', dest='thres', type=int, default=1, help='Number of differences allowed between CDR3 sequences')
    parser.add_argument('--only_v', action='store_true', help='Use only V idenity for clone assembly')
    parser.add_argument('--full_name', action='store_true', help='Retain full name of first V and J genes')
    parser.add_argument('--minimal', action='store_true', help='Work with and output only a minimal table')
    parser.add_argument('--short', action='store_true', help='Analysing short sequences')

    opts = parser.parse_args()

    return opts


def main():

    #argparse
    opts = parse_args()

    files = glob.glob(opts.in_file)
    # functional_cln, non_functional_cln = read_changeo_out(opts.in_file, plot=opts.plot)
    full_path = os.path.abspath(files[0])
    prefix = os.path.basename(files[0]).split('.')[0]

    if opts.out_dir == None:
        out_dir = os.path.dirname(full_path)
    else:
        out_dir = opts.out_dir

    # print('out_dir', out_dir, opts.out_dir)
    # files = glob.glob('/media/chovanec/My_Passport/Dan_VDJ-seq_cycles_Jan17/results/lane5279*')
    igblast_cln = read_changeo_out(files, out_dir, prefix, plot=opts.plot, retain_nam=opts.full_name, minimal=opts.minimal, short=opts.short)
    # igblast_cln = read_changeo_out(files, out_dir, prefix)

    clonotype_dict = make_bundle(igblast_cln, only_v=opts.only_v)

    ig_blast_asm = assemble_colonotype(igblast_cln, clonotype_dict, opts.thres)


    if opts.minimal:
        write_out(ig_blast_asm, out_dir + '/' + prefix + '_assembled_clones_min.tab')
    else:
        write_out(ig_blast_asm, out_dir + '/' + prefix + '_assembled_clones.tab')


if __name__ == "__main__":
    main()
