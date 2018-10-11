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

from collections import defaultdict, Counter
#Do not hard code matplotlib backend, use export MPLBACKEND=pdf instead if running on headless node
# import matplotlib
# matplotlib.use('pdf')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import pandas as pd
import numpy as np
import Levenshtein
from babrahamlinkon import general, igblast_wrapper, umi_correction
import os
import argparse
import glob
import tempfile
import shutil
import logging
import pyximport
from tqdm import tqdm
import json
from babrahamlinkon._dedup_umi import edit_distance
from babrahamlinkon.version import __version__
# %matplotlib inline
# %config InlineBackend.figure_format = 'svg'


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


def adj_list_adjacency(umis, counts, threshold=1):
    ''' identify all umis within the hamming distance threshold'''

    #should be 50% faster
    adj_list = defaultdict(list)
    for i,umi in enumerate(umis):
        a1 = adj_list[umi]
        for j in range(i+1,len(umis)):
            umi2 = umis[j] #dict_keys object doesn't support indexing
            # if Levenshtein.distance(umi, umi2) <= threshold:
            if len(umi) == len(umi2):
                if edit_distance(umi.encode('utf-8'), umi2.encode('utf-8')) <= threshold:
                    adj_list[umi].append(umi2)

    return adj_list



# def change_v_call(row):
#     count = 0
#     if row['V_SCORE'] <= 50:
#         count += 1
#         v_gene = row['SEQUENCE_ID'].split('_')[-3].upper() #would normally be the barcode!
#         #replace v gene call with bowtie alignment call
#
#         if len(v_gene) < 1:
#             return np.NaN
#     else:
#         v_gene = row['V_CALL']
#     return v_gene



def v_identity_igblast(V_end, J_end_fasta, custom_ref, thread_num, spe, aux, dj):
    '''
    :param V_fastq: original fastq
    :param fasta: deduplicate.py fasta output (J end)
    :param thread_num: number of threads to use
    :param spe: species (mmu, hsa)
    '''
    #retain full name
    subset_reads = defaultdict()
    for name, seq in general.fasta_iter(J_end_fasta):
        subset_reads[name.split('_')[0]] = name

    #make tmp directory with igblast run files
    tmp_dir = tempfile.mkdtemp()
    tmp_fmt = os.path.join(tmp_dir, "igblast.fmt7")

    if V_end.endswith('fastq'):
        #fastq to fasta
        fasta = ''
        reads_fasta = 0
        with general.file_open(V_end) as fq:
            for item in general.fastq_parse(fq):
                title = item[0]
                seq = item[1]

                try:
                    fasta += '>' + subset_reads[title.split(' ')[0][1:]] + '\n' + seq + '\n'
                    reads_fasta += 1
                except:
                    pass

        #need to write out fasta for changeo (find alternative which accepts stdin?)
        with open(tmp_dir + '/igblast.fasta', 'w') as fa_out:
            fa_out.write(fasta)

    elif V_end.endswith('fasta'):
        shutil.copy(V_end, tmp_dir + '/igblast.fasta')


    igblast_wrapper.run_igblast(tmp_dir + '/igblast.fasta', tmp_fmt, 10000, spe, thread_num, custom_ref, dj, aux_file=aux, additional_flags=['-num_alignments_V', '1'])
    igblast_wrapper.parse_igblast(tmp_fmt, tmp_dir + '/igblast.fasta', spe, custom_ref, dj)
    #make v_identity dict key=qname value=idenity

    #need to find the output of changeo
    tmp_tab = glob.glob(tmp_dir + '/*.tab')

    df = pd.read_table(tmp_tab[0], header=0)
    sub_df = df[['SEQUENCE_ID', 'V_CALL', 'V_SCORE']].copy()
    sub_df.rename(columns={'V_CALL': 'V_CALL_VEND', 'V_SCORE': 'V_SCORE_VEND'}, inplace=True)
    #Delete temporary files
    shutil.rmtree(tmp_dir)

    return sub_df


#
# v_identity_igblast('/media/chovanec/My_Passport/Old_vdj_seq_data/Deduplicated_all/A1BC/lane1_A2BC_TGACCA_L001_R1_val_1.fq.gz',
# '/media/chovanec/My_Passport/Old_vdj_seq_data/Deduplicated_all/A1BC/lane1_A2BC_TGACCA_L001_R1_val_1_dedup.fasta', 7, 'mmu', aux=None)
#
#
# def add_v_call(row):
#     V_END = row['SEQUENCE_ID'].split('_')[-3]
#     V_END_IDENTITY, V_END_SCORE = V_END.split('.')
#     return [V_END_IDENTITY, V_END_SCORE]

def read_changeo_out(tab_file, out, prefix, fasta, v_fastq=None, plot=False, retain_nam=False,
                     minimal=False, short=False, thread_num=1, spe='mmu', aux=None, custom_ref=False,
                     j_cutoff=35, v_cutoff=50, dj=False):
    '''
    :param retain_nam: keep full name of V and J calls
    '''

    igblast_out = pd.DataFrame()
    df_list = []
    for f in tab_file:
        df = pd.read_table(f, header=0)
        #add column with file identity
        # id_label = os.path.basename(f).split('.')[0]
        # df['file_ID'] = pd.Series(np.repeat(id_label, len(df.index)))
        df_list.append(df)
    igblast_out = pd.concat(df_list)
    igblast_out.reset_index(drop=True, inplace=True)
    print('Called reads', len(igblast_out.index))

    #Seperate out DJ calls into seperate dataframe
    igblast_dj = igblast_out[igblast_out['V_CALL'].str.contains('IGHVD', na=False)]
    # igblast_dj.reset_index(drop=True, inplace=True)
    igblast_out = igblast_out[~igblast_out['V_CALL'].str.contains('IGHVD', na=False)]
    # igblast_out.reset_index(drop=True, inplace=True)

    if short:
        # HWI-1KL136:214:D1MR5ACXX:5:1103:17395:138047_J4_Ighv13-2_GTGTCTAC_11

        # igblast_out_v = igblast_out.apply(change_v_call, axis=1)
        # igblast_out['V_CALL'] = igblast_out_v
        # v_iden, v_score = igblast_out.apply(add_v_call, axis=1)
        # igblast_out['V_END_CALL'] = v_iden
        # igblast_out['V_END_SCORE'] = v_score
        if v_fastq == None:
            raise Exception('Short option requires the V end fastq file')

        #run igblast on the v end
        v_end_calls = v_identity_igblast(v_fastq, fasta, custom_ref, thread_num, spe, aux, dj)
        # logging.info('igblast_out', len(igblast_out), 'v_end_calls', len(v_end_calls))
        #merge data fragments
        igblast_out_m = pd.merge(igblast_out, v_end_calls, how='left', on=['SEQUENCE_ID'])


    if plot:
        with PdfPages(out + '/' + prefix + '_score_plots.pdf') as pdf_out:
        #Plot V and J scores
            labels, values = zip(*Counter(igblast_out['V_SCORE']).items())
            non_dedup_values = tuple(l*v for l, v in zip(labels, values))

            plt.figure()
            plt.bar(labels, non_dedup_values)
            plt.title('V score')
            my_plot = plt.axvline(v_cutoff, linestyle='dashed', linewidth=2).get_figure()
            pdf_out.savefig(my_plot)

            labels, values = zip(*Counter(igblast_out['J_SCORE']).items())
            non_dedup_values = tuple(l*v for l, v in zip(labels, values))

            plt.figure()
            plt.bar(labels, non_dedup_values)
            plt.title('J score')
            my_plot = plt.axvline(j_cutoff, linestyle='dashed', linewidth=2).get_figure()
            pdf_out.savefig(my_plot)

            if short:
                labels, values = zip(*Counter(v_end_calls['V_SCORE_VEND']).items())
                non_dedup_values = tuple(l*v for l, v in zip(labels, values))

                plt.figure()
                plt.bar(labels, non_dedup_values)
                plt.title('V score v end')
                my_plot = plt.axvline(v_cutoff, linestyle='dashed', linewidth=2).get_figure()
                pdf_out.savefig(my_plot)


    if short:

        #drop low scoring V genes that don't have idenitity
        igblast_out_na = igblast_out_m.dropna(subset = ['V_CALL', 'V_CALL_VEND'], how='all')
        #drop low quality J calls
        igblast_out_hs = igblast_out_na[(igblast_out_na['J_SCORE'] > j_cutoff) & ((igblast_out_na['V_SCORE'] > v_cutoff) | (igblast_out_na['V_SCORE_VEND'] > v_cutoff))]
        low_score = len(igblast_out_m.index) - len(igblast_out_hs.index)
        #how many filtered out?
        logging.info('Low V and J score:' + str(low_score))
        print('Low V and J score:', low_score)

    else:
        #Filter out low quality scores
        igblast_out_hs = igblast_out[(igblast_out['V_SCORE'] > v_cutoff) & (igblast_out['J_SCORE'] > j_cutoff)]
        low_score = len(igblast_out.index) - len(igblast_out_hs.index)
        #how many filtered out?
        logging.info('Low V and J score:' + str(low_score))
        print('Low V and J score:', low_score)

    if dj:
        #DJ filtering
        #Drop DJ without a J calls, suggests misidentification of D
        igblast_dj_na = igblast_dj.dropna(subset = ['J_CALL'])
        if short:
            #merge V end calls with J end calls
            igblast_out_na_sub = igblast_out_na[['SEQUENCE_ID', 'V_CALL_VEND', 'V_SCORE_VEND']]
            igblast_dj_na_all = igblast_dj_na.set_index('SEQUENCE_ID').join(igblast_out_na_sub.set_index('SEQUENCE_ID'))
            #keep only those with high scores
            igblast_dj_out = igblast_dj_na_all[(igblast_dj_na_all['V_SCORE'] > v_cutoff) | (igblast_dj_na_all['V_SCORE_VEND'] > v_cutoff)]
        else:
            igblast_dj_out = igblast_dj_na[(igblast_dj_na['V_SCORE'] > v_cutoff)]

        dj_filt = len(igblast_dj.index)-len(igblast_dj_out.index)
        logging.info('Number of DJ reads filtered:' + str(dj_filt))
        print('Number of DJ reads filtered:', dj_filt)
    else:
        igblast_dj_out = []


    #Filter mutiple different V calls
    # igblast_out_hs['V_CALL'].str.split('(,|\*)').str[0]
    # igblast_out_hs['V_CALL'].str.split('(,|\*)').str[4]

    # pd.options.mode.chained_assignment = None #Turn of warning http://stackoverflow.com/questions/20625582/how-to-deal-with-settingwithcopywarning-in-pandas
    # igblast_out_hs.loc[:,'V_CALL'] = igblast_out_hs['V_CALL'].apply(ambigious_calls, args=(retain_nam,))
    # igblast_out_hs.loc[:,'J_CALL'] = igblast_out_hs['J_CALL'].apply(ambigious_calls, args=(retain_nam,))
    #
    # # Counter(igblast_out_hs['V_CALL'])
    # # Counter(igblast_out_hs['J_CALL'])
    #
    # #Remove rows with None call
    # igblast_out_hs = igblast_out_hs[igblast_out_hs['V_CALL'].notnull()]
    # igblast_out_hs = igblast_out_hs[igblast_out_hs['J_CALL'].notnull()]

    # len(igblast_out_hs)

    #Seperate into functional non-functional
    # functional = igblast_out_hs[igblast_out_hs['FUNCTIONAL'] == 'T']
    # non_functional = igblast_out_hs[igblast_out_hs['FUNCTIONAL'] == 'F']

    # len(functional)
    # len(non_functional)

    #If not cdr3 present drop record
    igblast_out_hs_cln = igblast_out_hs.dropna(subset = ['CDR3_IGBLAST_NT'])
    #how many have no CDR3?
    no_cdr3 = len(igblast_out_hs.index)-len(igblast_out_hs_cln.index)
    logging.info('Number of reads without CDR3:' + str(no_cdr3))
    print('Number of reads without CDR3:', no_cdr3)
    # functional_cln = functional.dropna(subset = ['CDR3_IMGT'])
    # non_functional_cln = non_functional.dropna(subset = ['CDR3_IMGT'])

    #only output a minimal table
    if minimal:
        igblast_out_hs_cln = igblast_out_hs_cln[['SEQUENCE_ID', 'SEQUENCE_INPUT', 'V_CALL', 'D_CALL', 'J_CALL', 'CDR3_IGBLAST_NT', 'file_ID']]
    # return (functional_cln, non_functional_cln)
    return (igblast_out_hs_cln, igblast_dj_out)

# functional_cln['SEQUENCE_ID']['HWI-M02293:218:000000000-AKGG1:1:1101:8139:9569_J3_CTGCTCCT_AGCGGA_5']
# functional_cln.SEQUENCE_ID[functional_cln.SEQUENCE_ID == 'HWI-M02293:218:000000000-AKGG1:1:1101:8139:9569_J3_CTGCTCCT_AGCGGA_5'].index.tolist()[0]

def make_bundle(pd_data_frame, only_v=False):
    '''Make dictionary of V-CRD3-J (bundle)
    '''
    clonotype_dict = defaultdict(lambda: defaultdict(dict))

    assert 'V_CALL' in pd_data_frame.columns and 'J_CALL' in pd_data_frame.columns \
    and 'CDR3_IGBLAST_NT' in pd_data_frame.columns and 'SEQUENCE_ID' in pd_data_frame.columns \
    and 'SEQUENCE_INPUT' in pd_data_frame.columns, 'Requried columns not in data frame'
    # print(pd_data_frame)
    # print(pd_data_frame['CDR3_IGBLAST_NT'])
    for line in pd_data_frame.index:
        if only_v:
            v_j = pd_data_frame['V_CALL'][line]
        else:
            v_j = pd_data_frame['V_CALL'][line] + '_' + pd_data_frame['J_CALL'][line]

        #use IgBlast CDR3
        cdr3 = pd_data_frame['CDR3_IGBLAST_NT'][line]

        #seperate group based on sequence length as well as V and J
        v_j_len = v_j + '_' + str(len(cdr3))

        try:
            clonotype_dict[v_j_len][cdr3]['qname'].append(pd_data_frame['SEQUENCE_ID'][line])
            clonotype_dict[v_j_len][cdr3]['read'].update([pd_data_frame['SEQUENCE_INPUT'][line]])
            clonotype_dict[v_j_len][cdr3]['count'] += 1
        except KeyError:
            clonotype_dict[v_j_len][cdr3]['qname'] = [pd_data_frame['SEQUENCE_ID'][line]]
            clonotype_dict[v_j_len][cdr3]['read'] = Counter([pd_data_frame['SEQUENCE_INPUT'][line]])
            clonotype_dict[v_j_len][cdr3]['count'] = 1

    return clonotype_dict

# make_bundle(functional_cln)

# Assemble clones allowing 1 mismatch in CDR3 (how many times does the same recombination take place?)
def assemble_colonotype(pd_data_frame, bundles, threshold):
    '''
    '''

    pd_data_frame['clonotype'] = ''

    bundle_count = 0
    for bundle in tqdm(bundles.values()):
        umis = bundle.keys()

        # len_umis = [len(x) for x in umis]
        # assert max(len_umis) == min(len_umis), ('not all umis are the same length(!):  %d - %d' % (min(len_umis), max(len_umis)))

        counts = {umi: bundle[umi]['count'] for umi in umis}

        adj_list = adj_list_adjacency(list(umis), counts, threshold)
        clusters = umi_correction.get_connected_components_adjacency(umis, adj_list, counts)

        #Assign clonotypes
        cluster_count = 0
        for cluster in clusters:
            for cdr3 in cluster:
                #write into original pandas table
                for qname in bundle[cdr3]['qname']:
                    row_loc = pd_data_frame.SEQUENCE_ID[pd_data_frame.SEQUENCE_ID == qname].index.tolist()[0]
                    # pd_data_frame.set_value(row_loc, 'clonotype', str(bundle_count) + '_' + str(cluster_count))
                    pd_data_frame.at[row_loc, 'clonotype'] = str(bundle_count) + '_' + str(cluster_count)
                    # pd_data_frame.iloc[:,('clonotype', row_loc)] = str(bundle_count) + '_' + str(cluster_count)
            cluster_count += 1
        bundle_count += 1

    return pd_data_frame

# assemble_colonotype(functional_cln, clonotype_dict)

#Write out pandas table
def write_out(pd_data_frame, out):
    pd_data_frame.to_csv(out, sep='\t')

# pd_data_frame.to_csv('/media/chovanec/My_Passport/Dan_VDJ-seq_cycles_new/J_merged_1c_Deduplicated_test/J_merged_1c_dedup.0.fasta_db-pass_assembled.tab', sep='\t')


def add_assembled_column(df, json_path):
    #add @ to read name to a new tmp column
    df['SEQUENCE_ID_tmp'] = '@' + df['SEQUENCE_ID'].astype(str)

    with open(json_path, 'r') as in_json:
        assembled_dict = json.load(in_json)

    #split str and take only first part (read name) and check if it in assembled on unassembled list
    df.loc[df['SEQUENCE_ID_tmp'].str.split('_').str.get(0).isin(assembled_dict['unassembled']), 'Status'] = 'Unassembled'
    df.loc[df['SEQUENCE_ID_tmp'].str.split('_').str.get(0).isin(assembled_dict['assembled']), 'Status'] = 'Assembled'
    #delete tmp column
    del df['SEQUENCE_ID_tmp']

    return(df)



#parser
def parse_args():
    parser = argparse.ArgumentParser(description='BabrahamLinkON Assemble Clones')

    parser.add_argument('--version', action='version', version='%(prog)s ' + __version__)

    sub = parser.add_subparsers(dest='action', description='Choose pipeline')

    sp1 = sub.add_parser('umi')
    sp2 = sub.add_parser('short')
    sp3 = sub.add_parser('assemble_only')
    # sp3 = sub.add_parser('short_anchor')
    # sp4 = sub.add_parser('no_anchor')

    # parser.add_argument('--tab_file', dest='in_file', type=str, required=True, help='Input tab file from changeo IgBlast MakeDb (or file wildcard)')
    for sp in [sp1, sp2]:

        sp.add_argument('-fa', '--fasta', dest='fasta', type=str, help='Input fasta file from deduplication.py')

        sp.add_argument('--plot', action='store_true', help='Plot V and J scores with cutoff')
        sp.add_argument('--full_name', action='store_true', help='Retain full name of first V and J genes')
        sp.add_argument('--minimal', action='store_true', help='Work with and output only a minimal table')

        sp.add_argument('--threads', dest='nthreads', default=1, type=int, help='Number of threads to use [1]')
        sp.add_argument('--species', dest='species', default='mmu', type=str, help='Which species (mmu hsa mmuk) [mmu]')
        sp.add_argument('--aux', dest='aux', type=str, default=None, help='aux file for igblast')
        sp.add_argument('--custom_ref', dest='custom_ref', action='store_true', help='Use AEC custom reference for igblast')
        sp.add_argument('--v_cutoff', dest='v_cutoff', default=50, type=int, help='IgBlast V_SCORE cutoff [>50]')
        sp.add_argument('--j_cutoff', dest='j_cutoff', default=35, type=int, help='IgBlast J_SCORE cutoff [>35]')

        sp.add_argument('--skip_assembly', dest='skip_assembly', action='store_true', help='Do not perform clone assembly into clonotypes')
        sp.add_argument('--call_dj', dest='call_dj', action='store_true', help='Call DJ recombination, else only VDJ will be called')

    sp2.add_argument('--json_path', dest='json_path', type=str, help='JSON file produced during precleaning that contains identity of assembled and unassembled reads')
    sp3.add_argument('-tsv', '--tsv_files', dest='tsv_files', type=str, metavar='x.tsv', nargs='+', help='Input tsv files from previous run of assemble_clones.py')

    for sp in [sp1, sp2, sp3]:
        sp.add_argument('--threshold', dest='thres', type=int, default=1, help='Number of differences allowed between CDR3 sequences [1]')
        sp.add_argument('--only_v', action='store_true', help='Use only V idenity and CDR3 for clone assembly')
        sp.add_argument('--out', dest='out_dir', type=str, help='Output directory, default: creates Deduplicated in main directory')

    sp2.add_argument('-v', '--v_end', dest='v_fastq', type=str, help='V end fastq or fasta file')
    # parser.add_argument('--short', action='store_true', help='Analysing short sequences')


    sp1.set_defaults(short=False, v_fastq=None, assemble_only=False)
    sp2.set_defaults(short=True, assemble_only=False)
    sp3.set_defaults(short=False, v_fastq=None, assemble_only=True, call_dj=False, skip_assembly=False, minimal=False)

    opts = parser.parse_args()

    return opts


def main():

    #argparse
    opts = parse_args()

    if opts.short and opts.json_path is None:
        mark_reads = False
        print('No JSON input provided, short reads will not be marked as assembled or unassembled')
    else:
        mark_reads = True

    if opts.assemble_only:
        #add file name column and then merge all input files
        df_list = []
        prefixes = set()
        out_dirs = set()
        for item in opts.tsv_files:
            full_path = os.path.abspath(item) #put in same folder as deduplication
            file_name = os.path.basename(item).split('.')[0]
            prefixes.add(file_name)
            ac_df = pd.read_csv(full_path, sep='\t', index_col=0)
            ac_df = ac_df.assign(file_name=pd.Series([file_name]*len(ac_df)).values)
            df_list.append(ac_df)

            if opts.out_dir == None:
                out_dirs.add(os.path.dirname(full_path))

        igblast_cln = pd.concat(df_list, ignore_index=True)

        #output file prefix, if more than 1 file supplied use merged output name
        if len(opts.tsv_files) > 1:
            prefix = 'Merged_' + str(len(prefixes))
        else:
            prefix = list(prefixes)[0]

        if opts.out_dir == None:
            out_dir = list(out_dirs)[0] #if there are multiple paths just pick first
        else:
            out_dir = opts.out_dir
            #create out dir if it doesn't exist
            try:
                os.mkdir(out_dir)
            except FileExistsError:
                pass

    else:
        # files = glob.glob(opts.in_file)
        # functional_cln, non_functional_cln = read_changeo_out(opts.in_file, plot=opts.plot)
        full_path = os.path.abspath(opts.fasta) #put in same folder as deduplication
        prefix = os.path.basename(opts.fasta).split('.')[0]

        if opts.out_dir == None:
            out_dir = os.path.dirname(full_path)
        else:
            out_dir = opts.out_dir
            #create out dir if it doesn't exist
            try:
                os.mkdir(out_dir)
            except FileExistsError:
                pass


        logging.basicConfig(level=logging.DEBUG, filename=out_dir +'/' + prefix + '_assembled_clones.log', filemode='a+',
                            format='%(asctime)-15s %(levelname)-8s %(message)s')


        #make tmp directory with igblast run files
        tmp_dir = tempfile.mkdtemp()
        tmp_fmt = os.path.join(tmp_dir, "igblast.fmt7")


        igblast_wrapper.run_igblast(opts.fasta, tmp_fmt, 10000, opts.species, opts.nthreads, opts.custom_ref, dj=opts.call_dj, aux_file=opts.aux)
        igblast_wrapper.parse_igblast(tmp_fmt, opts.fasta, opts.species, opts.custom_ref, opts.call_dj)

        #need to find the output of changeo
        tmp_tab = glob.glob(tmp_dir + '/*.tab')

        igblast_cln, igblast_dj = read_changeo_out(tmp_tab, out_dir, prefix, opts.fasta, v_fastq=opts.v_fastq, plot=opts.plot,
                                       retain_nam=opts.full_name, minimal=opts.minimal, short=opts.short, thread_num=opts.nthreads,
                                       spe=opts.species, aux=opts.aux, custom_ref=opts.custom_ref, j_cutoff=opts.j_cutoff, v_cutoff=opts.v_cutoff,
                                       dj=opts.call_dj)
        # igblast_cln = read_changeo_out(files, out_dir, prefix)
        out_reads_count = len(igblast_cln.index)
        logging.info('Out reads:' + str(out_reads_count))
        print('Out reads:', out_reads_count)

    if opts.call_dj:
        dj_count = len(igblast_dj.index)
        logging.info('DJ reads:' + str(dj_count))
        print('DJ reads:', dj_count)

        #for short reads add an additional column marking reads that have been assembled
        if mark_reads:
            igblast_dj = add_assembled_column(igblast_dj, opts.json_path)

        write_out(igblast_dj, out_dir + '/' + prefix + '_annotated_dj.tsv')


    if opts.skip_assembly:
        if mark_reads:
            igblast_cln = add_assembled_column(igblast_cln, opts.json_path)

        if opts.minimal:
            write_out(igblast_cln, out_dir + '/' + prefix + '_annotated_clones_min.tsv')
        else:
            write_out(igblast_cln, out_dir + '/' + prefix + '_annotated_clones.tsv')

    else:
        print('Assembling clones')
        clonotype_dict = make_bundle(igblast_cln, only_v=opts.only_v)

        ig_blast_asm = assemble_colonotype(igblast_cln, clonotype_dict, opts.thres)

        if mark_reads:
            ig_blast_asm = add_assembled_column(ig_blast_asm, opts.json_path)

        if opts.minimal:
            write_out(ig_blast_asm, out_dir + '/' + prefix + '_assembled_clones_min.tsv')
        else:
            write_out(ig_blast_asm, out_dir + '/' + prefix + '_assembled_clones.tsv')


    #delete tmp dir
    if not opts.assemble_only:
        shutil.rmtree(tmp_dir)

if __name__ == "__main__":
    main()
