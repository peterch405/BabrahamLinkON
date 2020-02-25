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

import argparse
from collections import defaultdict
import glob
from babrahamlinkon import general
from babrahamlinkon.version import __version__
import os


def parse_args():
    parser = argparse.ArgumentParser(description='BabrahamLinkON Split germline')

    parser.add_argument('--version', action='version', version='%(prog)s ' + __version__)

    parser.add_argument('--input_dir', dest='input', type=str, required=True, help='Input directory where results from preclean and deduplicate are')
    
    opts = parser.parse_args()

    return opts


def main():

    #argparse
    opts = parse_args()

    if os.path.isfile(opts.input):
        raise Exception('Input need to be a directory')
    elif os.path.isdir(opts.input):
        in_dir = opts.input

    
    deduplicated_fa = glob.glob(in_dir + '/*_Deduplicated/*dedup.fasta')[0]
    germline_fq = glob.glob(in_dir + '/*_preclean/*_germline_J')[0]
    print('Input Deduplicated fasta:', deduplicated_fa)
    print('Input Germline fastq:', germline_fq)

    germ_out_fa = os.path.dirname(deduplicated_fa) + '/germline.fasta'
    vdj_out_fa = os.path.dirname(deduplicated_fa) + '/vdj.fasta'

    separate_germline(deduplicated_fa, germline_fq, germ_out_fa, vdj_out_fa)



def separate_germline(split_fa, germline_fq, germ_out_fa, vdj_out_fq):
    '''Split germline reads from fastq based on read names
    '''

    germ_reads = set()
    germ_count = 0
    vdj_count = 0
    count = 0
    with general.file_open(germline_fq) as germ_fq:
        for name, seq, thrd, qual in general.fastq_parse(germ_fq):
            germ_reads.add(name.split(' ')[0].split('_')[0][1:])
    
    with open(germ_out_fa, 'w') as germ_out, \
    open(vdj_out_fq, 'w') as vdj_out:
        for name, seq in general.fasta_iter(split_fa):
            count += 1
            if name.split(' ')[0].split('_')[0] in germ_reads:
                germ_out.write('>' + name + '\n' + seq + '\n')
                germ_count += 1
            else:
                vdj_out.write('>' + name + '\n' + seq + '\n')
                vdj_count += 1

    print('Germline reads written:', germ_count)
    print('VDJ reads written:', vdj_count)
    print('Total reads processed:', count)




if __name__ == "__main__":
    main()

