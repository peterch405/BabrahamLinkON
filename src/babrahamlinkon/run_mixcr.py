#!/usr/bin/env python3

import argparse
from babrahamlinkon import mixcr_wrapper

def parse_args():
    parser = argparse.ArgumentParser(description='BabrahamLinkON Deduplicate')

    parser.add_argument('-v', '--V_dedup', dest='input_V', type=str, required=True, help='Input deduplicated fastq file(s) with V end sequences')
    parser.add_argument('-j', '--J_dedup', dest='input_J', type=str, required=True, help='Input deduplicated fastq file(s) with J end sequences')
    parser.add_argument('--out', dest='out_dir', type=str, help='Output directory, default: creates MiXCR in main directory')
    parser.add_argument('--cores', dest='nthreads', default=1, type=int, help='Number of cores to use, default: 1')
    parser.add_argument('--species', dest='species', default='mmu', type=str, help='Which species (mmu, hsa), default: mmu')

    opts = parser.parse_args()

    return opts


def main():

    #argparse
    opts = parse_args()

    mixcr_out = mixcr_wrapper.run_mixcr(opts.input_V, opts.input_J, nthreads=opts.nthreads, species=opts.species, out_dir=opts.out_dir)

    unprod_out = mixcr_wrapper.extract_unproductive(mixcr_out[0], mixcr_out[1]) #outputs same dir as run_mixcr

    mixcr_wrapper.plot_results(mixcr_out[0], mixcr_out[1], unprod_out)

    # mixcr_out = run_mixcr('/media/chovanec/My_Passport/Sync/BabrahamLinkON/tests/Deduplicated/lane3_WT_FrBC_1_GCCAAT_L003_R1_val_1_40k_V.fastq',
    #                       '/media/chovanec/My_Passport/Sync/BabrahamLinkON/tests/Deduplicated/lane3_WT_FrBC_1_GCCAAT_L003_R1_val_1_40k_J.fastq',
    #                       nthreads=8)
    #
    # unprod_out = extract_unproductive(mixcr_out[0], mixcr_out[1])
    #
    # plot_results(mixcr_out[0], mixcr_out[1], unprod_out)


if __name__ == "__main__":
    main()
