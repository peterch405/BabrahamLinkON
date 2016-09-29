#!/usr/bin/env python3

import os
import pysam
import glob
from pathlib import Path
import subprocess
import shutil
import re
from collections
import fileinput
import shlex
import argparse
from babrahamlinkon import general

import pyximport



class deduplicate_umi_tools: #if using UMI_tools
    '''Deduplicate using umi_tools
    :param file_directory: location of output from preclean
    '''

    def __init__(self, file_directory, samples=['']):
        self.file_directory = file_directory
        self.samples = samples
        self.V_files = []
        self.J_files = []
        self.tmp_dir = ''
        self.out_dir = ''
        self.out_dir_prefix = ''
        self.extract_original_run = False
        self.j_prefix = ''
        self.v_prefix = ''


    def initial_setup(self, out_dir=None, barcode_1='GACTCGT', barcode_2='CTGCTCCT', include_all_v=False):
        '''Create file lists and all required output directories
        :param out_dir: output directory
        :param prefix: a list of prefix(s) of files to be analysed (if multiple samples being run together)
        :param barcode_1:
        :param barcode_2:
        :param include_all_v: if identity of J doesn't matter in downstream analysis (only freq of V usage of interest)
        '''

        if len(self.V_files) > 0 or len(self.J_files) > 0:
            raise Exception('Looks like initial_setup() already run!')


        # if prefix == None:
        #     print('If preclean contains more than one sample, you need to specify a prefix!')
        # else:
        #     self.sample_names = prefix


        for name in glob.glob(self.file_directory + '/*J[0-9]_*'):
            if not name.endswith(('pdf','bam')): #in case --plot_QC used
                self.J_files.append(name)
        for name in glob.glob(self.file_directory + '/*V[0-9]_*'):
            if not name.endswith(('pdf','bam')):
                self.V_files.append(name)
            #If identity of J does not matter for downstream analysis
        if include_all_v:
            for name in glob.glob(self.file_directory + '/*unclear_J*'):
                self.J_files.append(name)
            for name in glob.glob(self.file_directory + '/*unclear_V*'):
                self.V_files.append(name)

        #In case files are in incorrect order
        self.V_files = sorted(self.V_files)
        self.J_files = sorted(self.J_files)

        #Check both V_files and J_files have same amount of files
        if len(self.V_files) == 0 or len(self.J_files) == 0:
            raise Exception('Files in directory not identified correctly.')

        #Check the prefix of all files matches
        V_prefix_lst = []
        J_prefix_lst = []


        for item in self.V_files:
            V_prefix_lst.append(re.split('(_V|_unclear)', os.path.basename(item))[0])

        for item in self.J_files:
            J_prefix_lst.append(re.split('(_J|_unclear)', os.path.basename(item))[0])

        #Needed for merge
        self.j_prefix = J_prefix_lst[0] #should be all the same
        self.v_prefix = V_prefix_lst[0]

        assert V_prefix_lst.count(V_prefix_lst[0]) == len(V_prefix_lst), 'Input directory has to contain only one sample!'
        assert J_prefix_lst.count(J_prefix_lst[0]) == len(J_prefix_lst), 'Input directory has to contain only one sample!'


        assert len(self.V_files) == len(self.J_files), 'Number of V and J files is not the same!'

        print('Number of V files:', len(self.V_files), 'Number of J files:', len(self.J_files))



        #Create directories
        dir_nam = Path(os.path.abspath(self.V_files[0])).parents[1] #1 dir up, create outside of preclean directory
        self.tmp_dir = str(dir_nam) + '/' + str(V_prefix_lst[0]) + '_UMI_tmp'
        if out_dir == None:
            self.out_dir = str(dir_nam) + '/' + str(V_prefix_lst[0]) + '_Deduplicated'
            self.out_dir_prefix = str(V_prefix_lst[0])
        else:
            self.out_dir = out_dir
            self.out_dir_prefix = os.path.basename(out_dir)

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



    def write_v_region_file(self, spe='mmu'):
        '''Required by samtools view -L command
        '''
        #write v_region file, required for samtools -L
        v_region = open(self.tmp_dir + '/' + 'DJ_region.txt', 'w')

        dj_r = general.species(spe).dj()
        v_region.write(dj_r[0] + '\t' + dj_r[1] + '\t' + dj_r[2]) #last D + end of igh
        v_region.close()



    def remove_duplicates(self, spe='mmu', barcode_1='GACTCGT', barcode_2='CTGCTCCT', verbose=True, write=True, plot=True, cores=1, dedup_full=False):
        '''Use umi_tools to deduplicate using UMI and start position
        :param spe: species/organism
        :param barcode_1:
        :param barcode_2:
        :param verbose: print stats
        :param write: keep temporary directory with intermediate files (e.g. BAM)
        :param plot: Plot before and after deduplication
        :param cores: number of threads to use
        '''

        try:
            subprocess.check_output(['umi_tools', '-h'])
            subprocess.check_output(['samtools', '--help'])
        except OSError:
            raise RuntimeError('umi_tools/samtools not found; put directory in $PATH\n')

        Before = []
        After = []


        for v_file in self.V_files:

            if verbose:
                print('Processing:', v_file)

            # Do all extractionss
            out_name = os.path.basename(v_file)

            if barcode_1 in v_file:
                extract_args = ['umi_tools', 'extract',
                                '-I', v_file,
                                '-S', self.tmp_dir + '/' + out_name + '_umi_ext',
                                '-p', 'NNNNNNXXXXXXX']
            elif barcode_2 in v_file:
                extract_args = ['umi_tools', 'extract',
                                '-I', v_file,
                                '-S', self.tmp_dir + '/' + out_name + '_umi_ext',
                                '-p', 'NNNNNNXXXXXXXX']
            else:
                raise Exception('No barcode found in file name:', v_file)

            subprocess.call(extract_args)

            #REVIEW: filter on MAPQ score?
            # map with bowtie?
            ref_index = general.species(spe).bowtie_index()
            igh = general.species(spe).igh() #exclude things mapping elsewhere in genome

            if barcode_1 in v_file:
                # bowtie2_wrapper.align_single(self.tmp_dir + '/' + out_name + '_umi_ext', nthreads=cores, region=igh,
                #                              trim5='7', score='L,0,-1', samtools_mapq=0, write=True, write_full=True) #Make more stringent?

                br1_bowtie2 = general.bowtie2()
                br1_bowtie2.align_single(fastq=self.tmp_dir + '/' + out_name + '_umi_ext',
                                         nthreads=cores, trim5='7', spe=spe, verbose=verbose)

                br1_bowtie2.pysam_out(algn=True)
                if dedup_full:
                    br1_bowtie2.write()
                else:
                    br1_bowtie2.write(region=igh)
                br1_bowtie2.del_tmp()

            elif barcode_2 in v_file:
                # bowtie2_wrapper. align_single(self.tmp_dir + '/' + out_name + '_umi_ext', nthreads=cores, region=igh,
                #                               trim5='8', score='L,0,-1', samtools_mapq=0, write=True, write_full=True) #Make more stringent?

                br2_bowtie2 = general.bowtie2()
                br2_bowtie2.align_single(fastq=self.tmp_dir + '/' + out_name + '_umi_ext',
                                         nthreads=cores, trim5='8', spe=spe, verbose=verbose)

                br2_bowtie2.pysam_out(algn=True)
                if dedup_full:
                    br2_bowtie2.write()
                else:
                    br2_bowtie2.write(region=igh)
                br2_bowtie2.del_tmp()


            if dedup_full:
                postfix = '_umi_ext_full.bam'
            else:
                postfix = '_umi_ext.bam'

            # create index for bam
            samtools_index = ['samtools', 'index', self.tmp_dir + '/' + out_name +postfix]
            subprocess.call(samtools_index)

            #List of files for ploting
            Before.append(self.tmp_dir + '/' + out_name + postfix)

            #print align stats
            if verbose:
                print('Pre deduplication:')
                print_bam_stats(self.tmp_dir + '/' + out_name + postfix)

            #run umi_tools dedup
            #umi_tools not fully optimized for python3
            # dedup_args = ['python2','/home/chovanec/Downloads/UMI-tools/build/lib.linux-x86_64-2.7/umi_tools/umi_tools.py', 'dedup',

            #FIXME: Don't print detailed stats!

            dedup_args = ['umi_tools', 'dedup',
                         '--method=directional-adjacency',
                         '--output-stats=' + self.tmp_dir + '/' + out_name,
                         '-I', self.tmp_dir + '/' + out_name + postfix,
                         '-S', self.tmp_dir + '/' + out_name + '_umi_ext_dedup.bam']

            subprocess.call(dedup_args)

            #TODO: Plot histogram of pre post dedup UMI counts

            #Need to re-sort bam
            #create index for bam

            samtools_sort = ['samtools', 'sort', '-O', 'bam', '-o', self.tmp_dir + '/' + out_name + '_umi_ext_dedup_sorted.bam',
                             '-T tmp', self.tmp_dir + '/' + out_name + '_umi_ext_dedup.bam']
            samtools_index = ['samtools', 'index', self.tmp_dir + '/' + out_name + '_umi_ext_dedup_sorted.bam']
            subprocess.call(samtools_sort)
            subprocess.call(samtools_index)

            #List of files for ploting
            After.append(self.tmp_dir + '/' + out_name + '_umi_ext_dedup_sorted.bam')

            #print align stats
            if verbose:
                print('Post deduplication:')
                print_bam_stats(self.tmp_dir + '/' + out_name + '_umi_ext_dedup_sorted.bam')


            #Test presence of DJ_region.txt
            if os.path.isfile(self.tmp_dir + '/' + 'DJ_region.txt'):
                samtools_region = ['samtools', 'view', '-f 0x10',
                                   self.tmp_dir + '/' + out_name + '_umi_ext_dedup_sorted.bam',
                                  '-L', self.tmp_dir + '/' + 'DJ_region.txt',
                                  '-o', self.tmp_dir + '/' + out_name + '_umi_ext_dedup_DJ.bam',
                                  '-U', self.tmp_dir + '/' + out_name + '_umi_ext_dedup_V.bam'] #V with everything else, should be excluded by bowtie region

                subprocess.call(samtools_region)
            else:
                raise Exception('Need to run write_v_region_file first')

            #bam to fastq
            samtools_bam2fq_V = ['samtools', 'bam2fq', self.tmp_dir + '/' + out_name + '_umi_ext_dedup_V.bam']
            samtools_bam2fq_DJ = ['samtools', 'bam2fq', self.tmp_dir + '/' + out_name + '_umi_ext_dedup_DJ.bam']

            with open(self.out_dir + '/' + out_name + '_umi_ext_dedup_V.fastq', 'wb') as out:
                subprocess.Popen(samtools_bam2fq_V, stdout=out)
            with open(self.out_dir + '/' + out_name + '_umi_ext_dedup_DJ.fastq', 'wb') as out:
                subprocess.Popen(samtools_bam2fq_DJ, stdout=out)

        #Plot results of deduplication
        if plot:
            #TODO: Make names more sensible! Limited space on side of plot
            bef = shlex.split(' '.join('Before_' + str(i+1) for i in range(int(len(Before)))))

            aft = shlex.split(' '.join('After_'+ str(i+1) for i in range(int(len(After)))))

            dj_list = general.species(spe).dj()
            dj_coord = dj_list[0] + ':' + dj_list[1] + '-' + dj_list[2]


            #location of general module (same location as R script)
            plot_igh_bef = ['Rscript', os.path.dirname(os.path.realpath(general.__file__)) +
            '/' + 'plot_igh.R','-o', self.out_dir + '/' + self.out_dir_prefix + '_before_dedup_coverage_DJ.pdf',
            '-n'] + bef + ['--genome', general.species(spe).genome(), '-r', dj_coord,
            '-b'] + Before

            plot_igh_aft = ['Rscript', os.path.dirname(os.path.realpath(general.__file__)) +
            '/' + 'plot_igh.R','-o', self.out_dir + '/' + self.out_dir_prefix + '_after_dedup_coverage_DJ.pdf',
            '-n'] + aft + ['--genome', general.species(spe).genome(), '-r', dj_coord,
            '-b'] + After

            subprocess.call(plot_igh_bef)
            subprocess.call(plot_igh_aft)

            plot_igh_bef = ['Rscript', os.path.dirname(os.path.realpath(general.__file__)) +
            '/' + 'plot_igh.R','-o', self.out_dir + '/' + self.out_dir_prefix + '_before_dedup_coverage_V.pdf',
            '-n'] + bef + ['--genome', general.species(spe).genome(), '-r', general.species(spe).v_region(),
            '-b'] + Before

            plot_igh_aft = ['Rscript', os.path.dirname(os.path.realpath(general.__file__)) +
            '/' + 'plot_igh.R','-o', self.out_dir + '/' + self.out_dir_prefix + '_after_dedup_coverage_V.pdf',
            '-n'] + aft + ['--genome', general.species(spe).genome(), '-r', general.species(spe).v_region(),
            '-b'] + After

            subprocess.call(plot_igh_bef)
            subprocess.call(plot_igh_aft)


        #Remove temp direcotry
        if write == False:
            shutil.rmtree(self.tmp_dir, ignore_errors=True) #ignore errors that directory not empty (doesn't remove direcotyr itself!!)

    #  extract equivalent
    def extract_original(self, merge=True):
        '''Writes out J equivalent of V reads from remove_duplicates
        :param merge: merge all the output files into one
        '''
        assert self.extract_original_run == False, 'extract_original has already been run (destructive function!)'
        self.extract_original_run = True

        #loop through all the files in directory and write V and J reads
        dir_files = glob.glob(self.out_dir + '/*.fastq') #files in directory
        #include a check to see if any J files present (only V files should be!)

        #Test if output directory contains J sequences
        if not re.split('(_J)', os.path.basename(dir_files[-1]))[0].endswith('.fastq'): #can't split V files
            raise Exception('J files already present in', self.out_dir)

        all_out_v = []
        all_out_j = []

        for name in dir_files:

            print('Processing:', name)

            base = re.split('(_umi)', os.path.basename(name))[0] #'lane3_WT_FrBC_1_GCCAAT_L003_R1_val_1_40k_V1_CTGCTCCT'
            postfix = ''.join(re.split('(_umi)', os.path.basename(name))[1:]) #'_umi_ext_dedup_DJ.fastq'

            idx = [i for i, s in enumerate(self.V_files) if base in s] #same idx applies for J_files

            out_v = self.out_dir + '/' + os.path.basename(self.V_files[idx[0]]) + postfix
            out_j = self.out_dir + '/' + os.path.basename(self.J_files[idx[0]]) + postfix

            if '_V.fastq' in postfix:
                all_out_v.append(out_v)
                all_out_j.append(out_j)

            fq_names = set()
            with general.file_open(name) as int_file:
                lines = int_file.read().splitlines()
                for item in general.fastq_parse(lines):
                    title = item[0]
                    fq_names.add(title.split('_')[0])



            #find base file in self.V_files open it, loop though it and write out stuff
            pos = [i for i, s in enumerate(self.V_files) if base in s]
            if len(pos) == 1:
                full_v_fq = self.V_files[pos[0]]
                full_j_fq = self.J_files[pos[0]]
            else:
                raise Exception('Multiple out files identified')

            # for item in fq_names:
            with open(out_v, 'w') as v_out:     #Can now overwrite int_file now
                with general.file_open(full_v_fq) as fq:
                    lines = fq.read().splitlines()
                    for item in general.fastq_parse(lines):
                        title = item[0]
                        seq = item[1]
                        thrd = item[2]
                        qual = item[3]

                        #@HWI-1KL136:214:D1MR5ACXX:5:1101:1821:2154 1:N:0:TGACCA remove 1:N:0:TGACCA
                        if title.split(' ')[0] in fq_names:
                            v_out.write(title + '\n' + seq + '\n' + thrd + '\n' + qual + '\n')

            with open(out_j, 'w') as j_out:
                with general.file_open(full_j_fq) as fq:
                    lines = fq.read().splitlines()
                    for item in general.fastq_parse(lines):
                        title = item[0]
                        seq = item[1]
                        thrd = item[2]
                        qual = item[3]

                        #@HWI-1KL136:214:D1MR5ACXX:5:1101:1821:2154 1:N:0:TGACCA remove 1:N:0:TGACCA
                        if title.split(' ')[0] in fq_names:
                            j_out.write(title + '\n' + seq + '\n' + thrd + '\n' + qual + '\n')


        if merge:
            #merge based on sample names!

            dir_files = glob.glob(self.out_dir + '/*.fastq') #new files in directory


            with open(self.out_dir + '/' + self.v_prefix + '_V.fastq', 'w') as out, fileinput.input(all_out_v) as fin:
                for line in fin:
                    out.write(line)

            with open(self.out_dir + '/' + self.j_prefix + '_J.fastq', 'w') as out, fileinput.input(all_out_j) as fin: #fin output file name
                for line in fin: #open all output files and write all lines into new file
                    out.write(line)


def print_bam_stats(bam):
    mapped = 0
    unmapped = 0
    idx = pysam.idxstats(bam)
    for line in idx.split('\n')[:-1]:
        rname, rlen, nm, nu = line.split('\t')
        mapped += int(nm)
        unmapped += int(nu)

    print('Mapped:', mapped, 'Unmapped:', unmapped)




class deduplicate:
    '''Deduplicate using J, V start and UMI
    '''

    def __init__(self, file_directory):
        self.file_directory = file_directory









def parse_args():
    parser = argparse.ArgumentParser(description='BabrahamLinkON Deduplicate')

    parser.add_argument('--input_dir', dest='in_dir', type=str, required=True, help='Input directory (created for/by preclean)')
    parser.add_argument('--cores', dest='nthreads', default=1, type=int, help='Number of cores to use, default: 1')
    parser.add_argument('--species', dest='species', default='mmu', type=str, help='Which species (mmu, hsa), default: mmu')
    parser.add_argument('--keep', action='store_true', help='Keep temporary files (good for troubleshooting)')
    parser.add_argument('--all_v', action='store_true', help='Include all V (use when identity of J does not matter)')
    parser.add_argument('--br1', dest='br1', default='GACTCGT', type=str, help='Default: GACTCGT')
    parser.add_argument('--br2', dest='br2', default='CTGCTCCT', type=str, help='Default: CTGCTCCT')
    parser.add_argument('--verbose', action='store_true', help='Print detailed progress')
    parser.add_argument('--plot', action='store_true', help='Plot V region before and after deduplication')
    parser.add_argument('--out', dest='out_dir', type=str, help='Output directory, default: creates Deduplicated in main directory')
    parser.add_argument('--skip-merge', action='store_false', help='Skip merging of deduplicated files into one (used as input for downstream analysis)')
    parser.add_argument('--dedup_full', action='store_true', help='Do not exclude deduplication of reads outside of VDJ region')


    opts = parser.parse_args()

    return opts



def main():

    # if sys.version_info.major != 3:
    #     raise RuntimeError('deduplicate requires Python 3')

    #argparse
    opts = parse_args()

    dedup = deduplicate_umi_tools(opts.in_dir) #initialise class
    dedup.initial_setup(out_dir=opts.out_dir, barcode_1=opts.br1, barcode_2=opts.br2, include_all_v=opts.all_v)
    dedup.write_v_region_file(spe=opts.species)
    dedup.remove_duplicates(spe=opts.species, barcode_1=opts.br1, barcode_2=opts.br2, verbose=opts.verbose, write=opts.keep,
    plot=opts.plot, cores=opts.nthreads, dedup_full=opts.dedup_full)
    dedup.extract_original(merge=opts.skip_merge)



    # dedup = deduplicate('/media/chovanec/My_Passport/Sync/BabrahamLinkON/tests/preclean/')
    # dedup.initial_setup(include_all_v=True)
    # dedup.write_v_region_file()
    # dedup.remove_duplicates()
    # dedup.extract_original(merge=True)

if __name__ == "__main__":
    main()
