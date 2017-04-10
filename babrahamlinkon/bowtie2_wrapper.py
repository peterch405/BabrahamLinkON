

import os
import subprocess
import pysam
import tempfile
import shutil
import shlex
from babrahamlinkon import presets
import tempfile
import logging
import re
import sys


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


    def align_single(self, fastq, nthreads, trim5='0', score='L,0,-1', flags=('--very-sensitive', '--no-unal'), #'--quiet'
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

        logger_bowtie = logging.getLogger('bowtie2.align')
        logger_bowtie.info('Subprocess: "' + ' '.join(bowtie2_args) + '"')


        # bowtie2_out = subprocess.Popen(bowtie2_args, stdout=subprocess.PIPE)

        try:
            bowtie2_out = subprocess.Popen(
                bowtie2_args,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
            )


            process_output, log =  bowtie2_out.communicate()

            #split into log and sam file, bowtie outputs both to stdout
            # log, sep, sam = re.split('(@)', process_output.decode('utf-8'), 1)

            samtools_out = subprocess.Popen(samtools_bam, stdin=subprocess.PIPE,
                                            stdout=open(self.tmp_prefix + '.bam', 'wb'))#.communicate() #avoid communicate
            #read bowtie2 output as stdin for samtools
            # samtools_out.stdin.write((sep+sam).encode('utf-8'))
            samtools_out.stdin.write(process_output)

            samtools_out.stdin.close()
            bowtie2_out.stdout.close()

            logger_bowtie.info(log.decode('utf-8'))

        except (OSError, subprocess.CalledProcessError) as exception:
            logger_bowtie.info('Exception occured: ' + str(exception))
            logger_bowtie.info('Bowtie2 subprocess failed')
            return False
        else:
            # no exception was raised
            logger_bowtie.info('Bowtie2 subprocess finished')



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
