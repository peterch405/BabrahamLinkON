import subprocess
import pysam
import os
import tempfile
import shutil
import shlex
from babrahamlinkon import general

# from general import species

def align_single(fastq, nthreads, region, trim5='0', score='L,0,-1', flags=('--very-sensitive', '--no-unal'),
                 spe='mmu', samtools_mapq=0, write=False, write_full=False, verbose=True, plot=False, plot_region='', out_dir=None):
    """Call bowtie2-align on single read data
    Adapted from https://gist.github.com/ArtPoon/62ccff9ef937d16c44af
    https://gist.github.com/ngcrawford/6074386
    http://stackoverflow.com/questions/8466926/using-python-subprocess-to-redirect-stdout-to-stdin
    :param fastq: full path to fastq to aligned
    :param nthreads: number of cores to Use
    :param trim: trim the barcode from 5' end
    :param score: score min to be used  by bowtie2
    :param flags: additional bowtie2 flags, remove reads that failed to align
    :param spe: species
    :param samtools_mapq: if true will filter with MAPQ 5
    :param write: keep aligned bam
    :param region: input for pysam functions
    :return: touple of (fetch pysam object, pileup pysam object)
    """

    # check bowtie2 and samtools is accessible
    try:
        subprocess.check_output(['bowtie2', '-h'])
        subprocess.check_output(['samtools', '--help'])
    except OSError:
        raise RuntimeError('bowtie2/samtools not found; put directory in $PATH\n')

    ref_index = general.species(spe).bowtie_index()

    # Create tmp files
    tmp_dir = tempfile.mkdtemp()
    tmp_prefix = os.path.join(tmp_dir, "bowtiepipe")

    if out_dir == None:
        dir_nam = os.path.dirname(fastq)
    else:
        dir_nam = out_dir

    prefix = os.path.basename(fastq).split('.')[0]

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
    samtools_out = subprocess.Popen(samtools_bam, stdin=bowtie2_out.stdout, stdout=open(tmp_prefix + '.bam', 'wb')).communicate() #avoid communicate

    #If bowtie returns error, throw exception
    # out, err = bowtie2_out.communicate()
    # if b'(ERR)' in err:
    #     raise Exception(err.decode("utf-8"))

    bowtie2_out.stdout.close()


    #Sort and index bam file
    samtools_sort = ['samtools', 'sort', '-O', 'bam', '-o', tmp_prefix + '_sorted.bam', '-T tmp', tmp_prefix + '.bam']
    samtools_index = ['samtools', 'index', tmp_prefix + '_sorted.bam']

    subprocess.call(samtools_sort)
    subprocess.call(samtools_index)


    # Read from tmp file
    sam_algn = pysam.AlignmentFile(tmp_prefix + '_sorted.bam', "rb")


    sam_fetch = sam_algn.fetch(region[0], region[1], region[2])
    # sam_pileup = sam_algn.pileup(region[0], region[1], region[2])

    if verbose:
        # print(err.decode('utf-8'))
        samtools_count = shlex.split('samtools view -F 0x904 -c ' + tmp_prefix + '_sorted.bam')
        read_count = subprocess.check_output(samtools_count)
        print('Number of reads that aligned and passed filter:', read_count.decode("utf-8"))

    if plot:

        print('Saving coverage plot to', dir_nam)

        #location of current scripts
        plot_igh_bef = ['Rscript', os.path.dirname(os.path.realpath(general.__file__)) +
        '/plot_igh.R','-o', dir_nam + '/' + prefix + '.pdf',
        '-n', 'Coverage','--genome', general.species(spe).genome(), '-r', plot_region,
        '-b', tmp_prefix + '_sorted.bam']

        subprocess.call(plot_igh_bef)


    #keep bam file
    if write:

        name_out = dir_nam + '/' + prefix + '.bam' #output same location as input
        print('Saving bam', prefix + '.bam', 'to', dir_nam)
        write_reads = pysam.AlignmentFile(name_out, 'wb', template=sam_algn)
        for read in sam_algn.fetch(region[0], region[1], region[2]):
            write_reads.write(read)

        write_reads.close()

    if write_full:

        name_out = dir_nam + '/' + prefix + '_full.bam'
        print('Saving bam', prefix + '_full.bam', 'to', dir_nam)
        write_reads = pysam.AlignmentFile(name_out, 'wb', template=sam_algn)
        for read in sam_algn.fetch():
            write_reads.write(read)
        write_reads.close()

    # Clean up the named pipe and associated temp directory
    # print(tmp_dir)
    shutil.rmtree(tmp_dir)

    return (sam_fetch)
