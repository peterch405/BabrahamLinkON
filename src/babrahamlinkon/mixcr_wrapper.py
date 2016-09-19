import os
import subprocess
from pathlib import Path
import re
import sys
from babrahamlinkon import general

#Should have run the importFromIMGT.sh beforehand
def run_mixcr(V_fastq, J_fastq, nthreads, species='mmu', out_dir=None):

    # check mixcr is accessible
    try:
        subprocess.check_output(['mixcr', '-v'])
    except OSError:
        raise RuntimeError('mixcr not found; put directory in $PATH\n')

    fp_V_fastq = os.path.abspath(V_fastq)
    fp_J_fastq = os.path.abspath(J_fastq)

    #Include only V prefix for simplicity
    prefix = re.split('(_V)', os.path.basename(fp_V_fastq))[0]
    # prefix_J = re.split('(_J)', os.path.basename(fp_J_fastq))[0]

    #Create directories
    #Both V and J should have same directory
    dir_nam = Path(os.path.dirname(fp_V_fastq)).parents[0]  #1 dir up, create outside of preclean directory

    if out_dir == None:
        out_dir = str(dir_nam) + '/' + prefix + '_MiXCR'


    try:
        os.mkdir(out_dir)
        print('Writing files to:', out_dir)
    except FileExistsError:
        print('Directory', out_dir, 'already exists')
        pass


    mixcr_align = ['mixcr', 'align',
        '--library', 'local',
        '-s', species,
        '--loci', 'IGH',
        '-f',
        '--threads', str(nthreads),
        '--save-description',
        '--report', out_dir + '/' + prefix + '_MiXCR_align.log',
        fp_V_fastq, fp_J_fastq,
        out_dir + '/' + prefix + '.vdjca']

    print('Aligning')
    subprocess.call(mixcr_align)

    #Since v1.8 MiXCR by default separates clones with equal clonal sequence and different V and J genes
    mixcr_assemble = ['mixcr', 'assemble',
                     '-f',
                     '-OcloneFactoryParameters.vParameters.featureToAlign=VRegion',
                     '--report', out_dir + '/' + prefix + '_MiXCR_assemble.log',
                     out_dir + '/' + prefix + '.vdjca',
                     out_dir + '/' + prefix + '.clns']

    print('Assembling')
    subprocess.call(mixcr_assemble)


    mixcr_export_clones = ['mixcr', 'exportClones',
                          '-f',
                          out_dir + '/' + prefix + '.clns',
                          out_dir + '/' + prefix + '_clones.txt']

    print('Exporting clones')
    subprocess.call(mixcr_export_clones)

    mixcr_export_unprod_clones = ['mixcr', 'exportClones',
                                   '-f',
                                   '--filter-out-of-frames',
                                   '--filter-stops',
                                    out_dir + '/' + prefix + '.clns',
                                    out_dir + '/' + prefix + '_prod_clones.txt']

    print('Exporting productive clones')
    subprocess.call(mixcr_export_unprod_clones)


    return [out_dir + '/' + prefix + '_clones.txt', out_dir + '/' + prefix + '_prod_clones.txt']



def extract_unproductive(clones, productive_clones):

    prefix = re.split('(_clones)', os.path.basename(clones))[0]

    dir_nam = os.path.dirname(clones)

    CDR3_seq = set()

    with open(productive_clones, 'r') as pc:
        header_prod = pc.readline()
        pos = (header_prod.split('\t')).index('Clonal sequence(s)')

        for line in pc:
            CDR3 = line.split()[pos]
            CDR3_seq.add(CDR3)

    out_file = dir_nam + '/' + prefix + '_unprod_clones.txt'

    with open(out_file, 'w') as out:
        with open(clones, 'r') as c:
            header = c.readline()
            out.write(header)

            count = 0
            for line in c:
                CDR3 = line.split()[pos]
                if CDR3 not in CDR3_seq:
                    out.write(line)
                    count+=1

    return out_file



def plot_results(clones, clones_prod, clones_unprod):
    '''Run rmd script to summarise results
    '''

    prefix = re.split('(_clones)', os.path.basename(clones))[0]

    out_dir = os.path.dirname(clones)

    #simplify!!
    command ='setwd(\'' + out_dir + '\'); thetitle=' + '\'' + prefix + '\'; vj_clones=\'' + clones + '\'; vj_filt_clones=' + '\'' + clones_prod + '\'' \
    '; vj_unprod_clones=' + '\'' + clones_unprod + '\'' + '; library(knitr); library(markdown); opts_knit$set(root.dir=\'' + \
    out_dir + '\'); knit2html(\'' + os.path.dirname(os.path.realpath(general.__file__)) + '/plot_v_usage.Rmd\', \'' + out_dir + '/' + prefix + '.html\')'

    # command ='thetitle=' + '\'' + prefix + '\'; vj_clones=\'' + clones + '\'; vj_filt_clones=' + '\'' + clones_prod + '\'' \
    # '; vj_unprod_clones=' + '\'' + clones_unprod + '\'' + '; library(knitr); library(markdown); opts_knit$set(root.dir=\'' + \
    # out_dir + '\'); knit(\'' + os.path.dirname(os.path.realpath(__file__)) + '/plot_v_usage.Rmd' + \
    # '\'); markdownToHTML(\'plot_v_usage.md\', \'' + out_dir + '/' + prefix + '.html\')'

    print(command)
    rmd_args = ['R', '-e', command]
    print(rmd_args)

    subprocess.call(rmd_args)
