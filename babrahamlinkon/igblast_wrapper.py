
import subprocess
import shlex
from babrahamlinkon.general import fasta_iter
from itertools import islice
from joblib import Parallel, delayed
import pkg_resources
import logging


def splice_fasta(fasta_path, chunk_size):
    '''Splice fasta into chunks
    '''

    fasta_chunk = ''
    chunk = 0
    for name, seq in fasta_iter(fasta_path):
        chunk += 1
        fasta_chunk = fasta_chunk + ('>' + name + '\n' + seq + '\n')
        if chunk_size == chunk:
            yield fasta_chunk
            chunk = 0
            fasta_chunk = ''
    if fasta_chunk: #if not empty
        yield fasta_chunk

#TODO:They fixed num_threads in version > 1.5, could test this later
def igblast_worker(fasta, spe, custom_ref, aux_file, additional_flags):

    #if fasta string - stdin else file in
    if fasta.startswith('>'):
        fasta_input = '-'
    else:
        fasta_input = fasta

    DATA_PATH = pkg_resources.resource_filename('babrahamlinkon', 'resources/IgBlast_database/')

    if spe == 'mmu':
        if custom_ref:
            cmd = ['igblastn',
                '-germline_db_V', DATA_PATH + 'Mus_musculus_IGHV_AEC',
                '-germline_db_D', DATA_PATH + 'Mus_musculus_IGHD_AEC',
                '-germline_db_J', DATA_PATH + 'Mus_musculus_IGHJ_AEC',
                '-auxiliary_data', aux_file,
                '-domain_system', 'imgt',
                '-ig_seqtype', 'Ig' ,
                '-organism', 'mouse',
                '-num_threads', '1',
                '-outfmt', '7 std qseq sseq btop',
                '-query', fasta_input, '-out', '-',
                ]
        else:
            cmd = ['igblastn',
                '-germline_db_V', DATA_PATH + 'Mus_musculus_IGHV',
                '-germline_db_D', DATA_PATH + 'Mus_musculus_IGHD',
                '-germline_db_J', DATA_PATH + 'Mus_musculus_IGHJ',
                '-auxiliary_data', aux_file,
                '-domain_system', 'imgt',
                '-ig_seqtype', 'Ig' ,
                '-organism', 'mouse',
                '-num_threads', '1',
                '-outfmt', '7 std qseq sseq btop',
                '-query', fasta_input, '-out', '-',
                ]
        # if additional_flags is not None:
        #     cmd = cmd + additional_flags
    elif spe == 'mmuk':
        cmd = ['igblastn',
            '-germline_db_V', DATA_PATH + 'Mus_musculus_IGKV',
            '-germline_db_J', DATA_PATH + 'Mus_musculus_IGKJ',
            '-auxiliary_data', aux_file,
            '-domain_system', 'imgt',
            '-ig_seqtype', 'Ig' ,
            '-organism', 'mouse',
            '-num_threads', '1',
            '-outfmt', '7 std qseq sseq btop',
            '-query', fasta_input, '-out', '-',
            ]
    elif spe == 'hsa':
        cmd = ['igblastn',
            '-germline_db_V', DATA_PATH + 'Homo_sapiens_IGHV',
            '-germline_db_D', DATA_PATH + 'Homo_sapiens_IGHD',
            '-germline_db_J', DATA_PATH + 'Homo_sapiens_IGHJ',
            '-auxiliary_data', aux_file,
            '-domain_system', 'imgt',
            '-ig_seqtype', 'Ig' ,
            '-organism', 'human',
            '-num_threads', '1',
            '-outfmt', '7 std qseq sseq btop',
            '-query', fasta_input, '-out', '-',
            ]

    if additional_flags is not None:
        cmd = cmd + additional_flags

    result = subprocess.check_output(cmd, input=fasta if fasta_input == '-' else None, universal_newlines=True)

    return result



def run_igblast(fasta, out_fmt, splice_size, spe, nprocs, custom_ref, aux_file=None, additional_flags=None):
    '''Run IgBlast in parallel
    :param fasta: fasta file path
    :param out_fmt: out path for fmt7 file
    :param splice_size: fasta file split size (number of lines)
    :param spe: species
    :param nprocs: number of processes to run in parallel
    :param additional_flags: a list of additional flags to pass to igblast
    '''

    DATA_PATH = pkg_resources.resource_filename('babrahamlinkon', 'resources/IgBlast_database/')
    #try locating the aux files for igblast
    if spe == 'mmu' and aux_file == None:
        #won't work on some systems
        # mouse_aux = subprocess.check_output(['locate', '-br', 'mouse_gl.aux$'], universal_newlines=True)
        # aux_file = mouse_aux.split('\n')[0]
        aux_file = DATA_PATH + 'optional_file/mouse_gl.aux'
    elif spe == 'hsa' and aux_file == None:
        # human_aux = subprocess.check_output(['locate', '-br', 'human_gl.aux$'], universal_newlines=True)
        # aux_file = human_aux.split('\n')[0]
        aux_file = DATA_PATH + 'optional_file/human_gl.aux'

    #returns a list of all the results
    results = Parallel(n_jobs=nprocs)(delayed(igblast_worker)(chunk, spe, custom_ref, aux_file, additional_flags) for chunk in splice_fasta(fasta, 10000))

    with open(out_fmt, 'w') as out_file:
        for item in range(len(results)):
            out_file.write(results[item])



def parse_igblast(fmt, fasta, spe, custom_ref):
    '''run changeo MakeDb.py
    '''
    DATA_PATH = pkg_resources.resource_filename('babrahamlinkon', 'resources/IgBlast_database/')

    if spe == 'mmu':
        if custom_ref:
            cmd = ['MakeDb.py', 'igblast', '-i', fmt, '-s', fasta,
            	'-r', DATA_PATH + 'Mus_musculus_IGH[JDV]_AEC.fasta',
            	'--regions', '--scores', '--cdr3', '--partial']
        else:
            cmd = ['MakeDb.py', 'igblast', '-i', fmt, '-s', fasta,
            	'-r', DATA_PATH + 'Mus_musculus_IGH[JDV].fasta',
            	'--regions', '--scores', '--cdr3', '--partial']
    elif spe == 'mmuk':
        cmd = ['MakeDb.py', 'igblast', '-i', fmt, '-s', fasta,
            '-r', DATA_PATH + 'Mus_musculus_IGK[JV].fasta',
            '--regions', '--scores', '--cdr3', '--partial']
    elif spe == 'hsa':
        cmd = ['MakeDb.py', 'igblast', '-i', fmt, '-s', fasta,
        	'-r', DATA_PATH + 'Homo_sapiens_IGH[JDV].fasta',
        	'--regions', '--scores', '--cdr3', '--partial']



    logger_igblast = logging.getLogger('igblast.changeo')
    logger_igblast.info('Subprocess: "' + ' '.join(cmd) + '"')

    try:
        db_tab = subprocess.Popen(
            cmd, universal_newlines=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
        )

        process_output, _ =  db_tab.communicate()

        logger_igblast.info(process_output)

    except (OSError, subprocess.CalledProcessError) as exception:
        logger_igblast.info('Exception occured: ' + str(exception))
        logger_igblast.info('Changeo subprocess failed')
        return False
    else:
        # no exception was raised
        logger_igblast.info('Changeo subprocess finished')


    # db_tab = subprocess.check_output(cmd, universal_newlines=True)
#
# my_fa = '/media/chovanec/My_Passport/Old_vdj_seq_data/lane5_TGACCA_WT_BC_L005_R1_val_1_40k_Deduplicated/lane5_TGACCA_WT_BC_L005_R1_val_1_40k_dedup.fasta'
# out_fmt = '/media/chovanec/My_Passport/Old_vdj_seq_data/lane5_TGACCA_WT_BC_L005_R1_val_1_40k_Deduplicated/lane5_TGACCA_WT_BC_L005_R1_val_1_40k_dedup.fmt7'
# out_tab = '/media/chovanec/My_Passport/Old_vdj_seq_data/lane5_TGACCA_WT_BC_L005_R1_val_1_40k_Deduplicated/lane5_TGACCA_WT_BC_L005_R1_val_1_40k_dedup'

# run_igblast(my_fa, out_fmt, 10, 'mmu', 2, additional_flags=['-num_alignments_V', '1'])

#
#
# DATA_PATH = pkg_resources.resource_filename('babrahamlinkon', 'resources/IgBlast_database/')
#
# with open(DATA_PATH, 'rb') as test:
#     for line in test:
#         print(line)
