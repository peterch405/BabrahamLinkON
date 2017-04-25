from collections import defaultdict
import gzip
import os
from babrahamlinkon import presets
from itertools import groupby



def file_open(filename):
    """
    Open as normal or as gzip
    Faster using zcat?
    """
    #does file exist?
    f = open(filename,'rb')
    if (f.read(2) == b'\x1f\x8b'): #compressed alsways start with these two bytes
        f.seek(0) #return to start of file
        return gzip.GzipFile(fileobj=f, mode='rb')
    else:
        f.seek(0)
        return f



def fastq_parse(fp):
    """
    Parse fastq file.
    """
    linecount = 0
    name, seq, thrd, qual = [None] * 4
    for line in fp:

        linecount += 1
        if linecount % 4 == 1:
            try:
                name = line.decode('UTF-8').rstrip()
            except AttributeError:
                name = line.rstrip()
            assert name.startswith('@'),\
                   "ERROR: The 1st line in fastq element does not start with '@'.\n\
                   Please check FastQ file near line number %s" % (linecount)
        elif linecount % 4 == 2:
            try:
                seq = line.decode('UTF-8').rstrip()
            except AttributeError:
                seq = line.rstrip()
        elif linecount % 4 == 3:
            try:
                thrd = line.decode('UTF-8').rstrip()
            except AttributeError:
                thrd = line.rstrip()
            assert thrd.startswith('+'),\
                   "ERROR: The 3st line in fastq element does not start with '+'.\n\
                   Please check FastQ file near line number %s" % (linecount)
        elif linecount % 4 == 0:
            try:
                qual = line.decode('UTF-8').rstrip()
            except AttributeError:
                qual = line.rstrip()
            assert len(seq) == len(qual),\
                    "ERROR: The length of Sequence and Quality aren't equal.\n\
                    Please check FastQ file near line number %s" % (linecount)

            yield name, seq, thrd, qual,
            name, seq, thrd, qual = [None] * 4




def fasta_parse(fp):
    '''
    Parse fasta file
    '''
    linecount = 0
    name, seq = [None] * 2
    for line in fp:

        linecount += 1
        if linecount % 2 == 1:

            try:
                name = line.decode('UTF-8').rstrip()
            except AttributeError:
                name = line.rstrip()

            assert name.startswith('>'),\
                   "ERROR: The 1st line in fasta element does not start with '>'.\n\
                   Please check FastQ file near line number %s" % (linecount)
        elif linecount % 2 == 0:
            try:
                seq = line.decode('UTF-8').rstrip()
            except AttributeError:
                seq = line.rstrip()

            yield name, seq
            name, seq = [None] * 2

#TODO: better fasta parser in case sequence is not all on the same line
def fasta_iter(fasta_name):
    '''
    given a fasta file. yield tuples of header, sequence
    https://drj11.wordpress.com/2010/02/22/python-getting-fasta-with-itertools-groupby/
    https://www.biostars.org/p/710/
    '''
    #read fasta from stdin
    if fasta_name.startswith('>'):
        file_opened = False
        fasta_in = fasta_name.splitlines()
    else:
        file_opened = True
        fasta_in = open(fasta_name, 'r')

    # with open(fasta_name) as fasta_in:
    # ditch the boolean (x[0]) and just keep the header or sequence since we know they alternate.
    faiter = (x[1] for x in groupby(fasta_in, lambda line: line[0] == '>'))
    for header in faiter:
        # drop the ">"
        header = header.__next__()[1:].strip() #first group
        # join all sequence lines to one.
        seq = "".join(s.strip() for s in faiter.__next__()) #second group
        yield header, seq

    if file_opened:
        fasta_in.close()


def fastq_to_fasta_iter(fastq):
    '''Convert fastq to fasta
    '''

    with file_open(fastq) as fq:
        for item in fastq_parse(fq):
            title = item[0]
            seq = item[1]

            yield '>' + title.split(' ')[0][1:] + '\n' + seq + '\n'




def trim_fastq(J_region, out_file, trim_size=100):
    '''Trim all reads to a set length
    :param J_region: J fastq
    :param out_file: path of output file
    :param trim_size: what size to trim all the reads down to (default: 100)
    :return: fastq with trimmed sequences
    '''
    with open(out_file, 'w') as out, file_open(J_region) as Jr:
            lines = Jr.read().splitlines()
            for item in fastq_parse(lines):
                title = item[0]
                seq = item[1]
                thrd = item[2]
                qual = item[3]

                if len(seq) <= trim_size:
                    out.write(title + '\n' + seq + '\n' + thrd + '\n' + qual + '\n')
                else:
                    seq = seq[:trim_size]
                    qual = qual[:trim_size]
                    out.write(title + '\n' + seq + '\n' + thrd + '\n' + qual + '\n')



def reverse_complement(seq):
    '''Return reverse complement of DNA sequence
    '''
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N':'N'}
    return ''.join([complement[base] for base in seq.upper()[::-1]])


def reverse_string(my_str):
    '''
    Reverse the qual of a sequence
    '''
    return ''.join(let for let in my_str.upper()[::-1])


def check_qual(umi_qual, q_score=30):
    for val in umi_qual:
        phred = ord(val)-33
        assert phred <= 41 and phred >= 0, 'Phred score out side of range 0-41'
        if phred < q_score:
            return True
        else:
            return False
