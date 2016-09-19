from collections import defaultdict, OrderedDict
import gzip
import re
import numpy as np
import pandas as pd
import os
import Levenshtein

try:
    from skbio.alignment import StripedSmithWaterman #put into setup.py install
except:
    print('ERROR: This script requires scikit-bio: pip install scikit-bio')
    sys.exit()


# import seaborn as sns
# import matplotlib.pyplot as plt


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


class fastqHolder:
    '''
    Store fastq split by gene (J gene)
    '''
    def __init__(self):
        self.gene_split = defaultdict(set) #split into seperate J's and germline
        self.preclean_v = defaultdict(set) #read has to be in both gene_split and preclean_v
        self.demultiplex = defaultdict(set)
        # self.germline = defaultdict(set)
        self.misprimed = defaultdict(lambda: dict())
        self.fastq_split = defaultdict(list)

    # def add_to_gene(self, gene, qname):
    #     '''Add to gene_split dictionary
    #     :param gene: name of key
    #     :param qname: fastq qname to be stored under key
    #     '''
    #     self.gene_split[gene].add(qname)

    def add_to_misprimed(self, gene, qname, cor_seq):
        '''Add to gene_split dictionary
        :param gene: name of key
        :param qname: fastq qname to be stored under key
        :param cor_seq: corrected misprimed sequence to be stored under key
        '''
        self.misprimed[gene][qname] = cor_seq

    def add_to_fastq(self, gene, qname, seq, third, qual):
        self.fastq_split[gene].append([qname, seq, third, qual])

    def write_to_file(self, gene, path):
        '''Write fastq_split to file
        '''
        with open(path, 'w') as out_file:
            for qname, seq, thrd, qual in self.fastq_split[gene]:
                out_file.write(qname + '\n' + seq + '\n' + thrd + '\n' + qual + '\n')

    def demultiplex_fastq(self, fastq, gene, path):
        '''
        :param fastq: fastq to subset
        :param gene: what to subset by (J gene in split_gene)
        :param path: write path
        :param misprimed: write out misprime corrected sequence
        '''

        with file_open(fastq) as fq:
            with open(path, 'w') as out_file:
                for qname, seq, thrd, qual in fastq_parse(fq):
                    if qname.split(' ')[0][1:] in self.demultiplex[gene]:
                        if self.misprimed: #empty dict evals to False
                            try:
                                seq = self.misprimed[gene][qname.split(' ')[0][1:]] #overwrite seq
                            except KeyError:
                                pass
                        out_file.write(qname + '\n' + seq + '\n' + thrd + '\n' + qual + '\n')

    def preclean_fastq(self, fastq, gene, path):
        '''
        :param fastq: fastq to subset
        :param gene: what to subset by (J gene in split_gene)
        :param path: write path
        :param misprimed: write out misprime corrected sequence
        '''

        with file_open(fastq) as fq:
            with open(path, 'w') as out_file:
                for qname, seq, thrd, qual in fastq_parse(fq):
                    if qname.split(' ')[0][1:] in self.gene_split[gene]:
                        if self.misprimed: #empty dict evals to False
                            try:
                                seq = self.misprimed[gene][qname.split(' ')[0][1:]] #overwrite seq
                            except KeyError:
                                pass
                        out_file.write(qname + '\n' + seq + '\n' + thrd + '\n' + qual + '\n')





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



#TODO: Find a more elegent way to store this information
class species:
    '''Which species to be analysed
    :return: Data specific to organism
    '''
    def __init__(self, name):
        self.name = name

    def J_seq(self):
        if self.name == 'mmu' or self.name == 'mouse' or self.name == 'mus musculus':
            return {'J1':'CCCTGTGCCCCAGACATCGAA',
                    'J2':'AGTGGTGCCTTGGCCCCAGTAG',
                    'J3':'ACCAGAGTCCCTTGGCCCCAGTAA',
                    'J4':'TGAGGTTCCTTGACCCCAGTAGTCCATA'}
        elif self.name == 'hsa' or self.name == 'human' or self.name == 'homo sapien':
            return {'J1':'GGTGCCCTGGCCCCAGTG',
                    'J2':'GGTGCCACGGCCCCAGAGA',
                    'J3':'ACCATTGTCCCTTGGCCCCAG',
                    'J4.1':'GACCAGGGTTCCTTGGCCCC',
                    'J4.3':'GACCAGGGTCCCTTGGCCCC',
                    'J5.1':'CAGGGTTCCTTGGCCCCAGG',
                    'J5.2':'CAGGGTTCCCTGGCCCCAGG',
                    'J6.1':'TGCCCCCAGACGTCCATACCGT',
                    'J6.2':'TGGCCCCAGACGTCCATACCGT',
                    'J6.3':'CCTTTGCCCCAGACGTCCATGTAGT',
                    'J6.4':'TTGCCCCAGACGTCCATACCGT'}
        else:
            print('Under construction, use mm for now')

    def igh(self):
        if self.name == 'mmu' or self.name == 'mouse' or self.name == 'mus musculus':
            return ['chr12', 113249830, 116015093] #IgH mm10
        elif self.name == 'hsa' or self.name == 'human' or self.name == 'homo sapien':
            return ['chr14', 105857867, 106890699] #IgH GRCh38.7 (VDJ only)
        else:
            print('Under construction, use mm for now')

    def bowtie_index(self):
        if not os.environ.get('BOWTIE2_INDEXES'):
            print("'echo 'export BOWTIE2_INDEXES=/path/to/bowtie2_indexes' >> ~/.bashrc \
                    source ~/.bashrc'")
            raise Exception('BOWTIE2_INDEXES enviromental variable not set!')

        if self.name == 'mmu' or self.name == 'mouse' or self.name == 'mus musculus':
            return 'mm10'
        elif self.name == 'hsa' or self.name == 'human' or self.name == 'homo sapien':
            return 'GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.bowtie_index'
        elif self.name == '129S1':
            return '129S1_SvImJ_GRCm38'
        else:
            print('Under construction, use mmu for now')

    ## Germline identification ##

    def germline(self):
        if self.name == 'mmu' or self.name == 'mouse' or self.name == 'mus musculus':
            return ['chr12', 113428237, 113430474] #J genes mm10
        elif self.name == 'hsa' or self.name == 'human' or self.name == 'homo sapien':
            return ['chr14', 105863049, 105865530] #J genes without J1P and IGHD7-27
        else:
            print('Under construction, use mmu for now')

    #TODO:Implement automatic import from biomart
    def J_location(self): #required to seperate germline into unique J's
        if self.name == 'mmu' or self.name == 'mouse' or self.name == 'mus musculus':
            return {'J1':{'start':113429781, 'end':113429833},
                    'J2':{'start':113429514, 'end':113429567},
                    'J3':{'start':113429085, 'end':113429132},
                    'J4':{'start':113428514, 'end':113428567}}
                    #mm10
        elif self.name == 'hsa' or self.name == 'human' or self.name == 'homo sapien':
            return {'J1':{'start':105865407, 'end':105865458},
                    'J2':{'start':105865199, 'end':105865250},
                    'J2P':{'start':105864793, 'end':105864852}, #psedogene, but still need to get rid of germline
                    'J3':{'start':105864587, 'end':105864635},
                    'J4':{'start':105864215, 'end':105864260},
                    'J5':{'start':105863814, 'end':105863862},
                    'J6':{'start':105863198, 'end':105863258}}
                    # GRCh38.p7
        else:
            print('Under construction, use mmu for now')

    ## Regions / genomes ##

    def dj(self):
        '''Location of DJ genes to remove DJ recombination
        '''
        if self.name == 'mmu' or self.name == 'mouse' or self.name == 'mus musculus':
            return ['chr12','113428112', '113484692']
        else:
            print('Under construction, use mmu for now')

    def genome(self):
        '''Which genome is being used
        '''
        if self.name == 'mmu' or self.name == 'mouse' or self.name == 'mus musculus':
            return 'mm10'
        elif self.name == 'hsa' or self.name == 'human' or self.name == 'homo sapien':
            return 'hg38'
        else:
            print('Under construction, use mmu for now')

    def v_region(self):
        '''Coordinates of V region
        '''
        if self.name == 'mmu' or self.name == 'mouse' or self.name == 'mus musculus':
            return 'chr12:113531809-116015193'
        else:
            print('Under construction, use mmu for now')

    ## Mispriming ##

    def J_for_mispriming(self): #ref created using these
        if self.name == 'mmu' or self.name == 'mouse' or self.name == 'mus musculus':
            return {'J1':'CCCTGTGCCCCAGACATC',
                    'J2':'AGTGGTGCCTTGGCCCCAGTAGTCAAA',
                    'J3':'ACCAGAGTCCCTTGGCCCCAGTAAGCAAA',
                    'J4':'TGAGGTTCCTTGACCCCAGTAGTCCAT'}
        else:
            print('Under construction, use mmu for now')

    def mispriming(self): #implement auto method to generate this dict
        if self.name == 'mmu' or self.name == 'mouse' or self.name == 'mus musculus':
            return {'J1':'ACATC',
                    'J2':'TCAAA',
                    'J3':'GCAAA',
                    'J4':'TCCAT'}
        else:
            print('Under construction, use mmu for now')

    def offset(self):
        if self.name == 'mmu' or self.name == 'mouse' or self.name == 'mus musculus':
            return {'J1':{'J2':2, 'J3':2, 'J4':2},
                    'J2':{'J1':8},
                    'J3':{'J1':8},
                    'J4':{'J1':8}}
        else:
            print('Under construction, use mmu for now')

    def replace(self):
        if self.name == 'mmu' or self.name == 'mouse' or self.name == 'mus musculus':
            return {'J1':'CCCTGTGCCCCAG',
                    'J2':'AGTGGTGCCTTGGCCCCAGTAG',
                    'J3':'ACCAGAGTCCCTTGGCCCCAGTAA',
                    'J4':'TGAGGTTCCTTGACCCCAGTAG'}
        else:
            print('Under construction, use mmu for now')






class SSW_align:
    '''Align reads using StripedSmithWaterman from scikit-bio
    http://scikit-bio.org/docs/latest/generated/skbio.alignment.StripedSmithWaterman.html#skbio.alignment.StripedSmithWaterman
    '''
    #align_coord = {'J1':0,'J2':27, 'J3':55, 'J4':85}
    def __init__(self):
        self.alignments = []
        self.ref_overlap = defaultdict(list) #created by reference function
        self.align_coord = defaultdict(int) #created by reference function

    def align(self, target, qry, score=0, misprim_cor=True, mismatch_allow=1, spe='mmu'):
        '''Align against reference and return which J sequence aligns to
        :param target: reference
        :param qry: sequence to align against reference
        :param score: score filter, min score to pass as an alignment
        :param misprim_cor: perform mispriming correction
        :param mismatch_allow: number of mismatches to allow in sequence beyond J primer
        :param spe: which organism
        :return: read identity (store alignment in alignments)
        '''
        #Query needs to be a single string
        assert isinstance(qry, str), 'Align query should not be a list'

        #for mm take only first 21 nts (smallest J)
        # j_size_ord = sorted(species(spe).J_seq().values(), key=len)
        # shortest_J = len(j_size_ord[0])
        #
        # off_size_ord = [list(item.values()) for item in species(spe).offset().values()]
        # values = [item for sublist in off_size_ord for item in sublist]
        # shortest_offset = min(values)

         #use few bp to do initial identification

        ini_len = 5
        restart = True
        while restart:
            restart = False
            J_identities = defaultdict(int)
            for J in species(spe).replace().keys():
                dist = Levenshtein.distance(qry[:ini_len], species(spe).J_seq()[J][:ini_len])
                J_identities[J] = dist #if incorrect will catch later
            #chose J with lowest score
            if J_identities and ini_len < 10: #if dict not empty and iteration is less then 10bp
                min_dist = min(J_identities.values())
                min_dist_key = [k for k, v in J_identities.items() if v == min_dist]
                if len(min_dist_key) > 1: #if two keys equal score increase size and repeat until only one hit obtained
                    restart = True
                    ini_len += 1
                else:
                    initial_identity = min_dist_key[0]
            else:
                return 'other' #no hits

        #while loop in case initial identity is wrong
        redo = True
        add = 0
        not_added = True #prevent infinit loop
        not_switched = True
        while redo:
            redo = False

             #Use initial identity to align x bp
            align_size = len(species(spe).J_for_mispriming()[initial_identity]) + add

            #allows 1 error per 5 bases
            qry_err = align_size-round(align_size/5)

            query_out = StripedSmithWaterman(qry[:align_size], match_score=1, mismatch_score=0,
                                                 gap_open_penalty=1, gap_extend_penalty=1,
                                                 score_filter=qry_err)
            #score_filter filter out low score alignments (lowers computational time)

            alignment = query_out(target)

            if alignment['target_begin'] != -1: #target_begin and cigar not returned if score low


                #Which J does alignment correspond to
                #Overlap
                try:
                    read_identity = self.ref_overlap[alignment['target_begin']]
                except KeyError:
                    print('Need to run reference function first!')


                #if read doesn't match initial identification
                if not_switched:
                    if read_identity != initial_identity:
                        redo = True
                        initial_identity = read_identity
                        not_switched = False
                        continue

                    elif misprim_cor:
                        self.alignments.append(alignment)


                        s = alignment['cigar']
                        matches = re.findall(r'(\d+)([A-Z]{1})', s)
                        cigar_dict = ([{'type':m[1], 'length':int(m[0])} for m in matches])

                        start = alignment['target_begin']
                        end = alignment['target_end_optimal']+1
                        q_start = alignment['query_begin']
                        q_end = alignment['query_end']+1

                        print_seq = '-'*start
                        poscounter = q_start

                        for cig in cigar_dict:
                            if cig['type'] == 'M':
                                #match just print seq
                                c_start = poscounter
                                poscounter += cig['length']
                                print_seq = print_seq + alignment['query_sequence'][c_start:poscounter]
                            elif cig['type'] == 'D':
                                print_seq = print_seq + '-'*cig['length']
                            elif cig['type'] == 'I' and q_end-poscounter < 5: #insertions may cause aligning sequence to be too short!
                                #don't skip insertions near the end
                                c_start = poscounter
                                poscounter += cig['length']  #skip over insertion
                                print_seq = print_seq + alignment['query_sequence'][c_start:poscounter]
                            elif cig['type'] == 'I':
                                poscounter += cig['length']  #skip over insertion

                        print_seq = print_seq + '-'*(len(alignment['target_sequence'])-end)


                        # alignment['query_sequence'][q_end:q_end+dist_from_end]

                        target_end = alignment['target_end_optimal']+1
                        ref_end = self.align_coord[read_identity]['end']

                        #if last characters don't align get from query
                        dist_from_end = self.align_coord[read_identity]['end']-end
                        last_let = alignment['query_sequence'][q_end:q_end+dist_from_end]


                        differences = defaultdict(int)
                        for key in species(spe).mispriming().keys():

                            try: #if in offset dict do ... else do ...
                                off = species(spe).offset()[read_identity][key]
                                seq_chk = qry[q_end-off:q_end-off+5]
                            except KeyError:
                                off = 0
                                seq_chk = print_seq[ref_end-5+off:ref_end-dist_from_end+off]+last_let

                            if len(seq_chk) >= 4:
                                differences[key] = (Levenshtein.distance(seq_chk, species(spe).mispriming()[key][:len(seq_chk)]))
                            elif not_added:
                                add = 5-len(seq_chk)
                                redo=True
                                not_added = False
                            else:
                                return 'unclear'

                        if redo: #redo with a longer read if too many insertions present
                            continue



                        min_val = min(differences.values())

                        if min_val <= 0: #allow only 0 mismatches by default in beyond J seq
                            #retrieve key of min val (J)
                            min_val_key = [k for k, v in differences.items() if v == min_val]

                            if len(min_val_key) > 1: #if there are multiple values with same score (should not happen)
                                return 'unclear'
                            else:
                                return (min_val_key[0], differences.values(), species(spe).replace()[min_val_key[0]] + qry[start:]) #identity of J and corrected sequence

                        else: #if more than allowed mismatches
                            return 'unclear'

                    else: #if not misprime correcting return identity of primer seq
                        return read_identity
                else:
                    return 'unclear'
        else: #doesn't match primer seq
            return 'other'




        #     if misprim_cor:
        #         try:
        #             corrections = species(spe).mispriming()[read_identity]
        #             align_end = alignment['target_end_optimal']
        #             repl = species(spe).replace()[read_identity]
        #             offset_vals = species(spe).offset()[read_identity]
        #
        #             #check for deletions and insertions
        #             matches = re.findall(r'(\d+)([A-Z]{1})', alignment['cigar'])
        #             cigar_dict = ([{'type':m[1], 'length':int(m[0])} for m in matches])
        #
        #             #Adjust start position depending on deletions and insertions
        #             add = 0
        #             subtract = 0
        #             for cig in cigar_dict:
        #                 if cig['type'] == 'D':
        #                     subtract = cig['length']
        #                 elif cig['type'] == 'I':
        #                     add = cig['length']
        #
        #
        #             #Differences between ref and qry 5bp beyond J
        #             differences = defaultdict(int)
        #
        #             for correct_key in corrections.keys():
        #                 len_of_misprim_seg = len(mispriming[read_identity][correct_key])
        #
        #                 start =  offset_vals[correct_key] + add - subtract
        #                 end = start + len_of_misprim_seg  #number of nucleotides after J primer
        #
        #                 differences[correct_key] = Levenshtein.distance(qry[start:end], corrections[correct_key])
        #
        #
        #             #which correct key has lowest difference
        #             min_val = min(differences.values())
        #
        #             if min_val <= mismatch_allow: #allow only 0 mismatches by default in beyond J seq
        #                 #retrieve key of min val (J)
        #                 min_val_key = [k for k, v in differences.items() if v == min_val]
        #
        #                 if len(min_val_key) > 1: #if there are multiple values with same score (should not happen)
        #                     return 'unclear'
        #                 else:
        #                     return (min_val_key[0], repl[min_val_key[0]] + qry[start:]) #identity of J and corrected sequence
        #             else: #if more than allowed mismatches
        #                 return 'unclear'
        #
        #         except KeyError:
        #             print('Identified J not in mispriming dictionary')
        #
        #
        #     return read_identity #if not misprime correcting return identity of primer seq
        # else:
        #     return 'other' #doesn't match primer seq


    def print_align(self, num=100, print_out=True):
        '''Print pretty alignemnt
        :param num: how many alignments to print
        '''
        #Target sequence should be the same for all (reference)
        seqs = []
        if print_out:
            print(self.alignments[0]['target_sequence'])

        seqs.append(self.alignments[0]['target_sequence'])
        count = 0

        for item in self.alignments:
            if count <= num: #print only num alignments
                s = item['cigar']
                matches = re.findall(r'(\d+)([A-Z]{1})', s)
                cigar_dict = ([{'type':m[1], 'length':int(m[0])} for m in matches])

                start = item['target_begin']
                end = item['target_end_optimal']+1
                q_start = item['query_begin']
                q_end = item['query_end']+1

                print_seq = '-'*start
                poscounter = q_start

                for cig in cigar_dict:
                    if cig['type'] == 'M':
                        #match just print seq
                        c_start = poscounter
                        poscounter += cig['length']
                        print_seq = print_seq + item['query_sequence'][c_start:poscounter]
                    elif cig['type'] == 'D':
                        #deletion print '-'
                        print_seq = print_seq + '-'*cig['length']
                    elif cig['type'] == 'I':
#                         print_seq = print_seq + '^'*cig['length']
                        poscounter += cig['length']  #skip over insertion

                print_seq = print_seq + '-'*(len(item['target_sequence'])-end)
                seqs.append(print_seq)
                if print_out:
                    print(print_seq)
                count += 1
            else:
                break

        return(seqs)



    def reference(self, spe='mmu'):
        '''Create reference from J sequences supplied (will seperate each J by 6 Ns)
        :param spe: which organism
        '''
        ref = ''

        for key in sorted(species(spe).J_for_mispriming().keys()):
            start = len(ref)
            ref = ref + species(spe).J_for_mispriming()[key] + 'N'*6
            end = len(ref)-6
            if start == 0:
                self.align_coord[key] = {'start':start, 'end':end}
            else:
                self.align_coord[key] = {'start':start, 'end':end} #0 based coordinates

            for val in range(start, end):
                self.ref_overlap[val] = key
        ref = ref[:-6]

        return ref



def trim_fastq(J_region, out_file, trim_size=100):
    '''Trim all reads to a set length
    :param J_region: J fastq
    :param out_file: path of output file
    :param trim_size: what size to trim all the reads down to (default: 100)
    :return: fastq with trimmed sequences
    '''
    with open(out_file, 'w') as out:
        with file_open(J_region) as Jr:
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


# def plot_pileup(pysam_pileup, pysam_fetch, plot_name, out_dir, spe='mm'):
#     '''Create a read pileup plot
#     :param pysam_pileup: pysam pileup object
#     :param plot_name: Name of saved image
#     :param spe: which organism
#     :return: .png
#     '''
#     read_dist = []
#     read_pos = []
#     j_pos = []
#
#     keys = species(spe).J_location().keys()
#     location_j = species(spe).J_location()
#
#
#     for x in pysam_pileup:
#         read_dist.append(x.nsegments)
#         read_pos.append(x.pos)
#
#     for pos in read_pos:
#         found_j = ''
#         for key in keys:
#             if pos >=(location_j[key]['start']) and pos <=(location_j[key]['end']):
#                     found_j = key
#         if not found_j: #is found_j empty?
#             j_pos.append(0)
#         else:
#             j_pos.append(int(re.sub(r"\D", "", found_j))) #only digit in string
#
#
#     d = np.asarray(read_dist)
#     df = pd.DataFrame({'Read pileup':read_dist, 'Position':read_pos, 'J genes':j_pos})
#
#     fig, (ax1, ax2, ax3) = plt.subplots(nrows=3, ncols=1)
#
#     df.plot.line(x="Position", y="Read pileup", ax=ax1)
#     df.plot.line(x="Position", y="J genes", yticks=(1,2,3,4), ax=ax2) #plot as gantt chart?
#
#     J_count = OrderedDict()
#     for key in sorted(keys, reverse=True):
#         J_count[key] = 0
#
#     count = 0
#
#     for y in pysam_fetch:
#         for key in keys:
#             if y.pos >=(location_j[key]['start']-50) and y.pos <=(location_j[key]['end']+100):
#                 J_count[key] += 1
#         count += 1
#
#     J_count['other'] = count - (sum(J_count.values()))
#
#
#     dict_df=pd.DataFrame.from_dict(J_count, orient='index')
#     dict_df.columns = ['Reads over Js']
#     dict_df.plot.bar(ax=ax3)
#     plt.tight_layout(rect=[0.5, 0, 2, 2])
#
#     fig.savefig(out_dir + '/' + plot_name + '.png', bbox_inches='tight')
