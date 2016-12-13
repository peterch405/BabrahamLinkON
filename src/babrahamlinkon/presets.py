
import os
#####################
###### Presets ######
#####################


class prs:
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
            # return ['chr14', 105857867, 106890699] #IgH GRCh38.7 (VDJ only)
            return ['chr14', 106324871, 107292532] #hg19
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
            # return 'GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.bowtie_index'
            return 'Homo_sapiens.GRCh37'
        else:
            print('Under construction, use mmu for now')

    ## Germline identification ##

    def germline(self):
        if self.name == 'mmu' or self.name == 'mouse' or self.name == 'mus musculus':
            return ['chr12', 113428237, 113430474] #J genes mm10
        elif self.name == 'hsa' or self.name == 'human' or self.name == 'homo sapien':
            # return ['chr14', 105863049, 105865530] #J genes without J1P and IGHD7-27 hg38
            return ['chr14',106329364,106331708]
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
                    'J2P':{'start':105864793, 'end':105864852}, #pseudogene, but still need to get rid of germline
                    'J3':{'start':105864587, 'end':105864635},
                    'J4':{'start':105864215, 'end':105864260},
                    'J5':{'start':105863814, 'end':105863862},
                    'J6':{'start':105863198, 'end':105863258}}
                    # GRCh38.p7
        else:
            print('Under construction, use mmu for now')

    ## Regions / genomes ##

    def dj(self):
        '''Location of DJ genes to remove DJ recombination (seqmonk)
        '''
        if self.name == 'mmu' or self.name == 'mouse' or self.name == 'mus musculus':
            return ['chr12','113428112', '113488696']
        elif self.name == 'hsa' or self.name == 'human' or self.name == 'homo sapien':
            # return ['chr14','105862474','105920591']
            return ['chr14','106327849','106392148'] #hg19
        else:
            print('Under construction, use mmu for now')

    def genome(self):
        '''Which genome is being used
        '''
        if self.name == 'mmu' or self.name == 'mouse' or self.name == 'mus musculus':
            return 'mm10'
        elif self.name == 'hsa' or self.name == 'human' or self.name == 'homo sapien':
            # return 'hg38'
            return 'hg19'
        else:
            print('Under construction, use mmu for now')

    def v_region(self):
        '''Coordinates of V region (seqmonk)
        '''
        if self.name == 'mmu' or self.name == 'mouse' or self.name == 'mus musculus':
            return 'chr12:113531809-116015193'
        elif self.name == 'hsa' or self.name == 'human' or self.name == 'homo sapien':
            # return 'chr14:105931103-106903920'
            return 'chr14:106394250-107295467' #hg19
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

    #REVIEW: reduce to 4 bp?
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
