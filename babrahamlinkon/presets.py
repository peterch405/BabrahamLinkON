
import os
#####################
###### Presets ######
#####################

#TODO: use https://docs.python.org/3/library/configparser.html ?

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
        elif self.name == 'mmuk' or self.name == 'mouse kappa' or self.name == 'mus musculus kappa':
            return {'J1':'TTTGATTTCCAGCTTGGTGCCTCC'[4:],
                    'J2':'TTTATTTCCAGCTTGGTCCCCCCT'[4:],
                    'J4':'CGTTTTATTTCCAACTTTGTCCCCGA'[4:],
                    'J5':'CAGCTCCAGCTTGGTCCCAGC'[4:]}
                    #first 4bp removed due to low quality resulting from low diversity
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
                    'J6.4':'TTGCCCCAGACGTCCATACCGT',
                    'J6.c':'TTCCCCCAGACGTCCATACCGT'}  #Misprimed sequence
        else:
            print('Only mmu, mmuk and hsa available for now')

    #needed for short pipeline where v genes identity additionally identified using bowtie (might remove in future)
    def ig(self):
        if self.name == 'mmu' or self.name == 'mouse' or self.name == 'mus musculus':
            return ['chr12', 113249830, 116015093] #IgH mm10
        if self.name == 'mmuk' or self.name == 'mouse kappa' or self.name == 'mus musculus kappa':
            return ['chr6', 67506809, 70747460] #IgK mm10
        elif self.name == 'hsa' or self.name == 'human' or self.name == 'homo sapien':
            # return ['chr14', 105857867, 106890699] #IgH GRCh38.7 (VDJ only)
            return ['chr14', 106324871, 107292532] #hg19
        else:
            print('Only mmu, mmuk and hsa available for now')


    def bowtie_index(self):
        if not os.environ.get('BOWTIE2_INDEXES'):
            print("'echo 'export BOWTIE2_INDEXES=/path/to/bowtie2_indexes' >> ~/.bashrc \
                    source ~/.bashrc'")
            raise Exception('BOWTIE2_INDEXES enviromental variable not set!')

        if not os.environ.get('BOWTIE2_REF'):
            if self.name == 'mmu' or self.name == 'mouse' or self.name == 'mus musculus':
                return 'mm10'
            elif self.name == 'mmuk' or self.name == 'mouse kappa' or self.name == 'mus musculus kappa':
                return 'mm10'
            elif self.name == 'hsa' or self.name == 'human' or self.name == 'homo sapien':
                # return 'GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.bowtie_index'
                return 'Homo_sapiens.GRCh37'
            else:
                print('Only mmu, mmuk and hsa available for now')
        else:
            return os.environ.get('BOWTIE2_REF')

    ## Germline identification ##

    def germline(self):
        if self.name == 'mmu' or self.name == 'mouse' or self.name == 'mus musculus':
            return ['chr12', 113428237, 113430474] #J genes mm10
        elif self.name == 'mmuk' or self.name == 'mouse kappa' or self.name == 'mus musculus kappa':
            return ['chr6', 70720525, 70725451] #J genes mm10
        elif self.name == 'hsa' or self.name == 'human' or self.name == 'homo sapien':
            # return ['chr14', 105863049, 105865530] #J genes without J1P and IGHD7-27 hg38
            return ['chr14',106329364,106331708]
        else:
            print('Only mmu, mmuk and hsa available for now')

    #TODO:Implement automatic import from biomart
    def J_location(self): #required to seperate germline into unique J's
        if self.name == 'mmu' or self.name == 'mouse' or self.name == 'mus musculus':
            return {'J1':{'start':113429781, 'end':113429833},
                    'J2':{'start':113429514, 'end':113429567},
                    'J3':{'start':113429085, 'end':113429132},
                    'J4':{'start':113428514, 'end':113428567}}
                    #mm10
        elif self.name == 'mmuk' or self.name == 'mouse kappa' or self.name == 'mus musculus kappa':
            return {'J1':{'start':70722562, 'end':70722599},
                    'J2':{'start':70722916, 'end':70722958},
                    # 'J3':{'start':70723223, 'end':70723258},
                    'J4':{'start':70723548, 'end':70723585},
                    'J5':{'start':70723886, 'end':70723927}}
                    #mm10 J3 pseudogene
        elif self.name == 'hsa' or self.name == 'human' or self.name == 'homo sapien':
            # return {'J1':{'start':105865407, 'end':105865458},
            #         'J2':{'start':105865199, 'end':105865250},
            #         'J2P':{'start':105864793, 'end':105864852}, #pseudogene, but still need to get rid of germline
            #         'J3':{'start':105864587, 'end':105864635},
            #         'J4':{'start':105864215, 'end':105864260},
            #         'J5':{'start':105863814, 'end':105863862},
            #         'J6':{'start':105863198, 'end':105863258}}
            #         # GRCh38.p7
            return {'J1':{'start':106331617, 'end':106331668},
                    'J2':{'start':106331409, 'end':106331460},
                    'J2P':{'start':106331003, 'end':106331062}, #pseudogene, but still need to get rid of germline
                    'J3':{'start':106330797, 'end':106330845},
                    'J4':{'start':106330425, 'end':106330470},
                    'J5':{'start':106330024, 'end':106330072},
                    'J6':{'start':106329408, 'end':106329468}}
                    # GRCh37
        else:
            print('Only mmu, mmuk and hsa available for now')

    ## Regions / genomes ##
    #
    # def dj(self):
    #     '''Location of DJ genes to remove DJ recombination (seqmonk)
    #     '''
    #     if self.name == 'mmu' or self.name == 'mouse' or self.name == 'mus musculus':
    #         return ['chr12','113428112', '113488696']
    #     elif self.name == 'hsa' or self.name == 'human' or self.name == 'homo sapien':
    #         # return ['chr14','105862474','105920591']
    #         return ['chr14','106327849','106392148'] #hg19
    #     else:
    #         print('Only mmu, mmuk and hsa available for now')

    def genome(self):
        '''Which genome is being used
        '''
        if self.name == 'mmu' or self.name == 'mouse' or self.name == 'mus musculus':
            return 'mm10'
        elif self.name == 'mmuk' or self.name == 'mouse kappa' or self.name == 'mus musculus kappa':
            return 'mm10'
        elif self.name == 'hsa' or self.name == 'human' or self.name == 'homo sapien':
            # return 'hg38'
            return 'hg19'
        else:
            print('Only mmu, mmuk and hsa available for now')

    # def v_region(self):
    #     '''Coordinates of V region (seqmonk)
    #     '''
    #     if self.name == 'mmu' or self.name == 'mouse' or self.name == 'mus musculus':
    #         return 'chr12:113531809-116015193'
    #     elif self.name == 'hsa' or self.name == 'human' or self.name == 'homo sapien':
    #         # return 'chr14:105931103-106903920'
    #         return 'chr14:106394250-107295467' #hg19
    #     else:
    #         print('Only mmu, mmuk and hsa available for now')

    ## Mispriming ##

    def J_for_mispriming(self): #ref created using these
        if self.name == 'mmu' or self.name == 'mouse' or self.name == 'mus musculus':
            return {'J1':'CCCTGTGCCCCAGACATC',
                    'J2':'AGTGGTGCCTTGGCCCCAGTAGTCAAA',
                    'J3':'ACCAGAGTCCCTTGGCCCCAGTAAGCAAA',
                    'J4':'TGAGGTTCCTTGACCCCAGTAGTCCAT'}
        if self.name == 'mmuk' or self.name == 'mouse kappa' or self.name == 'mus musculus kappa':
            return {'J1':'TTTGATTTCCAGCTTGGTGCCTCCACCGA'[4:],
                    'J2':'TTTATTTCCAGCTTGGTCCCCCCTCCGA'[4:],
                    'J4':'CGTTTTATTTCCAACTTTGTCCCCGAGCCGA'[4:],
                    'J5':'CAGCTCCAGCTTGGTCCCAGCA'[4:]}
        elif self.name == 'hsa' or self.name == 'human' or self.name == 'homo sapien':
            return {'J1':'GGTGCCCTGGCCCCAGTGCTGGAA',
                    'J2':'GGTGCCACGGCCCCAGAGATCGAA',
                    'J3':'ACCATTGTCCCTTGGCCCCAGATATCAA',
                    'J4.1':'GACCAGGGTTCCTTGGCCCCAGTAG',
                    'J4.3':'GACCAGGGTCCCTTGGCCCCAGTAG',
                    'J4.1c':'GACCAGGGTTCCTTGGCCCCAGGAGT',
                    'J4.3c':'GACCAGGGTCCCTTGGCCCCAGGAGT',
                    'J5.1':'CAGGGTTCCTTGGCCCCAGGAGTCG',
                    'J5.2':'CAGGGTTCCCTGGCCCCAGGGGTCG',
                    'J6.1':'TGCCCCCAGACGTCCATACCGTAGTAGA',
                    'J6.2':'CCTTTGCCCCAGACGTCCATGTAGTAGTAGA',
                    'J6.3':'TGGCCCCAGACGTCCATACCGTAGTAGA',
                    'J6.4':'TTGCCCCAGACGTCCATACCGTAGTAGA',
                    'J6.c':'TTCCCCCAGACGTCCATACCGTA',
                    'J6.c2':'GGTGCCCTGGCCCCAGACGTCCATACCGTA'} #Misprimed sequence (with c)
        else:
            print('Only mmu, mmuk and hsa available for now')

    #REVIEW: reduce to 4 bp?
    def mispriming(self): #implement auto method to generate this dict
        if self.name == 'mmu' or self.name == 'mouse' or self.name == 'mus musculus':
            return {'J1':'ACATC',
                    'J2':'TCAAA',
                    'J3':'GCAAA',
                    'J4':'TCCAT'}
        elif self.name == 'mmuk' or self.name == 'mouse kappa' or self.name == 'mus musculus kappa':
            return {'J1':'CCTCC',
                    'J2':'TCCGA',
                    'J4':'GCCGA',
                    'J5':'CAGCA',
                    'J5.c':'ACGTG'}
        elif self.name == 'hsa' or self.name == 'human' or self.name == 'homo sapien':
            return {'J1':'TGGAA',
                    'J2':'CGAAG',
                    'J3':'ATCAA',
                    'J4.1':'AGTAG',
                    'J4.3':'GGAGT',
                    'J5.1':'AGTCG',
                    'J5.2':'GGTCG',
                    'J6.1':'GTAGT'}
        else:
            print('Only mmu, mmuk and hsa available for now')

    def offset(self):
        if self.name == 'mmu' or self.name == 'mouse' or self.name == 'mus musculus':
            return {'J1':{'J2':[2], 'J3':[2], 'J4':[2]},
                    'J2':{'J1':[8]},
                    'J3':{'J1':[8]},
                    'J4':{'J1':[8]}}
        elif self.name == 'mmuk' or self.name == 'mouse kappa' or self.name == 'mus musculus kappa':
            return {'J1':{'J1':[10], 'J4':[-2], 'J5':[9]},
                    'J2':{'J2':[5], 'J5':[4]},
                    'J4':{},
                    'J5':{'J2':[1], 'J4':[2], 'J5':[5]},
                    'J5.c':{}}
        elif self.name == 'hsa' or self.name == 'human' or self.name == 'homo sapien':
            return {'J1':{'J2':[4], 'J3':[6], 'J4.1':[1], 'J5.1':[7], 'J5.2':[7]},
                    'J2':{'J2':[4], 'J3':[6], 'J4.1':[1], 'J5.1':[7], 'J5.2':[7], 'J6.1':[-3]},
                    'J3':{'J1':[4], 'J2':[3], 'J4.1':[0], 'J5.1':[6], 'J5.2':[6], 'J6.1':[-4]},
                    'J4.1':{'J1':[0], 'J2':[-1, -2], 'J3':[1], 'J5.1':[2], 'J5.2':[2], 'J6.1':[-8]},
                    'J4.3':{'J1':[0], 'J2':[-1, -2], 'J3':[1], 'J5.1':[2], 'J5.2':[2], 'J6.1':[-8]},
                    'J4.1c':{'J1':[1], 'J2':[0, -1], 'J3':[1], 'J5.1':[3], 'J5.2':[3], 'J6.1':[-7]},
                    'J4.3c':{'J1':[1], 'J2':[0, -1], 'J3':[1], 'J5.1':[3], 'J5.2':[3], 'J6.1':[-7]},
                    'J5.1':{'J1':[3], 'J2':[2], 'J3':[4], 'J4.1':[8], 'J4.3':[7], 'J6.1':[-5]},
                    'J5.2':{'J1':[3], 'J2':[2], 'J3':[4], 'J4.1':[8], 'J4.3':[7], 'J6.1':[-5]},
                    'J6.c':{'J6.1':[0]},
                    'J6.c2':{'J6.1':[0]},
                    'J6.3':{},
                    'J6.4':{}} #TODO:check multiple human germline
        else:
            print('Only mmu, mmuk and hsa available for now')

    def replace(self):
        if self.name == 'mmu' or self.name == 'mouse' or self.name == 'mus musculus':
            return {'J1':'CCCTGTGCCCCAG',
                    'J2':'AGTGGTGCCTTGGCCCCAGTAG',
                    'J3':'ACCAGAGTCCCTTGGCCCCAGTAA',
                    'J4':'TGAGGTTCCTTGACCCCAGTAG'}
        elif self.name == 'mmuk' or self.name == 'mouse kappa' or self.name == 'mus musculus kappa':
            return {'J1':'TTTGATTTCCAGCTTGGTGCCTCC',
                    'J2':'TTTATTTCCAGCTTGGTCCCCCC',
                    'J4':'CGTTTTATTTCCAACTTTGTCCCCGA',
                    'J5':'CAGCTCCAGCTTGGTCC',
                    'J5.c':'CAGCTCCAGCTTGGTCC'}
        elif self.name == 'hsa' or self.name == 'human' or self.name == 'homo sapien':
            return {'J1':'GGTGCCCTGGCCCCAGTGC',
                    'J2':'GGTGCCACGGCCCCAGAGA',
                    'J3':'ACCATTGTCCCTTGGCCCCAGAT',
                    'J4.1':'GACCAGGGTTCCTTGGCCCC',
                    'J4.3':'GACCAGGGTCCCTTGGCCCCA',
                    'J5.1':'CAGGGTTCCTTGGCCCCAGGA',
                    'J5.2':'CAGGGTTCCCTGGCCCCAGGG',
                    'J6.1':'TGCCCCCAGACGTCCATACCGT',
                    'J6.2':'CCTTTGCCCCAGACGTCCATGTAGT',
                    'J6.3':'TGGCCCCAGACGTCCATACCGT',
                    'J6.4':'TTGCCCCAGACGTCCATACCGT'}
        else:
            print('Only mmu, mmuk and hsa available for now')
