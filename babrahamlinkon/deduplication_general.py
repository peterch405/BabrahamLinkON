
import itertools
from collections import Counter, defaultdict
import subprocess
from babrahamlinkon import general
import operator
import re
from babrahamlinkon._dedup_umi import edit_distance

def qual_highest(quals_list_of_lists):
    '''if two identical sequences, take highest quality'''
    max_qual = ''
    for pos in zip(*quals_list_of_lists):
        phred = [ord(q)-33 for q in pos]
        max_qual += chr(max(phred)+33)

    return max_qual



def consensus_qual(seqs, qual_dict):
    '''
    :param list_of_lists: string split into individual letters
                          [['C', 'A', 'C', 'A', 'T', 'A', 'T', 'A'],
                          ['G', 'A', 'T', 'A', 'T', 'A', 'T', 'A'],
                          ['G', 'A', 'T', 'A', 'T', 'A', 'T', 'A']]
    If base ambigous return N
    Else return most frequent base
    '''

    quals = [qual_dict.get(seq) for seq in seqs]

    # print(seqs, quals)
    #make a list of lists for each position of seq and qual

    lofls_seqs = [list(item) for item in seqs]
    lof_quals = list(itertools.chain(*quals)) #unlist nested list
    lofls_quals = [list(item) for item in lof_quals]
    # print('lofls_quals', lofls_quals)

    pos_quals = []
    for pos in zip(*lofls_quals):
        pos_quals.append(pos)
    # print('pos_quals', pos_quals)

    bp_pos = 0
    consensus_seq = ''
    consensus_q = ''
    for pos in zip(*lofls_seqs):
        base_counts = Counter(pos)
        most_freq = base_counts.most_common(2)
        if len(most_freq) > 1:
            if most_freq[0][1] == most_freq[1][1]:
                #get positions of the top two bases
                frst = [i for i, x in enumerate(pos) if x == most_freq[0][0]]
                scnd = [i for i, x in enumerate(pos) if x == most_freq[1][0]]
                #get the max score for each of them
                frst_max = max([ord(pos_quals[bp_pos][i])-33 for i in frst])
                scnd_max = max([ord(pos_quals[bp_pos][i])-33 for i in scnd])

                #take nt with highest phred score, if both have same score return N
                if frst_max == scnd_max:
                    consensus_seq += 'N' #ambigous base
                    consensus_q += chr(frst_max+33)
                elif frst_max > scnd_max:
                    consensus_seq += most_freq[0][0]
                    consensus_q += chr(frst_max+33)
                else:
                    consensus_seq += most_freq[1][0]
                    consensus_q += chr(scnd_max+33)
            else:
                consensus_seq += most_freq[0][0]
                position_of_cons = [i for i, x in enumerate(pos) if x == most_freq[0][0]]
                phred_max = max([ord(pos_quals[bp_pos][i])-33 for i in position_of_cons])
                consensus_q += chr(phred_max+33)

        else:
            consensus_seq += most_freq[0][0]
            position_of_cons = [i for i, x in enumerate(pos) if x == most_freq[0][0]]
            #TODO: do something more sophisticated than max qual
            phred_max = max([ord(pos_quals[bp_pos][i])-33 for i in position_of_cons])
            consensus_q += chr(phred_max+33)

        bp_pos += 1

    return (consensus_seq, consensus_q)


# seqs = ['ATCGATA', 'CTGGATA']
#
# qual_dict = {'ATCGATA':'@><@JAA', 'CTGGATA':'KA""@>>'}



def consensus(list_of_lists):
    '''
    :param list_of_lists: string split into individual letters
                          [['C', 'A', 'C', 'A', 'T', 'A', 'T', 'A'],
                          ['G', 'A', 'T', 'A', 'T', 'A', 'T', 'A'],
                          ['G', 'A', 'T', 'A', 'T', 'A', 'T', 'A']]
    If base ambigous return N
    Else return most frequent base
    '''

    consensus_seq = ''
    # count = 0
    for pos in zip(*list_of_lists):
        base_counts = Counter(pos)
        most_freq = base_counts.most_common(2)
        if len(most_freq) > 1:
            if most_freq[0][1] == most_freq[1][1]:
                consensus_seq += 'N' #ambigous base
            else:
                consensus_seq += most_freq[0][0]
        else:
            consensus_seq += most_freq[0][0]

    return consensus_seq



def consensus_unequal(list_of_lists):
    '''
    :param list_of_lists: string split into individual letters
                          [['C', 'A', 'C', 'A', 'T', 'A', 'T', 'A'],
                          ['G', 'A', 'T', 'A', 'T', 'A', 'T', 'A'],
                          ['G', 'A', 'T', 'A', 'T', 'A', 'T', 'A']]
    If base ambigous return N
    Else return most frequent base
    '''

    consensus_seq = ''
    # count = 0
    for pos in itertools.zip_longest(*list_of_lists):
        base_counts = Counter(pos)
        most_freq = base_counts.most_common(2)
        if len(most_freq) > 1:
            if most_freq[0][1] == most_freq[1][1]:
                consensus_seq += 'N' #ambigous base
            elif most_freq[0][0] == None:
                #skip position
                continue
            else:
                consensus_seq += most_freq[0][0]
        elif most_freq[0][0] == None:
            continue
        else:
            consensus_seq += most_freq[0][0]

    return consensus_seq



def kalign_msa(seq_counter_dict, qual_dict=None): #, umi=None
    '''Multiple sequence alignment for read loss analysis
    :param dict seq_counter_dict: dict of umi:Counter(sequences) object with sequences
    :param qual_dict: perserve fq quality of aligned reads
    '''
    #Only aligning single copy of a duplicate sequence and then subsequently
    #multiplying the output alignment by Counter

    #Convert Counter into fasta input for kalign
    #Need to sort beacuse of Instability in progressive multiple sequence alignment algorithms?
    seq_fasta = ''
    seq_qual = defaultdict(list)
    out_qual = defaultdict(list)
    count = 0
    reads = 0
    for umi, cntr in sorted(seq_counter_dict.items(), key=operator.itemgetter(0)):
        for seq, freq in sorted(cntr.items(), key=lambda x:x[0]):
            seq_fasta = seq_fasta + '>' + str(count) + '_' + umi + '_' + str(freq) + '\n' + seq + '\n'
            if qual_dict:
                seq_qual['>' + str(count) + '_' + umi + '_' + str(freq)].append(qual_dict.get(seq))
            count += 1  #how many different reads
            reads += freq #how many total reads

    if qual_dict:
        #collapse qual for same reads
        for k,v in seq_qual.items():
            if len(v) > 1:
                qual_lofls = [list(item) for item in v]
                seq_qual[k] = [qual_highest(qual_lofls)]

    #Can't get consensus from single sequence, return single sequence
    if count == 1:

        assert len(list(seq_counter_dict.values())) == 1, 'Not a single sequence'
        seq_out= re.sub('\n', '', seq_fasta)
        seq_out = list(filter(None, re.split('(>\d+_[A-Z]+_\d+)', seq_out)))
        if qual_dict:
            out_seq_qual = seq_qual.get(seq_out[0])[0]
            out_qual[seq_out[1]].extend(out_seq_qual)

        return (seq_out, out_qual)

    #Multiple sequence alignment
    #http://msa.sbc.su.se/cgi-bin/msa.cgi
    kalign_cmd = ['kalign', '-f', 'fasta'] #'-s', '80.0', '-e', '3.0', '-t', '3.0', '-m', '0.0'


    p = subprocess.Popen(kalign_cmd, stdout=subprocess.PIPE, stdin=subprocess.PIPE, stderr=subprocess.STDOUT)
    kalign_stdout = p.communicate(input=seq_fasta.encode('utf-8'))[0]

    #parse kalign output and derive consensus seq
    head, sep, tail = kalign_stdout.partition(b'>') #remove head kalign intro

    alignment = sep.decode() + tail.decode()
    alignment = re.sub('\n', '', alignment)
    alignment = list(filter(None, re.split('(>\d+_[A-Z]+_\d+)', alignment)))

    # assert len(alignment) == count*2, 'Alignment output fewer reads'

    #fill spaces in quality if space in alignment introduced -

    for item in general.fasta_parse(alignment):
        qname = item[0]
        seq = item[1]
        if qual_dict:
            qual = seq_qual.get(qname)[0]

            new_qual = add_at_hypens(seq, qual[0], '#')
            out_qual[seq].append(new_qual)

    return (alignment, out_qual)


def add_at_hypens(hypen_str, str_add, sym_to_add):
    '''add_at_hypens('ATATGAGATA-A', '>>>>>>>>>>>', '#')
    produces '>>>>>>>>>>#>'
    '''
    fill_pos = [i for i, x in enumerate(hypen_str) if x == '-']
    count = 0
    added = 0
    new_str = ''
    for p in range(len(str_add)+len(fill_pos)):
        if count in fill_pos:
            new_str += sym_to_add
            added += 1
        else:
            new_str += str_add[p-added]
        count += 1
    return new_str

# add_at_hypens('ATATGAGATA-A', 'AAAAAAAAAAA', '#')



def get_best_higher_counts(cluster, counts):
    ''' return the UMI with the highest counts'''
    count = 0

    if len(cluster) == 1:
        return list(cluster)[0]
    else:
        sorted_nodes = sorted(cluster, key=lambda x: counts[x], reverse=True)
        return sorted_nodes[0]



#TODO fill N's with highest quality base
def fill_Ns(N_positions, seq_diffs, qual_dict):
    '''
    :param N_positions: which positions in seq contain N's
    :param seq_diffs: how many differences from consensus
    :param qual_dict: quality dictionary with seq:qual
    '''
    replace_nts = defaultdict()
    replace_qual = defaultdict()

    qual_dict

    #get best sequence (least diffs from consensus)
    seq_with_least_diff = sorted(seq_diffs, key=seq_diffs.get, reverse=False)

    for pos in N_positions:
        filling = True
        seq_to_use = 0
        while filling:

            if seq_with_least_diff[seq_to_use][pos] != '-':
                replace_nts[pos] = seq_with_least_diff[seq_to_use][pos]
                qual = qual_dict.get(seq_with_least_diff[seq_to_use], '')
                if qual != '':
                    replace_qual[pos] = qual[0][pos]
                filling = False
            else:
                #if run out of sequences to correct the N with (should not happen)
                if len(seq_with_least_diff) <= seq_to_use:
                    filling = False
                    replace_nts[pos] = '-'
                    replace_qual[pos] = '#'
                seq_to_use += 1

    return replace_nts, replace_qual




def read_loss(seq_counter_dict, qual_dict, differences=5, no_msa=False, short=False,
              cons_no_qual=False, j_trim=25, with_N=False): #, umi=None
    '''Read loss analysis
    :param seq_counter_dict: {'AATATA':Counter(['AGACCACTGCTCCTAGAAGCACAGAAGTG',
                                                'AGACCACTGCTCCTAGAAGCACAGAAGTG'])}
                             or fasta from kalign_msa
    :param qual_dict: {'AGACCACTGCTCCTAGAAGCACAGAAGTG':['IIIIIIIIIIIIIIIIIIIIIIIIIIIII']}
    :param differences: number of differences from consensus allowed
    :param no_msa: the input seq_counter_dict is not output from kalign_msa
    :param short: are sequences short? will do additional trimming to shortest sequences
    :param cons_no_qual: don't do quality consensus
    :param j_trim: trim off J primers so in case of mispriming it won't result in mismatches
                   (gives misprimed sequences a chance to correct J)
    '''

    good = 0
    total = 0
    diffs_from_cons = []
    seq_diffs = defaultdict()

    if no_msa:

        list_of_lists = []
        for umi, seqs in seq_counter_dict.items():
            for seq in seqs.elements():
                list_of_lists.append(list(seq))

        # cons_seq = consensus(seq_list_lsts)
        cons_no_qual = True
        cons_seq = consensus_unequal(list_of_lists)
        cons_qual = ''

        #Skip poor consensus (i.e. if building consensus from only two different seqs)
        if cons_seq.count('N') > differences: #5
            return(0, cons_seq, cons_qual, '0')
        if len(cons_seq) < 25:
            return(0, cons_seq, cons_qual, '0')

        #get positions of N's
        if 'N' in cons_seq:
            N_positions = [i for i in range(0,len(cons_seq)) if cons_seq[i] == 'N']

        #umi:Counter(seqs)
        for umi, seqs in seq_counter_dict.items():
            for seq in seqs.elements():
                total += 1

                if short: #trim sequence to shortest one
                    min_len = min(len(seq),len(cons_seq))
                    cons_seq_el = cons_seq[:min_len]
                    seq = seq[:min_len]

                #Need to pad seq if length unequal!
                elif len(seq) != len(cons_seq):
                    len_diff = len(cons_seq) - len(seq)
                    if len_diff < 0: #pad cons_seq_1
                        cons_seq_el = cons_seq + '-'*abs(len_diff)
                    elif len_diff > 0:
                        seq = seq + '-'*len_diff
                        cons_seq_el = cons_seq
                else:
                    cons_seq_el = cons_seq

                assert len(seq) == len(cons_seq_el), 'Length of sequences into hamming distance unequal'

                consensus_diff = edit_distance(seq[j_trim:].encode('utf-8'), cons_seq_el[j_trim:].encode('utf-8')) #how many mismatches present
                diffs_from_cons.append(str(consensus_diff))

                seq_diffs[seq] = consensus_diff

                if consensus_diff <= differences:
                    good += 1

    else: #performed msa (default)

        seq_dict = defaultdict(list)

        for item in general.fasta_parse(seq_counter_dict):
            qname = item[0]
            seq = item[1]

            freq = int(qname.split('_')[-1]) #freq saved in name
            for i in range(freq):
                seq_dict[qname+str(i)] = seq #each seq has unique key

        #all values in list of lists

        if cons_no_qual:
            seq_list_lsts = [list(item) for item in seq_dict.values()]
            cons_seq = consensus(seq_list_lsts)
            cons_qual = ''
        else:
            cons_seq, cons_qual = consensus_qual(list(seq_dict.values()), qual_dict)

        #Skip poor consensus (i.e. if building consensus from only two different seqs)
        if cons_seq.count('N') > differences: #5
            return(0, cons_seq, cons_qual, '0')
        if len(cons_seq) < 25:
            return(0, cons_seq, cons_qual, '0')

        #get positions of N's

        if 'N' in cons_seq:
            N_positions = [i for i in range(0,len(cons_seq)) if cons_seq[i] == 'N']

        if short:
            seq_lengths = []
            for qname, seq in general.fasta_parse(seq_counter_dict):
                seq_lengths.append(len(seq))


        for qname, seq in general.fasta_parse(seq_counter_dict):

            freq = int(qname.split('_')[-1])
            total += 1

            if short: #trim sequence to shortest one
                min_len = min(seq_lengths)
                assert len(seq[:min_len]) == len(cons_seq[:min_len]), \
                'Length of sequences into hamming distance unequal'
            else:
                assert len(seq) == len(cons_seq), \
                'Length of sequences into hamming distance unequal'
                min_len = len(cons_seq)

            #how many mismatches present
            consensus_diff = edit_distance(seq[j_trim:min_len].encode('utf-8'),
                                           cons_seq[j_trim:min_len].encode('utf-8'))
            diffs_from_cons.append(','.join([str(consensus_diff)]*freq))
            #best_seq is the seq with fewest differences from consensus seq
            seq_diffs[seq] = consensus_diff

            if consensus_diff <= differences:
                good += 1

            # test = False
            # if 'AGTAATAAGCCGTCCT' in qname:
            #     print(seq, cons_seq, min_len, j_trim, consensus_diff)
            #     test = True



    #replace N's in cons_seq by values in best_seq
    if cons_no_qual:

        if with_N or 'N' not in cons_seq:
            assert len(cons_seq) > 0, "consensus str is empty"
            return (good/total, cons_seq, '', diffs_from_cons)
        else:
            replace_nts, replace_qual = fill_Ns(N_positions, seq_diffs, qual_dict)
            cons_str = list(cons_seq)
            for k, v in replace_nts.items():
                cons_str[k] = v

            new_str = ''.join(cons_str)
            #0 bad 1 good
            assert len(new_str) == len(cons_seq), "new str is not same length as cons seq"
            if test:
                print(new_str)
            return (good/total, new_str, '', diffs_from_cons)



    else:

        if with_N or 'N' not in cons_seq:
            assert len(cons_seq) > 0, "consensus str is empty"
            return (good/total, cons_seq, cons_qual, diffs_from_cons)
        else:
            #0 bad 1 good
            replace_nts, replace_qual = fill_Ns(N_positions, seq_diffs, qual_dict)
            cons_str = list(cons_seq)
            qual_str = list(cons_qual)
            for k, v in replace_nts.items():
                cons_str[k] = v
            for k, v in replace_qual.items():
                qual_str[k] = v

            new_str = ''.join(cons_str)
            new_qual = ''.join(qual_str)

            # assert new_str != '', "new str is empty"
            assert len(new_str) == len(cons_seq), "new str is not same length as cons seq"
            assert len(new_str) == len(new_qual), "new str is not same length as new qual"
            return (good/total, new_str, new_qual, diffs_from_cons)



# import pickle

#
# def rcs_check(func):
#     def func_wrapper(bundle, clusters, counts, stats, mismtch, gt_threshold, msa, short):
#         count = len(glob.glob('/media/chovanec/My_Passport/test/ratio_file*'))
#         with open('/media/chovanec/My_Passport/test/ratio_file' + str(count) + '.pkl', 'wb') as out_pkl:
#             reads, consensus, final_umis, umi_counts, low_gt, corrected, low_gt_corrected, gt_list = func(bundle, clusters, counts, stats, mismtch, gt_threshold, msa=False, short=False)
#             pickle.dump(gt_list, out_pkl)
#         return reads, consensus, final_umis, umi_counts, low_gt, corrected, low_gt_corrected
#     return func_wrapper
#
# @rcs_check
def reduce_clusters_single(bundle, clusters, counts, mismtch, gt_threshold,
                           qual_dict, j_trim, no_msa=False, short=False,
                           cons_no_qual=False, with_N=False):
    ''' collapse clusters down to the UMI which accounts for the cluster
    and return the list of final UMIs using consensus sequence'''

    reads = []
    consensus_seqs = []
    consensus_quals = []
    final_umis = []
    umi_counts = []
    low_gt = 0
    corrected = 0
    low_gt_corrected = 0

    cons_diffs = defaultdict()
    cons_algn = defaultdict()

    for cluster in clusters:
        # total += 1
        umi_in_cluster = len(cluster)
        #Consensus for read loss filter
        out_dict = {umi:bundle[umi]['seq'] for umi in cluster}


        assert len(out_dict) > 0, 'No sequence from umi'

        if umi_in_cluster > 1: #contains umi with 1 error
            corrected += 1

        if no_msa:
            alignment = ''
            gt_ratio, consensus_seq, consensus_qual, diffs_from_cons = read_loss(out_dict, qual_dict, differences=mismtch,
                                                                                short=short, no_msa=no_msa, cons_no_qual=cons_no_qual,
                                                                                j_trim=j_trim, with_N=with_N) #umi=umi_cons
        else:
            alignment, new_qual_dict = kalign_msa(out_dict, qual_dict)
            gt_ratio, consensus_seq, consensus_qual, diffs_from_cons = read_loss(alignment, new_qual_dict, differences=mismtch,
                                                                                 no_msa=no_msa, cons_no_qual=cons_no_qual,
                                                                                 j_trim=j_trim, with_N=with_N) #umi=umi_cons
        #keep record of the distance between sequence and consensus per UMI bases (clustered UMIs seperated by ,)
        if not isinstance(diffs_from_cons, int):
            cons_algn[','.join(cluster)] = ','.join(x for x in alignment)
            cons_diffs[','.join(cluster)] = ','.join(x for x in diffs_from_cons)
        else:
            cons_algn[','.join(cluster)] = diffs_from_cons #key = each cluster umis
            cons_diffs[','.join(cluster)] = diffs_from_cons # values = fasta of differences

        # gt_list.append(gt_ratio)
        if gt_ratio >= gt_threshold:
            #Parent umi = highest count umi which account for the cluster
            parent_umi = get_best_higher_counts(cluster, counts)
            reads.append(bundle[parent_umi]['read'])
            consensus_seqs.append(consensus_seq.replace('-', '')) #remove padding or indels from msa
            consensus_quals.append(consensus_qual.replace('#', '')) #should not have any # qual as preclean removed them

            final_umis.append(parent_umi)
            #Number of UMI's in the cluster (how many have been collapsed)
            umi_counts.append(sum([counts[x] for x in cluster]))
        else:
            low_gt += 1

            if umi_in_cluster > 1: #contains umi with 1 error and low ratio
                low_gt_corrected += 1

    assert len(set(reads)) == len(reads), 'Not all reads unique!'

    #list of reads, final umi's used, list of umi counts within clusters
    return reads, consensus_seqs, consensus_quals, final_umis, umi_counts, low_gt, corrected, low_gt_corrected, cons_diffs, cons_algn#, gt_list
