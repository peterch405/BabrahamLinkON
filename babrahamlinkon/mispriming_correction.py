
# Copyright: (C) 2017 Peter Chovanec <peter.chovanec@babraham.ac.uk>
# Copyright: (C) 2017, BabrahamLinkON
# GNU General Public License v3.0+ (see LICENSE or https://www.gnu.org/licenses/gpl-3.0.txt)


from collections import defaultdict
import Levenshtein
from babrahamlinkon import presets
import re

try:
    from skbio.alignment import StripedSmithWaterman #put into setup.py install
except:
    print('ERROR: This script requires scikit-bio: pip install scikit-bio')
    sys.exit()



################################################################################
###### J alignment for mispriming ######
################################################################################

class SSW_align:
    '''Align reads using StripedSmithWaterman from scikit-bio
    http://scikit-bio.org/docs/latest/generated/skbio.alignment.StripedSmithWaterman.html#skbio.alignment.StripedSmithWaterman
    '''
    #align_coord = {'J1':0,'J2':27, 'J3':55, 'J4':85}
    def __init__(self):
        self.alignments = []
        # self.print_full = []
        self.ref_overlap = defaultdict(list) #created by reference function
        self.align_coord = defaultdict(int) #created by reference function

    def align(self, target, qry, score=0, no_misprim_cor=False, mismatch_allow=1, spe='mmu', quick_align=False):
        '''Align against reference and return which J sequence aligns to
        :param target: reference
        :param qry: sequence to align against reference
        :param score: score filter, min score to pass as an alignment
        :param no_misprim_cor: do not perform mispriming correction
        :param mismatch_allow: number of mismatches to allow in sequence beyond J primer
        :param spe: which organism
        :param quick_align: skip SSW and mispriming and return initial identification using only levenshtein distance on first 5bps
        :return: read identity (store alignment in alignments)
        '''
        #Query needs to be a single string
        assert isinstance(qry, str), 'Align query should not be a list'

        #use 5 bp to do initial identification
        ini_len = 5
        restart = True
        while restart:
            restart = False
            J_identities = defaultdict(int)
            for J in presets.prs(spe).J_seq().keys():
                dist = Levenshtein.distance(qry[:ini_len], presets.prs(spe).J_seq()[J][:ini_len])
                J_identities[J] = dist #if incorrect will catch later
            #chose J with lowest score
            if J_identities and ini_len <= 10: #if dict not empty and iteration is less then 10bp
                min_dist = min(J_identities.values())
                min_dist_key = [k for k, v in J_identities.items() if v == min_dist]
                if len(min_dist_key) > 1: #if two keys equal score increase size and repeat until only one hit obtained
                    restart = True
                    ini_len += 1
                else:
                    initial_identity = min_dist_key[0]
            else:
                #if all the same keys are haplotypes, take one of them
                min_dist = min(J_identities.values())
                min_dist_key = [k for k, v in J_identities.items() if v == min_dist]
                keys = set()
                for key in min_dist_key:
                    keys.add(key.split('.')[0])
                if len(keys) > 1:
                    return 'other' #no hits
                else:
                    initial_identity = min_dist_key[0]

        if quick_align: #in case only interested in V genes skipp SSW and mispriming
            return [initial_identity]

        #while loop in case initial identity is wrong
        redo = True
        add = 0
        not_added = True #prevent infinit loop
        not_switched = True
        while redo:
            redo = False

            #Use initial identity to align x bp
            align_size = len(presets.prs(spe).J_for_mispriming()[initial_identity]) + add

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
                if read_identity != initial_identity:
                    if not_switched:
                        redo = True
                        initial_identity = read_identity
                        not_switched = False
                        continue
                    else: #initial and after alignment identity didn't match for second time
                        return 'other'
                else:

                    #If mispriming include original read identity in header; output: (before, after, correct_seq)
                    if no_misprim_cor: #if not misprime correcting return identity of primer seq
                        return [read_identity]
                    else:
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

                        #get sequence beyond reference/alignment
                        #know where reference ends so can go 5bp back even if the end did not align / add the full query sequcence back
                        print_seq_full = print_seq + qry[poscounter:] #+ '-'*(len(alignment['target_sequence'])-end-len(qry[poscounter:]))
                        # CCCTGTGCCCCAGACATCNNNNNNAGTGGTGCCTTGGCCCCAGTAGTCAAANNNNNNACCAGAGTCCCTTGGCCCCAGTAAGCAAANNNNNNTGAGGTTCCTTGACCCCAGTAGTCCAT
                        # CCCCTGTGCCCCAGACATCGAAGTACCAGTA raw read
                        # CCTTGTGCCCCAGACATCGAAGTACCAGTA-----------------------------------------------------------------------------------------
                        # print_seq = print_seq + '-'*(len(alignment['target_sequence'])-end)
                        # CCTTGTGCCCCAGACAT------------------------------------------------------------------------------------------------------
                        # self.print_full = print_seq_full

                        ref_end = self.align_coord[read_identity]['end'] #End of J in the reference


                        differences = defaultdict(int)
                        primer_end = defaultdict(int)
                        for key in presets.prs(spe).mispriming().keys():
                            seq_chk_lst = []
                            try: #if in offset dict do ... else do ...
                                #Special case for J1 (rare) in mouse
                                #TODO: don't hard code the len of the sequence to check
                                offst = presets.prs(spe).offset()[read_identity][key]
                                for off in offst:
                                    seq_chk_lst.append(print_seq_full[ref_end-off:ref_end-off+5])
                                    primer_end[key] = ref_end-off
                                # read_offset = True
                            except KeyError: #if don't need to offset using offset dict then just use last 5bp
                                #reference end - 5 : reference end - distance from actual end of reference + last letters that have been cut off
                                seq_chk_lst.append(print_seq_full[ref_end-5:ref_end])
                                primer_end[key] = ref_end-5
                                # read_offset = False
                                # pass

                            diffs = []
                            for seq_chk in seq_chk_lst:
                                diffs.append(Levenshtein.distance(seq_chk, presets.prs(spe).mispriming()[key][:len(seq_chk)]))
                            #Pick smallest diff
                            differences[key] = min(diffs)

                            # if len(seq_chk) >= 4:
                            #     # print('seq_chk', seq_chk)
                            #     differences[key] = (Levenshtein.distance(seq_chk, presets.prs(spe).mispriming()[key][:len(seq_chk)]))
                            # #REVIEW: Can skip this now?
                            # elif not_added: #if sequence beyond primer is shorter than 4bp redo adding x bp to initial alignment
                            #     add = 5-len(seq_chk)
                            #     redo=True
                            #     not_added = False
                            # else: #shouldn't return anything here
                            #     return 'unclear-' + read_identity + '-bpa'

                        if redo: #redo with a longer read if too many insertions present
                            continue

                        min_val = min(differences.values())

                        if min_val <= 0: #allow only 0 mismatches by default in beyond J seq
                            #retrieve key of min val (J)
                            min_val_key = [k for k, v in differences.items() if v == min_val]

                            if len(min_val_key) > 1: #if there are multiple values with same score (shouldn't happen)
                                return 'unclear-' + read_identity + '-mh'
                            else:
                                replace_seq = presets.prs(spe).replace()[min_val_key[0]]

                                corrected_seq = replace_seq + print_seq_full[primer_end[min_val_key[0]]:]
                                if '-' in corrected_seq: #issue for kappa with long stretch of C's
                                    # print('Possible deletion in read', corrected_seq, 'print_seq_full', print_seq_full, 'qry', qry)
                                    corrected_seq = corrected_seq.replace('-', '')
                                    # print(corrected_seq)
                                assert '-' not in corrected_seq, 'Unknown space in corrected sequence'

                                #Will correct qual when writing out fastq (fastqHolder), just trim start or add 'I'
                                # print(read_identity, min_val_key[0], corrected_seq)
                                return [read_identity, min_val_key[0], corrected_seq] #identity of J and corrected sequence

                        else: #if more than allowed mismatches (chewing back?)
                            return 'unclear-' + read_identity + '-chw'
                            # #return 15bp beyond primer for lowest score
                            # min_val_key = [k for k, v in differences.items() if v == min_val]
                            # if len(min_val_key) > 1: #if there are multiple values with same score
                            #     return 'unclear-' + read_identity + '-mh'
                            # else:
                            #     return 'unclear-' + read_identity + '-' + print_seq_full[primer_end[min_val_key[0]]:primer_end[min_val_key[0]]+15]

            else: #doesn't match primer seq within acceptable paramenters
                return 'other'


    def print_align(self, num=100, print_out=True):
        '''Print pretty alignemnt
        :param num: how many alignments to print
        '''
        #Target sequence should be the same for all (reference)
        seqs = []
        if self.alignments:
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



    def reference(self, spe='mmu', quick_align=False):
        '''Create reference from J sequences supplied (will seperate each J by 6 Ns)
        :param spe: which organism
        '''

        ref = ''

        if quick_align:
            return ref

        for key in sorted(presets.prs(spe).J_for_mispriming().keys()):
            start = len(ref)
            ref = ref + presets.prs(spe).J_for_mispriming()[key] + 'N'*6
            end = len(ref)-6
            if start == 0:
                self.align_coord[key] = {'start':start, 'end':end}
            else:
                self.align_coord[key] = {'start':start, 'end':end} #0 based coordinates

            for val in range(start, end):
                self.ref_overlap[val] = key
        ref = ref[:-6]

        return ref
