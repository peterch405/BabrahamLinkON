

# Copyright: (C) 2017 Peter Chovanec <peter.chovanec@babraham.ac.uk>
# Copyright: (C) 2017, BabrahamLinkON
# GNU General Public License v3.0+ (see LICENSE or https://www.gnu.org/licenses/gpl-3.0.txt)


from weblogolib import *
from babrahamlinkon import general
from collections import Counter
from io import StringIO

umis = []

def umi_seq_logo(fastq, out_eps):
    with general.file_open(fastq) as in_fq:
        for item in general.fastq_parse(in_fq):
            qname = item[0]
            seq = item[1]
            thrd = item[2]
            qual = item[3]

            #remove duplicates
            umis.append(qname.split(' ')[0].split('_')[-1])

    deduped_umis = Counter(umis)

    #get UMI in fasta format
    fa = ''
    count = 0
    for umi in deduped_umis.elements():
        fa += '>' + str(count) + '\n' + umi + '\n'
        count += 1


    # fin = open('cap.fa')
    # with open('/media/chovanec/My_Passport/Sync/BabrahamLinkON/tests/logo.fa', 'w') as out_fa:
    #     out_fa.write(fa)

    fa_out = StringIO()
    fa_out.write(fa)

    seqs = read_seq_data(fa_out)

    data = LogoData.from_seqs(seqs)
    options = LogoOptions()
    options.title = "UMI base usage"
    options.unit_name = "probability"
    format = LogoFormat(data, options)
    eps = eps_formatter(data, format)
    # pdf = pdf_formatter(data, format)

    fa_out.close()

    with open(out_eps, 'w') as eps_out:
        eps_out.write(eps.decode())
