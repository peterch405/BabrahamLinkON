
'''
BabrahamLinkON: Analysis pipeline for VDJ-seq
Copyright (C) 2017  Peter Chovanec

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <https://www.gnu.org/licenses/>.
'''


from weblogo import *
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
