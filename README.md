
# BabrahamLinkON: Analysis pipeline for VDJ-seq


Babrahamlinkon is a tool for the analysis of immunoglobulin receptor
sequences from NGS data generated using the DNA VDJ-seq assay.

#### Relevant publications:
-----------------------------

Chovanec, P., Bolland, D.J., Matheson, L.S., Wood, A.L., Corcoran, A.E. (2018). Unbiased quantification of
immunoglobulin diversity at the DNA level with VDJ-seq. [Nat. Protoc. 13, 1232–1252.](https://doi.org/10.1038/nprot.2018.021)

Matheson, L.S., Bolland, D.J., Chovanec, P., Krueger, F., Andrews, S., Koohy, H., and Corcoran, A. (2017). Local
chromatin features including PU.1 and IKAROS binding and H3K4 methylation shape the repertoire of
immunoglobulin kappa genes chosen for V(D)J recombination. [Front. Immunol. 8, 1550.](https://doi.org/10.3389/fimmu.2017.01550)

Bolland, D.J., Koohy, H., Wood, A.L., Matheson, L.S., Krueger, F., Stubbington, M.J.T., Baizan-Edge, A., Chovanec, P.,
Stubbs, B.A., Tabbada, K., Andrews, S.R., Spivakov, M., Corcoran, A.E. (2016). Two Mutually Exclusive Local
Chromatin States Drive Efficient V(D)J Recombination. [Cell Rep. 15, 2475–2487.](https://doi.org/10.1016/j.celrep.2016.05.020)


### Installation

Babrahamlinkon is only compatible with Python 3.


## Pre-requisites

### Software:


With [bioconda](https://bioconda.github.io/) (recommended) or follow tool specific instructions available on their website:


[IgBlast 1.7.0](https://www.ncbi.nlm.nih.gov/igblast/faq.html#standalone)

```bash
  conda install igblast
```

[Samtools](http://samtools.sourceforge.net/)

```bash
  conda install samtools
```

[Bowtie 2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml)

```bash
  conda install bowtie2
```

[Kalign2](http://msa.sbc.su.se)

Ubuntu install:

```bash
  sudo apt-get install kalign
```

[Pear](http://www.exelixis-lab.org/web/software/pear)


### Python modules:

BabrahamLinkON is dependent on:
 * numpy>=1.11.0,
 * pandas>=0.18.1,
 * scikit-bio>=0.5.0,
 * python-Levenshtein>=0.12.0,
 * pysam>=0.9.1.3,
 * joblib>=0.9.3,
 * changeo>=0.3.7,
 * tqdm>=4.13.0,
 * weblogo>=3.7.9.


Installation time with all dependencies: ~5 minutes

### Enviroment variables:

```bash
  export BOWTIE2_INDEXES='/path/to/bowtie2/indexes'
  export BOWTIE2_REF='Basename_of_reference'
```

If running in cluster enviroment:

```bash

  #Home directory
  export home='/path/to/working/directory'
  #Folder for all the log/output files
  export log_folder=${home}/logs

  #matplotlib backend for headless nodes
  export MPLBACKEND=pdf

  #specify tmp dir (needed for nodes as they don't have much memory)
  export TMPDIR='/state/partition1'
```

## Setup

I would recommend installing BabrahamLinkON within its own virtual enviroment:

```bash
conda env create -f environment.yaml
conda activate babrahamlinkon
```

To install Babrahamlinkon straight from the git repository:

```bash

  git clone https://github.com/peterch405/BabrahamLinkON
  cd BabrahamLinkON
  pip install .
```

## Basic usage for data with Unique Molecular Identifiers (UMI's)


### Precleaning


```bash
  preclean.py umi -v <v_end.fastq> -j <j_end_fastq> --species <mmu or hsa or mmuk> --threads <int> --umi_len <int>
```

### Deduplication

```bash
  deduplicate.py umi --input_dir <preclean output directory> --stats --threads <int>
```


### Annotation and clone assembly

```bash
 assemble_clones.py umi -fa <fasta from deduplication> --full_name --threads <int> --species <mmu or hsa or mmuk>
```

### Running partis


Partis expects sequences to be input in the VDJ direction. BabrahamLinkON returns reads in the JDV orientation.
To make the fasta/q compatible with partis, simply run:

```bash
  deduplicate.py reverse_complement --input <fasta/q file or directory of files>
```

If providing a fastq, use the `--fq` flag.


## Test dataset


A small dataset can be found in the test folder. This can be used to test your installation:

```bash
 . run_test
```

The expected output is in `expected_test_output` folder

Run time for test data on a i7-4790 running on all 8 threads: ~9 minutes
