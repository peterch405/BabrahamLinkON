=============================================
BabrahamLinkON: Analysis pipeline for VDJ-seq
=============================================

Babrahamlinkon is a tool for the analysis of immunoglobulin receptor
sequences from NGS data generated using the DNA VDJ-seq assay.


------------
Installation
------------

Babrahamlinkon is only compatible with Python 3.

Pre-requisites
===============

Software:
------------------

With bioconda (recommended):
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

`IgBlast 1.7.0 <https://www.ncbi.nlm.nih.gov/igblast/faq.html#standalone>`_

.. code:: bash

  conda install igblast

`Samtools <http://samtools.sourceforge.net/>`_

.. code:: bash

  conda install samtools

`Bowtie 2 <http://bowtie-bio.sourceforge.net/bowtie2/index.shtml>`_

.. code::

  conda install bowtie2


`Kalign2 <http://msa.sbc.su.se>`_

Ubuntu install:

.. code:: bash

  sudo apt-get install kalign

`Pear <http://www.exelixis-lab.org/web/software/pear>`_

Python modules:
---------------

BabrahamLinkON is dependent on:
 * numpy>=1.11.0,
 * pandas>=0.18.1,
 * scikit-bio>=0.5.0,
 * python-Levenshtein>=0.12.0,
 * pysam>=0.9.1.3,
 * joblib>=0.9.3,
 * changeo>=0.3.7,
 * tqdm>=4.13.0.


Enviroment variables:
------------------------------

.. code:: bash

  export BOWTIE2_INDEXES='/path/to/bowtie2/indexes'
  export BOWTIE2_REF='Basename_of_reference'

If running in cluster enviroment:

.. code:: bash

  #Home directory
  export home='/path/to/working/directory'
  #Folder for all the log/output files
  export log_folder=${home}/logs

  #matplotlib backend for headless nodes
  export MPLBACKEND=pdf

  #specify tmp dir (needed for nodes as they don't have much memory)
  export TMPDIR='/state/partition1'


Setup
=====

To install Babrahamlinkon straight from the git repository:

.. code:: bash

  git clone https://github.com/peterch405/BabrahamLinkON
  cd BabrahamLinkON
  pip install .

Basic usage for data with Unique Molecular Identifiers (UMI's)
==============================================================

Precleaning
-----------

.. code:: bash

  preclean.py umi -v <v_end.fastq> -j <j_end_fastq> --species <mmu or hsa or mmuk> --threads <int> --umi_len <int>

Deduplication
-------------

.. code:: bash

  deduplicate.py umi --input_dir <preclean output directory> --stats --threads <int>

Annotation and clone assembly
-----------------------------

.. code:: bash

 assemble_clones.py umi -fa <fasta from deduplication> --full_name --threads <int> --species <mmu or hsa or mmuk>



Running partis
--------------

Partis expects sequences to be input in the VDJ direction. BabrahamLinkON returns reads in the JDV orientation.
To make the fasta/q with partis, simply run:

.. code:: bash

  deduplicate.py reverse_complement --input <fasta/q file or directory of files>

If providing a fastq, use the `--fq` flag.
