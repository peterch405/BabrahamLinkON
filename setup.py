# https://github.com/pypa/sampleproject/blob/master/setup.py
# https://packaging.python.org/distributing/#requirements-for-packaging-and-distributing

import sys
from setuptools import setup, find_packages
import subprocess
import re

if sys.version_info.major != 3:
    raise RuntimeError('BabrahamLinkON requires Python 3')

try:
    import Cython
except ImportError:
    raise ImportError(
        "BabrahamLinkON requires cython to be installed before running setup.py (pip install cython)")

# check kalign is accessible
try:
    subprocess.check_output(['kalign', '-q'])
except OSError:
    print('kalign not found, some options won\'t work\n')

# check bowtie2 and samtools is accessible
try:
    subprocess.check_output(['bowtie2', '-h'])
    subprocess.check_output(['samtools', '--help'])
except OSError:
    raise RuntimeError('bowtie2/samtools not found; put directory in $PATH\n')

try:
    subprocess.check_output(['MakeDb.py', '-h'])
except OSError:
    print('MakeDb.py from Changeo not found. Is changeo installed? (pip install changeo, \
    https://bitbucket.org/kleinstein/changeo/downloads/)\nSome options won\'t work')


#check igblast and version_info
try:
    version = subprocess.check_output(['igblastn', '-version'])
    version = version.decode('utf-8').split('\n')[1]
    major, minor, micro = re.split(',| ', version)[3].split('.')
    if int(major) < 1 and int(minor) < 5:
        raise Exception('IgBlast version 1.5.0 or higher required')
except OSError:
    print('IgBlast not found. Some scripts won\'t work')


from Cython.Build import cythonize



setup(name='BabrahamLinkON',
      version='0.1',
      description='BabrahamLinkON pipeline for preprocessing and analysing VDJ-seq data',
      author='Peter Chovanec',
      author_email='peter.chovanec@babraham.ac.uk',
      packages=['babrahamlinkon'],
      package_dir={'babrahamlinkon': 'babrahamlinkon'},
    #   include_package_data=True,
    #   package_data={'babrahamlinkon': ['plot_igh.R', 'plot_v_usage.Rmd']},
      package_data={'babrahamlinkon': ['resources/IgBlast_database/Mus*', 'resources/IgBlast_database/Homo*',
                                       'resources/IgBlast_database/optional_file/*',
                                       'resources/igh_genes_GRCm38.p4',
                                       'resources/igh_genes_GRCh37']},
      install_requires=[
          'numpy>=1.11.0',
          'pandas>=0.18.1',
          'scikit-bio>=0.5.0',
          'python-Levenshtein>=0.12.0',
          'pysam>=0.9.1.3',
          'joblib>=0.9.3'],
      scripts=['babrahamlinkon/deduplicate.py',
               'babrahamlinkon/preclean.py',
               'babrahamlinkon/assemble_clones.py'],
      classifiers = [
		'Environment :: Console',
		'Intended Audience :: Science/Research',
		'License :: ',
		'Natural Language :: English',
        'Operating System :: OS Independent',
		'Programming Language :: Python :: 3',
		'Topic :: Scientific/Engineering :: Bio-Informatics'],
      ext_modules=cythonize('babrahamlinkon/_dedup_umi.pyx')
)
