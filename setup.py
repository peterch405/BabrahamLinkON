# https://github.com/pypa/sampleproject/blob/master/setup.py
# https://packaging.python.org/distributing/#requirements-for-packaging-and-distributing

import sys
from setuptools import setup, find_packages


if sys.version_info.major != 3:
    raise RuntimeError('BabrahamLinkON requires Python 3')


setup(name='BabrahamLinkON',
      version='0.1',
      description='BabrahamLinkON pipeline for preprocessing and analysing VDJ-seq data',
      author='Peter Chovanec',
      author_email='peter.chovanec@babraham.ac.uk',
      package_dir={'': 'src'},
      packages=['babrahamlinkon'],
      package_data={'babrahamlinkon': ['plot_igh.R', 'plot_v_usage.Rmd']},
      install_requires=[
          'numpy>=1.11.0',
          'pandas>=0.18.1',
          'scikit-bio>=0.5.0',
          'python-Levenshtein>=0.12.0',
          'pysam>=0.9.1.3'],
      scripts=['src/babrahamlinkon/deduplicate.py',
               'src/babrahamlinkon/preclean.py',
               'src/babrahamlinkon/run_mixcr.py',
               'src/babrahamlinkon/germline_mispriming.py'],
)
