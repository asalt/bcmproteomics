import os
import re
from setuptools import setup, find_packages

_version_re = re.compile(r'__version__\s+=\s+(.*)') # from Armin Ronacher

with open('bcmproteomics/__init__.py', 'r') as f:
    version = _version_re.search(f.read()).group(1).strip("'")

def read(fname):
    return open(os.path.join(os.path.dirname(__file__), fname)).read()

setup(
    name='bcmproteomics',
    version=version,
    author='Alex Saltzman',
    author_email='saltzman@bcm.edu',
    description='A library for conveniently grabbing data from BCM_Proteomics iSPEC',
    license='MIT',
    url='https://github.com/asalt/bcmproteomics',
    long_description = read('README.md'),
    packages=find_packages(),
    # install_requires=['pyodbc'],
    keywords=['bioinformatics', 'data wrangling'],
    classifiers=[
        'Development Status :: 3 - Alpha',
        'Topic :: Utilities',
        'License :: OSI Approved :: MIT License',
        'Programming Language :: Python :: 3.4',
        'Programming Language :: Python :: 3.5',
    ],
    package_data={
        #'bcmproteomics': ['training_data/*.txt'],
    },
)
