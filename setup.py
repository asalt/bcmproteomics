import os
from setuptools import setup, find_packages

def read(fname):
    return open(os.path.join(os.path.dirname(__file__), fname)).read()

setup(
    name = 'bcmproteomics',
    version = '0.0.1',
    author = 'Alex Saltzman',
    author_email = 'saltzman@bcm.edu',
    description = 'A library for conveniently grabbing data from BCM_Proteomics iSPEC',
    license = 'MIT',
    url = 'http://www.epicome.org/',
    long_description = read('README.md'),
    packages=find_packages(),
    install_requires=['pandas', 'pyodbc'],
    keywords=['bioinformatics', 'data wrangling'],
    classifiers=[
        'Development Status :: 3 - Alpha',
        'Topic :: Utilities',
        'License :: OSI Approved :: MIT License',
        'Programming Language :: Python :: 3.4'        
    ],    
)


