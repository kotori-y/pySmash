# pySmash: Python package and software for structural alerts generation and application

[![Travis (.com)](https://www.travis-ci.com/kotori-y/pySmash.svg?branch=master)](https://www.travis-ci.com/github/kotori-y/pySmash) ![Codecov](https://img.shields.io/codecov/c/github/kotori-y/pySmash) [![GitHub last commit](https://img.shields.io/github/last-commit/kotori-y/pySmash)](https://github.com/kotori-y/pySmash/commits/master) [![MIT License](https://img.shields.io/badge/license-MIT-black)](https://anaconda.org/kotori_y/scopy) [![Kouhai](https://img.shields.io/badge/contributor-Ziyi-%23B3D0BE)](https://github.com/Yangziyi1997) [![DOI](https://img.shields.io/badge/doi-Briefings%20in%20Bioinformatics-informational)](https://doi.org/10.1093/bib/bbab017)

<div align=center>
    <img src='tutorial/image/pysmash.png' width='200'>
</div>

## Overview

*PySmash*, a user-friendly and integrated tool which offers convenience and advantage for the automatic structural alerts derivation and application is presented. The current version of *PySmash* provides both software and Python package, which achieves both flexibility and ease of pipeline integration. Three kinds of substructure derivation algorithm, including circular, path-based and function group-based algorithms are provided, where researchers can customize own parameters about substructure size, accuracy and coverage, statistical significance and server deployment. Besides, *PySmash* also provides prediction function for derived substructure application. 

## Installation

### Install RDKit

```
>>> conda install -c conda-forge rdkit
```

### Intall pySmash

#### Source

```
>>> git clone git@github.com:kotori-y/pysmash.git && cd pysmash
>>> [sudo] python setup.py install
```

#### PyPI [![PyPI](https://img.shields.io/pypi/v/smasher?style=flat-square)](https://pypi.org/project/smasher/)

```shell
>>> pip install smasher
```

### Software

You can also download the executable program, which need not any dependencies, from https://github.com/kotori-y/pySmash/releases/latest/download/pysmash.tar.gz

## Quick Start

running test

```powershell
# testing circular fragments generation
>>> python ./tests/Carc/testCircular.py

# testing path-based fragments generation
>>> python ./tests/Carc/testPath.py

# testing function group-based fragments generation
>>> python ./tests/Carc/testFg.py
```

## Documentation

(1)The online version of the documentation is available here: https://pysmash.iamkotori.com/<br>(2)Quick start with package: https://pysmash.iamkotori.com/smash_tutorial.html<br>(3)Quick start with software: https://pysmash.iamkotori.com/pysmash_gui.html<br>(3)Application examples(pipelines): https://pysmash.iamkotori.com/Application.html

## Cite us

Yang, Zi-Yi, et al. "PySmash: Python package and individual executable program for representative substructure generation and application." *Briefings in Bioinformatics* (2021).

```
@article{10.1093/bib/bbab017,
    author = {Yang, Zi-Yi and Yang, Zhi-Jiang and Zhao, Yue and Yin, Ming-Zhu and Lu, Ai-Ping and Chen, Xiang and Liu, Shao and Hou, Ting-Jun and Cao, Dong-Sheng},
    title = "{PySmash: Python package and individual executable program for representative substructure generation and application}",
    journal = {Briefings in Bioinformatics},
    year = {2021},
    month = {03},
    abstract = "{Substructure screening is widely applied to evaluate the molecular potency and ADMET properties of compounds in drug discovery pipelines, and it can also be used to interpret QSAR models for the design of new compounds with desirable physicochemical and biological properties. With the continuous accumulation of more experimental data, data-driven computational systems which can derive representative substructures from large chemical libraries attract more attention. Therefore, the development of an integrated and convenient tool to generate and implement representative substructures is urgently needed.In this study, PySmash, a user-friendly and powerful tool to generate different types of representative substructures, was developed. The current version of PySmash provides both a Python package and an individual executable program, which achieves ease of operation and pipeline integration. Three types of substructure generation algorithms, including circular, path-based and functional group-based algorithms, are provided. Users can conveniently customize their own requirements for substructure size, accuracy and coverage, statistical significance and parallel computation during execution. Besides, PySmash provides the function for external data screening.PySmash, a user-friendly and integrated tool for the automatic generation and implementation of representative substructures, is presented. Three screening examples, including toxicophore derivation, privileged motif detection and the integration of substructures with machine learning (ML) models, are provided to illustrate the utility of PySmash in safety profile evaluation, therapeutic activity exploration and molecular optimization, respectively. Its executable program and Python package are available at https://github.com/kotori-y/pySmash.}",
    issn = {1477-4054},
    doi = {10.1093/bib/bbab017},
    url = {https://doi.org/10.1093/bib/bbab017},
    note = {bbab017},
    eprint = {https://academic.oup.com/bib/advance-article-pdf/doi/10.1093/bib/bbab017/36622735/bbab017.pdf},
}
```

## Contact

If you have questions or suggestions, please contact: kotori@cbdd.me,and oriental-cds@163.com.<br>Please see the file LICENSE for details about the "MIT" license which covers this software and its associated data and documents.

## Acknowledgement

*Thanks to my colleague, [Ziyi](https://github.com/Yangziyi1997), for helping me to complete the writing of document and article.*