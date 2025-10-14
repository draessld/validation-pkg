# validation-pkg

## Overview

The **Validation Package** is designed to check the integrity, format, and encoding of genomic data files used in all follow-up pipelines. It supports multiple file types (FASTA, FASTQ, BAM, GTF, GFF, GBK, BED) and can handle compressed files (gzip, bzip2 and tgz). The module is extensible, allowing for custom validators for each file type. Filetypes to proccess are **Genome, Reads & Feature**

---

## Directory Structure
```
validation_pkg
├── config_guide.md
├── LICENSE
├── README.md
├── requirements.txt
├── setup.py
├── tests
│   ├── fixtures
│   │   ├── ID_testcase
│   │   └── ...
│   ├── __init__.py
│   ├── test_coordinator.py
│   ├── test_exceptions.py
│   └── test_file_processing.py
└── validation_pkg
    ├── __init__.py
    ├── coordinator.py
    ├── exceptions.py
    ├── statistics.py
    ├── editors
    │   └── __init__.py
    ├── parsers
    │   ├──  __init__.py
    │   ├── fasta_parser.py
    │   └── genbank_parser.py
    ├── utils
    │   └── file_handler.py
    └── validators
        ├── __init__.py
        ├── feature_validator.py
        ├── genome_validator.py
        └── read_validator.py
```
---

##  Main Components

### Configuration handling

### Exceptions handling

### Coding/Encoding handling
