# Validation Package

A comprehensive validation package for genomic data files used in bioinformatics pipelines. Validates and standardizes genome files, sequencing reads, and feature annotations with detailed error reporting and quality checks.

---

## Features

- **Multi-format support:** FASTA, FASTQ, BAM, GTF, GFF, GBK, BED
- **Compression handling:** Automatic detection and processing of gzip (.gz) and bzip2 (.bz2) files
- **Quality checks:** Sequence validation, duplicate ID detection, invalid character checks, empty sequence detection
- **Format conversion:** GenBank to FASTA, BAM to FASTQ, BED to GFF
- **Plasmid handling:** Automatic plasmid splitting for bacterial genomes
- **Flexible API:** Multiple usage patterns from simple to advanced
- **Configuration-driven:** JSON-based configuration for batch processing
- **Detailed reporting:** Structured validation reports with statistics
- **CLI support:** Command-line interface for automation

---

## Quick Start

### Installation

```bash
# Install from source
cd /path/to/validation_pkg
pip install -e .

# Or use directly
python -m validation_pkg validate config.json
```

### Simple Usage

```python
from validation_pkg import ValidationCoordinator

# Validate all files in config
coordinator = ValidationCoordinator("config.json")
report = coordinator.validate_all()

if report.passed:
    print("✓ All validations passed!")
else:
    print(f"✗ Found {report.errors} error(s)")
```

### Configuration Example

Create a `config.json` file:

```json
{
  "output_dir": "./output",
  "ref_genome_filename": {
    "filename": "reference.fasta"
  },
  "mod_genome_filename": {
    "filename": "modified.fasta.gz"
  },
  "reads": [
    {
      "filename": "sample1.fastq.gz",
      "ngs_type": "illumina"
    }
  ],
  "ref_feature_filename": {
    "filename": "annotations.gff"
  }
}
```

See [USAGE.md](USAGE.md) for detailed usage patterns and examples.

---

## Usage Patterns

The package supports multiple usage patterns to fit different needs:

### Pattern 1: Simple Workflow (Recommended)

Use `ValidationCoordinator` for automatic workflow orchestration:

```python
from validation_pkg import ValidationCoordinator

coordinator = ValidationCoordinator("config.json")
report = coordinator.validate_all()
print(report.summary())
```

### Pattern 2: Selective Validation

Validate only specific file types:

```python
coordinator = ValidationCoordinator("config.json")

# Validate only genomes
genome_report = coordinator.validate_genomes()

# Validate only reads
reads_report = coordinator.validate_reads()

# Validate only features
features_report = coordinator.validate_features()
```

### Pattern 3: Custom Settings (Functional API)

Use custom settings for each file type:

```python
from validation_pkg import (
    ConfigManager,
    GenomeValidator,
    ReadValidator,
    FeatureValidator,
    validate_genome,
    validate_reads,
    validate_feature
)

# Load configuration
config = ConfigManager.load("config.json")

# Create custom settings
ref_settings = GenomeValidator.Settings()
ref_settings = ref_settings.update(
    plasmid_split=True,
    coding_type='gz',
    output_filename_suffix='ref',
    replace_id_with='chr',
    min_sequence_length=500
)

reads_settings = ReadValidator.Settings()
reads_settings = reads_settings.update(
    coding_type='gz',
    check_invalid_chars=True
)

# Validate with custom settings
if config.ref_genome:
    stats = validate_genome(config.ref_genome, config.output_dir, ref_settings)
    print(f"Validated {stats['total_sequences']} sequences")

if config.reads:
    stats_list = validate_reads(config.reads, config.output_dir, reads_settings)
```

See [simple_usage.py](simple_usage.py) and [example_usage.py](example_usage.py) for complete examples.

---

## API Reference

### Main Classes

#### `ConfigManager`

Load and parse configuration files:

```python
config = ConfigManager.load("config.json")
print(config.output_dir)
print(config.ref_genome.filename)
print(len(config.reads))
```

#### `ValidationCoordinator`

Orchestrate complete validation workflows:

```python
coordinator = ValidationCoordinator("config.json")
report = coordinator.validate_all()
```

#### Validator Classes

Each validator has a `Settings` class for customization:

- `GenomeValidator` - Validate genome files (FASTA, GenBank)
- `ReadValidator` - Validate sequencing reads (FASTQ, BAM)
- `FeatureValidator` - Validate feature annotations (GFF, GTF, BED)

**Example Settings:**

```python
# Genome settings
genome_settings = GenomeValidator.Settings()
genome_settings = genome_settings.update(
    allow_empty_sequences=False,
    allow_empty_id=False,
    warn_n_sequences=2,
    plasmid_split=True,
    replace_id_with='chr',
    min_sequence_length=100,
    coding_type='gz',
    output_filename_suffix='ref'
)

# Read settings
read_settings = ReadValidator.Settings()
read_settings = read_settings.update(
    check_invalid_chars=True,
    allow_empty_id=False,
    allow_duplicate_ids=False,
    keep_bam=False,
    coding_type='gz'
)

# Feature settings
feature_settings = FeatureValidator.Settings()
feature_settings = feature_settings.update(
    sort_by_position=True,
    check_coordinates=True,
    allow_zero_length=False,
    coding_type='gz'
)
```

### Functional API

Simplified functions for direct validation:

- `validate_genome(genome_config, output_dir, settings=None)` - Validate a genome file
- `validate_read(read_config, output_dir, settings=None)` - Validate a single read file
- `validate_reads(read_configs, output_dir, settings=None)` - Validate multiple read files
- `validate_feature(feature_config, output_dir, settings=None)` - Validate a feature file

---

## Command Line Interface

```bash
# Validate all files
python -m validation_pkg validate config.json

# Validate only genomes
python -m validation_pkg validate config.json --only genomes

# Validate with verbose output
python -m validation_pkg validate config.json --verbose

# Specify output directory
python -m validation_pkg validate config.json --output-dir ./custom_output

# Write logs to file
python -m validation_pkg validate config.json --log-file validation.log
```

---

## Supported File Formats

### Genome Files
- **FASTA:** `.fasta`, `.fa`, `.fna`
- **GenBank:** `.gb`, `.gbk`, `.genbank`
- **Compression:** `.gz`, `.bz2`

### Read Files
- **FASTQ:** `.fastq`, `.fq`
- **BAM:** `.bam`
- **Compression:** `.gz`, `.bz2`

### Feature Files
- **GFF:** `.gff`, `.gff3`, `.gtf`
- **BED:** `.bed`
- **Compression:** `.gz`, `.bz2`

---

## Validation Workflow

1. **Load configuration** - Parse JSON config and detect file formats
2. **Validate genomes** - Check reference and modified genome files
3. **Validate features** - Validate feature annotations (optional)
4. **Validate reads** - Process all sequencing read files
5. **Generate report** - Create detailed validation report with statistics

---

## Error Handling

The package uses a hierarchical exception system:

- `ValidationError` - Base exception for all validation errors
  - `GenomeValidationError` - Genome-specific errors
    - `FastaFormatError`
    - `GenBankFormatError`
  - `ReadValidationError` - Read-specific errors
    - `FastqFormatError`
    - `BamFormatError`
  - `FeatureValidationError` - Feature-specific errors
  - `ConfigurationError` - Config file errors
  - `FileNotFoundError` - Missing file errors

All errors are logged with structured details for debugging.

---

## Directory Structure

```
validation_pkg/
├── README.md
├── USAGE.md                      # Detailed usage guide
├── config_guide.md               # Configuration documentation
├── setup.py                      # Package installation
├── requirements.txt              # Dependencies
├── simple_usage.py               # Simple example script
├── example_usage.py              # Comprehensive examples
├── validation_pkg/               # Main package
│   ├── __init__.py              # Public API
│   ├── coordinator.py           # Workflow orchestration
│   ├── config_manager.py        # Configuration handling
│   ├── logger.py                # Structured logging
│   ├── exceptions.py            # Exception hierarchy
│   ├── validators/              # Validator modules
│   │   ├── genome_validator.py
│   │   ├── read_validator.py
│   │   └── feature_validator.py
│   ├── parsers/                 # File parsers
│   │   ├── fasta_parser.py
│   │   └── genbank_parser.py
│   └── utils/                   # Utility functions
│       └── file_handler.py
└── tests/                       # Test suite
    ├── fixtures/                # Test data
    └── test_*.py                # Test modules
```

---

## Tips and Best Practices

### 1. Always use Settings.update()

Settings classes use an immutable pattern. Always assign the result of `update()`:

```python
# ✓ CORRECT
settings = settings.update(coding_type='gz')

# ✗ WRONG - changes are lost!
settings.update(coding_type='gz')
```

### 2. Check Config Before Validation

Always check if config fields exist before using them:

```python
if config.ref_genome:
    validate_genome(config.ref_genome, config.output_dir, settings)

if config.reads:
    validate_reads(config.reads, config.output_dir, settings)
```

### 3. Use Appropriate Compression

Choose compression based on your needs:
- `None` - No compression (fast, large files)
- `'gz'` - gzip compression (balanced)
- `'bz2'` - bzip2 compression (slower, smaller files)

### 4. Error Handling

Always wrap validations in try-except blocks:

```python
try:
    stats = validate_genome(config.ref_genome, config.output_dir, settings)
    print(f"✓ Success: {stats['total_sequences']} sequences")
except Exception as e:
    print(f"✗ Validation failed: {e}")
```

---

## Documentation

- **Configuration Guide:** [config_guide.md](config_guide.md)
- **Usage Guide:** [USAGE.md](USAGE.md)
- **Example Scripts:** [simple_usage.py](simple_usage.py), [example_usage.py](example_usage.py)

---

## Requirements

- Python >= 3.7
- Standard library only (no external dependencies for core functionality)
- Optional: `samtools` for BAM file processing

---

## Version Information

- **Version:** 0.1.0
- **Author:** Dominika Bohuslavová
- **License:** MIT
- **Repository:** https://github.com/draessld/validation-pkg

---

## License

MIT License - See LICENSE file for details
