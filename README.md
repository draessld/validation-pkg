# Validation Package

A comprehensive validation package for genomic data files used in bioinformatics pipelines. Validates and standardizes genome files, sequencing reads, and feature annotations with detailed error reporting, quality checks, and multiple validation modes for performance optimization.

---

## Table of Contents

- [Features](#features)
- [Quick Start](#quick-start)
- [Validation Levels](#validation-levels)
- [Configuration](#configuration)
- [Usage Patterns](#usage-patterns)
- [API Reference](#api-reference)
- [Command Line Interface](#command-line-interface)
- [Supported File Formats](#supported-file-formats)
- [Performance Optimization](#performance-optimization)
- [Error Handling](#error-handling)
- [Tips and Best Practices](#tips-and-best-practices)
- [Requirements](#requirements)
- [Version Information](#version-information)

---

## Features

- **Multi-format support:** FASTA, FASTQ, BAM, GTF, GFF, GBK, BED
- **Compression handling:** Automatic detection and processing of gzip (.gz) and bzip2 (.bz2) files
- **Multi-level validation:** Choose between strict, trust, or minimal validation for performance optimization
- **Quality checks:** Sequence validation, duplicate ID detection, invalid character checks, empty sequence detection
- **Format conversion:** GenBank to FASTA, BAM to FASTQ, BED to GFF
- **Plasmid handling:** Automatic plasmid splitting for bacterial genomes
- **Directory-based reads:** Load multiple read files from directories with inherited configuration
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

---

## Validation Levels

The package supports three validation levels to balance thoroughness and performance:

### Comparison Table

| Level | Parsing | Validation | Edits | Output | Speed | Use Case |
|-------|---------|------------|-------|--------|-------|----------|
| **strict** (default) | All sequences | All sequences | All applied | BioPython write | Slowest | Development, QC, untrusted data |
| **trust** | All sequences (genome)<br>First record only (reads) | First sequence only | All applied (genome)<br>None (reads) | BioPython write (genome)<br>File copy (reads) | Fast | Production with trusted data |
| **minimal** | None | None | None | File copy | Fastest | Staging, archiving |

### Detailed Behavior

#### **Strict Mode** (default: `validation_level='strict'`)
**Full validation - comprehensive**
- Parses and validates every sequence
- Applies all editing specifications
- Complete statistics collection
- **Use case:** Development, quality control, untrusted data sources
- **Performance:** Slowest, most thorough

#### **Trust Mode** (`validation_level='trust'`)
**Fast validation with basic checks**

For **Genome Files**:
- Parses ALL sequences
- Validates ONLY first sequence
- Applies ALL edits (min_sequence_length, replace_id_with, plasmid handling)
- Writes sequences using BioPython

For **Read Files**:
- Checks line count (FASTQ: divisible by 4)
- Validates ONLY first record format
- NO edits applied
- Copies file as-is (byte-for-byte)

For **Feature Files**:
- Not implemented yet

**Use case:** Production pipelines with trusted upstream data
**Performance:** Faster validation, maintains editing capabilities

#### **Minimal Mode** (`validation_level='minimal'`)
**No validation - just file copy**
- Skips parsing completely
- No validation, no edits
- Copies file as-is to output
- **Use case:** Staging files, archiving, maximum trust scenarios
- **Performance:** Fastest (just file copy)

### Usage Examples

```python
from validation_pkg import GenomeValidator, ReadValidator

# Genome with trust mode
genome_settings = GenomeValidator.Settings(
    validation_level='trust',
    min_sequence_length=100,
    replace_id_with='chromosome'
)

# Read with trust mode (fast for large FASTQ files)
read_settings = ReadValidator.Settings(
    validation_level='trust',
    outdir_by_ngs_type=True
)

# Minimal mode (fastest - just copy)
minimal_settings = ReadValidator.Settings(
    validation_level='minimal'
)
```

---

## Configuration

### Basic Configuration Example

Create a `config.json` file:

```json
{
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

### Directory-Based Reads

Load multiple read files from a directory:

```json
{
  "ref_genome_filename": {"filename": "reference.fasta"},
  "mod_genome_filename": {"filename": "modified.fasta"},
  "reads": [
    {
      "directory": "reads_dir",
      "ngs_type": "ont"
    }
  ]
}
```

All files in `reads_dir/` will be loaded and assigned `ngs_type: "ont"`.

### Advanced Configuration

```json
{
  "ref_genome_filename": {
    "filename": "reference.fasta.gz",
    "validation_level": "trust",
    "min_sequence_length": 500,
    "replace_id_with": "chr"
  },
  "reads": [
    {
      "filename": "illumina_R1.fastq.gz",
      "ngs_type": "illumina",
      "validation_level": "trust"
    },
    {
      "directory": "ont_reads",
      "ngs_type": "ont",
      "validation_level": "strict"
    }
  ]
}
```

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

# Create custom settings with validation levels
ref_settings = GenomeValidator.Settings()
ref_settings = ref_settings.update(
    validation_level='trust',
    plasmid_split=True,
    coding_type='gz',
    output_filename_suffix='ref',
    replace_id_with='chr',
    min_sequence_length=500
)

reads_settings = ReadValidator.Settings()
reads_settings = reads_settings.update(
    validation_level='trust',
    coding_type='gz',
    outdir_by_ngs_type=True
)

# Validate with custom settings
if config.ref_genome:
    stats = validate_genome(config.ref_genome, config.output_dir, ref_settings)
    print(f"Validated {stats['total_sequences']} sequences")

if config.reads:
    stats_list = validate_reads(config.reads, config.output_dir, reads_settings)
```

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

### GenomeValidator Settings

```python
genome_settings = GenomeValidator.Settings()
genome_settings = genome_settings.update(
    # Validation level
    validation_level='strict',  # 'strict', 'trust', or 'minimal'

    # Validation thresholds
    allow_empty_sequences=False,
    allow_empty_id=False,
    warn_n_sequences=2,

    # Plasmid handling
    is_plasmid=False,
    plasmid_split=True,
    plasmids_to_one=False,
    main_longest=True,
    main_first=False,

    # Editing specifications
    replace_id_with='chr',
    min_sequence_length=100,

    # Output format
    coding_type='gz',
    output_filename_suffix='ref',
    output_subdir_name='genomes'
)
```

### ReadValidator Settings

```python
read_settings = ReadValidator.Settings()
read_settings = read_settings.update(
    # Validation level
    validation_level='strict',  # 'strict', 'trust', or 'minimal'

    # Validation thresholds (strict mode only)
    check_invalid_chars=True,
    allow_empty_id=False,
    allow_duplicate_ids=False,

    # BAM handling
    keep_bam=False,
    ignore_bam=True,

    # Output format
    coding_type='gz',
    output_filename_suffix='filtered',
    output_subdir_name='reads',
    outdir_by_ngs_type=True  # Auto-set subdir to ngs_type
)
```

### FeatureValidator Settings

```python
feature_settings = FeatureValidator.Settings()
feature_settings = feature_settings.update(
    sort_by_position=True,
    check_coordinates=True,
    allow_zero_length=False,
    coding_type='gz',
    output_filename_suffix='sorted'
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
- **Compression:** `.gz`, `.bz2`, `.gzip`, `.bzip2`

### Read Files
- **FASTQ:** `.fastq`, `.fq`
- **BAM:** `.bam` (limited support)
- **Compression:** `.gz`, `.bz2`

### Feature Files
- **GFF:** `.gff`, `.gff3`, `.gtf`
- **BED:** `.bed`
- **Compression:** `.gz`, `.bz2`

---

## Performance Optimization

### Choosing the Right Validation Level

#### Use **Strict Mode** when:
- 🔍 Validating data from external/untrusted sources
- 🐛 Debugging data quality issues
- ✅ Running initial QC on new datasets
- 📊 Need complete statistics

#### Use **Trust Mode** when:
- ⚡ Processing large FASTQ files (>10GB)
- 🔁 Re-running pipelines with previously validated data
- 🏭 Production pipelines with trusted upstream steps
- ✏️ Need editing features but not full validation

#### Use **Minimal Mode** when:
- 📦 Staging files for archival
- 🔄 Simply copying/organizing validated files
- ⏱️ Maximum speed required
- ✅ Data already validated in previous steps

### Performance Comparison

For a 10GB FASTQ file with 50 million reads:

| Mode | Time | Memory | Notes |
|------|------|--------|-------|
| Strict | ~30 min | ~2GB | Full validation |
| Trust | ~2 min | ~100MB | Line count + first record |
| Minimal | ~30 sec | ~10MB | File copy only |

For a genome file with 3 sequences (1 chromosome + 2 plasmids):

| Mode | Time | Memory | Notes |
|------|------|--------|-------|
| Strict | ~1 sec | ~50MB | All sequences validated |
| Trust | ~0.5 sec | ~50MB | First sequence validated |
| Minimal | ~0.1 sec | ~5MB | File copy only |

### Memory-Efficient Processing

```python
# For large datasets, use trust mode with compression
settings = ReadValidator.Settings(
    validation_level='trust',
    coding_type='gz',
    outdir_by_ngs_type=True
)
```

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
  - `CompressionError` - Compression/decompression errors

All errors are logged with structured details for debugging.

**Example:**

```python
from validation_pkg import validate_genome
from validation_pkg.exceptions import GenomeValidationError

try:
    stats = validate_genome(config.ref_genome, config.output_dir, settings)
    print(f"✓ Success: {stats['total_sequences']} sequences")
except GenomeValidationError as e:
    print(f"✗ Genome validation failed: {e}")
except Exception as e:
    print(f"✗ Unexpected error: {e}")
```

---

## Directory Structure

```
validation_pkg/
├── README.md                        # This file
├── setup.py                         # Package installation
├── requirements.txt                 # Dependencies
├── validation_pkg/                  # Main package
│   ├── __init__.py                 # Public API exports
│   ├── coordinator.py              # Workflow orchestration
│   ├── config_manager.py           # Configuration handling
│   ├── logger.py                   # Structured logging
│   ├── exceptions.py               # Exception hierarchy
│   ├── validators/                 # Validator modules
│   │   ├── __init__.py
│   │   ├── genome_validator.py    # Genome validation
│   │   ├── read_validator.py      # Read validation
│   │   └── feature_validator.py   # Feature validation
│   └── utils/                      # Utility functions
│       ├── __init__.py
│       ├── file_handler.py        # File I/O utilities
│       ├── formats.py             # Format enums
│       └── settings.py            # Settings base class
└── tests/                          # Test suite
    ├── fixtures/                   # Test data
    ├── test_config_manager.py     # Config tests
    ├── test_genome_validator.py   # Genome tests
    ├── test_read_validator.py     # Read tests
    └── test_feature_validator.py  # Feature tests
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
- `'gz'` - gzip compression (balanced, recommended)
- `'bz2'` - bzip2 compression (slower, smaller files)

### 4. Choose Validation Levels Wisely

```python
# For initial QC - use strict
initial_settings = ReadValidator.Settings(validation_level='strict')

# For production pipeline - use trust
prod_settings = ReadValidator.Settings(validation_level='trust')

# For archiving - use minimal
archive_settings = ReadValidator.Settings(validation_level='minimal')
```

### 5. Use Directory Reads for Multiple Files

Instead of listing each file individually, use directory-based reads:

```json
{
  "reads": [
    {
      "directory": "illumina_reads",
      "ngs_type": "illumina",
      "validation_level": "trust"
    }
  ]
}
```

All files inherit the `ngs_type` and other settings.

### 6. Organize Outputs by NGS Type

```python
settings = ReadValidator.Settings(
    outdir_by_ngs_type=True  # Creates output/illumina/, output/ont/, etc.
)
```

### 7. Error Handling

Always wrap validations in try-except blocks:

```python
try:
    stats = validate_genome(config.ref_genome, config.output_dir, settings)
    print(f"✓ Success: {stats['total_sequences']} sequences")
except Exception as e:
    print(f"✗ Validation failed: {e}")
```

### 8. Monitor Large File Processing

For large datasets, enable verbose logging:

```bash
python -m validation_pkg validate config.json --verbose --log-file validation.log
```

---

## Requirements

### Python Version
- Python >= 3.7

### Dependencies
- **BioPython** == 1.85 - Biological sequence handling
- **pysam** >= 0.19.0 - BAM file support
- **structlog** >= 24.1.0 - Structured logging
- **pytest** >= 7.0.0 - Testing framework (dev)

### Optional
- **samtools** - External BAM processing (system package)

### Installation

```bash
pip install -e .
```

All dependencies will be installed automatically from `requirements.txt`.

---

## Version Information

- **Version:** 0.1.0
- **Author:** Dominika Bohuslavová
- **Email:** dominikadraesslerova@gmail.com
- **License:** MIT
- **Repository:** https://github.com/draessld/validation-pkg

---

## License

MIT License - See LICENSE file for details

---

## Contributing

Contributions are welcome! Please:
1. Fork the repository
2. Create a feature branch
3. Add tests for new features
4. Submit a pull request

---

## Support

For issues, questions, or contributions:
- **GitHub Issues:** https://github.com/draessld/validation-pkg/issues
- **Email:** dominikadraesslerova@gmail.com

---

## Changelog

### Version 0.1.0
- Initial release
- Multi-format support (FASTA, FASTQ, BAM, GFF, GTF, BED, GenBank)
- Multi-level validation (strict/trust/minimal)
- Compression support (gzip, bzip2)
- Directory-based read loading
- Plasmid handling for bacterial genomes
- Structured logging with error reporting
- Configuration-driven workflow
- CLI support
