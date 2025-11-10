# Validation Package

A comprehensive validation package for genomic data files used in bioinformatics pipelines. Validates and standardizes genome files, sequencing reads, and feature annotations with detailed error reporting, quality checks, and multiple validation modes for performance optimization.

---

## Table of Contents

- [Documentation](#documentation)
- [Features](#features)
- [Quick Start](#quick-start)
- [Validation Levels](#validation-levels)
- [Configuration](#configuration)
- [Usage Patterns](#usage-patterns)
- [API Reference](#api-reference)
- [Supported File Formats](#supported-file-formats)
- [Performance Optimization](#performance-optimization)
- [Error Handling](#error-handling)
- [Tips and Best Practices](#tips-and-best-practices)
- [Requirements](#requirements)
- [Version Information](#version-information)

---

## Documentation

### Getting Started
- **[Quick Start](#quick-start)** - Get up and running in 5 minutes
- **[examples/README.md](examples/README.md)** - runnable examples with explanations
- **[examples/QUICK_REFERENCE.md](examples/QUICK_REFERENCE.md)** - Syntax cheat sheet

### Core Documentation
- **[docs/CONFIG_GUIDE.md](docs/CONFIG_GUIDE.md)** - Configuration file reference and examples
- **[docs/LOGGING_GUIDE.md](docs/LOGGING_GUIDE.md)** - Logging system and debugging
- **[docs/VALIDATOR_DESIGN.md](docs/VALIDATOR_DESIGN.md)** - Validator architecture and design patterns
- **[docs/PERFORMANCE_GUIDE.md](docs/PERFORMANCE_GUIDE.md)** - Performance optimization techniques

### For Developers
- **[CLAUDE.md](CLAUDE.md)** - AI assistant guidance and architectural decisions
- **[tests/README.md](tests/README.md)** - Testing guide, fixtures, and debugging
- **[examples/INDEX.md](examples/INDEX.md)** - Examples navigation guide

---

## Features

- **Multi-format support:** FASTA, FASTQ, BAM, GTF, GFF, GBK, BED
- **Compression handling:** Automatic detection and processing of gzip (.gz) and bzip2 (.bz2) files
- **Multi-level validation:** Choose between strict, trust, or minimal validation for performance optimization
- **Quality checks:** Sequence validation, duplicate ID detection, invalid character checks, empty sequence detection
- **Format conversion:** GenBank to FASTA, BAM to FASTQ, BED to GFF
- **Plasmid handling:** Automatic plasmid splitting for bacterial genomes
- **Directory-based reads:** Load multiple read files from directories with inherited configuration
- **Parallel processing:** Multi-file and multi-threaded compression for performance
- **Flexible API:** Functional API and direct validator access
- **Configuration-driven:** JSON-based configuration for batch processing
- **Detailed reporting:** Structured validation reports with statistics

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

### Optional Dependencies (Recommended)

For significantly faster compression/decompression of large genomic files, install parallel compression tools:

```bash
# Ubuntu/Debian
sudo apt-get install pigz pbzip2

# macOS (via Homebrew)
brew install pigz pbzip2

# From source (if package manager not available)
# pigz: https://zlib.net/pigz/
# pbzip2: https://launchpad.net/pbzip2
```

**Performance improvement with parallel tools:**
- **pigz** (parallel gzip): 3-4x faster than standard gzip
- **pbzip2** (parallel bzip2): 2-3x faster than standard bzip2
- Automatically detected and used when available
- Graceful fallback to standard gzip/bzip2 if not installed

The package works perfectly fine without these tools, but installing them is highly recommended for production workloads with large files (>1GB).

### Simple Usage

```python
from validation_pkg import ConfigManager, validate_genome, validate_reads, validate_feature

# Load configuration
config = ConfigManager.load("config.json")

# Validate genome files
if config.ref_genome:
    validate_genome(config.ref_genome, config.output_dir)
if config.mod_genome:
    validate_genome(config.mod_genome, config.output_dir)

# Validate all read files
if config.reads:
    validate_reads(config.reads, config.output_dir)

# Validate feature files
if config.ref_feature:
    validate_feature(config.ref_feature, config.output_dir)

print("✓ All validations complete!")
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
- **Use case:** Development, untrusted data sources
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

### Advanced Configuration with Global Options ⭐ NEW

**New in v0.2.0:** Use the `options` field to set global defaults for all files:

```json
{
  "ref_genome_filename": {
    "filename": "reference.fasta.gz",
    "min_sequence_length": 500,
    "replace_id_with": "chr"
  },
  "reads": [
    {
      "filename": "illumina_R1.fastq.gz",
      "ngs_type": "illumina"
    },
    {
      "directory": "ont_reads",
      "ngs_type": "ont",
      "validation_level": "strict"
    }
  ],
  "ref_feature_filename": {
    "filename": "annotations.gff",
    "sort_by_position": true
  },
  "options": {
    "threads": 8,
    "validation_level": "trust",
    "logging_level": "INFO"
  }
}
```

**Global Options:**
- `threads`: Number of parallel threads (default: 8, warns if > CPU cores)
- `validation_level`: `"strict"`, `"trust"`, or `"minimal"` (default: `"strict"`)
- `logging_level`: `"DEBUG"`, `"INFO"`, `"WARNING"`, `"ERROR"`, or `"CRITICAL"` (default: `"INFO"`)

**How the 4-layer settings system works:**

1. **Layer 1 (Defaults):** Built-in defaults from Settings dataclass
2. **Layer 2 (Global options):** `options` field applies to ALL files
   - Only `threads`, `validation_level`, and `logging_level` allowed
   - Prevents mistakes (e.g., can't set `plasmid_split` globally)
3. **Layer 3 (File-level):** Per-file settings override global
   - Any Settings field can be specified
   - WARNING logged when overriding global option
4. **Layer 4 (User code):** Settings in Python completely replace all config settings

**Example result:**
- `ref_genome`: Uses `threads=8`, `validation_level='trust'` (from global), plus file-specific settings
- `illumina_R1`: Uses `threads=8`, `validation_level='trust'` (from global)
- `ont_reads`: Uses `threads=8` (global), `validation_level='strict'` (file overrides global with WARNING)
- `ref_feature`: Uses `threads=8`, `validation_level='trust'` (from global), plus file-specific settings

**Available settings:**
- **Global options** (in `options` field): ONLY `threads` and `validation_level`
- **File-level settings** (per-file): ANY validator setting
  - **Genome files**: `validation_level`, `plasmid_split`, `min_sequence_length`, `replace_id_with`, etc.
  - **Read files**: `validation_level`, `check_invalid_chars`, `allow_duplicate_ids`, `keep_bam`, etc.
  - **Feature files**: `validation_level`, `sort_by_position`, `check_coordinates`, etc.

**Best practices:**
- Use global `options` for common settings (threads, validation_level)
- Override per-file only when needed
- Avoid passing Settings in Python code unless you need full control (it ignores config!)

See [CONFIG_GUIDE.md](docs/CONFIG_GUIDE.md) for complete documentation.

---

## Usage Patterns

The package supports multiple usage patterns to fit different needs:

### Pattern 1: Simple Workflow (Recommended)

Use the functional API for straightforward validation:

```python
from validation_pkg import ConfigManager, validate_genome, validate_reads, validate_feature

# Load configuration
config = ConfigManager.load("config.json")

# Validate all genome files
if config.ref_genome:
    validate_genome(config.ref_genome, config.output_dir)
if config.mod_genome:
    validate_genome(config.mod_genome, config.output_dir)

# Validate all read files
if config.reads:
    validate_reads(config.reads, config.output_dir)

# Validate feature files
if config.ref_feature:
    validate_feature(config.ref_feature, config.output_dir)
```

### Pattern 2: Custom Settings (Functional API)

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

#### Functional API

High-level validation functions:

```python
from validation_pkg import validate_genome, validate_read, validate_reads, validate_feature

# Single genome file
validate_genome(genome_config, output_dir, settings)

# Single read file
validate_read(read_config, output_dir, settings)

# Multiple read files (with optional parallel processing)
validate_reads(read_configs, output_dir, settings)

# Single feature file
validate_feature(feature_config, output_dir, settings)
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

### Thread Optimization ⭐ NEW

The package automatically detects CPU cores and warns if you specify too many threads:

```json
{
  "options": {
    "threads": 16
  }
}
```

If your system has only 4 cores, you'll see:
```
WARNING: Requested 16 threads but system only has 4 CPU cores.
         Performance may degrade due to context switching overhead.
         Consider using threads ≤ 4 for optimal performance.
```

**Best practices:**
- Default threads: 8 (good for most systems)
- Set `threads ≤ CPU cores` for best performance
- Use `threads=1` for single-file processing on low-end systems
- Use `threads > 1` only for multi-file batches or large FASTQ files in strict mode

**Compression tool selection:**
- `threads=1`: Uses gzip/bzip2 (better performance than pigz/pbzip2 with 1 thread)
- `threads > 1`: Uses pigz/pbzip2 if available (3-4x faster)

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
