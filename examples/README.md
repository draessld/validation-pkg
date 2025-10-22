# Validation Package Examples

This directory contains comprehensive examples demonstrating how to use the `validation_pkg` package for validating genomic data files.

## Overview

The validation package provides tools for validating and processing:
- **Genome files**: FASTA (.fasta, .fa, .fna), GenBank (.gb, .gbk)
- **Read files**: FASTQ (.fastq, .fq), BAM (.bam)
- **Feature files**: GFF (.gff, .gff3, .gtf), BED (.bed)

All formats support compression: gzip (.gz), bzip2 (.bz2), or uncompressed.

## Directory Structure

```
examples/
├── README.md                       # This file
├── data/                          # Test data files
│   ├── config.json               # Sample configuration file
│   ├── genome.*                  # Genome files in various formats
│   ├── read.*                    # Read files in various formats
│   └── feature.*                 # Feature files in various formats
├── 01_simple_validation.py       # Basic usage with Functional API
├── 02_selective_validation.py    # Validate specific file types
├── 03_custom_genome_settings.py  # Advanced genome validation
├── 04_custom_read_settings.py    # Advanced read validation
├── 05_custom_feature_settings.py # Advanced feature validation
├── 06_validation_levels.py       # Understanding strict/trust/minimal modes
├── 07_direct_instantiation.py    # Creating validators directly
├── 08_error_handling.py          # Exception handling patterns
├── 09_settings_patterns.py       # Working with immutable settings
├── 10_complete_pipeline.py       # Production-ready pipeline example
└── 11_directory_based_reads.py   # Directory-based read validation
```

## Test Data Files

The `data/` directory contains sample files in all supported formats:

### Genome Files
- `genome.fasta` - Uncompressed FASTA
- `genome.fa.gz` - Gzip compressed FASTA
- `genome.gb.bz2` - Bzip2 compressed GenBank
- `genome.gbk.bz2` - Bzip2 compressed GenBank (alternative extension)

### Read Files
- `read.fastq` - Uncompressed FASTQ (Illumina reads)
- `read.fastq.gz` - Gzip compressed FASTQ
- `read.fq` - Uncompressed FASTQ (ONT reads)
- `read.fq.gz` - Gzip compressed FASTQ

### Feature Files
- `feature.gff` - Uncompressed GFF3
- `feature.gff3.gz` - Gzip compressed GFF3
- `feature.gtf.bz2` - Bzip2 compressed GTF
- `feature.bed` - Uncompressed BED
- `feature.bed.gz` - Gzip compressed BED

### Configuration
- `config.json` - Sample configuration file referencing the test data

## Running the Examples

All examples are self-contained and can be run directly from the examples directory:

```bash
cd examples

# Run individual examples
python 01_simple_validation.py
python 02_selective_validation.py
# ... etc
```

Or make them executable and run directly:

```bash
chmod +x *.py
./01_simple_validation.py
```

Each example will:
1. Use test data from the `data/` directory
2. Create an `output/` subdirectory for results
3. Print progress and results to the console

## Example Descriptions

### 1. Simple Validation (`01_simple_validation.py`)
**Best for**: First-time users, quick validation of all files

The simplest way to use the package:
- Load a configuration file
- Validate all specified files
- Get a comprehensive report

```python
from validation_pkg import Functional API

coordinator = Functional API("config.json")
report = coordinator.validate_all()
```

### 2. Selective Validation (`02_selective_validation.py`)
**Best for**: Processing specific file types separately

Validate only the file types you need:
- Validate only genomes
- Validate only reads
- Validate only features

```python
coordinator = Functional API("config.json")
genome_report = coordinator.validate_genomes()
reads_report = coordinator.validate_reads()
features_report = coordinator.validate_features()
```

### 3. Custom Genome Settings (`03_custom_genome_settings.py`)
**Best for**: Advanced genome processing with custom requirements

Demonstrates:
- Validation levels (strict, trust, minimal)
- Sequence filtering by length
- ID prefix replacement
- Plasmid handling (split/merge)
- Output compression control

```python
settings = GenomeValidator.Settings().update(
    validation_level='strict',
    min_sequence_length=500,
    replace_id_with='chr',
    plasmid_split=True,
    coding_type='gz'
)
```

### 4. Custom Read Settings (`04_custom_read_settings.py`)
**Best for**: Processing sequencing data with specific requirements

Demonstrates:
- Different settings for different NGS types (Illumina, ONT, PacBio)
- Automatic organization by NGS type
- Quality checks (invalid characters, duplicate IDs)
- Compression format conversion
- Fast vs thorough validation

```python
settings = ReadValidator.Settings().update(
    validation_level='trust',
    outdir_by_ngs_type=True,
    check_invalid_chars=True,
    coding_type='gz'
)
```

### 5. Custom Feature Settings (`05_custom_feature_settings.py`)
**Best for**: Processing annotation files

Demonstrates:
- Position-based sorting
- Coordinate validation
- Sequence name replacement
- Format handling (GFF, GTF, BED)

```python
settings = FeatureValidator.Settings().update(
    validation_level='strict',
    sort_by_position=True,
    replace_id_with='chromosome',
    coding_type='gz'
)
```

### 6. Validation Levels (`06_validation_levels.py`)
**Best for**: Understanding performance trade-offs

Compares the three validation levels:

| Level | Speed | Validation | Use Case |
|-------|-------|------------|----------|
| **STRICT** | Slowest | Complete | Development, QC, untrusted data |
| **TRUST** | 10-15x faster | Selective | Production, trusted data |
| **MINIMAL** | Instant | None | Archiving, staging |

```python
# Strict: Full validation
strict_settings = GenomeValidator.Settings().update(
    validation_level='strict'
)

# Trust: Fast validation (10-15x faster)
trust_settings = GenomeValidator.Settings().update(
    validation_level='trust'
)

# Minimal: No validation (instant copy)
minimal_settings = GenomeValidator.Settings().update(
    validation_level='minimal'
)
```

### 7. Direct Instantiation (`07_direct_instantiation.py`)
**Best for**: Programmatic configuration, custom pipelines

Shows how to create validators without config files:
- Manually creating Config objects
- Direct validator instantiation
- Programmatic batch processing

```python
from validation_pkg.config_manager import GenomeConfig
from validation_pkg.utils.formats import GenomeFormat, CodingType

config = GenomeConfig(
    filename="genome.fasta",
    filepath=Path("genome.fasta"),
    coding_type=CodingType.NONE,
    detected_format=GenomeFormat.FASTA,
    extra={}
)

validator = GenomeValidator(config, output_dir, settings)
validator.validate()
```

### 8. Error Handling (`08_error_handling.py`)
**Best for**: Building robust, production-ready pipelines

Demonstrates:
- Exception hierarchy
- Catching specific errors
- Validation report checking
- Graceful error recovery
- Best practices

```python
from validation_pkg.exceptions import (
    ValidationError,
    GenomeValidationError,
    ConfigurationError
)

try:
    coordinator = Functional API("config.json")
    report = coordinator.validate_all()
    if not report.passed:
        # Handle validation failures
        pass
except GenomeValidationError as e:
    # Handle genome-specific errors
    pass
except ConfigurationError as e:
    # Handle config errors
    pass
```

### 9. Settings Patterns (`09_settings_patterns.py`)
**Best for**: Understanding the settings system

Critical for effective use of the package:
- Immutable settings pattern
- Creating and updating settings
- Converting to/from dictionaries
- Reusable settings templates
- Settings for different validators

```python
# IMPORTANT: Settings are immutable!
settings = GenomeValidator.Settings()

# CORRECT: Assign the result
settings = settings.update(min_sequence_length=500)

# WRONG: Changes are lost
settings.update(min_sequence_length=500)  # Don't do this!

# Chaining updates
settings = settings.update(
    validation_level='strict'
).update(
    coding_type='gz'
)
```

### 10. Complete Pipeline (`10_complete_pipeline.py`)
**Best for**: Production deployment template

A complete, production-ready pipeline combining all concepts:
- Configuration loading
- Custom settings for each file type
- Comprehensive error handling
- Progress reporting
- Organized output structure
- Timing and performance metrics

This serves as a template for building your own pipelines.

### 11. Directory-Based Reads (`11_directory_based_reads.py`)
**Best for**: Batch processing multiple read files from the same source

Demonstrates the directory option for reads:
- Using "directory" instead of individual "filename" entries
- Automatic discovery of all files in a directory
- Applying same settings to all files in directory
- Different settings for different sequencing platforms

This is ideal when you have many read files from the same sequencing run.

```python
# Config with directory option
{
  "reads": [
    {
      "directory": "illumina_reads/",
      "ngs_type": "illumina",
      "validation_level": "trust"
    },
    {
      "directory": "ont_reads/",
      "ngs_type": "ont",
      "validation_level": "strict"
    }
  ]
}
```

## Quick Start Guide

### 1. Install the package
```bash
pip install -e /path/to/validation_pkg
```

### 2. Create a configuration file

```json
{
  "ref_genome_filename": {
    "filename": "reference.fasta.gz"
  },
  "reads": [
    {
      "filename": "sample1.fastq.gz",
      "ngs_type": "illumina"
    }
  ],
  "ref_feature_filename": {
    "filename": "annotations.gff.gz"
  }
}
```

### 3. Run validation

```python
from validation_pkg import Functional API

coordinator = Functional API("config.json")
report = coordinator.validate_all()

if report.passed:
    print(f"Success! Validated {len(report.validated_files)} files")
else:
    print(f"Failed: {report.errors} errors")
```

## Common Use Cases

### Use Case 1: Quality Control for New Data
```python
# Use STRICT mode with thorough checks
settings = GenomeValidator.Settings().update(
    validation_level='strict',
    min_sequence_length=100,
    allow_empty_sequences=False,
    coding_type='gz'
)
```

### Use Case 2: Fast Production Pipeline
```python
# Use TRUST mode for speed
settings = ReadValidator.Settings().update(
    validation_level='trust',
    outdir_by_ngs_type=True,
    coding_type='gz'
)
```

### Use Case 3: Format Conversion
```python
# Convert between compression formats
settings = ReadValidator.Settings().update(
    validation_level='trust',
    coding_type='bz2'  # Convert to bzip2
)
```

### Use Case 4: Data Archiving
```python
# Minimal mode for fast copying
settings = GenomeValidator.Settings().update(
    validation_level='minimal'
)
```

## Configuration File Format

### Minimal Configuration
```json
{
  "ref_genome_filename": {
    "filename": "genome.fasta"
  }
}
```

### Full Configuration with Options
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
      "directory": "ont_reads/",
      "ngs_type": "ont",
      "validation_level": "strict"
    }
  ],
  "ref_feature_filename": {
    "filename": "annotations.gff.gz",
    "sort_by_position": true
  }
}
```

## Validation Levels Decision Tree

```
Do you trust the data source?
├─ No / New data
│  └─ Use STRICT mode
│     - Full validation
│     - Catches all errors
│     - Best for development/QC
│
├─ Yes / Previously validated
│  └─ Use TRUST mode
│     - 10-15x faster
│     - Selective validation
│     - Best for production
│
└─ Archiving only (no validation needed)
   └─ Use MINIMAL mode
      - Instant (file copy)
      - No validation
      - Best for staging/archiving
```

## Performance Tips

1. **Use TRUST mode for large files** (>10GB FASTQ)
   - 10-15x faster than STRICT
   - Still validates first record/sequence

2. **Enable `outdir_by_ngs_type`** for reads
   - Automatically organizes by sequencing platform
   - Easier downstream processing

3. **Use compression** (coding_type='gz')
   - Saves disk space
   - Standard in bioinformatics pipelines

4. **Batch process similar files** together
   - Reuse settings objects
   - Process by file type

## Troubleshooting

### Problem: "Config file not found"
**Solution**: Check the path to your config.json file
```python
config_path = Path("config.json")
if not config_path.exists():
    print(f"Config not found at: {config_path.absolute()}")
```

### Problem: "Validation failed" but no clear error
**Solution**: Check the validation report details
```python
report = coordinator.validate_all()
if not report.passed:
    print(report.summary())
    for failed in report.failed_files:
        print(f"File: {failed['filename']}")
        print(f"Error: {failed['error']}")
```

### Problem: Settings changes not applied
**Solution**: Remember settings are immutable!
```python
# WRONG
settings.update(coding_type='gz')  # Changes are lost!

# CORRECT
settings = settings.update(coding_type='gz')  # Assign result
```

### Problem: Minimal mode not working
**Solution**: Input and output compression must match
```python
# If input is .gz, output must be .gz (or None)
settings = GenomeValidator.Settings().update(
    validation_level='minimal',
    coding_type='gz'  # Must match input
)
```

## Additional Resources

- **Package Documentation**: See the main README.md in the package root
- **API Reference**: Check docstrings in the source code
- **Test Suite**: See `tests/` directory for more examples
- **Issue Tracker**: Report bugs or request features on GitHub

## Getting Help

If you have questions or run into issues:

1. Check these examples for similar use cases
2. Review the main package documentation
3. Look at the test suite for edge cases
4. Check the docstrings in the source code
5. Open an issue on GitHub

## License

This examples directory is part of the validation_pkg package.
See the main package LICENSE file for details.
