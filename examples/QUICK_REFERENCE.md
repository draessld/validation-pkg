# Validation Package Quick Reference

## Installation

```bash
pip install -e /path/to/validation_pkg
```

## Supported File Formats

| Type | Formats | Compression |
|------|---------|-------------|
| **Genome** | FASTA (.fasta, .fa, .fna), GenBank (.gb, .gbk) | .gz, .bz2, or none |
| **Read** | FASTQ (.fastq, .fq), BAM (.bam) | .gz, .bz2, or none |
| **Feature** | GFF (.gff, .gff3, .gtf), BED (.bed) | .gz, .bz2, or none |

## Basic Usage

### Simple Validation (Recommended)

```python
from validation_pkg import Functional API

coordinator = Functional API("config.json")
report = coordinator.validate_all()

if report.passed:
    print(f"Success! {len(report.validated_files)} files validated")
```

### Selective Validation

```python
coordinator = Functional API("config.json")

genome_report = coordinator.validate_genomes()
reads_report = coordinator.validate_reads()
features_report = coordinator.validate_features()
```

### Custom Settings

```python
from validation_pkg import GenomeValidator, validate_genome

settings = GenomeValidator.Settings().update(
    validation_level='trust',
    min_sequence_length=500,
    coding_type='gz'
)

validate_genome(config.ref_genome, output_dir, settings)
```

## Configuration File

### Minimal Example

```json
{
  "ref_genome_filename": {
    "filename": "genome.fasta"
  },
  "reads": [
    {
      "filename": "reads.fastq.gz",
      "ngs_type": "illumina"
    }
  ]
}
```

### With Custom Settings

```json
{
  "ref_genome_filename": {
    "filename": "genome.fasta.gz",
    "validation_level": "trust",
    "min_sequence_length": 500,
    "replace_id_with": "chr",
    "coding_type": "gz",
    "threads": 4
  },
  "reads": [
    {
      "filename": "sample1.fastq.gz",
      "ngs_type": "illumina",
      "validation_level": "trust",
      "outdir_by_ngs_type": true
    }
  ],
  "ref_feature_filename": {
    "filename": "features.gff.gz",
    "sort_by_position": true,
    "coding_type": "gz"
  },
  "options": {
    "threads": 8
  }
}
```

### With Directory-Based Reads

```json
{
  "ref_genome_filename": {
    "filename": "genome.fasta.gz"
  },
  "reads": [
    {
      "directory": "illumina_reads/",
      "ngs_type": "illumina",
      "validation_level": "trust"
    },
    {
      "directory": "ont_reads/",
      "ngs_type": "ont",
      "validation_level": "strict",
      "check_invalid_chars": true
    }
  ]
}
```

**Note**: When using "directory", all files in that directory will be automatically discovered and validated with the specified settings.

## Validation Levels

| Level | Speed | Validation | Use Case |
|-------|-------|------------|----------|
| `strict` | Slowest | Complete | Development, QC, untrusted data |
| `trust` | 10-15x faster | First record only | Production, trusted data |
| `minimal` | Instant | None (copy only) | Archiving, staging |

```python
# Strict - full validation
settings = GenomeValidator.Settings().update(
    validation_level='strict'
)

# Trust - fast validation (recommended for production)
settings = GenomeValidator.Settings().update(
    validation_level='trust'
)

# Minimal - no validation (archiving only)
settings = GenomeValidator.Settings().update(
    validation_level='minimal'
)
```

## Common Settings

### Genome Validator Settings

```python
settings = GenomeValidator.Settings().update(
    validation_level='strict',          # 'strict', 'trust', or 'minimal'
    min_sequence_length=100,            # Filter short sequences
    replace_id_with='chr',              # Add prefix to sequence IDs
    plasmid_split=True,                 # Split plasmids to separate files
    plasmids_to_one=False,              # Merge all plasmids to one file
    is_plasmid=False,                   # Treat all as plasmids
    coding_type='gz',                   # Output compression: 'gz', 'bz2', or None
    output_filename_suffix='validated', # Add suffix to output filename
    output_subdir_name='genomes',       # Output subdirectory
    threads=8,                          # Unified threads parameter (automatic split)
    max_workers=None,                   # File-level parallelization (power users)
    compression_threads=None            # Compression-level threads (power users)
)
```

### Read Validator Settings

```python
settings = ReadValidator.Settings().update(
    validation_level='trust',           # 'strict', 'trust', or 'minimal'
    check_invalid_chars=True,           # Check for invalid nucleotides
    allow_empty_id=False,               # Allow reads without IDs
    allow_duplicate_ids=True,           # Allow duplicate IDs
    keep_bam=True,                      # Keep original BAM file
    ignore_bam=True,                    # Skip BAM files (don't convert)
    outdir_by_ngs_type=True,            # Organize by NGS type (illumina/ont/pacbio)
    coding_type='gz',                   # Output compression: 'gz', 'bz2', or None
    output_filename_suffix='validated', # Add suffix to output filename
    threads=8,                          # Unified threads parameter (automatic split)
    max_workers=None,                   # File-level parallelization (power users)
    compression_threads=None            # Compression-level threads (power users)
)
```

### Feature Validator Settings

```python
settings = FeatureValidator.Settings().update(
    validation_level='strict',          # 'strict', 'trust', or 'minimal'
    sort_by_position=True,              # Sort features by genomic position
    check_coordinates=True,             # Validate feature coordinates
    allow_zero_length=False,            # Allow zero-length features
    replace_id_with='chr1',             # Replace sequence name field
    coding_type='gz',                   # Output compression: 'gz', 'bz2', or None
    output_filename_suffix='validated', # Add suffix to output filename
    threads=8,                          # Unified threads parameter (automatic split)
    max_workers=None,                   # File-level parallelization (power users)
    compression_threads=None            # Compression-level threads (power users)
)
```

## Settings Pattern (IMPORTANT!)

**Settings are IMMUTABLE** - always assign the result of `update()`:

```python
# CORRECT ✓
settings = GenomeValidator.Settings()
settings = settings.update(min_sequence_length=500)

# WRONG ✗ - changes are lost!
settings = GenomeValidator.Settings()
settings.update(min_sequence_length=500)  # Does nothing!
```

### Chaining Updates

```python
settings = GenomeValidator.Settings().update(
    validation_level='trust',
    min_sequence_length=500
).update(
    coding_type='gz',
    replace_id_with='chr'
)
```

### Converting to/from Dictionary

```python
# To dictionary
settings_dict = settings.to_dict()

# From dictionary
settings = GenomeValidator.Settings.from_dict(settings_dict)
```

## Error Handling

```python
from validation_pkg import Functional API
from validation_pkg.exceptions import (
    ValidationError,
    GenomeValidationError,
    ConfigurationError
)

try:
    coordinator = Functional API("config.json")
    report = coordinator.validate_all()

    if not report.passed:
        print(f"Errors: {report.errors}")
        for failed in report.failed_files:
            print(f"  {failed['filename']}: {failed['error']}")

except ConfigurationError as e:
    print(f"Config error: {e}")
except GenomeValidationError as e:
    print(f"Genome validation error: {e}")
except ValidationError as e:
    print(f"Validation error: {e}")
```

## Direct Instantiation

```python
from pathlib import Path
from validation_pkg.validators.genome_validator import GenomeValidator
from validation_pkg.config_manager import GenomeConfig
from validation_pkg.utils.formats import GenomeFormat, CodingType

# Create config manually
config = GenomeConfig(
    filename="genome.fasta.gz",
    filepath=Path("data/genome.fasta.gz"),
    detected_format=GenomeFormat.FASTA,
    coding_type=CodingType.GZIP,
    extra={}
)

# Create settings
settings = GenomeValidator.Settings().update(
    validation_level='strict',
    min_sequence_length=500
)

# Create and run validator
validator = GenomeValidator(config, output_dir, settings)
validator.validate()
```

## Common Use Cases

### 1. Quality Control (New Data)
```python
settings = GenomeValidator.Settings().update(
    validation_level='strict',
    min_sequence_length=100,
    allow_empty_sequences=False,
    coding_type='gz'
)
```

### 2. Production Pipeline (Trusted Data)
```python
settings = ReadValidator.Settings().update(
    validation_level='trust',
    outdir_by_ngs_type=True,
    coding_type='gz'
)
```

### 3. Format Conversion
```python
# Convert to bzip2 compression
settings = ReadValidator.Settings().update(
    validation_level='trust',
    coding_type='bz2'
)
```

### 4. Feature Sorting and Cleanup
```python
settings = FeatureValidator.Settings().update(
    validation_level='strict',
    sort_by_position=True,
    replace_id_with='chromosome',
    coding_type='gz'
)
```

### 5. Fast Archiving
```python
settings = GenomeValidator.Settings().update(
    validation_level='minimal'
    # Note: input and output compression must match
)
```

### 6. Parallel Processing (Multiple Files)
```python
from validation_pkg import validate_reads

# RECOMMENDED: Unified threads parameter (automatic splitting)
settings = ReadValidator.Settings().update(
    threads=8,  # Automatically splits based on file count
    validation_level='trust'
)
validate_reads(config.reads, output_dir, settings)
# 4 files → 4 workers × 2 compression threads each
# 1 file → 1 worker × 8 compression threads
# 8+ files → 8 workers × 1 compression thread each

# ADVANCED: Manual control (power users)
settings = ReadValidator.Settings().update(
    max_workers=4,  # Process 4 files concurrently
    compression_threads=2,  # Each uses 2 threads for compression
    validation_level='trust'
)
```

## Parallel Processing

### Unified Threads Parameter (Recommended)

The `threads` parameter automatically splits between file-level and compression-level parallelization:

```python
# Simple approach - automatic optimization
settings = ReadValidator.Settings().update(
    threads=8,  # Auto-splits intelligently
    validation_level='trust'
)
validate_reads(config.reads, output_dir, settings)
```

**How it works:**
- 1 file: All threads for compression (1 worker × 8 threads)
- 4 files: Balanced split (4 workers × 2 threads each)
- 8+ files: Prioritize file parallelization (8 workers × 1 thread each)

### Manual Control (Power Users)

You can still set `max_workers` and `compression_threads` separately:

```python
settings = ReadValidator.Settings().update(
    max_workers=4,  # File-level parallelization
    compression_threads=2,  # Per-file compression threads
    validation_level='trust'
)
```

### Config-Level Threads

Specify threads in config.json:

```json
{
  "ref_genome_filename": {
    "filename": "genome.fasta",
    "threads": 4
  },
  "options": {
    "threads": 8
  }
}
```

## Validation Report

```python
report = coordinator.validate_all()

# Check success
if report.passed:
    print("All validations passed")

# Get statistics
print(f"Files validated: {len(report.validated_files)}")
print(f"Errors: {report.errors}")
print(f"Warnings: {report.warnings}")

# Get failed files
for failed in report.failed_files:
    print(f"Failed: {failed['filename']}")
    print(f"  Error: {failed['error']}")

# Print summary
print(report.summary())
```

## NGS Types

For read files, specify the sequencing platform:

- `"illumina"` - Illumina short reads
- `"ont"` - Oxford Nanopore long reads
- `"pacbio"` - PacBio long reads

```json
{
  "reads": [
    {
      "filename": "illumina_reads.fastq.gz",
      "ngs_type": "illumina"
    },
    {
      "filename": "ont_reads.fastq.gz",
      "ngs_type": "ont"
    }
  ]
}
```

## Output Organization

### Default Output
```
output_dir/
├── validated_genome.fasta.gz
├── validated_reads.fastq.gz
└── validated_features.gff.gz
```

### With `output_subdir_name`
```
output_dir/
└── genomes/
    └── validated_genome.fasta.gz
```

### With `outdir_by_ngs_type` (reads only)
```
output_dir/
├── illumina/
│   └── sample1.fastq.gz
├── ont/
│   └── sample2.fastq.gz
└── pacbio/
    └── sample3.fastq.gz
```

### With `plasmid_split` (genomes only)
```
output_dir/
├── chromosome1.fasta.gz
├── plasmid1.fasta.gz
└── plasmid2.fasta.gz
```

## Performance Tips

1. **Use `trust` mode for large files** - 10-15x faster than `strict`
2. **Use parallel processing for multiple files** - Set `threads=8` for automatic optimization (35-40x faster when combined with trust mode)
3. **Enable compression** (`coding_type='gz'`) - saves disk space, pigz/pbzip2 auto-detected for 3-4x compression speedup
4. **Use `outdir_by_ngs_type`** for reads - organizes by platform
5. **Batch similar files** - reuse settings objects
6. **Use config-level settings** - specify validation settings in config.json for production
7. **Use `minimal` mode** only for archiving (no validation)

## Decision Tree

```
What do you need to do?

├─ Validate new/untrusted data
│  └─ Use validation_level='strict'
│
├─ Fast processing of trusted data
│  └─ Use validation_level='trust'
│
├─ Archive/stage files (no validation)
│  └─ Use validation_level='minimal'
│
├─ Convert compression format
│  └─ Use validation_level='trust' + coding_type='gz'/'bz2'
│
├─ Filter sequences by length
│  └─ Use min_sequence_length=N
│
├─ Split plasmids
│  └─ Use plasmid_split=True
│
├─ Organize reads by platform
│  └─ Use outdir_by_ngs_type=True
│
└─ Sort features by position
   └─ Use sort_by_position=True
```

## Exception Hierarchy

```
ValidationError (base)
├── ConfigurationError
├── FileNotFoundError
├── FileFormatError
│   ├── FastaFormatError
│   ├── GenBankFormatError
│   ├── FastqFormatError
│   ├── BamFormatError
│   ├── GffFormatError
│   └── BedFormatError
├── CompressionError
├── GenomeValidationError
├── ReadValidationError
├── FeatureValidationError
└── InterFileValidationError
```

## Functional API

```python
from validation_pkg import (
    validate_genome,
    validate_genomes,  # Parallel genome processing
    validate_read,
    validate_reads,  # Parallel read processing
    validate_feature,
    validate_features_list  # Parallel feature processing
)

# Single genome
validate_genome(genome_config, output_dir, settings)

# Multiple genomes (parallel)
validate_genomes(genome_configs, output_dir, settings)

# Single read
validate_read(read_config, output_dir, settings)

# Multiple reads (parallel)
validate_reads(read_configs, output_dir, settings)

# Single feature
validate_feature(feature_config, output_dir, settings)

# Multiple features (parallel)
validate_features_list(feature_configs, output_dir, settings)
```

## Examples Directory

See the `examples/` directory for complete working examples:

- `01_simple_validation.py` - Basic usage
- `02_selective_validation.py` - Validate by file type
- `03_custom_genome_settings.py` - Genome-specific options
- `04_custom_read_settings.py` - Read-specific options
- `05_custom_feature_settings.py` - Feature-specific options
- `06_validation_levels.py` - Performance comparison
- `07_direct_instantiation.py` - Manual config creation
- `08_error_handling.py` - Exception handling
- `09_settings_patterns.py` - Settings system details
- `10_complete_pipeline.py` - Production-ready template
- `11_directory_based_reads.py` - Directory-based read validation
- `12_config_validation_levels.py` - Config-level settings
- `13_parallel_processing.py` - Parallel file processing

## Getting Help

1. Check the examples directory for similar use cases
2. Read the main README.md for detailed documentation
3. Check docstrings in the source code
4. Review the test suite for edge cases
5. Open an issue on GitHub
