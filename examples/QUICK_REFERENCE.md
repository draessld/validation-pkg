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
    "coding_type": "gz"
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
    output_subdir_name='genomes'        # Output subdirectory
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
    output_filename_suffix='validated'  # Add suffix to output filename
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
    output_filename_suffix='validated'  # Add suffix to output filename
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
2. **Enable compression** (`coding_type='gz'`) - saves disk space
3. **Use `outdir_by_ngs_type`** for reads - organizes by platform
4. **Batch similar files** - reuse settings objects
5. **Use `minimal` mode** only for archiving (no validation)

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
    validate_read,
    validate_reads,
    validate_feature
)

# Single genome
validate_genome(genome_config, output_dir, settings)

# Single read
validate_read(read_config, output_dir, settings)

# Multiple reads
validate_reads(read_configs, output_dir, settings)

# Single feature
validate_feature(feature_config, output_dir, settings)
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

## Getting Help

1. Check the examples directory for similar use cases
2. Read the main README.md for detailed documentation
3. Check docstrings in the source code
4. Review the test suite for edge cases
5. Open an issue on GitHub
