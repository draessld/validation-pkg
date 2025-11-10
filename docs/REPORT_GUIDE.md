# Validation Report Guide

## Overview

The validation report system (`validation_pkg.report`) provides comprehensive reporting of validation results, including input settings, output metadata, statistics, and inter-file validation results. Reports can be generated in both human-readable text format and machine-parseable JSON format.

## Table of Contents

- [Quick Start](#quick-start)
- [Report Components](#report-components)
- [Using ValidationReport](#using-validationreport)
- [Report Formats](#report-formats)
- [Displayed Information](#displayed-information)
- [Examples](#examples)
- [API Reference](#api-reference)

---

## Quick Start

```python
from pathlib import Path
from validation_pkg import ConfigManager, validate_genome, validate_reads
from validation_pkg.report import ValidationReport

# Create report
report = ValidationReport(Path("validation_report.txt"))

# Run validations and add to report
config = ConfigManager.load("config.json")

# Validate genome
genome_result = validate_genome(config.ref_genome, config.output_dir)
report.write(
    genome_result,
    file_type="genome",
    input_file=config.ref_genome.filename,
    elapsed_time=2.5
)

# Validate reads
read_results = validate_reads(config.reads, config.output_dir)
report.write(
    read_results,
    file_type="read",
    input_file=config.reads[0].filename
)

# Generate reports
report.flush(format="text")  # Human-readable
report.flush(format="json")  # Machine-parseable
```

---

## Report Components

A validation report contains three main sections:

### 1. Summary Section
- Overall validation status (PASSED/FAILED)
- File counts by type (genomes, reads, features)
- Inter-file validation counts (passed/failed)
- Total execution time

### 2. File Validation Results
For each validated file:
- **Input information**: Original filename, validator type, settings used
- **Output information**: Output file path, validation level
- **Metadata**: Type-specific information and statistics
  - **Genomes**: Sequence count, IDs, lengths, GC content
  - **Reads**: Read count, NGS type, pairing info, N50, lengths, total bases
  - **Features**: Feature count, types, sequence coverage
- **Timing**: Elapsed time for validation

### 3. Inter-file Validation Results
- Validation type (genome×genome, read×read, feature×genome)
- Status (PASSED/FAILED)
- Errors and warnings
- Metadata about cross-file checks

---

## Using ValidationReport

### Creating a Report

```python
from pathlib import Path
from validation_pkg.report import ValidationReport

# Create report (directory will be created if it doesn't exist)
report = ValidationReport(Path("reports/validation_report.txt"))
```

### Adding File Validation Results

**Single file result:**
```python
result = validate_genome(config.ref_genome, output_dir, settings)
report.write(
    result,
    file_type="genome",              # "genome", "read", or "feature"
    input_file="genome.fasta",       # Original filename
    settings=settings,                # Optional: settings object
    elapsed_time=12.5                 # Optional: validation time in seconds
)
```

**Multiple files (e.g., multiple reads):**
```python
results = validate_reads(config.reads, output_dir)
# Results is a list - each item will be recorded separately
report.write(
    results,
    file_type="read",
    input_file="reads.fastq"  # Will be inferred per-file if not provided
)
```

### Adding Inter-file Validation Results

```python
# Genome inter-file validation
from validation_pkg.validators.interfile_genome import genomexgenome_validation

interfile_result = genomexgenome_validation(
    ref_result,
    mod_result,
    settings
)

report.write(interfile_result, file_type="genomexgenome")
```

**Inter-file validation types:**
- `"genomexgenome"` - Compare reference vs modified genomes
- `"readxread"` - Check paired-end read consistency
- `"featurexgenome"` - Verify feature coordinates against genome

### Generating Report Files

```python
# Generate human-readable text report
report.flush(format="text")  # Creates: validation_report.txt

# Generate machine-parseable JSON report
report.flush(format="json")  # Creates: validation_report.json
```

---

## Report Formats

### Text Format

Human-readable format with clear sections and formatting:

```
====================================================================================================
                                     VALIDATION PIPELINE REPORT
====================================================================================================
Generated: 2025-11-10 21:30:45
Total Duration: 15.32s

SUMMARY
====================================================================================================
Overall Status: ✓ PASSED

Files Processed: 3
  Genomes:  2
  Reads:    1
  Features: 0

Inter-file Validations: 1
  Passed:  1
  Failed:  0

====================================================================================================
                                      FILE VALIDATION RESULTS
====================================================================================================

[1] GENOME: ref_genome.fasta.gz
----------------------------------------------------------------------------------------------------
Input File:  ref_genome.fasta
Output File: /output/genomes/ref_genome.fasta.gz
Output Metadata:
  Sequences: 3
    IDs: chr1, chr2, chr3
  Total Length: 5,000,000 bp

[2] READ: reads_R1.fastq.gz
----------------------------------------------------------------------------------------------------
Input File:  reads.fastq.gz
Output File: /output/reads/reads_R1.fastq.gz
Output Metadata:
  Reads: 1,000,000
  NGS Type: illumina
  Paired-End: R1 (base: reads)
  Total Bases: 150,000,000 bp
  Mean Read Length: 150.0 bp
  N50: 155 bp
  Longest Read: 200 bp
  Shortest Read: 100 bp

====================================================================================================
                                    INTER-FILE VALIDATION RESULTS
====================================================================================================

[1] GENOME×GENOME
----------------------------------------------------------------------------------------------------
Status: ✓ PASSED
Warnings:
  - Minor length difference detected (0.01%)
```

### JSON Format

Structured format for programmatic processing:

```json
{
  "report_metadata": {
    "generated": "2025-11-10T21:30:45.123456",
    "duration_seconds": 15.32,
    "version": "1.0.0"
  },
  "summary": {
    "overall_status": "PASSED",
    "files_processed": 3,
    "genomes": 2,
    "reads": 1,
    "features": 0,
    "interfile_validations": 1,
    "interfile_passed": 1,
    "interfile_failed": 0
  },
  "file_validations": [
    {
      "input_file": "ref_genome.fasta",
      "validator_type": "genome",
      "input_settings": {
        "validation_level": "strict",
        "plasmid_split": true,
        "min_sequence_length": 1000
      },
      "output_file": "/output/genomes/ref_genome.fasta.gz",
      "output_metadata": {
        "num_sequences": 3,
        "sequence_ids": ["chr1", "chr2", "chr3"],
        "total_length": 5000000
      },
      "elapsed_time": 12.5,
      "timestamp": "2025-11-10T21:30:32.123456"
    }
  ],
  "inter_file_validations": [
    {
      "validation_type": "genomexgenome",
      "status": "PASSED",
      "passed": true,
      "errors": [],
      "warnings": ["Minor length difference detected (0.01%)"],
      "metadata": {
        "reference_sequences": 3,
        "modified_sequences": 3
      },
      "timestamp": "2025-11-10T21:30:45.123456"
    }
  ]
}
```

---

## Displayed Information

### Genome Validation Results

**Metadata displayed:**
- Number of sequences
- Sequence IDs (first 5 if >5 total)
- Total genome length
- Validation level used

**Example:**
```
Output Metadata:
  Sequences: 3
    IDs: chr1, chr2, chr3
  Total Length: 5,000,000 bp
```

### Read Validation Results

**Metadata displayed:**
- Total read count
- NGS type (illumina, ont, pacbio)
- Paired-end information (if detected)
  - Base name
  - Read number (R1/R2)
- **Statistics** (NEW):
  - Total bases
  - Mean read length
  - N50
  - Longest read length
  - Shortest read length
- Validation level used

**Example:**
```
Output Metadata:
  Reads: 1,000,000
  NGS Type: illumina
  Paired-End: R1 (base: reads)
  Total Bases: 150,000,000 bp
  Mean Read Length: 150.0 bp
  N50: 155 bp
  Longest Read: 200 bp
  Shortest Read: 100 bp
```

**Note**: Statistics are only available in `strict` validation mode. In `trust` or `minimal` mode, only basic information (read count, NGS type, pairing) is shown.

### Feature Validation Results

**Metadata displayed:**
- Number of features
- Feature types (gene, mRNA, exon, CDS, etc.)
- Sequence IDs covered
- Validation level used

**Example:**
```
Output Metadata:
  Features: 5,000
    Types: gene, mRNA, exon, CDS
    Sequences: chr1, chr2, chr3
```

### Inter-file Validation Results

**Information displayed:**
- Validation type
- Status (PASSED/FAILED)
- List of errors (if any)
- List of warnings (if any)
- Validation metadata

**Example:**
```
[1] GENOME×GENOME
----------------------------------------------------------------------------------------------------
Status: ✓ PASSED
Warnings:
  - Sequence length mismatch: chr1 (reference: 2,000,000 bp, modified: 2,000,100 bp)
  - Minor length difference detected (0.005%)
```

---

## Examples

### Example 1: Basic Validation Report

```python
from pathlib import Path
from validation_pkg import ConfigManager, validate_genome
from validation_pkg.report import ValidationReport

# Load config
config = ConfigManager.load("config.json")

# Create report
report = ValidationReport(Path("reports/validation.txt"))

# Validate and record
result = validate_genome(config.ref_genome, config.output_dir)
report.write(result, file_type="genome", input_file="genome.fasta")

# Generate report
report.flush(format="text")
print(f"Report saved to: {report.report_path}")
```

### Example 2: Complete Pipeline Report

```python
from pathlib import Path
from validation_pkg import (
    ConfigManager,
    validate_genome,
    validate_reads,
    validate_feature
)
from validation_pkg.report import ValidationReport
from validation_pkg.validators.interfile_genome import genomexgenome_validation
import time

# Load config
config = ConfigManager.load("config.json")

# Create report
report = ValidationReport(Path("reports/full_pipeline.txt"))

# Validate reference genome
start = time.time()
ref_result = validate_genome(config.ref_genome, config.output_dir)
ref_time = time.time() - start
report.write(ref_result, file_type="genome",
            input_file=config.ref_genome.filename,
            elapsed_time=ref_time)

# Validate modified genome
start = time.time()
mod_result = validate_genome(config.mod_genome, config.output_dir)
mod_time = time.time() - start
report.write(mod_result, file_type="genome",
            input_file=config.mod_genome.filename,
            elapsed_time=mod_time)

# Validate reads
start = time.time()
read_results = validate_reads(config.reads, config.output_dir)
read_time = time.time() - start
report.write(read_results, file_type="read", elapsed_time=read_time)

# Inter-file validation: genome comparison
interfile_result = genomexgenome_validation(ref_result, mod_result)
report.write(interfile_result, file_type="genomexgenome")

# Generate both formats
report.flush(format="text")
report.flush(format="json")

print(f"✓ Text report: {report.report_path}")
print(f"✓ JSON report: {report.report_path.with_suffix('.json')}")
```

### Example 3: Custom Settings in Report

```python
from pathlib import Path
from validation_pkg import ConfigManager, validate_genome, GenomeValidator
from validation_pkg.report import ValidationReport

config = ConfigManager.load("config.json")
report = ValidationReport(Path("reports/custom_settings.txt"))

# Custom settings
custom_settings = GenomeValidator.Settings(
    validation_level='trust',
    plasmid_split=True,
    min_sequence_length=500
)

# Validate with custom settings
result = validate_genome(config.ref_genome, config.output_dir, custom_settings)

# Record with settings for report
report.write(
    result,
    file_type="genome",
    input_file=config.ref_genome.filename,
    settings=custom_settings  # Settings will be serialized in report
)

report.flush(format="text")
```

---

## API Reference

### ValidationReport Class

#### `__init__(report_path: Path)`

Create a new validation report.

**Parameters:**
- `report_path` (Path): Path where the report will be written (directory created if needed)

**Example:**
```python
report = ValidationReport(Path("reports/validation.txt"))
```

---

#### `write(result, file_type, input_file=None, settings=None, elapsed_time=None)`

Add a validation result to the report.

**Parameters:**
- `result` (dict | List[dict]): Result from validator or list of results
- `file_type` (str): Type of result
  - File validations: `"genome"`, `"read"`, `"feature"`
  - Inter-file validations: `"genomexgenome"`, `"readxread"`, `"featurexgenome"`
- `input_file` (str, optional): Original input filename (inferred from output if not provided)
- `settings` (Settings object, optional): Settings object used for validation
- `elapsed_time` (float, optional): Validation time in seconds

**Example:**
```python
report.write(
    genome_result,
    file_type="genome",
    input_file="genome.fasta",
    settings=my_settings,
    elapsed_time=12.5
)
```

---

#### `flush(format="text")`

Generate and write the final report to file.

**Parameters:**
- `format` (str): Report format - `"text"` (default) or `"json"`

**Raises:**
- `ValueError`: If format is not "text" or "json"

**Example:**
```python
report.flush(format="text")  # Creates .txt file
report.flush(format="json")  # Creates .json file
```

---

### FileValidationRecord Dataclass

Represents a single file's validation result.

**Fields:**
- `input_file` (str): Original input filename
- `validator_type` (str): "genome", "read", or "feature"
- `input_settings` (dict | None): Serialized settings object
- `output_file` (str | None): Output file path
- `output_metadata` (dict): Validation metadata and statistics
- `elapsed_time` (float | None): Time taken for validation
- `timestamp` (str): ISO 8601 timestamp

---

### InterFileValidationRecord Dataclass

Represents an inter-file validation check.

**Fields:**
- `validation_type` (str): "genomexgenome", "readxread", or "featurexgenome"
- `status` (str): "PASSED" or "FAILED"
- `passed` (bool): True if validation passed
- `errors` (List[str]): List of error messages
- `warnings` (List[str]): List of warning messages
- `metadata` (dict): Validation-specific metadata
- `timestamp` (str): ISO 8601 timestamp

---

## Best Practices

### 1. Create Reports Early
Create the report object at the start of your pipeline:
```python
report = ValidationReport(Path("reports/validation.txt"))
```

### 2. Record Timing Information
Track validation times for performance analysis:
```python
import time
start = time.time()
result = validate_genome(config, output_dir)
elapsed = time.time() - start
report.write(result, file_type="genome", elapsed_time=elapsed)
```

### 3. Include Settings
Record settings used for reproducibility:
```python
report.write(
    result,
    file_type="genome",
    settings=custom_settings  # Important for reproducibility!
)
```

### 4. Generate Both Formats
Create both text (human) and JSON (machine) formats:
```python
report.flush(format="text")  # For humans
report.flush(format="json")  # For automation
```

### 5. Organize Reports by Date
Use date-based directory structure:
```python
from datetime import datetime
date_str = datetime.now().strftime("%Y%m%d_%H%M%S")
report = ValidationReport(Path(f"reports/{date_str}/validation.txt"))
```

---

## Troubleshooting

### Report shows "unknown" for input files
**Problem**: Input filename not specified and couldn't be inferred.

**Solution**: Always provide `input_file` parameter:
```python
report.write(result, file_type="genome", input_file="genome.fasta")
```

### Missing statistics in read reports
**Problem**: Statistics only calculated in `strict` validation mode.

**Solution**: Use `validation_level='strict'` for complete statistics:
```python
settings = ReadValidator.Settings(validation_level='strict')
result = validate_read(config, output_dir, settings)
```

### Settings not showing in report
**Problem**: Settings parameter not provided to `write()`.

**Solution**: Pass the settings object:
```python
report.write(result, file_type="genome", settings=my_settings)
```

### All read files show same input filename
**Problem**: This was a bug in earlier versions (fixed in v0.3.0).

**Solution**: Upgrade to latest version or provide `input_file=None` to auto-infer per-file.

---

## See Also

- [CONFIG_GUIDE.md](CONFIG_GUIDE.md) - Configuration file format
- [LOGGING_GUIDE.md](LOGGING_GUIDE.md) - Logging system documentation
- [VALIDATOR_DESIGN.md](VALIDATOR_DESIGN.md) - Validator architecture
- Main README - API reference and examples
