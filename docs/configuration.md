# validation_pkg

A comprehensive Python package for validating and processing bioinformatics files including genomes, features, and sequencing reads.

## Table of Contents

- [Overview](#overview)
- [Installation](#installation)
- [Quick Start](#quick-start)
- [Configuration](#configuration)
- [Architecture](#architecture)
- [Error Handling](#error-handling)
- [Logging](#logging)
- [Testing](#testing)
- [Development](#development)

---

## Overview

`validation_pkg` is designed to validate bioinformatics data files and ensure consistency across multiple file types. The package handles:

- **Genome files**: FASTA and GenBank formats
- **Feature files**: BED, GFF, and GTF formats
- **Read files**: FASTQ format (Illumina, Oxford Nanopore, PacBio)
- **Compressed files**: Transparent handling of `.gz`, `.bz2` and `tar.gz` compression
- **Inter-file validation**: Ensuring consistency between genomes and features

### Key Features

✅ **Configuration-driven**: All inputs specified in a single JSON configuration file  
✅ **Comprehensive validation**: Format checking, consistency verification, and quality control  
✅ **Detailed logging**: Multiple logging levels with file and console output  
✅ **Structured error handling**: Clear, actionable error messages  
✅ **Automated testing**: Extensive test suite with fixture-based integration tests  
✅ **Compression support**: Automatic detection and handling of compressed files  

---

## Installation

### For Users

```bash
pip install validation_pkg
```

### For Development

```bash
# Clone the repository
git clone https://github.com/yourusername/validation_pkg.git
cd validation_pkg

# Install in editable mode
pip install -e .

---

## Quick Start

### 1. Create a Configuration File

Create a `config.json` file specifying your input files:

```json
{
  "ref_genome_filename": {"filename": "data/reference.fasta"},
  "mod_genome_filename": {"filename": "data/modified.fasta"},
  "reads": [
    {"filename": "data/reads_R1.fastq.gz", "ngs_type": "illumina"},
    {"filename": "data/reads_R2.fastq.gz", "ngs_type": "illumina"}
  ],
  "ref_feature_filename": {"filename": "data/features.gff"}
}
```

### 2. Run Validation

```bash
python main.py config.json
```

### 3. Check Results

The package will:
- Validate all input files
- Generate a detailed log file in `logs/`
- Create a validation report in `logs/`
- Exit with code 0 (success) or 1 (failure)

---

## Configuration

### Configuration File Structure

The configuration file (`config.json`) defines all inputs and processing options.

#### Required Fields

| Field | Type | Description |
|-------|------|-------------|
| `ref_genome_filename` | GenomeConfig | Reference genome (FASTA or GenBank) |
| `mod_genome_filename` | GenomeConfig | Modified genome (FASTA or GenBank) |
| `reads` | List[ReadConfig] | At least one read file or directory |

#### Optional Fields

| Field | Type | Description |
|-------|------|-------------|
| `ref_plasmid_filename` | GenomeConfig | Reference plasmid (FASTA or GenBank) |
| `mod_plasmid_filename` | GenomeConfig | Modified plasmid (FASTA or GenBank) |
| `ref_feature_filename` | FeatureConfig | Reference features (BED, GFF, or GTF) |
| `mod_feature_filename` | FeatureConfig | Modified features (BED, GFF, or GTF) |
| `options` | Object | Additional processing options |

#### GenomeConfig Object

```json
{
  "filename": "path/to/genome.fasta",
  "output_coding": "gz"  // "none", "gz", "bzip", or "tgz"
}
```

#### ReadConfig Object

Option 1: Single file
```json
{
  "filename": "path/to/reads.fastq",
  "ngs_type": "illumina",  // "illumina", "ont", or "pacbio"
  "output_coding": "gz"
}
```

Option 2: Directory of files
```json
{
  "directory": "path/to/reads/",
  "ngs_type": "illumina"
}
```

#### FeatureConfig Object

```json
{
  "filename": "path/to/features.gff",
  "output_coding": "none"
}
```

### Supported File Formats

**Genomes**: `.fa`, `.fasta`, `.fna`, `.gb`, `.gbk`, `.genbank`  
**Features**: `.bed`, `.gff`, `.gff3`, `.gtf`  
**Reads**: `.fastq`, `.fq`  
**Compression**: Any of the above with `.gz` or `.bz2` suffix

### Example Configurations

**Minimal configuration:**
```json
{
  "ref_genome_filename": {"filename": "ref.fasta"},
  "mod_genome_filename": {"filename": "mod.fasta"},
  "reads": [
    {"filename": "reads.fastq", "ngs_type": "illumina"}
  ]
}
```

**Full configuration with all options:**
```json
{
  "ref_genome_filename": {"filename": "ref.gbk", "output_coding": "gz"},
  "mod_genome_filename": {"filename": "mod.fasta.gz"},
  "ref_plasmid_filename": {"filename": "plasmid_ref.gbk"},
  "mod_plasmid_filename": {"filename": "plasmid_mod.fasta"},
  "reads": [
    {"filename": "illumina_R1.fastq.gz", "ngs_type": "illumina", "output_coding": "gz"},
    {"filename": "illumina_R2.fastq.gz", "ngs_type": "illumina"},
    {"filename": "nanopore.fastq", "ngs_type": "ont"}
  ],
  "ref_feature_filename": {"filename": "features_ref.gff3"},
  "mod_feature_filename": {"filename": "features_mod.bed"},
  "options": {
    "threads": 8,
    "verbose": true
  }
}
```

---

## Architecture

### Package Structure

```
validation_pkg/
├── main.py                    # Entry point script
├── coordinator.py             # Configuration loader and validator
├── logger.py                  # Logging system
├── exceptions.py              # Custom exception classes
├── utils/
│   └── file_handler.py       # File I/O and compression handling
└── validators/               # (In development)
    ├── genome_validator.py   # Genome file validation
    ├── feature_validator.py  # Feature file validation
    └── read_validator.py     # Read file validation
```

### Processing Pipeline

The validation pipeline consists of five main stages:

```
┌─────────────────────────────────────────────────────────────┐
│ 1. Configuration Loading                                    │
│    - Parse JSON configuration                               │
│    - Validate required fields                               │
│    - Check file existence                                   │
└─────────────────────────────────────────────────────────────┘
                            ↓
┌─────────────────────────────────────────────────────────────┐
│ 2. Genome Validation                                        │
│    - Detect file format (FASTA/GenBank)                     │
│    - Validate format correctness                            │
│    - Handle compression (if needed)                         │
│    - Collect statistics                                     │
└─────────────────────────────────────────────────────────────┘
                            ↓
┌─────────────────────────────────────────────────────────────┐
│ 3. Feature Validation                                       │
│    - Detect file format (BED/GFF/GTF)                       │
│    - Validate format correctness                            │
│    - Check feature coordinates                              │
│    - Collect statistics                                     │
└─────────────────────────────────────────────────────────────┘
                            ↓
┌─────────────────────────────────────────────────────────────┐
│ 4. Read Validation                                          │
│    - Validate FASTQ format                                  │
│    - Check quality scores                                   │
│    - Verify paired-end consistency                          │
│    - Collect statistics                                     │
└─────────────────────────────────────────────────────────────┘
                            ↓
┌─────────────────────────────────────────────────────────────┐
│ 5. Inter-File Validation                                    │
│    - Check feature coordinates vs genome length             │
│    - Verify sequence ID consistency                         │
│    - Validate cross-references                              │
└─────────────────────────────────────────────────────────────┘
                            ↓
                    ┌───────────────┐
                    │ Final Report  │
                    └───────────────┘
```

### Component Descriptions

#### ConfigManager

**Purpose**: Configuration management and validation entry point

**Responsibilities**:
- Parse and validate JSON configuration
- Resolve file paths relative to config directory
- Check file existence before processing
- Instantiate validators

**Key Classes**:
- `Config`: Data holder for configuration
- `ConfigManager`: Configuration loader and validator
- `GenomeConfig`, `ReadConfig`, `FeatureConfig`: Type-safe configuration objects

#### Logger

**Purpose**: Centralized logging and issue tracking

**Features**:
- Colored console output for readability
- Detailed debug log file
- Structured validation report
- Issue categorization (ERROR, WARNING, INFO)
- Summary statistics

**Outputs**:
- Console: User-friendly colored messages
- Log file: Complete debug information with timestamps
- Report file: Summary of all validation issues

#### Exception System

**Purpose**: Clear, structured error handling

**Exception Hierarchy**:
```
ValidationError (base)
├── ConfigurationError
├── FileNotFoundError
├── FileFormatError
│   ├── FastaFormatError
│   ├── GenBankFormatError
│   ├── BedFormatError
│   ├── GffFormatError
│   └── FastqFormatError
├── CompressionError
├── GenomeValidationError
├── FeatureValidationError
├── ReadValidationError
└── InterFileValidationError
```

---

## Error Handling

### Error Categories

The package uses custom exceptions to provide clear, actionable error messages:

#### 1. Configuration Errors

**Cause**: Issues in the configuration file  
**Examples**:
- Missing required fields
- Invalid values (e.g., unsupported `ngs_type`)
- Malformed JSON
- Referenced files don't exist

**Error Message Example**:
```
ConfigurationError: Missing required field: ref_genome_filename
```

#### 2. File Format Errors

**Cause**: Files exist but have invalid format  
**Examples**:
- Invalid FASTA headers
- Malformed GenBank records
- Wrong number of BED columns
- Invalid FASTQ quality scores

**Error Message Example**:
```
FastaFormatError: Invalid FASTA format at line 42: expected '>' but found 'A'
```

#### 3. Inter-File Validation Errors

**Cause**: Inconsistencies between files  
**Examples**:
- Feature coordinates exceed genome length
- Sequence IDs don't match
- Missing referenced sequences

**Error Message Example**:
```
InterFileValidationError: Feature 'gene1' at position 5000 exceeds genome length of 4500
```

### How Errors Are Reported

When an error occurs:

1. **Logged to console** with appropriate level (ERROR/WARNING)
2. **Added to validation report** with full details
3. **Recorded in debug log** with timestamp and context
4. **Exception raised** with clear message
5. **Exit code set** to 1 (failure)

### Example Error Flow

```python
try:
    config = ConfigManager.load("config.json")
except ConfigurationError as e:
    # User sees: "ConfigurationError: Missing required field: reads"
    # Log shows: Full stack trace and config file location
    # Report shows: Structured issue with details
    pass
```

---

## Logging

### Logging Levels

The package uses standard Python logging levels:

| Level | Console | Log File | Purpose |
|-------|---------|----------|---------|
| DEBUG | ❌ | ✅ | Detailed diagnostic information |
| INFO | ✅ | ✅ | Progress updates and confirmations |
| WARNING | ✅ | ✅ | Non-critical issues that don't prevent validation |
| ERROR | ✅ | ✅ | Critical issues that cause validation to fail |

### Log Outputs

#### 1. Console Output

**Purpose**: Real-time feedback to the user  
**Format**: Colored, user-friendly messages  
**Content**: Progress updates, warnings, and errors

**Example**:
```
INFO - [1/4] Loading configuration...
INFO - ✓ Configuration loaded successfully
ERROR - [genome] Invalid FASTA format at line 42
```

#### 2. Debug Log File

**Purpose**: Complete diagnostic information  
**Location**: `logs/validation_YYYYMMDD_HHMMSS.log`  
**Format**: Timestamped, detailed entries  
**Content**: All messages including DEBUG level

**Example**:
```
2025-10-13 14:30:15 - validation_pkg - INFO - load:125 - Loading configuration from: config.json
2025-10-13 14:30:15 - validation_pkg - DEBUG - load:130 - Reading configuration file...
2025-10-13 14:30:15 - validation_pkg - ERROR - validate:250 - Invalid FASTA format
```

#### 3. Validation Report

**Purpose**: Summary of all validation issues  
**Location**: `logs/report_YYYYMMDD_HHMMSS.txt`  
**Format**: Structured, human-readable report  
**Content**: Categorized issues with details

**Example**:
```
================================================================================
GMO VALIDATION REPORT
================================================================================
Generated: 2025-10-13 14:30:20

SUMMARY
--------------------------------------------------------------------------------
Total Issues: 2
  Errors:   1
  Warnings: 1

DETAILS
--------------------------------------------------------------------------------

1. [ERROR] genome
   Invalid FASTA format
   - file: ref.fasta
   - line: 42
   - expected: Sequence line
   - found: Invalid character 'X'

2. [WARNING] feature
   Feature outside genome bounds
   - feature_id: gene1
   - position: 5000
   - genome_length: 4500

================================================================================
```

### Using the Logger

#### In Your Code

```python
from logger import get_logger

logger = get_logger()

# Simple logging
logger.info("Processing genome file...")
logger.warning("Unusual GC content detected")
logger.error("Validation failed")

# Structured issue tracking
logger.add_validation_issue(
    level='ERROR',
    category='genome',
    message='Invalid FASTA format',
    details={
        'file': 'ref.fasta',
        'line': 42,
        'reason': 'Unexpected character'
    }
)
```

#### Command Line Options

```bash
# Standard output
python main.py config.json

# Verbose output (includes DEBUG messages)
python main.py config.json --verbose

# Custom log directory
python main.py config.json --log-dir results/

# Console only (no log files)
python main.py config.json --no-log-file --no-report
```

---

## Testing

### Test Structure

```
tests/
├── test_coordinator.py      # Unit tests for configuration
├── test_logger.py           # Unit tests for logging
├── test_cases.py            # Integration tests
├── fixtures/                # Test data
│   ├── valid_minimal/
│   ├── valid_full/
│   ├── invalid_*/
│   └── ...
└── create_fixtures.sh       # Script to generate test data
```

### Test Categories

#### 1. Unit Tests

**Purpose**: Test individual components in isolation  
**Files**: `test_coordinator.py`, `test_logger.py`  
**Approach**: Use temporary files and mock data

**Run unit tests**:
```bash
pytest tests/test_coordinator.py -v
pytest tests/test_logger.py -v
```

#### 2. Integration Tests

**Purpose**: Test complete workflow with real data  
**File**: `test_cases.py`  
**Approach**: Use fixture directories with actual config and data files

**Run integration tests**:
```bash
pytest tests/test_cases.py -v -s
```

### Test Fixtures

Test fixtures are organized by expected behavior:

- **`valid_*`**: Cases that should pass validation
- **`invalid_*`**: Cases that should fail with specific errors

**Example fixtures**:
- `valid_minimal`: Minimal valid configuration
- `valid_full`: Full configuration with all optional fields
- `valid_compressed_inputs`: Tests compressed file handling
- `invalid_missing_ref_genome`: Missing required field
- `invalid_bad_ngs_type`: Invalid NGS type value

### Creating Test Fixtures

```bash
# Generate all test fixtures
./create_fixtures.sh

# Or use the Python script
python create_fixtures.py
```

### Running All Tests

```bash
# Run all tests
pytest tests/ -v

# Run with coverage
pytest tests/ --cov=validation_pkg --cov-report=html

# Run specific test
pytest tests/test_cases.py::test_all_cases -v -s
```

### Test Naming Convention

When creating new test fixtures:

1. **Valid cases**: Name starts with `valid_`
   - Example: `valid_paired_end_reads`
   
2. **Invalid cases**: Name starts with `invalid_` or `error_`
   - Example: `invalid_missing_genome`

3. **Each fixture directory must contain**:
   - `config.json`: Configuration file
   - All referenced data files

---

## Development

### Setting Up Development Environment

```bash
# Clone repository
git clone https://github.com/yourusername/validation_pkg.git
cd validation_pkg

# Create virtual environment
python -m venv venv
source venv/bin/activate  # On Windows: venv\Scripts\activate

# Install in editable mode with dev dependencies
pip install -e .[dev]

# Create test fixtures
./create_fixtures.sh

# Run tests
pytest tests/ -v
```

### Development Workflow

1. **Create feature branch**
   ```bash
   git checkout -b feature/new-validator
   ```

2. **Make changes and test locally**
   ```bash
   # Run tests
   pytest tests/ -v
   
   # Test with your own data
   python main.py my_config.json --verbose
   ```

3. **Add tests for new features**
   - Create new test fixtures if needed
   - Add unit tests for new components
   - Update integration tests

4. **Push and create Pull Request**
   ```bash
   git push origin feature/new-validator
   ```

5. **GitHub Actions runs tests automatically**

### Adding New Validators

To add a new validator (e.g., for a new file type):

1. **Create validator class** in `validators/`
2. **Use logger** for progress and issues
3. **Raise appropriate exceptions** from `exceptions.py`
4. **Add unit tests** in `tests/`
5. **Add test fixtures** in `tests/fixtures/`
6. **Update documentation**

**Example template**:
```python
from logger import get_logger
from exceptions import ValidationError

class NewValidator:
    def __init__(self, config):
        self.logger = get_logger()
        self.config = config
    
    def validate(self):
        self.logger.info("Validating new file type...")
        # Your validation logic here
        self.logger.info("✓ Validation complete")
```

### Code Style

The project follows PEP 8 style guidelines:
- Use 4 spaces for indentation
- Maximum line length: 100 characters
- Use type hints where appropriate
- Document all public functions and classes

### Contributing

1. Fork the repository
2. Create a feature branch
3. Make your changes
4. Add tests
5. Ensure all tests pass
6. Submit a Pull Request

---

## Command Line Reference

### Basic Usage

```bash
python main.py <config.json>
```

### Options

| Option | Short | Description |
|--------|-------|-------------|
| `--verbose` | `-v` | Enable verbose output (DEBUG level) |
| `--log-dir DIR` | | Directory for log files (default: `logs/`) |
| `--no-report` | | Don't generate validation report |
| `--no-log-file` | | Console output only, no log file |
| `--help` | `-h` | Show help message |

### Examples

```bash
# Standard validation
python main.py config.json

# Verbose output with custom log directory
python main.py config.json --verbose --log-dir results/validation/

# Quick check (console only)
python main.py config.json --no-log-file --no-report

# Test with fixture
python main.py tests/fixtures/valid_minimal/config.json
```

---

## Troubleshooting

### Common Issues

**Issue**: `ConfigurationError: Missing required field`  
**Solution**: Check that your `config.json` includes all required fields: `ref_genome_filename`, `mod_genome_filename`, and `reads`.

**Issue**: `FileNotFoundError: File not found`  
**Solution**: Verify that all file paths in `config.json` are relative to the config file location and that the files exist.

**Issue**: `ValidationError: Invalid output_coding`  
**Solution**: Use only supported values: `"none"`, `"gz"`, `"bzip"`, or `"tgz"`.

**Issue**: Tests fail with "Fixtures directory not found"  
**Solution**: Run `./create_fixtures.sh` to generate test fixtures.

### Getting Help

- Check the validation report in `logs/` for detailed error information
- Run with `--verbose` flag for debug output
- Review the debug log file for complete diagnostic information
- Open an issue on GitHub with your config file and error message

---

## License

[Your License Here]

## Citation

If you use this package in your research, please cite:

```
[Citation information to be added]
```

## Contact

- **Repository**: https://github.com/yourusername/validation_pkg
- **Issues**: https://github.com/yourusername/validation_pkg/issues
- **Email**: your.email@example.com

---

**Last Updated**: October 2025  
**Version**: 0.1.0 (Early Development)