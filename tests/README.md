# Validation Package Test Suite

Complete testing documentation for the validation package.

## Overview

This test suite provides comprehensive coverage with **156 tests** across unit and integration testing:
- **Unit Tests**: 126 tests
- **Integration Tests**: 30 tests (using 20 fixture directories)
- **Test Coverage**: ConfigManager, GenomeValidator, ReadValidator, FeatureValidator, Logger

---

## Table of Contents

1. [Test Types](#test-types)
2. [Unit Tests](#unit-tests)
3. [Integration Tests](#integration-tests)
4. [Running Tests](#running-tests)
5. [Creating New Tests](#creating-new-tests)
6. [Fixture Management](#fixture-management)
7. [Quick Reference](#quick-reference)

---

## Test Types

### Unit Tests
- **Purpose**: Test individual methods and classes in isolation
- **Naming**: `test_{module_name}.py`
- **Location**: `tests/`
- **Example**: `test_config_manager.py` tests the ConfigManager class

### Integration Tests
- **Purpose**: Test complex pipelines across the entire package
- **Naming**: Fixture directories (`config_test_*` or `test_case_*`)
- **Location**: `tests/fixtures/`
- **Runner**: `test_integration.py`
- **Example**: Load config → Run validator → Verify output files

---

## Unit Tests

Unit tests directly test methods and classes within individual modules.

### Test Files (126 tests total)

| Test File | Tests | Module Tested | Description |
|-----------|-------|---------------|-------------|
| **test_config_manager.py** | 45 | config_manager.py | Configuration loading, validation, path resolution |
| **test_genome_validator.py** | 24 | validators/genome_validator.py | Genome parsing, validation, plasmid splitting |
| **test_read_validator.py** | 23 | validators/read_validator.py | Read file parsing, quality control, NGS types |
| **test_feature_validator.py** | 24 | validators/feature_validator.py | Feature file parsing, format conversion |
| **test_logger.py** | 13 | logger.py | Logging, issue tracking, report generation |

### Unit Test Categories

#### test_config_manager.py (45 tests)
```
TestConfigManager (26 tests)
├── Configuration loading (minimal, full)
├── Missing field validation
├── File existence checking
├── Directory-based reads
├── Compression detection
└── Format detection

TestConfigManagerUtilities (15 tests)
├── Compression type detection
├── File format detection
└── Config value parsing

TestConfigManagerOutputDirectory (4 tests)
├── Output directory creation
├── Path structure validation
└── Config object consistency
```

#### test_genome_validator.py (24 tests)
```
TestGenomeValidatorInitialization (2 tests)
TestGenomeValidatorParsing (4 tests)
TestGenomeValidatorCompression (2 tests)
TestGenomeValidatorValidation (1 test)
TestGenomeValidatorEditing (2 tests)
TestGenomeValidatorStatistics (1 test)
TestGenomeValidatorOutput (6 tests)
TestGenomeValidatorPlasmidSplit (6 tests)
```

#### test_read_validator.py (23 tests)
```
TestReadValidatorInitialization (2 tests)
TestReadValidatorParsing (3 tests)
TestReadValidatorCompression (2 tests)
TestReadValidatorValidation (4 tests)
TestReadValidatorNGSTypes (3 tests)
TestReadValidatorOutput (5 tests)
TestReadValidatorStatistics (1 test)
TestReadValidatorEdgeCases (2 tests)
TestReadValidatorBAMHandling (1 test)
```

#### test_feature_validator.py (24 tests)
```
TestFeatureValidatorInitialization (2 tests)
TestFeatureValidatorParsing (4 tests)
TestFeatureValidatorCompression (2 tests)
TestFeatureValidatorValidation (4 tests)
TestFeatureValidatorEditing (2 tests)
TestFeatureValidatorOutput (5 tests)
TestFeatureValidatorStatistics (1 test)
TestFeatureValidatorFormatConversion (1 test)
TestFeatureValidatorEdgeCases (3 tests)
```

#### test_logger.py (13 tests)
```
TestLogger (13 tests)
├── Logger initialization and singleton
├── Validation issue tracking
├── Report generation
├── Summary statistics
└── File/directory creation
```

### Running Unit Tests

```bash
# Run all unit tests
pytest tests/test_*.py -v

# Run specific unit test file
pytest tests/test_config_manager.py -v

# Run specific test class
pytest tests/test_config_manager.py::TestConfigManager -v

# Run specific test
pytest tests/test_config_manager.py::TestConfigManager::test_load_minimal_valid_config -v

# Run with pattern matching
pytest tests/ -k "compression" -v
pytest tests/ -k "validation" -v

# Run unit tests only (exclude integration)
pytest tests/ --ignore=tests/test_integration.py --ignore=tests/test_config_manager_integration.py --ignore=tests/test_genome_validator_integration.py -v
```

---

## Integration Tests

Integration tests verify complex workflows across the package using realistic test fixtures.

### Test Structure (30 tests, 20 fixtures)

```
tests/
├── test_integration.py                          # Master runner (runs all integration tests)
├── test_config_manager_integration.py          # Config integration tests (17 tests)
├── test_genome_validator_integration.py        # Genome integration tests (13 tests)
└── fixtures/
    ├── config_test_01_minimal_valid/           # Config fixture 1
    │   ├── data/
    │   │   ├── config.json
    │   │   ├── ref_genome.fasta
    │   │   ├── mod_genome.fasta
    │   │   └── reads.fastq
    │   └── description.txt
    ├── config_test_02_full_config/             # Config fixture 2
    ├── ...                                      # (10 config fixtures total)
    ├── test_case_01_simple_fasta/              # Validator fixture 1
    │   ├── genome.fasta
    │   └── description.txt
    ├── test_case_02_multi_sequence/            # Validator fixture 2
    └── ...                                      # (10 validator fixtures total)
```

### Fixture Naming Convention

#### Config Manager Fixtures: `config_test_##_description`
```
config_test_01_minimal_valid
config_test_02_full_config
config_test_03_directory_reads
config_test_04_compressed
config_test_05_multiple_ngs_types
config_test_06_missing_config          (negative test)
config_test_07_missing_required_field  (negative test)
config_test_08_missing_file            (negative test)
config_test_09_invalid_json            (negative test)
config_test_10_invalid_ngs_type        (negative test)
```

#### Validator Fixtures: `test_case_##_description`
```
test_case_01_simple_fasta
test_case_02_multi_sequence
test_case_03_compressed_gzip
test_case_04_compressed_bzip2
test_case_05_genbank
test_case_06_mixed_lengths
test_case_07_high_gc
test_case_08_id_replacement
test_case_09_empty_sequence            (negative test)
test_case_10_complex_plasmids
```

### Integration Test Coverage

#### Config Manager Integration (17 tests)
Tests complete configuration loading workflow:
1. Valid configurations (minimal, full, compressed)
2. Directory-based reads discovery
3. Multiple NGS types
4. Error handling (missing files, invalid JSON, bad values)
5. Path resolution
6. Format detection
7. Output directory setup

#### Genome Validator Integration (13 tests)
Tests complete validation pipeline:
1. File format parsing (FASTA, GenBank)
2. Compression handling (gzip, bzip2)
3. Sequence filtering by length
4. Plasmid splitting with multiple plasmids
5. Statistics collection (GC content, lengths)
6. ID replacement
7. Output generation (compression, subdirectories, suffixes)
8. Error handling (empty sequences)

### Running Integration Tests

```bash
# Run all integration tests (RECOMMENDED)
pytest tests/test_integration.py -v

# Run only config integration tests
pytest tests/test_integration.py -k config -v

# Run only validator integration tests
pytest tests/test_integration.py -k genome -v

# Run specific integration test
pytest tests/test_integration.py::TestConfigIntegration::test_config_01_minimal_valid -v

# Run all integration test files separately
pytest tests/test_config_manager_integration.py tests/test_genome_validator_integration.py -v

# Verbose with detailed output
pytest tests/test_integration.py -vv --tb=short
```

---

## Running Tests

### Run Everything

```bash
# Run ALL tests (unit + integration)
pytest tests/ -v

# Run with coverage report
pytest tests/ --cov=validation_pkg --cov-report=html

# Quick run (stop on first failure)
pytest tests/ -x

# Parallel execution (if pytest-xdist installed)
pytest tests/ -n auto
```

### Run by Category

```bash
# Unit tests only
pytest tests/test_config_manager.py tests/test_genome_validator.py tests/test_read_validator.py tests/test_feature_validator.py tests/test_logger.py -v

# Integration tests only
pytest tests/test_integration.py -v

# Specific module
pytest tests/test_config_manager.py -v
```

### Run by Pattern

```bash
# All compression-related tests
pytest tests/ -k "compression" -v

# All validation tests
pytest tests/ -k "validation" -v

# All output tests
pytest tests/ -k "output" -v

# Specific test case number
pytest tests/ -k "test_case_01" -v
```

### Debug Options

```bash
# Show print statements
pytest tests/ -s

# Show local variables on failure
pytest tests/ -l

# Very verbose
pytest tests/ -vv

# Drop into debugger on failure
pytest tests/ --pdb

# Stop on first failure
pytest tests/ -x
```

---

## Creating New Tests

### Creating Unit Tests

1. **Create test file**: `tests/test_{module_name}.py`
2. **Import module**: Import the module/class to test
3. **Create test class**: Group related tests
4. **Write tests**: Name tests as `test_{what_it_tests}`
5. **Use fixtures**: For setup/teardown

**Example:**
```python
import pytest
from validation_pkg.my_module import MyClass

class TestMyClass:
    """Tests for MyClass."""

    @pytest.fixture
    def my_object(self):
        """Create test object."""
        return MyClass()

    def test_method_returns_expected_value(self, my_object):
        """Test that method returns expected value."""
        result = my_object.my_method()
        assert result == "expected"

    def test_method_raises_error_on_invalid_input(self, my_object):
        """Test error handling."""
        with pytest.raises(ValueError):
            my_object.my_method(invalid_input)
```

### Creating Integration Tests

#### 1. Create Fixture Directory

**For config tests:**
```bash
mkdir tests/fixtures/config_test_##_description
mkdir tests/fixtures/config_test_##_description/data
```

**For validator tests:**
```bash
mkdir tests/fixtures/test_case_##_description
```

#### 2. Add Test Files

Create realistic input files in the fixture directory.

#### 3. Add Description

Create `description.txt`:
```
Test Case ##: Brief Title
- Description of what this tests
- Expected behavior
- Special considerations
```

#### 4. Add Test in Integration Test File

Add test method to `test_config_manager_integration.py` or `test_genome_validator_integration.py`:

```python
def test_config_##_description(self, fixtures_dir):
    """Test Case ##: Brief description."""
    test_case = fixtures_dir / "config_test_##_description" / "data"
    config_file = test_case / "config.json"

    config = ConfigManager.load(str(config_file))

    # Assertions
    assert config.ref_genome is not None
    # ... more assertions ...
```

#### 5. Run Tests

```bash
pytest tests/test_integration.py -k "config_##" -v
```

### Regenerating Fixtures

If you modify fixture generation scripts:

```bash
cd tests/fixtures

# Regenerate config fixtures
./generate_config_fixtures.sh

# Regenerate validator fixtures
./generate_fixtures.sh

# Or use the utility script
python test_cases.py --regenerate
```

---

## Fixture Management

### Fixture Generation Scripts

Located in `tests/fixtures/`:

1. **generate_config_fixtures.sh**
   - Generates 10 config test fixtures
   - Creates config.json and associated files
   - Includes positive and negative test cases

2. **generate_fixtures.sh**
   - Generates 10 validator test fixtures
   - Creates genome files in various formats
   - Includes compression variants

### Regenerating Fixtures

```bash
# From tests/fixtures directory
./generate_config_fixtures.sh
./generate_fixtures.sh

# Or using utility script
python test_cases.py --regenerate
```

### Fixture Structure

Each fixture directory should contain:

- **description.txt**: What the test case does
- **Input files**: All necessary data files
- **data/**: Subdirectory for config tests (optional)

### Adding New Fixtures

1. Edit the generation script (`generate_config_fixtures.sh` or `generate_fixtures.sh`)
2. Add new fixture case following the pattern
3. Run generation script
4. Add corresponding test in integration test file
5. Verify test passes

---

## Quick Reference

### Test Statistics

```
Total Tests: 156
├── Unit Tests: 126
│   ├── test_config_manager.py: 45
│   ├── test_genome_validator.py: 24
│   ├── test_read_validator.py: 23
│   ├── test_feature_validator.py: 24
│   └── test_logger.py: 13
└── Integration Tests: 30
    ├── Config integration: 17
    └── Genome integration: 13

Test Fixtures: 20
├── Config fixtures: 10
└── Validator fixtures: 10
```

### Essential Commands

```bash
# Run all tests
pytest tests/ -v

# Run unit tests only
pytest tests/test_*.py --ignore=tests/test_integration.py -v

# Run integration tests only
pytest tests/test_integration.py -v

# Run specific module
pytest tests/test_config_manager.py -v

# Run with pattern
pytest tests/ -k "compression" -v

# Quick run (stop on first failure)
pytest tests/ -x

# Debug mode
pytest tests/ --pdb -x

# Coverage report
pytest tests/ --cov=validation_pkg --cov-report=html
```

### File Naming Conventions

| Type | Pattern | Example |
|------|---------|---------|
| Unit test file | `test_{module}.py` | `test_config_manager.py` |
| Integration runner | `test_integration.py` | (single file) |
| Config fixture | `config_test_##_description` | `config_test_01_minimal_valid` |
| Validator fixture | `test_case_##_description` | `test_case_01_simple_fasta` |

### Best Practices

1. **Unit tests**: Test one thing per test method
2. **Integration tests**: Test complete workflows
3. **Fixtures**: Use realistic data
4. **Naming**: Be descriptive and follow conventions
5. **Documentation**: Update this README when adding tests
6. **Run often**: Run tests before committing
7. **Coverage**: Aim for high code coverage
8. **Maintenance**: Keep fixtures up to date

---

## Troubleshooting

### Tests Failing

1. **Check dependencies**: `pip install -r requirements.txt`
2. **Regenerate fixtures**: `cd tests/fixtures && ./generate_config_fixtures.sh && ./generate_fixtures.sh`
3. **Check file permissions**: Ensure scripts are executable
4. **Run verbose**: `pytest tests/ -vv` for detailed output

### Missing Fixtures

```bash
cd tests/fixtures
./generate_config_fixtures.sh
./generate_fixtures.sh
```

### Permission Errors

```bash
chmod +x tests/fixtures/*.sh
```

### Import Errors

```bash
# Run from project root
cd /path/to/validation_pkg
pytest tests/ -v
```

---

## Contributing

When adding new functionality:

1. **Add unit tests** for new methods/classes
2. **Add integration tests** if testing workflows
3. **Update this README** with new test information
4. **Ensure all tests pass** before committing
5. **Maintain >90% coverage** when possible

---

**Last Updated**: 2025-10-15
**Test Count**: 156 (100% pass rate)
**Maintained by**: Validation Package Team
