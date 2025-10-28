# Validator Design and Architecture

This document explains the internal architecture and design patterns used in the validation package. It's intended for developers who want to understand, extend, or contribute to the codebase.

## Table of Contents

- [Design Philosophy](#design-philosophy)
- [Core Design Pattern](#core-design-pattern)
- [Validation Architecture](#validation-architecture)
- [Settings Pattern](#settings-pattern)
- [Exception Hierarchy](#exception-hierarchy)
- [Format Detection System](#format-detection-system)
- [Extending Validators](#extending-validators)

---

## Design Philosophy

The validation package is built on several key principles:

1. **Configuration-driven**: All file inputs specified in a single JSON configuration
2. **Immutable settings**: Settings use copy-on-update pattern for safety
3. **Multi-level validation**: Balance between thoroughness and performance
4. **Clear error handling**: Structured exceptions with actionable messages
5. **Functional API**: Simple high-level functions for common tasks
6. **Extensibility**: Easy to add new file types and validation rules

---

## Core Design Pattern

### Validator + Settings Architecture

Each file type follows a consistent **validator pattern**:

```
┌─────────────────────────────────────────┐
│         Validator Class                 │
│  (GenomeValidator, ReadValidator, etc.) │
│                                         │
│  ┌───────────────────────────────────┐ │
│  │  Settings (nested dataclass)      │ │
│  │  - validation_level               │ │
│  │  - file-specific parameters       │ │
│  │  - immutable update pattern       │ │
│  └───────────────────────────────────┘ │
│                                         │
│  Methods:                               │
│  - __init__(config, settings)          │
│  - validate() → output_path            │
│  - _validate_strict()                  │
│  - _validate_trust()                   │
│  - _validate_minimal()                 │
└─────────────────────────────────────────┘
         ↓
┌─────────────────────────────────────────┐
│      Functional API Wrapper             │
│  validate_genome(config, output, ...)   │
│  validate_read(config, output, ...)     │
│  validate_feature(config, output, ...)  │
└─────────────────────────────────────────┘
```

### Three-Layer API

1. **Functional API** (recommended): Simple functions for common tasks
2. **Validator Classes**: Direct access for customization
3. **ConfigManager**: Configuration loading and path resolution

**Example:**

```python
# Layer 1: Functional API (simplest)
from validation_pkg import validate_genome, ConfigManager

config = ConfigManager.load("config.json")
output_path = validate_genome(config.ref_genome, config.output_dir)

# Layer 2: Validator Classes (custom settings)
from validation_pkg import GenomeValidator

settings = GenomeValidator.Settings()
settings = settings.update(validation_level='trust', plasmid_split=True)
validator = GenomeValidator(config.ref_genome, config.output_dir, settings)
output_path = validator.validate()

# Layer 3: ConfigManager only (load config, manual validation)
config = ConfigManager.load("config.json")
# ... manual validation logic
```

---

## Validation Architecture

### Validation Levels

The package implements a **three-tier validation system** to balance thoroughness vs. performance:

| Level | Parsing | Validation | Edits Applied | Performance | Use Case |
|-------|---------|------------|---------------|-------------|----------|
| **strict** | Full | All records | Yes | Slowest (baseline) | Production, critical data |
| **trust** | Partial* | First record only | Partial** | 10-15x faster | Pre-validated data |
| **minimal** | None | None | No (file copy) | Fastest | Known-good data |

\* **Genomes**: Parse all sequences, validate only first
\*\* **Reads**: No parsing, file copy only; **Genomes**: All edits applied

### Validation Level Implementation

Each validator's `validate()` method uses a switch statement:

```python
def validate(self) -> str:
    """Main validation entry point."""
    level = self.settings.validation_level

    if level == 'minimal':
        return self._validate_minimal()
    elif level == 'trust':
        return self._validate_trust()
    else:  # strict (default)
        return self._validate_strict()
```

**Key insight**: The validation level determines the entire code path, not just validation depth. This allows for dramatic performance improvements in `minimal` and `trust` modes.

### Validator Structure

All validators share a common structure:

```python
class GenomeValidator:
    """Validates and processes genome files."""

    @dataclass
    class Settings(BaseSettings):
        """Validator-specific settings (immutable)."""
        validation_level: str = 'strict'
        min_sequence_length: int = 0
        plasmid_split: bool = False
        # ... other settings

        def __post_init__(self):
            """Validate settings after initialization."""
            # Field validation logic

    def __init__(self, config: GenomeConfig, output_dir: str,
                 settings: Optional[Settings] = None):
        """Initialize validator with config and settings."""
        self.config = config
        self.output_dir = output_dir
        self.settings = settings or GenomeValidator.Settings()
        self.logger = get_logger()

    def validate(self) -> str:
        """Main validation entry point."""
        # Validation level switch
        # Returns output file path

    def _validate_strict(self) -> str:
        """Full validation with all edits."""
        # BioPython parsing
        # Validate all records
        # Apply all edits

    def _validate_trust(self) -> str:
        """Fast validation, first record only."""
        # Partial parsing
        # Validate first record
        # Apply edits (genomes) or copy (reads)

    def _validate_minimal(self) -> str:
        """No validation, just file copy."""
        # Direct file copy with compression handling
```

---

## Settings Pattern

### BaseSettings Class

All validator settings inherit from `BaseSettings` (in `utils/settings.py`):

```python
from dataclasses import dataclass, field, fields
from typing import Dict, Any

@dataclass
class BaseSettings:
    """Base class for immutable settings with update pattern."""

    def update(self, **kwargs) -> 'BaseSettings':
        """Return new instance with updated fields (immutable pattern)."""
        # Validate field names
        valid_fields = {f.name for f in fields(self)}
        invalid = set(kwargs.keys()) - valid_fields
        if invalid:
            raise ValueError(f"Invalid settings: {invalid}")

        # Create new instance with updated values
        current_values = {f.name: getattr(self, f.name) for f in fields(self)}
        current_values.update(kwargs)
        return self.__class__(**current_values)

    def to_dict(self) -> Dict[str, Any]:
        """Serialize to dictionary."""
        return {f.name: getattr(self, f.name) for f in fields(self)}

    @classmethod
    def from_dict(cls, data: Dict[str, Any]) -> 'BaseSettings':
        """Deserialize from dictionary."""
        valid_fields = {f.name for f in fields(cls)}
        filtered_data = {k: v for k, v in data.items() if k in valid_fields}
        return cls(**filtered_data)
```

### Immutable Update Pattern

**CRITICAL**: Settings are immutable. Always assign the result of `update()`:

```python
# ✅ CORRECT
settings = GenomeValidator.Settings()
settings = settings.update(validation_level='trust')  # Assign result!

# ❌ WRONG - Changes are lost!
settings = GenomeValidator.Settings()
settings.update(validation_level='trust')  # Doesn't modify original!
```

**Why immutable?**
- Prevents accidental modification during validation
- Allows safe sharing across parallel processes
- Makes settings changes explicit and trackable
- Enables rollback by keeping previous instances

### Field Validation

Settings validate themselves in `__post_init__()`:

```python
@dataclass
class Settings(BaseSettings):
    validation_level: str = 'strict'
    min_sequence_length: int = 0

    def __post_init__(self):
        """Validate settings after initialization."""
        # Validation level check
        valid_levels = ['strict', 'trust', 'minimal']
        if self.validation_level not in valid_levels:
            raise ValueError(
                f"validation_level must be one of {valid_levels}, "
                f"got '{self.validation_level}'"
            )

        # Min sequence length check
        if self.min_sequence_length < 0:
            raise ValueError(
                f"min_sequence_length must be >= 0, "
                f"got {self.min_sequence_length}"
            )
```

---

## Exception Hierarchy

### Exception Structure

The package uses a hierarchical exception system for clear error handling:

```
ValidationError (base)
├── ConfigurationError         # Config file issues
├── FileNotFoundError          # Missing files
├── FileFormatError            # Format parsing errors
│   ├── FastaFormatError
│   ├── GenBankFormatError
│   ├── GffFormatError
│   ├── GtfFormatError
│   ├── BedFormatError
│   └── FastqFormatError
├── CompressionError           # Compression/decompression errors
├── GenomeValidationError      # Genome-specific validation errors
├── ReadValidationError        # Read-specific validation errors
└── FeatureValidationError     # Feature-specific validation errors
```

### Exception Design Principles

1. **Catch at the right level**: Catch specific exceptions for targeted handling
2. **Actionable messages**: Include context (file, line, expected vs. actual)
3. **Preserve stack trace**: Use `raise ... from e` to chain exceptions
4. **Log before raising**: Use logger to record error context

**Example:**

```python
from exceptions import FastaFormatError

try:
    # Parse FASTA file
    for i, record in enumerate(SeqIO.parse(file_path, 'fasta')):
        if not record.id:
            raise FastaFormatError(
                f"Record {i} missing sequence ID in {file_path}"
            )
except Exception as e:
    self.logger.error(f"FASTA parsing failed: {e}")
    raise GenomeValidationError(
        f"Failed to validate {file_path}"
    ) from e
```

### Error Handling Best Practices

**In validators:**
- Catch specific exceptions (e.g., `FastaFormatError`)
- Add context before re-raising
- Log errors with `logger.error()`
- Track issues with `logger.add_validation_issue()`

**In user code:**
- Catch `ValidationError` to handle all package errors
- Catch specific exceptions for targeted recovery
- Check validation report for detailed issue list

---

## Format Detection System

### Two-Stage Detection

Format detection happens in two stages within `ConfigManager`:

**Stage 1: Compression Detection**
```python
def _detect_compression_type(filename: str) -> CodingType:
    """Detect compression from last extension."""
    ext = Path(filename).suffix.lower()

    if ext in ['.gz', '.gzip']:
        return CodingType.GZIP
    elif ext in ['.bz2', '.bzip2']:
        return CodingType.BZIP
    else:
        return CodingType.NONE
```

**Stage 2: File Format Detection**
```python
def _detect_file_format(filename: str) -> GenomeFormat:
    """Detect file format from first extension (after removing compression)."""
    # Remove compression extension first
    base = get_base_filename(filename)  # genome.fasta.gz → genome.fasta
    ext = Path(base).suffix.lower()

    # Use format enum's _missing_() for flexible matching
    return GenomeFormat(ext)
```

### Format Enums with Flexible Matching

Format enums use `_missing_()` for flexible extension handling:

```python
from enum import Enum

class GenomeFormat(str, Enum):
    """Genome file format enum with flexible extension matching."""
    FASTA = 'fasta'
    GENBANK = 'genbank'

    @classmethod
    def _missing_(cls, value):
        """Handle flexible extension matching."""
        value_lower = str(value).lower().lstrip('.')

        # FASTA variations
        if value_lower in ['fa', 'fasta', 'fna']:
            return cls.FASTA

        # GenBank variations
        if value_lower in ['gb', 'gbk', 'genbank']:
            return cls.GENBANK

        # Not found
        raise ValueError(f"Unknown genome format: {value}")
```

**Important**: Format enums only extract the **last** extension. For compressed files:

```python
# ❌ WRONG - Will fail (tries to parse .gz as format)
format = GenomeFormat('genome.fasta.gz')

# ✅ CORRECT - Strip compression first
from utils.file_handler import get_base_filename
base = get_base_filename('genome.fasta.gz')  # → 'genome.fasta'
format = GenomeFormat(base)  # → GenomeFormat.FASTA
```

### BioPython Format Mapping

Format enums provide BioPython format strings:

```python
class GenomeFormat(str, Enum):
    FASTA = 'fasta'
    GENBANK = 'genbank'

    def to_biopython(self) -> str:
        """Get BioPython format string."""
        return self.value  # 'fasta' or 'genbank'

# Usage
format = GenomeFormat.FASTA
SeqIO.parse(file_handle, format.to_biopython())
```

---

## Extending Validators

### Adding a New Validator

To add support for a new file type (e.g., VCF files):

**Step 1: Create Format Enum** (in `utils/formats.py`)

```python
class VcfFormat(str, Enum):
    """VCF file format enum."""
    VCF = 'vcf'
    BCF = 'bcf'

    @classmethod
    def _missing_(cls, value):
        value_lower = str(value).lower().lstrip('.')
        if value_lower == 'vcf':
            return cls.VCF
        elif value_lower == 'bcf':
            return cls.BCF
        raise ValueError(f"Unknown VCF format: {value}")
```

**Step 2: Create Config Class** (in `config_manager.py`)

```python
@dataclass
class VcfConfig:
    """Configuration for VCF files."""
    filename: str
    file_format: VcfFormat
    coding_type: CodingType
    settings_dict: Dict[str, Any] = field(default_factory=dict)
    extra: Dict[str, Any] = field(default_factory=dict)
```

**Step 3: Create Validator Class** (in `validators/vcf_validator.py`)

```python
from dataclasses import dataclass
from utils.settings import BaseSettings
from logger import get_logger

class VcfValidator:
    """Validates and processes VCF files."""

    @dataclass
    class Settings(BaseSettings):
        """VCF validator settings."""
        validation_level: str = 'strict'
        check_alt_alleles: bool = True
        min_quality: float = 0.0

        def __post_init__(self):
            """Validate settings."""
            valid_levels = ['strict', 'trust', 'minimal']
            if self.validation_level not in valid_levels:
                raise ValueError(
                    f"validation_level must be one of {valid_levels}"
                )
            if self.min_quality < 0:
                raise ValueError("min_quality must be >= 0")

    def __init__(self, config: VcfConfig, output_dir: str,
                 settings: Optional[Settings] = None):
        """Initialize validator."""
        self.config = config
        self.output_dir = output_dir
        self.settings = settings or VcfValidator.Settings()
        self.logger = get_logger()

    def validate(self) -> str:
        """Main validation entry point."""
        level = self.settings.validation_level

        if level == 'minimal':
            return self._validate_minimal()
        elif level == 'trust':
            return self._validate_trust()
        else:
            return self._validate_strict()

    def _validate_strict(self) -> str:
        """Full validation."""
        self.logger.info(f"Validating VCF: {self.config.filename}")
        # Full validation logic
        return output_path

    def _validate_trust(self) -> str:
        """Fast validation."""
        # Trust mode logic
        return output_path

    def _validate_minimal(self) -> str:
        """Minimal validation (copy only)."""
        # File copy with compression handling
        return output_path
```

**Step 4: Create Functional API** (in `__init__.py`)

```python
def validate_vcf(
    vcf_config: VcfConfig,
    output_dir: str,
    settings: Optional[VcfValidator.Settings] = None
) -> str:
    """Validate a VCF file (functional API)."""
    validator = VcfValidator(vcf_config, output_dir, settings)
    return validator.validate()
```

**Step 5: Add Tests** (in `tests/test_vcf_validator.py`)

```python
import pytest
from validation_pkg.validators.vcf_validator import VcfValidator

class TestVcfValidator:
    """Tests for VCF validator."""

    def test_settings_immutable(self):
        """Test settings immutable pattern."""
        settings = VcfValidator.Settings()
        new_settings = settings.update(min_quality=20.0)

        assert settings.min_quality == 0.0  # Original unchanged
        assert new_settings.min_quality == 20.0  # New instance updated

    def test_validation_levels(self):
        """Test all validation levels work."""
        # Test strict, trust, minimal modes
        pass
```

### Validator Checklist

When creating a new validator, ensure:

- [ ] Settings inherit from `BaseSettings`
- [ ] Settings validate in `__post_init__()`
- [ ] All three validation levels implemented
- [ ] Format enum created with `_missing_()` method
- [ ] Config dataclass created
- [ ] Functional API wrapper added to `__init__.py`
- [ ] Comprehensive tests added
- [ ] Documentation updated
- [ ] Examples added to `examples/` directory

---

## Best Practices

### For Validator Developers

1. **Always use immutable settings pattern**
   ```python
   settings = settings.update(param=value)  # ✅
   settings.update(param=value)             # ❌
   ```

2. **Validate settings in `__post_init__()`**
   - Check value ranges
   - Validate enum values
   - Provide clear error messages

3. **Implement all three validation levels**
   - `strict`: Full validation (default)
   - `trust`: Fast validation for known-good data
   - `minimal`: File copy only

4. **Use structured logging**
   ```python
   self.logger.info("Starting validation...")
   self.logger.add_validation_issue(
       level='ERROR',
       category='vcf',
       message='Invalid variant',
       details={'line': 42, 'reason': 'Missing ALT'}
   )
   ```

5. **Handle compression transparently**
   - Use `file_handler` utilities
   - Support all compression types
   - Use parallel compression tools when available

6. **Write comprehensive tests**
   - Test all validation levels
   - Test settings immutability
   - Test error cases
   - Test compression handling

### For Package Users

1. **Use the functional API for simplicity**
   ```python
   from validation_pkg import validate_genome
   output = validate_genome(config.ref_genome, output_dir)
   ```

2. **Use validators directly for customization**
   ```python
   settings = GenomeValidator.Settings()
   settings = settings.update(validation_level='trust')
   validator = GenomeValidator(config, output_dir, settings)
   ```

3. **Always assign update() results**
   ```python
   settings = settings.update(param=value)  # Don't forget assignment!
   ```

4. **Catch ValidationError for all package errors**
   ```python
   try:
       output = validate_genome(config, output_dir)
   except ValidationError as e:
       print(f"Validation failed: {e}")
   ```

---

## See Also

- [CONFIG_GUIDE.md](CONFIG_GUIDE.md) - Configuration file reference
- [PERFORMANCE_GUIDE.md](PERFORMANCE_GUIDE.md) - Performance optimization
- [LOGGING_GUIDE.md](LOGGING_GUIDE.md) - Logging system documentation
- [README.md](../README.md) - Main package documentation
- [tests/README.md](../tests/README.md) - Testing guide
