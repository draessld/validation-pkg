# Settings Reference

This document provides a comprehensive reference for all settings across the validation_pkg validators, including their categorization, usage patterns, and future refactoring plans.

## Table of Contents
- [Overview](#overview)
- [BaseSettings Pattern](#basesettings-pattern)
- [GenomeValidator Settings](#genomevalidator-settings)
- [ReadValidator Settings](#readvalidator-settings)
- [FeatureValidator Settings](#featurevalidator-settings)
- [Common Patterns](#common-patterns)
- [Settings Categories](#settings-categories)
- [Future Refactoring Plans](#future-refactoring-plans)

---

## Overview

All validators in validation_pkg use a consistent Settings pattern based on `BaseSettings`. Settings are:

- **Immutable**: Use `settings.update(param=value)` to create a new instance
- **Type-safe**: Dataclass with type hints
- **Validated**: Field validation in `__post_init__()`
- **Serializable**: `to_dict()` and `from_dict()` methods

### Quick Example

```python
from validation_pkg import GenomeValidator

# Create settings
settings = GenomeValidator.Settings()

# Update settings (returns NEW instance)
settings = settings.update(
    plasmid_split=True,
    coding_type='gz',
    output_filename_suffix='ref'
)

# WRONG - changes are lost!
settings.update(plasmid_split=True)  # Returns new instance, not assigned

# Serialize/deserialize
settings_dict = settings.to_dict()
settings_copy = GenomeValidator.Settings.from_dict(settings_dict)
```

---

## BaseSettings Pattern

**Location**: `validation_pkg/utils/settings.py`

All validator Settings classes inherit from `BaseSettings`, which provides:

### Methods

| Method | Description | Returns |
|--------|-------------|---------|
| `copy()` | Create deep copy of settings | New Settings instance |
| `update(**kwargs)` | Create new instance with updated fields | New Settings instance |
| `to_dict(include_unset=True)` | Serialize to dictionary | Dict[str, Any] |
| `from_dict(data)` | Deserialize from dictionary | Settings instance |
| `__str__()` | Pretty-print settings | str |
| `__repr__()` | Debug representation | str |

### Key Features

1. **Immutable Update Pattern**
   ```python
   # CORRECT
   settings = settings.update(param=value)

   # WRONG
   settings.update(param=value)  # Changes lost!
   ```

2. **Field Name Validation**
   ```python
   settings.update(invalid_field=True)  # Raises ValueError
   ```

3. **Deep Copy Support**
   ```python
   settings_copy = settings.copy()
   ```

4. **Dictionary Serialization**
   ```python
   data = settings.to_dict()
   settings = Settings.from_dict(data)
   ```

---

## GenomeValidator Settings

**Location**: `validation_pkg/validators/genome_validator.py`

### Complete Field Reference

| Field | Type | Default | Category | Description |
|-------|------|---------|----------|-------------|
| `allow_empty_sequences` | `bool` | `False` | Validation | Allow SeqRecord with empty sequence |
| `allow_empty_id` | `bool` | `False` | Validation | Allow SeqRecord with empty ID |
| `warn_n_sequences` | `int` | `2` | Validation | Warn if sequence count exceeds this |
| `is_plasmid` | `bool` | `False` | Edit | Treat all sequences as plasmids |
| `plasmid_split` | `bool` | `False` | Edit | Split plasmids into separate files |
| `plasmids_to_one` | `bool` | `False` | Edit | Merge all plasmids into one file |
| `main_longest` | `bool` | `True` | Edit | Select longest sequence as main chromosome |
| `main_first` | `bool` | `False` | Edit | Select first sequence as main chromosome |
| `replace_id_with` | `Optional[str]` | `None` | Edit | Prefix to add to sequence IDs (e.g., "chr") |
| `min_sequence_length` | `int` | `100` | Edit | Minimum sequence length to keep |
| `coding_type` | `Optional[str]` | `None` | Output | Output compression ('gz', 'bz2', None) |
| `output_filename_suffix` | `Optional[str]` | `None` | Output | Suffix for output filename (e.g., 'ref', 'mod') |
| `output_subdir_name` | `Optional[str]` | `None` | Output | Subdirectory for output files |

### Field Dependencies

**Mutually Exclusive Pairs** (validated in `__post_init__`):

1. **`plasmid_split` ⚠️ `plasmids_to_one`**
   - Cannot use both simultaneously
   - Raises `ValueError` if both are `True`

2. **`main_longest` ⚠️ `main_first`**
   - Cannot use both simultaneously
   - Raises `ValueError` if both are `True`

### Usage Examples

#### Basic Genome Validation
```python
settings = GenomeValidator.Settings()
settings = settings.update(
    min_sequence_length=500,
    coding_type='gz'
)
```

#### Plasmid Splitting
```python
settings = GenomeValidator.Settings()
settings = settings.update(
    plasmid_split=True,
    main_longest=True,
    min_sequence_length=1000
)
```

#### ID Replacement
```python
settings = GenomeValidator.Settings()
settings = settings.update(
    replace_id_with='chr',  # sequence1 → chr1
    output_filename_suffix='ref'
)
```

---

## ReadValidator Settings

**Location**: `validation_pkg/validators/read_validator.py`

### Complete Field Reference

| Field | Type | Default | Category | Description |
|-------|------|---------|----------|-------------|
| `check_invalid_chars` | `bool` | `False` | Validation | Validate sequences contain only ATCGN |
| `allow_empty_id` | `bool` | `False` | Validation | Allow sequences without IDs |
| `allow_duplicate_ids` | `bool` | `True` | Validation | Allow duplicate sequence IDs |
| `keep_bam` | `bool` | `True` | Edit | Keep original BAM when converting to FASTQ |
| `ignore_bam` | `bool` | `True` | Edit | Skip BAM files entirely |
| `coding_type` | `Optional[str]` | `None` | Output | Output compression ('gz', 'bz2', None) |
| `output_filename_suffix` | `Optional[str]` | `None` | Output | Suffix for output filename |
| `output_subdir_name` | `Optional[str]` | `None` | Output | Subdirectory for output files |
| `outdir_by_ngs_type` | `bool` | `False` | Output | Auto-set subdirectory to NGS type |

### Field Dependencies

**`outdir_by_ngs_type` overrides `output_subdir_name`**:
- If `outdir_by_ngs_type=True`, subdirectory is set to NGS type (illumina/ont/pacbio)
- Manual `output_subdir_name` is ignored with a warning

### Usage Examples

#### Strict Read Validation
```python
settings = ReadValidator.Settings()
settings = settings.update(
    check_invalid_chars=True,
    allow_duplicate_ids=False,
    coding_type='gz'
)
```

#### Organize by NGS Type
```python
settings = ReadValidator.Settings()
settings = settings.update(
    outdir_by_ngs_type=True,  # Creates illumina/, ont/, pacbio/ subdirs
    coding_type='gz'
)
```

#### BAM Handling
```python
# Skip BAM files
settings = ReadValidator.Settings()
settings = settings.update(ignore_bam=True)

# Convert BAM to FASTQ, keep original
settings = ReadValidator.Settings()
settings = settings.update(
    ignore_bam=False,
    keep_bam=True
)
```

---

## FeatureValidator Settings

**Location**: `validation_pkg/validators/feature_validator.py`

### Complete Field Reference

| Field | Type | Default | Category | Description |
|-------|------|---------|----------|-------------|
| `sort_by_position` | `bool` | `True` | Edit | Sort features by genomic position |
| `check_coordinates` | `bool` | `True` | Validation | Validate that start < end |
| `replace_id_with` | `Optional[str]` | `None` | Edit | Replace seqname/chromosome field |
| `coding_type` | `Optional[str]` | `None` | Output | Output compression ('gz', 'bz2', None) |
| `output_filename_suffix` | `Optional[str]` | `None` | Output | Suffix for output filename |
| `output_subdir_name` | `Optional[str]` | `None` | Output | Subdirectory for output files |

### Field Dependencies

**No inter-dependencies** in FeatureValidator.Settings.

### Usage Examples

#### Basic Feature Validation
```python
settings = FeatureValidator.Settings()
settings = settings.update(
    sort_by_position=True,
    check_coordinates=True
)
```

#### Sequence ID Replacement
```python
settings = FeatureValidator.Settings()
settings = settings.update(
    replace_id_with='chr',  # sequence1 → chr1 in seqname column
    output_filename_suffix='ref'
)
```

#### Skip Coordinate Checking (Trusted Data)
```python
settings = FeatureValidator.Settings()
settings = settings.update(
    check_coordinates=False,  # Trust that start < end
    sort_by_position=True
)
```

---

## Common Patterns

### Settings Shared Across ALL Validators

#### OutputSettings (3 fields)

These fields appear in **every** validator:

| Field | Type | Default | Description |
|-------|------|---------|-------------|
| `coding_type` | `Optional[str]` | `None` | Output compression format |
| `output_filename_suffix` | `Optional[str]` | `None` | Filename suffix |
| `output_subdir_name` | `Optional[str]` | `None` | Output subdirectory |

**coding_type values**:
- `None` - No compression (raw format)
- `'gz'` - gzip compression (.gz)
- `'bz2'` - bzip2 compression (.bz2)

**Example**:
```python
settings = settings.update(
    coding_type='gz',
    output_filename_suffix='validated',
    output_subdir_name='genomes'
)
# Output: valid/genomes/genome_validated.fasta.gz
```

#### Shared Validation Settings

**`allow_empty_id`** - Shared by GenomeValidator and ReadValidator:
- Genome: Allow SeqRecord with empty ID field
- Read: Allow FASTQ sequences without IDs

### Validator-Specific Patterns

#### Plasmid Handling (GenomeValidator only)

Bacterial genomes often have chromosome + plasmids:

```python
# Split plasmids into separate files
settings = settings.update(
    plasmid_split=True,
    main_longest=True,       # Longest = main chromosome
    min_sequence_length=1000 # Filter short contigs
)
# Output: genome_chr.fasta, genome_plasmid1.fasta, genome_plasmid2.fasta

# Merge all plasmids into one file
settings = settings.update(
    plasmids_to_one=True,
    main_longest=True
)
# Output: genome_chr.fasta, genome_plasmid.fasta
```

#### Character Validation (ReadValidator only)

Strict nucleotide validation for reads:

```python
settings = settings.update(
    check_invalid_chars=True  # Only allows ATCGN
)
# Raises ReadValidationError if invalid chars found
```

#### Position Sorting (FeatureValidator only)

GFF/GTF features can be sorted by position:

```python
settings = settings.update(
    sort_by_position=True  # Sort by seqname, start, end
)
```

---

## Settings Categories

All settings can be categorized into three groups:

### 1. ValidationSettings

Control **what validation checks** to perform:

**GenomeValidator**:
- `allow_empty_sequences` - Allow empty sequences
- `allow_empty_id` - Allow empty IDs
- `warn_n_sequences` - Warning threshold for sequence count

**ReadValidator**:
- `check_invalid_chars` - Validate nucleotide characters
- `allow_empty_id` - Allow empty IDs
- `allow_duplicate_ids` - Allow duplicate IDs

**FeatureValidator**:
- `check_coordinates` - Validate start < end

### 2. EditSettings

Control **file modifications/transformations**:

**GenomeValidator**:
- `is_plasmid` - Treat as plasmid
- `plasmid_split` - Split plasmids to files
- `plasmids_to_one` - Merge plasmids
- `main_longest` - Select longest as main
- `main_first` - Select first as main
- `replace_id_with` - ID prefix replacement
- `min_sequence_length` - Length filter

**ReadValidator**:
- `keep_bam` - Keep original BAM
- `ignore_bam` - Skip BAM files

**FeatureValidator**:
- `sort_by_position` - Sort by position
- `replace_id_with` - Seqname replacement

### 3. OutputSettings

Control **output file naming/format**:

**All Validators**:
- `coding_type` - Compression format
- `output_filename_suffix` - Filename suffix
- `output_subdir_name` - Subdirectory

**ReadValidator only**:
- `outdir_by_ngs_type` - Auto-subdirectory by NGS type

---

## Future Refactoring Plans

### Current Structure (Flat)

All settings in a single flat dataclass:

```python
@dataclass
class GenomeValidator.Settings(BaseSettings):
    # Mix of validation, edit, and output settings
    allow_empty_sequences: bool = False
    plasmid_split: bool = False
    coding_type: Optional[str] = None
    # ... all mixed together
```

**Issues**:
- Hard to see categorization
- No separation of concerns
- Validation logic spread across `__post_init__`

### Proposed Structure (Nested)

Separate Settings classes for each category:

```python
@dataclass
class GenomeValidator:
    @dataclass
    class ValidationSettings(BaseSettings):
        """Settings for validation checks."""
        allow_empty_sequences: bool = False
        allow_empty_id: bool = False
        warn_n_sequences: int = 2

    @dataclass
    class EditSettings(BaseSettings):
        """Settings for file transformations."""
        is_plasmid: bool = False
        plasmid_split: bool = False
        plasmids_to_one: bool = False
        main_longest: bool = True
        main_first: bool = False
        replace_id_with: Optional[str] = None
        min_sequence_length: int = 100

        def __post_init__(self):
            # Validation specific to edit settings
            if self.plasmid_split and self.plasmids_to_one:
                raise ValueError("...")

    @dataclass
    class OutputSettings(BaseSettings):
        """Settings for output files."""
        coding_type: Optional[str] = None
        output_filename_suffix: Optional[str] = None
        output_subdir_name: Optional[str] = None

    @dataclass
    class Settings(BaseSettings):
        """Combined settings for GenomeValidator."""
        validation: ValidationSettings = field(default_factory=ValidationSettings)
        edit: EditSettings = field(default_factory=EditSettings)
        output: OutputSettings = field(default_factory=OutputSettings)

        @classmethod
        def from_flat(cls, **kwargs):
            """Migrate from old flat settings."""
            # Backward compatibility helper
```

### Usage After Refactoring

```python
# New nested structure
settings = GenomeValidator.Settings()
settings = settings.update(
    validation=GenomeValidator.ValidationSettings(allow_empty_id=True),
    edit=GenomeValidator.EditSettings(plasmid_split=True),
    output=GenomeValidator.OutputSettings(coding_type='gz')
)

# Access: settings.validation.allow_empty_id
# Access: settings.edit.plasmid_split
# Access: settings.output.coding_type

# Backward compatibility
settings = GenomeValidator.Settings.from_flat(
    allow_empty_id=True,
    plasmid_split=True,
    coding_type='gz'
)
```

### Benefits of Nested Structure

✅ **Clear separation of concerns** - Validation vs edit vs output
✅ **Self-documenting** - Code structure reflects design intent
✅ **Type safety** - IDEs understand nested structure
✅ **Isolated validation** - Each Settings class has its own `__post_init__`
✅ **Extensible** - Easy to add new categories or fields
✅ **Better organization** - Scales as settings grow

### Migration Strategy

1. **Phase 1**: Create new nested structure alongside old flat structure
2. **Phase 2**: Add deprecation warnings to old flat access
3. **Phase 3**: Provide automatic migration helper (`from_flat()`)
4. **Phase 4**: Update all examples and documentation
5. **Phase 5**: Remove deprecated flat structure (major version bump)

### Timeline

- **Current**: Flat structure (all validators)
- **Planned**: Nested structure (after inter-file validation is stable)
- **Target**: v2.0.0 release (breaking change)

---

## See Also

- [VALIDATOR_DESIGN.md](VALIDATOR_DESIGN.md) - Overall validator architecture
- [CONFIG_GUIDE.md](CONFIG_GUIDE.md) - Configuration file specification
- [README.md](../README.md) - API reference and usage examples
- [BaseSettings source](../validation_pkg/utils/settings.py) - Implementation details
