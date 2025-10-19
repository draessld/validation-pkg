# Validation Package - Usage Guide

This guide demonstrates different ways to use the validation package, from simple to advanced.

## Table of Contents

1. [Quick Start](#quick-start)
2. [Usage Patterns](#usage-patterns)
3. [Example Scripts](#example-scripts)
4. [API Reference](#api-reference)

---

## Quick Start

### Installation

```bash
# Install the package (if distributed)
pip install validation_pkg

# Or use it locally
cd /path/to/validation_pkg
python3 -m validation_pkg validate config.json
```

### Simplest Usage

```python
from validation_pkg import ValidationCoordinator

# Validate all files in config
coordinator = ValidationCoordinator("config.json")
report = coordinator.validate_all()

print(report.summary())
```

---

## Usage Patterns

The package supports **5 different usage patterns** depending on your needs:

### Pattern 1: Simple Workflow (Recommended for most cases)

**Use when:** You want to validate all files with default settings.

```python
from validation_pkg import ValidationCoordinator

coordinator = ValidationCoordinator("config.json")
report = coordinator.validate_all()

if report.passed:
    print("✓ All validations passed!")
else:
    print(f"✗ Found {report.errors} error(s)")
```

### Pattern 2: Selective Validation

**Use when:** You want to validate only specific file types.

```python
from validation_pkg import ValidationCoordinator

coordinator = ValidationCoordinator("config.json")

# Validate only genomes
genome_report = coordinator.validate_genomes()

# Validate only reads
reads_report = coordinator.validate_reads()

# Validate only features
features_report = coordinator.validate_features()
```

### Pattern 3: Custom Settings with Functional API

**Use when:** You need different settings for different files.

**This is the pattern you requested!**

```python
import sys
from validation_pkg import (
    ConfigManager,
    GenomeValidator,
    ReadValidator,
    FeatureValidator,
    validate_genome,
    validate_reads,
    validate_feature
)

# 1. Read and validate config
config_path = sys.argv[1]
config = ConfigManager.load(config_path)

# 2. Edit settings
ref_genome_settings = GenomeValidator.Settings()
ref_genome_settings = ref_genome_settings.update(
    plasmid_split=True,
    coding_type='gz',
    output_filename_suffix='ref',
    replace_id_with='chr',
    min_sequence_length=500
)

mod_genome_settings = GenomeValidator.Settings()
mod_genome_settings = mod_genome_settings.update(
    plasmid_split=True,
    coding_type='gz',
    output_filename_suffix='mod'
)

reads_settings = ReadValidator.Settings()
reads_settings = reads_settings.update(
    coding_type='gz',
    check_invalid_chars=True
)

feature_settings = FeatureValidator.Settings()
feature_settings = feature_settings.update(
    coding_type='gz',
    sort_by_position=True
)

# 3. Run validation
if config.ref_genome:
    stats = validate_genome(config.ref_genome, config.output_dir, ref_genome_settings)
    print(f"✓ Validated {stats['total_sequences']} sequences")

if config.mod_genome:
    stats = validate_genome(config.mod_genome, config.output_dir, mod_genome_settings)

if config.reads:
    stats_list = validate_reads(config.reads, config.output_dir, reads_settings)

if config.ref_feature:
    stats = validate_feature(config.ref_feature, config.output_dir, feature_settings)
```

### Pattern 4: Direct Validator Usage (Maximum Control)

**Use when:** You need complete control over the validation workflow.

```python
from pathlib import Path
from validation_pkg import ConfigManager, GenomeValidator

config = ConfigManager.load("config.json")
output_dir = Path(config.output_dir)

# Create custom settings
settings = GenomeValidator.Settings()
settings = settings.update(
    plasmid_split=True,
    coding_type='gz',
    output_filename_suffix='validated'
)

# Instantiate and run validator
validator = GenomeValidator(config.ref_genome, output_dir, settings)
validator.validate()

# Get detailed statistics
stats = validator.get_statistics()
print(f"Total sequences: {stats['total_sequences']}")
print(f"Total length: {stats['total_length']} bp")
print(f"Longest sequence: {stats['longest_sequence']} bp")
```

### Pattern 5: Settings from Dictionary

**Use when:** You're loading settings from a config file or database.

```python
from validation_pkg import ConfigManager, GenomeValidator, validate_genome

config = ConfigManager.load("config.json")

# Define settings as dictionary
settings_dict = {
    'plasmid_split': True,
    'coding_type': 'gz',
    'output_filename_suffix': 'validated',
    'replace_id_with': 'seq',
    'min_sequence_length': 100
}

# Create settings from dictionary
settings = GenomeValidator.Settings.from_dict(settings_dict)

# Validate
stats = validate_genome(config.ref_genome, config.output_dir, settings)
```
