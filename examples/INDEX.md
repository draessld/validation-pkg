# Examples Index

Quick navigation guide for the validation package examples.

## Getting Started

**New users start here:**
1. Read [QUICK_REFERENCE.md](QUICK_REFERENCE.md) - Essential syntax and patterns
2. Run [01_simple_validation.py](01_simple_validation.py) - Your first validation
3. Read [README.md](README.md) - Comprehensive guide

## Examples by Use Case

### I want to validate files quickly
→ [01_simple_validation.py](01_simple_validation.py) - Simplest approach

### I only need to validate genomes/reads/features
→ [02_selective_validation.py](02_selective_validation.py) - Selective validation

### I need to process genome files with custom requirements
→ [03_custom_genome_settings.py](03_custom_genome_settings.py) - Genome options

### I need to process sequencing reads with quality checks
→ [04_custom_read_settings.py](04_custom_read_settings.py) - Read options

### I need to sort and validate annotation files
→ [05_custom_feature_settings.py](05_custom_feature_settings.py) - Feature options

### I want to understand performance trade-offs
→ [06_validation_levels.py](06_validation_levels.py) - Strict vs Trust vs Minimal

### I need to integrate validation into my code
→ [07_direct_instantiation.py](07_direct_instantiation.py) - Programmatic usage

### I need robust error handling for production
→ [08_error_handling.py](08_error_handling.py) - Exception handling

### I'm confused about how settings work
→ [09_settings_patterns.py](09_settings_patterns.py) - Settings deep dive

### I need a complete pipeline template
→ [10_complete_pipeline.py](10_complete_pipeline.py) - Production template

### I have many read files in directories that need validation
→ [11_directory_based_reads.py](11_directory_based_reads.py) - Directory-based processing

### I want to configure validation in config.json instead of code
→ [12_config_validation_levels.py](12_config_validation_levels.py) - Config-level settings

### I need to process many files quickly with parallel processing
→ [13_parallel_processing.py](13_parallel_processing.py) - Parallel processing

## Examples by Experience Level

### Beginner
- [01_simple_validation.py](01_simple_validation.py)
- [02_selective_validation.py](02_selective_validation.py)
- [QUICK_REFERENCE.md](QUICK_REFERENCE.md)

### Intermediate
- [03_custom_genome_settings.py](03_custom_genome_settings.py)
- [04_custom_read_settings.py](04_custom_read_settings.py)
- [05_custom_feature_settings.py](05_custom_feature_settings.py)
- [06_validation_levels.py](06_validation_levels.py)
- [12_config_validation_levels.py](12_config_validation_levels.py)

### Advanced
- [07_direct_instantiation.py](07_direct_instantiation.py)
- [08_error_handling.py](08_error_handling.py)
- [09_settings_patterns.py](09_settings_patterns.py)
- [10_complete_pipeline.py](10_complete_pipeline.py)
- [11_directory_based_reads.py](11_directory_based_reads.py)
- [13_parallel_processing.py](13_parallel_processing.py)

## Examples by File Type

### Genome Files (FASTA, GenBank)
- [03_custom_genome_settings.py](03_custom_genome_settings.py)

### Read Files (FASTQ, BAM)
- [04_custom_read_settings.py](04_custom_read_settings.py)

### Feature Files (GFF, BED, GTF)
- [05_custom_feature_settings.py](05_custom_feature_settings.py)

### All File Types
- [01_simple_validation.py](01_simple_validation.py)
- [02_selective_validation.py](02_selective_validation.py)
- [10_complete_pipeline.py](10_complete_pipeline.py)

## Examples by Topic

### Configuration
- [01_simple_validation.py](01_simple_validation.py) - Config file usage
- [07_direct_instantiation.py](07_direct_instantiation.py) - Manual config
- [12_config_validation_levels.py](12_config_validation_levels.py) - Config-level settings

### Validation Modes
- [06_validation_levels.py](06_validation_levels.py) - Strict/Trust/Minimal

### Settings Management
- [09_settings_patterns.py](09_settings_patterns.py) - Immutable pattern
- [03_custom_genome_settings.py](03_custom_genome_settings.py) - Custom settings
- [12_config_validation_levels.py](12_config_validation_levels.py) - Config-based settings

### Error Handling
- [08_error_handling.py](08_error_handling.py) - Exceptions and reports

### Performance & Optimization
- [06_validation_levels.py](06_validation_levels.py) - Strict/Trust/Minimal
- [13_parallel_processing.py](13_parallel_processing.py) - Parallel processing

### Production Pipelines
- [10_complete_pipeline.py](10_complete_pipeline.py) - Complete workflow
- [11_directory_based_reads.py](11_directory_based_reads.py) - Batch processing
- [13_parallel_processing.py](13_parallel_processing.py) - Parallel processing

## Test Data

All examples use test data from the [data/](data/) directory:

**Genome files:**
- `genome.fasta` (uncompressed FASTA)
- `genome.fa.gz` (gzip FASTA)
- `genome.gb.bz2` (bzip2 GenBank)
- `genome.gbk.bz2` (bzip2 GenBank)

**Read files:**
- `read.fastq` (uncompressed Illumina)
- `read.fastq.gz` (gzip Illumina)
- `read.fq` (uncompressed ONT)
- `read.fq.gz` (gzip ONT)

**Feature files:**
- `feature.gff` (uncompressed GFF3)
- `feature.gff3.gz` (gzip GFF3)
- `feature.gtf.bz2` (bzip2 GTF)
- `feature.bed` (uncompressed BED)
- `feature.bed.gz` (gzip BED)

**Configuration:**
- `config.json` (sample configuration)

## Running Examples

All examples can be run directly:

```bash
cd examples
python 01_simple_validation.py
```

Or make executable:

```bash
chmod +x *.py
./01_simple_validation.py
```

Each example creates output in `examples/output/<example_name>/`

## Documentation

- **[README.md](README.md)** - Comprehensive guide with detailed explanations
- **[QUICK_REFERENCE.md](QUICK_REFERENCE.md)** - Syntax cheat sheet
- **[INDEX.md](INDEX.md)** - This file - navigation guide

## Common Workflows

### First Time Setup
1. Read [QUICK_REFERENCE.md](QUICK_REFERENCE.md)
2. Create a config.json file
3. Run [01_simple_validation.py](01_simple_validation.py)

### Quality Control Workflow
1. Use [03_custom_genome_settings.py](03_custom_genome_settings.py) with `validation_level='strict'`
2. Enable all quality checks
3. Review validation reports

### Production Pipeline
1. Use [10_complete_pipeline.py](10_complete_pipeline.py) as template
2. Set `validation_level='trust'` for speed
3. Add error handling from [08_error_handling.py](08_error_handling.py)

### Format Conversion
1. Use [04_custom_read_settings.py](04_custom_read_settings.py)
2. Set `validation_level='trust'`
3. Set `coding_type='gz'` or `'bz2'`

## Learning Path

**Day 1: Basics**
1. [QUICK_REFERENCE.md](QUICK_REFERENCE.md)
2. [01_simple_validation.py](01_simple_validation.py)
3. [02_selective_validation.py](02_selective_validation.py)

**Day 2: Customization**
4. [03_custom_genome_settings.py](03_custom_genome_settings.py)
5. [04_custom_read_settings.py](04_custom_read_settings.py)
6. [05_custom_feature_settings.py](05_custom_feature_settings.py)

**Day 3: Advanced Features**
7. [06_validation_levels.py](06_validation_levels.py)
8. [09_settings_patterns.py](09_settings_patterns.py)

**Day 4: Production**
9. [08_error_handling.py](08_error_handling.py)
10. [10_complete_pipeline.py](10_complete_pipeline.py)

**Day 5: Optimization (Optional)**
11. [12_config_validation_levels.py](12_config_validation_levels.py)
12. [13_parallel_processing.py](13_parallel_processing.py)

## Quick Answers

**Q: How do I validate all my files?**
A: Use [01_simple_validation.py](01_simple_validation.py)

**Q: How do I make validation faster?**
A: Use [06_validation_levels.py](06_validation_levels.py), set `validation_level='trust'`. For multiple files, use [13_parallel_processing.py](13_parallel_processing.py) with `threads=8`.

**Q: How do I handle errors?**
A: Use [08_error_handling.py](08_error_handling.py)

**Q: Settings aren't working?**
A: Read [09_settings_patterns.py](09_settings_patterns.py) - settings are immutable!

**Q: Need a production template?**
A: Use [10_complete_pipeline.py](10_complete_pipeline.py)

**Q: How do I configure validation in config.json?**
A: Use [12_config_validation_levels.py](12_config_validation_levels.py) - specify settings in config file

**Q: How do I process many files faster?**
A: Use [13_parallel_processing.py](13_parallel_processing.py) - parallel processing with `threads=8`

## Support

- Check examples for similar use cases
- Read [README.md](README.md) for detailed info
- See main package documentation
- Open GitHub issue for bugs
