# Configuration File Guide

This document defines the **JSON** configuration file that `ConfigManager` uses to locate and process your bioinformatics files (genomes, reads, features).

## Table of Contents

- [Quick Overview](#quick-overview)
- [Configuration Structure](#configuration-structure)
- [File Format Support](#file-format-support)
- [Compression Support](#compression-support)
- [Configuration Objects](#configuration-objects)
- [Validator Settings in Config](#validator-settings-in-config)
- [Examples](#examples)
- [Validation Rules](#validation-rules)

---

## Quick Overview

- **Filename:** `config.json` (standard)
- **Location:** Same directory as your input files (paths are relative to config location)
- **Format:** JSON (UTF-8)
- **Purpose:** Specify reference/modified genomes, reads, and optional plasmids & features
- **Created by:** Users or analysis pipelines
- **Processed by:** `ConfigManager.load("config.json")`

---

## Configuration Structure

### Top-Level Keys

| Key | Type | Required | Description |
|---|---|:---:|---|
| `ref_genome_filename` | GenomeConfig | ✅ | Reference genome (FASTA or GenBank) |
| `mod_genome_filename` | GenomeConfig | ❌ | Modified genome (FASTA or GenBank) |
| `ref_plasmid_filename` | GenomeConfig | ❌ | Reference plasmid (FASTA or GenBank) |
| `mod_plasmid_filename` | GenomeConfig | ❌ | Modified plasmid (FASTA or GenBank) |
| `reads` | List[ReadConfig] | ❌ | Read files (FASTQ/BAM), minimum one if provided |
| `ref_feature_filename` | FeatureConfig | ❌ | Features for reference genome (BED, GFF, GTF) |
| `mod_feature_filename` | FeatureConfig | ❌ | Features for modified genome (BED, GFF, GTF) |
| `options` | dict | ❌ | Additional options (e.g., `{"threads": 8}`) |

**Important:**
- All file paths must be **relative to the config file directory**
- At minimum, provide `ref_genome_filename` or include `reads`
- Genome validation requires at least one genome file

---

## File Format Support

### Genome and Plasmid Files

| Format | Extensions | Notes |
|--------|-----------|-------|
| FASTA | `.fa`, `.fasta`, `.fna` | Most common genome format |
| GenBank | `.gb`, `.gbk`, `.genbank` | Includes annotations |

### Feature Files

| Format | Extensions | Notes |
|--------|-----------|-------|
| GFF | `.gff`, `.gff3` | General Feature Format |
| GTF | `.gtf` | Gene Transfer Format |
| BED | `.bed` | Browser Extensible Data |

### Read Files

| Format | Extensions | Notes |
|--------|-----------|-------|
| FASTQ | `.fastq`, `.fq` | Quality scores included |
| BAM | `.bam` | Binary alignment format (converted to FASTQ) |

---

## Compression Support

All file types support transparent compression. The package automatically detects and handles compressed files.

### Supported Compression Formats

| Type | Extensions | Detection | Performance |
|------|-----------|-----------|-------------|
| **gzip** | `.gz`, `.gzip` | Automatic | Fast with `pigz` (parallel) |
| **bzip2** | `.bz2`, `.bzip2` | Automatic | Fast with `pbzip2` (parallel) |

**Performance tip:** Install `pigz` and `pbzip2` for 3-4x faster compression/decompression:
```bash
# Ubuntu/Debian
sudo apt-get install pigz pbzip2

# macOS
brew install pigz pbzip2
```

### Compression Examples

```json
{
  "ref_genome_filename": {"filename": "genome.fasta.gz"},
  "reads": [
    {"filename": "reads.fastq.bz2", "ngs_type": "illumina"}
  ]
}

---

## Configuration Objects

### GenomeConfig Object

Specifies a genome or plasmid file.

**Structure:**
```json
{
  "filename": "path/to/genome.fasta"
}
```

**Optional fields (validator settings):**
```json
{
  "filename": "path/to/genome.fasta.gz",
  "validation_level": "trust",
  "min_sequence_length": 500,
  "plasmid_split": true
}
```

**Example:**
```json
"ref_genome_filename": {"filename": "data/reference.gbk"}
```

### ReadConfig Object

Specifies sequencing read files. Supports individual files or directories.

#### Option 1: Individual Files

```json
"reads": [
  {
    "filename": "samples/reads_R1.fastq.gz",
    "ngs_type": "illumina"
  },
  {
    "filename": "samples/reads_R2.fastq.gz",
    "ngs_type": "illumina"
  }
]
```

#### Option 2: Directory of Files

All files in the directory inherit the same settings:

```json
"reads": [
  {
    "directory": "ont_reads/",
    "ngs_type": "ont",
    "validation_level": "trust"
  }
]
```

#### NGS Type Values

| Value | Description | Typical Read Length |
|-------|-------------|-------------------|
| `"illumina"` | Illumina short reads (SE or PE) | 50-300 bp |
| `"ont"` | Oxford Nanopore long reads | 1-100+ kb |
| `"pacbio"` | PacBio long reads | 10-50+ kb |

**Note:** BAM files are automatically detected and converted to FASTQ. Default NGS type for BAM: `pacbio`.

### FeatureConfig Object

Specifies feature annotation files (BED, GFF, GTF).

**Structure:**
```json
{
  "filename": "path/to/features.gff"
}
```

**Optional fields (validator settings):**
```json
{
  "filename": "path/to/features.gff",
  "validation_level": "strict",
  "sort_by_position": true,
  "check_coordinates": true
}
```

**Example:**
```json
"ref_feature_filename": {"filename": "annotations/reference.gff"}
```

---

## Validator Settings in Config

**NEW:** You can specify validator-specific settings directly in your config file. These settings customize validation behavior without modifying code.

### Supported Settings

#### Genome Validator Settings
- `validation_level`: `"strict"`, `"trust"`, or `"minimal"`
- `min_sequence_length`: Minimum sequence length (int)
- `plasmid_split`: Split plasmids into separate files (bool)
- `main_longest`: Select longest sequence as main (bool)
- `threads`: Number of threads for compression (int)

#### Read Validator Settings
- `validation_level`: `"strict"`, `"trust"`, or `"minimal"`
- `check_invalid_chars`: Check for invalid characters (bool)
- `allow_duplicate_ids`: Allow duplicate sequence IDs (bool)
- `ignore_bam`: Skip BAM files (bool)
- `threads`: Number of threads for compression (int)

#### Feature Validator Settings
- `validation_level`: `"strict"`, `"trust"`, or `"minimal"`
- `sort_by_position`: Sort features by position (bool)
- `check_coordinates`: Validate feature coordinates (bool)
- `threads`: Number of threads for compression (int)

### Example with Settings

```json
{
  "ref_genome_filename": {
    "filename": "genome.fasta",
    "validation_level": "trust",
    "min_sequence_length": 500
  },
  "reads": [
    {
      "filename": "reads.fastq.gz",
      "ngs_type": "illumina",
      "validation_level": "minimal",
      "check_invalid_chars": false
    }
  ],
  "ref_feature_filename": {
    "filename": "features.gff",
    "sort_by_position": true
  },
  "options": {
    "threads": 8
  }
}
```

**Settings precedence:** Settings in your Python code override settings in config.json.

---

## Examples

### Minimal Configuration

Simplest valid configuration with required fields only:

```json
{
  "ref_genome_filename": {"filename": "ref.fasta"},
  "reads": [
    {"filename": "reads.fastq", "ngs_type": "illumina"}
  ]
}
```

### Full Configuration

Complete example with all optional fields and settings:

```json
{
  "ref_genome_filename": {
    "filename": "ref.gbk",
    "validation_level": "strict",
    "plasmid_split": true
  },
  "mod_genome_filename": {
    "filename": "mod.fasta.gz",
    "validation_level": "trust"
  },
  "ref_plasmid_filename": {"filename": "plasmid_ref.gbk"},
  "mod_plasmid_filename": {"filename": "plasmid_mod.fasta"},
  "reads": [
    {
      "filename": "illumina_R1.fastq.gz",
      "ngs_type": "illumina",
      "validation_level": "trust"
    },
    {
      "filename": "illumina_R2.fastq.gz",
      "ngs_type": "illumina",
      "validation_level": "trust"
    },
    {
      "directory": "ont_reads/",
      "ngs_type": "ont",
      "validation_level": "minimal"
    }
  ],
  "ref_feature_filename": {
    "filename": "features_ref.gff3",
    "sort_by_position": true
  },
  "mod_feature_filename": {"filename": "features_mod.bed"},
  "options": {
    "threads": 8
  }
}
```

### Directory-Based Reads

Process all files in a directory with the same settings:

```json
{
  "ref_genome_filename": {"filename": "genome.fasta"},
  "reads": [
    {
      "directory": "illumina_reads/",
      "ngs_type": "illumina",
      "validation_level": "trust"
    },
    {
      "directory": "ont_reads/",
      "ngs_type": "ont",
      "validation_level": "minimal"
    }
  ]
}
```

---

## Validation Rules

The package validates your configuration before processing:

1. **Required fields:** At least `ref_genome_filename` OR `reads` must be present
2. **File existence:** All specified files must exist and be accessible
3. **Path resolution:** All paths are resolved relative to config directory
4. **Format detection:** File extensions must match supported formats
5. **NGS type:** Read files must specify a valid `ngs_type`
6. **Settings validation:** Validator settings must use valid field names and values

**Common validation errors:**
- Missing required fields
- File not found
- Invalid NGS type
- Invalid validator setting name
- Unsupported file format
