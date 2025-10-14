# Configuration File — Specification & Guide

This document defines the **JSON** configuration that `coordinator.py` consumes to understand your inputs (what’s reference vs. modified; genome vs. features vs. reads) and to validate them robustly.

---

## Quick overview

- **Format:** JSON (UTF-8)
- **Purpose:** Tell where to find the reference/modified genomes, reads, and optional plasmids & features.
- **Who creates it:** User, Modification Applier
- **Who reads it:** `Coordinator.load("config.json")` (and maybe EFSA)

---

## Top-level keys

| Key | Type | Required | Description |
|---|---|:---:|---|
| `ref_genome_filename` | GenomeConfig | ✅ | Path to **reference genome** (FASTA or GenBank). |
| `mod_genome_filename` | GenomeConfig | ✅ | Path to **modified genome** (FASTA or GenBank). |
| `ref_plasmid_filename` | GenomeConfig | ❌ | Path to **reference plasmid** (FASTA or GenBank). |
| `mod_plasmid_filename` | GenomeConfig | ❌ | Path to **modified plasmid** (FASTA or GenBank). |
| `reads` | ReadConfig | ✅ | **At least one** read input. See “Reads structure” below. |
| `ref_feature_filename` | FeatureConfig | ❌ | Features for the **reference genome** (BED, GFF or GTF). |
| `mod_feature_filename` | FeatureConfig | ❌ | Features for the **modified genome** (BED, GFF or GTF). |
| `options` | object | ❌ | Free-form extra knobs (e.g., `{"threads": 8}`) - not implemented yet. |

- **GenomeConfig**


> ℹ️ **Paths** must be relative to the **config directory** where your configuration file is.

---

## Supported file formats

- **Genomes & plasmids:** `.fa`, `.fasta`, `.fna` → FASTA; `.gb`, `.gbk`, `.genbank` → GenBank  
- **Features:** `.gff`, `.gff3` → GFF; `.gtf` → GTF; `.bed` → BED
- **Reads:** `.fastq`, `.fq`, `.bam` → FASTQ

---

## Output compression (coding)

All file-type or directory-type entries may include an optional `output_coding` field to control how the results will be stored.  
If omitted, `"none"` (uncompressed) is assumed.

Allowed values:

- `"gz"` — gzip compression (`.gz`)
- `"bzip"` — bzip2 compression (`.bz2`)
- `"tgz"` — tar + gzip archive (`.tgz`)
- `"none"` — uncompressed (default when the `output_coding` is not included)

---
## Genome Config
Each entry is **an object**: include `filename`. May also include `output_coding` field.


Example:
```json
"ref_genome_filename": {"filename":"data/reference.gbk", "output_coding": "gz"}
```

---
## Read Config

#### a) file-type non-empty list
Each entry is **an object**: include `filename` and `ngs_type`.
```json
"reads": [
  { "filename": "samples/reads_R1.fastq.gz", "ngs_type": "illumina" },
  { "filename": "samples/reads_R2.fastq.gz", "ngs_type": "illumina" }
]
```

####  b) directory-type path

`reads` must be a **non-empty list** or **directory path**. Each entry is **an object**: include `filename` and `ngs_type`.
```json
"reads": [
  { "directory": "samples/", "ngs_type": "illumina" },
]
```

####  Every record in Reads config may have output coding requirements
```json
"reads": [
  { "filename": "samples/reads_R1.fastq.gz", "ngs_type": "illumina", "output_coding": "gz" },
  { "filename": "samples/reads_R2.fastq.gz", "ngs_type": "illumina" }
]
```

#### `ngs_type` allowed values

- `"illumina"` — short reads (SE or PE)
- `"ont"` — Oxford Nanopore
- `"pacbio"` — PacBio

---
## Feature Config
Each entry is **an object**: include `filename`. May also include `output_coding` field.


Example:
```json
"ref_feature_filename": {"filename":"data/reference.gbk", "output_coding": "gz"}
```

---
## Validation rules

1. `ref_genome_filename` and `mod_genome_filename` must exist.
2. `reads` must be non-empty.
3. All paths must exist or be accessible (but it will be checked).
4. File extensions determine format.
