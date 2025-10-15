#!/bin/bash
# Script to generate test fixtures for config_manager integration tests
# Creates realistic test cases with config files and genome/read files

set -e  # Exit on error

SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
cd "$SCRIPT_DIR"

echo "Generating config_manager test fixtures..."

# Clean up old config fixtures
rm -rf config_test_*

# ============================================================================
# Config Test Case 1: Valid minimal configuration
# ============================================================================
echo "Creating config_test_01_minimal_valid..."
mkdir -p config_test_01_minimal_valid/data

# Create genome files
cat > config_test_01_minimal_valid/data/ref_genome.fasta << 'EOF'
>chr1
ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG
GCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTA
EOF

cat > config_test_01_minimal_valid/data/mod_genome.fasta << 'EOF'
>chr1
ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG
GCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTA
TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT
EOF

# Create read file
cat > config_test_01_minimal_valid/data/reads.fastq << 'EOF'
@read1
ATCGATCGATCGATCGATCG
+
IIIIIIIIIIIIIIIIIIII
@read2
GCTAGCTAGCTAGCTAGCTA
+
IIIIIIIIIIIIIIIIIIII
EOF

# Create config file
cat > config_test_01_minimal_valid/data/config.json << 'EOF'
{
  "ref_genome_filename": "ref_genome.fasta",
  "mod_genome_filename": "mod_genome.fasta",
  "reads": [
    {
      "filename": "reads.fastq",
      "ngs_type": "illumina"
    }
  ]
}
EOF

cat > config_test_01_minimal_valid/description.txt << 'EOF'
Config Test Case 1: Valid minimal configuration
- All required fields present
- Single read file
- All files exist and are valid
- Expected: Should load successfully
EOF

# ============================================================================
# Config Test Case 2: Full configuration with all optional fields
# ============================================================================
echo "Creating config_test_02_full_config..."
mkdir -p config_test_02_full_config/data

# Create all required and optional files
cat > config_test_02_full_config/data/ref_genome.fasta << 'EOF'
>chr1
ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG
EOF

cat > config_test_02_full_config/data/mod_genome.fasta << 'EOF'
>chr1
ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG
EOF

cat > config_test_02_full_config/data/ref_plasmid.gbk << 'EOF'
LOCUS       PLASMID_REF             50 bp    DNA     circular     01-JAN-2024
DEFINITION  Reference plasmid
ACCESSION   PLAS001
VERSION     PLAS001.1
KEYWORDS    .
SOURCE      Test organism
  ORGANISM  Test organism
            Bacteria.
ORIGIN
        1 atcgatcgat cgatcgatcg atcgatcgat cgatcgatcg atcgatcgat
//
EOF

cat > config_test_02_full_config/data/mod_plasmid.gbk << 'EOF'
LOCUS       PLASMID_MOD             50 bp    DNA     circular     01-JAN-2024
DEFINITION  Modified plasmid
ACCESSION   PLAS002
VERSION     PLAS002.1
KEYWORDS    .
SOURCE      Test organism
  ORGANISM  Test organism
            Bacteria.
ORIGIN
        1 atcgatcgat cgatcgatcg atcgatcgat cgatcgatcg atcgatcgat
//
EOF

cat > config_test_02_full_config/data/reads_R1.fastq << 'EOF'
@read1
ATCGATCGATCGATCGATCG
+
IIIIIIIIIIIIIIIIIIII
EOF

cat > config_test_02_full_config/data/reads_R2.fastq << 'EOF'
@read2
GCTAGCTAGCTAGCTAGCTA
+
IIIIIIIIIIIIIIIIIIII
EOF

cat > config_test_02_full_config/data/ref_features.gff << 'EOF'
##gff-version 3
chr1	.	gene	100	200	.	+	.	ID=gene1;Name=testGene
EOF

cat > config_test_02_full_config/data/mod_features.gff << 'EOF'
##gff-version 3
chr1	.	gene	100	200	.	+	.	ID=gene1;Name=testGene
chr1	.	gene	300	400	.	+	.	ID=gene2;Name=newGene
EOF

# Create full config
cat > config_test_02_full_config/data/config.json << 'EOF'
{
  "ref_genome_filename": "ref_genome.fasta",
  "mod_genome_filename": "mod_genome.fasta",
  "ref_plasmid_filename": "ref_plasmid.gbk",
  "mod_plasmid_filename": "mod_plasmid.gbk",
  "ref_feature_filename": "ref_features.gff",
  "mod_feature_filename": "mod_features.gff",
  "reads": [
    {
      "filename": "reads_R1.fastq",
      "ngs_type": "illumina"
    },
    {
      "filename": "reads_R2.fastq",
      "ngs_type": "illumina"
    }
  ],
  "options": {
    "threads": 8,
    "verbose": true
  }
}
EOF

cat > config_test_02_full_config/description.txt << 'EOF'
Config Test Case 2: Full configuration with all optional fields
- All required and optional fields present
- Multiple read files
- Plasmids and features included
- Options specified
- Expected: Should load successfully with all fields populated
EOF

# ============================================================================
# Config Test Case 3: Directory-based reads
# ============================================================================
echo "Creating config_test_03_directory_reads..."
mkdir -p config_test_03_directory_reads/data/reads_dir

# Create genome files
cat > config_test_03_directory_reads/data/ref_genome.fasta << 'EOF'
>chr1
ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG
EOF

cat > config_test_03_directory_reads/data/mod_genome.fasta << 'EOF'
>chr1
ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG
EOF

# Create multiple read files in directory
for i in {1..3}; do
cat > config_test_03_directory_reads/data/reads_dir/sample${i}.fastq << EOF
@read${i}
ATCGATCGATCGATCGATCG
+
IIIIIIIIIIIIIIIIIIII
EOF
done

# Create config with directory
cat > config_test_03_directory_reads/data/config.json << 'EOF'
{
  "ref_genome_filename": "ref_genome.fasta",
  "mod_genome_filename": "mod_genome.fasta",
  "reads": [
    {
      "directory": "reads_dir",
      "ngs_type": "ont"
    }
  ]
}
EOF

cat > config_test_03_directory_reads/description.txt << 'EOF'
Config Test Case 3: Directory-based reads
- Reads specified as directory instead of individual files
- Multiple read files in directory
- Expected: Should load all files from directory
EOF

# ============================================================================
# Config Test Case 4: Compressed files
# ============================================================================
echo "Creating config_test_04_compressed..."
mkdir -p config_test_04_compressed/data

# Create and compress genome files
cat > config_test_04_compressed/data/ref_genome.fasta << 'EOF'
>chr1
ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG
EOF
gzip -c config_test_04_compressed/data/ref_genome.fasta > config_test_04_compressed/data/ref_genome.fasta.gz
rm config_test_04_compressed/data/ref_genome.fasta

cat > config_test_04_compressed/data/mod_genome.fasta << 'EOF'
>chr1
ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG
EOF
bzip2 -c config_test_04_compressed/data/mod_genome.fasta > config_test_04_compressed/data/mod_genome.fasta.bz2
rm config_test_04_compressed/data/mod_genome.fasta

# Create compressed read file
cat > config_test_04_compressed/data/reads.fastq << 'EOF'
@read1
ATCGATCGATCGATCGATCG
+
IIIIIIIIIIIIIIIIIIII
EOF
gzip -c config_test_04_compressed/data/reads.fastq > config_test_04_compressed/data/reads.fastq.gz
rm config_test_04_compressed/data/reads.fastq

# Create config
cat > config_test_04_compressed/data/config.json << 'EOF'
{
  "ref_genome_filename": "ref_genome.fasta.gz",
  "mod_genome_filename": "mod_genome.fasta.bz2",
  "reads": [
    {
      "filename": "reads.fastq.gz",
      "ngs_type": "illumina"
    }
  ]
}
EOF

cat > config_test_04_compressed/description.txt << 'EOF'
Config Test Case 4: Compressed files
- Gzip compressed reference genome
- Bzip2 compressed modified genome
- Gzip compressed reads
- Expected: Should detect compression types and load successfully
EOF

# ============================================================================
# Config Test Case 5: Multiple NGS types
# ============================================================================
echo "Creating config_test_05_multiple_ngs_types..."
mkdir -p config_test_05_multiple_ngs_types/data

cat > config_test_05_multiple_ngs_types/data/ref_genome.fasta << 'EOF'
>chr1
ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG
EOF

cat > config_test_05_multiple_ngs_types/data/mod_genome.fasta << 'EOF'
>chr1
ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG
EOF

cat > config_test_05_multiple_ngs_types/data/illumina.fastq << 'EOF'
@illumina_read
ATCGATCGATCGATCGATCG
+
IIIIIIIIIIIIIIIIIIII
EOF

cat > config_test_05_multiple_ngs_types/data/ont.fastq << 'EOF'
@ont_read
ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG
+
IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
EOF

cat > config_test_05_multiple_ngs_types/data/pacbio.fastq << 'EOF'
@pacbio_read
ATCGATCGATCGATCGATCGATCGATCGATCG
+
IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
EOF

cat > config_test_05_multiple_ngs_types/data/config.json << 'EOF'
{
  "ref_genome_filename": "ref_genome.fasta",
  "mod_genome_filename": "mod_genome.fasta",
  "reads": [
    {
      "filename": "illumina.fastq",
      "ngs_type": "illumina"
    },
    {
      "filename": "ont.fastq",
      "ngs_type": "ont"
    },
    {
      "filename": "pacbio.fastq",
      "ngs_type": "pacbio"
    }
  ]
}
EOF

cat > config_test_05_multiple_ngs_types/description.txt << 'EOF'
Config Test Case 5: Multiple NGS types
- Three different NGS types: illumina, ont, pacbio
- Tests handling of different sequencing platforms
- Expected: Should load all read types correctly
EOF

# ============================================================================
# Config Test Case 6: Missing config file (should fail)
# ============================================================================
echo "Creating config_test_06_missing_config..."
mkdir -p config_test_06_missing_config/data

cat > config_test_06_missing_config/description.txt << 'EOF'
Config Test Case 6: Missing config file
- No config.json file present
- Expected: Should raise FileNotFoundError
EOF

# ============================================================================
# Config Test Case 7: Missing required field (should fail)
# ============================================================================
echo "Creating config_test_07_missing_required_field..."
mkdir -p config_test_07_missing_required_field/data

cat > config_test_07_missing_required_field/data/ref_genome.fasta << 'EOF'
>chr1
ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG
EOF

cat > config_test_07_missing_required_field/data/mod_genome.fasta << 'EOF'
>chr1
ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG
EOF

# Config missing 'reads' field
cat > config_test_07_missing_required_field/data/config.json << 'EOF'
{
  "ref_genome_filename": "ref_genome.fasta",
  "mod_genome_filename": "mod_genome.fasta"
}
EOF

cat > config_test_07_missing_required_field/description.txt << 'EOF'
Config Test Case 7: Missing required field
- Config missing 'reads' field
- Expected: Should raise ConfigurationError about missing 'reads'
EOF

# ============================================================================
# Config Test Case 8: Missing genome file (should fail)
# ============================================================================
echo "Creating config_test_08_missing_file..."
mkdir -p config_test_08_missing_file/data

cat > config_test_08_missing_file/data/mod_genome.fasta << 'EOF'
>chr1
ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG
EOF

cat > config_test_08_missing_file/data/reads.fastq << 'EOF'
@read1
ATCGATCGATCGATCGATCG
+
IIIIIIIIIIIIIIIIIIII
EOF

# Config references non-existent file
cat > config_test_08_missing_file/data/config.json << 'EOF'
{
  "ref_genome_filename": "ref_genome.fasta",
  "mod_genome_filename": "mod_genome.fasta",
  "reads": [
    {
      "filename": "reads.fastq",
      "ngs_type": "illumina"
    }
  ]
}
EOF

cat > config_test_08_missing_file/description.txt << 'EOF'
Config Test Case 8: Missing genome file
- Config references ref_genome.fasta which doesn't exist
- Expected: Should raise FileNotFoundError
EOF

# ============================================================================
# Config Test Case 9: Invalid JSON (should fail)
# ============================================================================
echo "Creating config_test_09_invalid_json..."
mkdir -p config_test_09_invalid_json/data

# Create invalid JSON
cat > config_test_09_invalid_json/data/config.json << 'EOF'
{
  "ref_genome_filename": "ref_genome.fasta",
  "mod_genome_filename": "mod_genome.fasta",
  "reads": [
    {
      "filename": "reads.fastq"
      "ngs_type": "illumina"
    }
  ]
}
EOF

cat > config_test_09_invalid_json/description.txt << 'EOF'
Config Test Case 9: Invalid JSON
- Malformed JSON (missing comma)
- Expected: Should raise JSONDecodeError or ConfigurationError
EOF

# ============================================================================
# Config Test Case 10: Invalid NGS type (should fail)
# ============================================================================
echo "Creating config_test_10_invalid_ngs_type..."
mkdir -p config_test_10_invalid_ngs_type/data

cat > config_test_10_invalid_ngs_type/data/ref_genome.fasta << 'EOF'
>chr1
ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG
EOF

cat > config_test_10_invalid_ngs_type/data/mod_genome.fasta << 'EOF'
>chr1
ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG
EOF

cat > config_test_10_invalid_ngs_type/data/reads.fastq << 'EOF'
@read1
ATCGATCGATCGATCGATCG
+
IIIIIIIIIIIIIIIIIIII
EOF

# Config with invalid ngs_type
cat > config_test_10_invalid_ngs_type/data/config.json << 'EOF'
{
  "ref_genome_filename": "ref_genome.fasta",
  "mod_genome_filename": "mod_genome.fasta",
  "reads": [
    {
      "filename": "reads.fastq",
      "ngs_type": "454"
    }
  ]
}
EOF

cat > config_test_10_invalid_ngs_type/description.txt << 'EOF'
Config Test Case 10: Invalid NGS type
- Config specifies invalid ngs_type '454'
- Valid types are: illumina, ont, pacbio
- Expected: Should raise ValueError about invalid ngs_type
EOF

# No README needed - see tests/README.md for documentation

echo ""
echo "✓ Config manager test fixtures generated successfully!"
echo ""
echo "Generated test cases:"
ls -d config_test_* | while read dir; do
    echo "  - $dir"
done
echo ""
echo "To use these fixtures, run:"
echo "  pytest tests/test_integration.py -v"
echo ""
echo "For documentation, see: tests/README.md"
