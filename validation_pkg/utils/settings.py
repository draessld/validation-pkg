"""
Global settings and configuration for validators.

This file contains default settings for all validators including:
- Output directories
- Compression settings
- Editing specifications
- Validation thresholds
"""

from pathlib import Path
from typing import Dict, Any, Optional


# ===== GLOBALS =====
output_dir_name = "output"


# ===== Validator settings =====
class ValidatorSettings:
    """
    Centralized settings for all validators.
    
    Can be customized by creating an instance and modifying attributes,
    or by loading from a settings file.
    """
    
    def __init__(self):
        # ===== OUTPUT SETTINGS =====
        self.output_base_dir = "output"  # Base output directory
        self.output_subdirs = {
            'genomes': 'genomes',
            'plasmids': 'plasmids',
            'features': 'features',
            'reads': 'reads'
        }
        
        # ===== COMPRESSION SETTINGS =====
        # Override output compression (None means use config file setting)
        self.force_output_compression = None  # None, 'none', 'gz', 'bzip', 'tgz'
        
        # ===== GENOME VALIDATOR SETTINGS =====
        self.genome_settings = {
            # Editing specifications, default
            'edits': {
                'plasmid_split':True,
                'sequence_prefix':None,
                # 'remove_short_sequences': True,
                # 'min_sequence_length': 100,
                # 'convert_to_uppercase': True,
                # 'check_duplicate_ids': True,
                # 'remove_ambiguous_nucleotides': False,
                # 'ambiguous_threshold': 0.05,  # Max 5% ambiguous (N) allowed
                # 'trim_sequences': False,
                # 'trim_start': 0,
                # 'trim_end': 0,
                # 'filter_by_gc': False,
                # 'min_gc_content': 20.0,
                # 'max_gc_content': 80.0,
            },
            
            # Validation thresholds
            'validation': {
                'warn_gc_content_low': 20.0,
                'warn_gc_content_high': 80.0,
                'warn_sequence_length_low': 1000,
                'warn_sequence_length_high': 10000000,
                'allow_empty_sequences': False,
                'allow_duplicate_ids': True,  # Just warn, don't fail
            },
            
            # Output format
            'output': {
                'always_fasta': True,  # Always convert to FASTA
                'output_subdir_name':None,
                'output_filename_suffix':None,
                'line_width': 80,  # Characters per line in FASTA
            }
        }
        
        # ===== FEATURE VALIDATOR SETTINGS =====
        self.feature_settings = {
            'edits': {
                'sort_by_position': True,
                'remove_overlapping': False,
                'merge_adjacent': False,
                'filter_by_type': False,
                'allowed_types': ['gene', 'CDS', 'exon'],
            },
            
            'validation': {
                'check_coordinates': True,
                'allow_negative_strand': True,
                'check_feature_types': False,
                'required_attributes': [],
            },
            
            'output': {
                'preferred_format': 'gff3',  # 'gff3', 'gtf', or 'bed'
            }
        }
        
        # ===== READ VALIDATOR SETTINGS =====
        self.read_settings = {
            'validation': {
                'check_quality_scores': True,
                'min_quality_score': 20,
                'min_read_length': 50,
                'max_read_length': 1000000,
                'check_paired_end': True,
            },
            
            'edits': {
                'trim_quality': False,
                'quality_threshold': 20,
                'trim_adapters': False,
            }
        }
        
        # ===== INTER-FILE VALIDATION SETTINGS =====
        self.inter_file_settings = {
            'check_feature_bounds': True,  # Features within genome bounds
            'check_sequence_ids': True,    # Sequence IDs match
            'allow_missing_sequences': False,  # Features can reference missing seqs
        }