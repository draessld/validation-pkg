"""
Tests for unified threads parameter and smart thread splitting.

This module tests the new unified `threads` parameter that automatically
splits between file-level parallelization (max_workers) and compression-level
parallelization (compression_threads).
"""

import pytest
from validation_pkg.utils.file_handler import calculate_thread_distribution
from validation_pkg.validators.genome_validator import GenomeValidator
from validation_pkg.validators.read_validator import ReadValidator
from validation_pkg.validators.feature_validator import FeatureValidator


class TestCalculateThreadDistribution:
    """Tests for the calculate_thread_distribution function."""

    def test_single_file_all_threads_to_compression(self):
        """Test that with 1 file, all threads go to compression."""
        max_workers, compression_threads = calculate_thread_distribution(8, 1)
        assert max_workers == 1
        assert compression_threads == 8

    def test_two_files_balanced_split(self):
        """Test balanced split with 2 files."""
        max_workers, compression_threads = calculate_thread_distribution(16, 2)
        assert max_workers == 2
        assert compression_threads == 8

    def test_four_files_balanced_split(self):
        """Test balanced split with 4 files."""
        max_workers, compression_threads = calculate_thread_distribution(8, 4)
        assert max_workers == 4
        assert compression_threads == 2

    def test_many_files_prioritize_file_parallelization(self):
        """Test that with many files, prioritize file parallelization."""
        max_workers, compression_threads = calculate_thread_distribution(8, 20)
        assert max_workers == 8
        assert compression_threads == 1

    def test_more_files_than_threads(self):
        """Test when there are more files than threads."""
        max_workers, compression_threads = calculate_thread_distribution(4, 8)
        assert max_workers == 4
        assert compression_threads == 1

    def test_minimal_threads(self):
        """Test with minimal thread count."""
        max_workers, compression_threads = calculate_thread_distribution(2, 2)
        assert max_workers == 2
        assert compression_threads == 1

    def test_compression_threads_always_at_least_one(self):
        """Test that compression_threads is always at least 1."""
        # Even with many files
        max_workers, compression_threads = calculate_thread_distribution(4, 100)
        assert compression_threads >= 1

    def test_max_workers_never_exceeds_num_files(self):
        """Test that max_workers never exceeds number of files."""
        max_workers, compression_threads = calculate_thread_distribution(100, 4)
        assert max_workers <= 4

    def test_three_files_balanced(self):
        """Test 3 files with 12 threads."""
        max_workers, compression_threads = calculate_thread_distribution(12, 3)
        assert max_workers == 3
        assert compression_threads == 4


class TestUnifiedThreadsParameter:
    """Tests for the unified threads parameter in Settings classes."""

    def test_genome_settings_has_threads_parameter(self):
        """Test that GenomeValidator.Settings has threads parameter."""
        settings = GenomeValidator.Settings()
        assert hasattr(settings, 'threads')
        assert settings.threads is None  # Default value

    def test_read_settings_has_threads_parameter(self):
        """Test that ReadValidator.Settings has threads parameter."""
        settings = ReadValidator.Settings()
        assert hasattr(settings, 'threads')
        assert settings.threads is None

    def test_feature_settings_has_threads_parameter(self):
        """Test that FeatureValidator.Settings has threads parameter."""
        settings = FeatureValidator.Settings()
        assert hasattr(settings, 'threads')
        assert settings.threads is None

    def test_threads_parameter_can_be_set(self):
        """Test that threads parameter can be set."""
        settings = GenomeValidator.Settings(threads=8)
        assert settings.threads == 8

    def test_threads_parameter_in_update(self):
        """Test that threads parameter works with update()."""
        settings = GenomeValidator.Settings()
        updated = settings.update(threads=4)
        assert updated.threads == 4
        assert settings.threads is None  # Original unchanged

    def test_threads_parameter_in_to_dict(self):
        """Test that threads parameter is included in to_dict()."""
        settings = ReadValidator.Settings(threads=8)
        settings_dict = settings.to_dict()
        assert 'threads' in settings_dict
        assert settings_dict['threads'] == 8

    def test_threads_parameter_in_from_dict(self):
        """Test that threads parameter can be loaded from dict."""
        settings_dict = {'threads': 6, 'validation_level': 'trust'}
        settings = ReadValidator.Settings.from_dict(settings_dict)
        assert settings.threads == 6
        assert settings.validation_level == 'trust'


class TestThreadsPrecedence:
    """Tests for precedence of threads vs explicit max_workers/compression_threads."""

    def test_explicit_max_workers_overrides_threads(self):
        """Test that explicit max_workers takes precedence over threads."""
        settings = GenomeValidator.Settings(threads=8, max_workers=2)
        # threads should not be used if max_workers is explicit
        assert settings.max_workers == 2
        assert settings.threads == 8  # Both are stored

    def test_explicit_compression_threads_overrides_threads(self):
        """Test that explicit compression_threads takes precedence."""
        settings = ReadValidator.Settings(threads=8, compression_threads=4)
        assert settings.compression_threads == 4
        assert settings.threads == 8  # Both are stored

    def test_both_explicit_override_threads(self):
        """Test that both explicit settings override threads."""
        settings = FeatureValidator.Settings(
            threads=8,
            max_workers=3,
            compression_threads=2
        )
        assert settings.max_workers == 3
        assert settings.compression_threads == 2
        assert settings.threads == 8  # Still stored but won't be used


class TestThreadSplittingEdgeCases:
    """Tests for edge cases in thread splitting."""

    def test_threads_equal_to_files(self):
        """Test when threads equals number of files."""
        max_workers, compression_threads = calculate_thread_distribution(4, 4)
        assert max_workers == 4
        assert compression_threads == 1

    def test_one_thread_multiple_files(self):
        """Test with only 1 thread and multiple files."""
        max_workers, compression_threads = calculate_thread_distribution(1, 5)
        # Should process files sequentially with 1 thread
        assert max_workers == 1
        assert compression_threads == 1

    def test_large_thread_count(self):
        """Test with large number of threads."""
        max_workers, compression_threads = calculate_thread_distribution(64, 8)
        # Should distribute reasonably
        assert max_workers > 0
        assert compression_threads > 0
        # Product should be close to but not exceed total threads
        assert max_workers * compression_threads <= 64

    def test_five_files_eight_threads(self):
        """Test 5 files with 8 threads (sqrt heuristic)."""
        max_workers, compression_threads = calculate_thread_distribution(8, 5)
        # Should use balanced approach
        assert max_workers >= 2
        assert compression_threads >= 1


class TestBackwardCompatibility:
    """Tests to ensure backward compatibility."""

    def test_old_code_without_threads_still_works(self):
        """Test that code not using threads parameter still works."""
        settings = GenomeValidator.Settings(
            max_workers=4,
            compression_threads=2
        )
        assert settings.max_workers == 4
        assert settings.compression_threads == 2
        assert settings.threads is None

    def test_none_threads_doesnt_break(self):
        """Test that threads=None doesn't cause issues."""
        settings = ReadValidator.Settings(threads=None)
        assert settings.threads is None

    def test_default_settings_unchanged(self):
        """Test that default Settings() behavior is unchanged."""
        genome_settings = GenomeValidator.Settings()
        read_settings = ReadValidator.Settings()
        feature_settings = FeatureValidator.Settings()

        # All should have None for parallelization params
        assert genome_settings.threads is None
        assert genome_settings.max_workers is None
        assert genome_settings.compression_threads is None

        assert read_settings.threads is None
        assert feature_settings.threads is None


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
