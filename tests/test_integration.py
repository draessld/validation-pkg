#!/usr/bin/env python3
"""
Integration Tests Runner

This module runs all integration tests using fixtures from the fixtures/ directory.
Integration tests verify complex pipelines across the package, such as:
- Loading config files
- Running validators
- Checking output files

Fixture Types:
- config_test_*: ConfigManager integration tests
- test_case_*: Validator integration tests (genome, read, feature)

Usage:
    pytest tests/test_integration.py              # Run all integration tests
    pytest tests/test_integration.py -k config    # Run only config tests
    pytest tests/test_integration.py -k genome    # Run only genome tests
    pytest tests/test_integration.py -v           # Verbose output
"""

import pytest
from pathlib import Path

# Import existing integration test classes
from tests.test_config_manager_integration import TestConfigManagerIntegration
from tests.test_genome_validator_integration import TestGenomeValidatorIntegration


class TestIntegrationSuite:
    """
    Master integration test suite that runs all fixture-based tests.

    This class serves as an entry point for all integration tests.
    Individual test implementations are in:
    - test_config_manager_integration.py
    - test_genome_validator_integration.py
    """

    @pytest.fixture
    def fixtures_dir(self):
        """Get path to fixtures directory."""
        return Path(__file__).parent / "fixtures"

    def test_fixtures_exist(self, fixtures_dir):
        """Verify that fixture directories exist."""
        assert fixtures_dir.exists(), "Fixtures directory not found"

        # Check for config fixtures
        config_fixtures = list(fixtures_dir.glob("config_test_*"))
        assert len(config_fixtures) > 0, "No config test fixtures found"

        # Check for validator fixtures
        validator_fixtures = list(fixtures_dir.glob("test_case_*"))
        assert len(validator_fixtures) > 0, "No validator test fixtures found"

        print(f"\nFound {len(config_fixtures)} config fixtures")
        print(f"Found {len(validator_fixtures)} validator fixtures")

    def test_config_fixtures_structure(self, fixtures_dir):
        """Verify config fixtures have required structure."""
        config_fixtures = list(fixtures_dir.glob("config_test_*"))

        for fixture in config_fixtures:
            # Each config fixture should have a data directory
            data_dir = fixture / "data"
            assert data_dir.exists(), f"{fixture.name} missing data/ directory"

            # Should have a description
            desc_file = fixture / "description.txt"
            assert desc_file.exists(), f"{fixture.name} missing description.txt"

    def test_validator_fixtures_structure(self, fixtures_dir):
        """Verify validator fixtures have required structure."""
        validator_fixtures = list(fixtures_dir.glob("test_case_*"))

        for fixture in validator_fixtures:
            # Should have a description
            desc_file = fixture / "description.txt"
            assert desc_file.exists(), f"{fixture.name} missing description.txt"

            # Should have at least one input file
            input_files = list(fixture.glob("*.*"))
            # Filter out description.txt
            input_files = [f for f in input_files if f.name != "description.txt"]
            assert len(input_files) > 0, f"{fixture.name} has no input files"


# Include all tests from config manager integration
class TestConfigIntegration(TestConfigManagerIntegration):
    """Config Manager integration tests from test_config_manager_integration.py"""
    pass


# Include all tests from genome validator integration
class TestGenomeIntegration(TestGenomeValidatorIntegration):
    """Genome Validator integration tests from test_genome_validator_integration.py"""
    pass


class TestFeatureValidatorIntegration:
    """
    Integration tests for FeatureValidator covering real-world scenarios.
    Tests full validation pipeline from input files to output files.
    """

    @pytest.fixture
    def temp_dir(self):
        """Create temporary directory for test files."""
        import tempfile
        with tempfile.TemporaryDirectory() as tmpdir:
            yield Path(tmpdir)

    @pytest.fixture
    def output_dir(self, temp_dir):
        """Create output directory."""
        out_dir = temp_dir / "output"
        out_dir.mkdir(parents=True, exist_ok=True)
        return out_dir

    def test_gff_complete_workflow(self, temp_dir, output_dir):
        """
        Test complete GFF validation workflow:
        - Parse GFF file
        - Validate coordinates
        - Sort by position
        - Collect statistics
        - Write output
        """
        from validation_pkg.config_manager import FeatureConfig
        from validation_pkg.validators.feature_validator import FeatureValidator
        from validation_pkg.utils.formats import FeatureFormat, CodingType

        # Create test GFF file
        gff_file = temp_dir / "test.gff"
        with open(gff_file, "w") as f:
            f.write("##gff-version 3\n")
            f.write("chr2\ttest\tgene\t1000\t2000\t.\t+\t.\tID=gene2;Name=TestGene2\n")
            f.write("chr1\ttest\tgene\t500\t1500\t.\t-\t.\tID=gene3;Name=TestGene3\n")
            f.write("chr1\ttest\tgene\t100\t500\t.\t+\t.\tID=gene1;Name=TestGene1\n")
            f.write("chr1\ttest\tCDS\t200\t400\t.\t+\t0\tID=cds1;Parent=gene1\n")

        # Configure validator
        feature_config = FeatureConfig(
            filename="test.gff",
            filepath=gff_file,
            coding_type=CodingType.NONE,
            detected_format=FeatureFormat.GFF
        )

        settings = FeatureValidator.Settings(
            sort_by_position=True,
            check_coordinates=True,
            allow_zero_length=False
        )

        validator = FeatureValidator(feature_config, output_dir, settings)
        validator.validate()

        # Verify results
        assert len(validator.features) == 4
        # Check sorting: chr1:100, chr1:200, chr1:500, chr2:1000
        assert validator.features[0].seqname == "chr1"
        assert validator.features[0].start == 100
        assert validator.features[3].seqname == "chr2"

        # Check statistics
        stats = validator.get_statistics()
        assert stats['num_features'] == 4
        assert stats['num_by_type']['gene'] == 3
        assert stats['num_by_type']['CDS'] == 1
        assert stats['num_by_strand']['+'] == 3
        assert stats['num_by_strand']['-'] == 1

        # Check output file exists
        output_file = output_dir / "test.gff"
        assert output_file.exists()

        # Verify output content
        with open(output_file, "r") as f:
            content = f.read()
            assert "##gff-version 3" in content
            assert "ID=gene1" in content

    def test_bed_to_gff_conversion_integration(self, temp_dir, output_dir):
        """
        Test BED to GFF conversion integration:
        - Parse BED file (0-based)
        - Convert to GFF (1-based)
        - Write GFF output
        """
        from validation_pkg.config_manager import FeatureConfig
        from validation_pkg.validators.feature_validator import FeatureValidator
        from validation_pkg.utils.formats import FeatureFormat, CodingType

        # Create test BED file (0-based coordinates)
        bed_file = temp_dir / "test.bed"
        with open(bed_file, "w") as f:
            f.write("chr1\t0\t1000\tfeature1\t100\t+\n")
            f.write("chr1\t1500\t2500\tfeature2\t200\t-\n")
            f.write("chr2\t500\t1500\tfeature3\t150\t+\n")

        feature_config = FeatureConfig(
            filename="test.bed",
            filepath=bed_file,
            coding_type=CodingType.NONE,
            detected_format=FeatureFormat.BED
        )

        validator = FeatureValidator(feature_config, output_dir)
        validator.validate()

        # Verify conversion to 1-based
        assert validator.features[0].start == 1  # 0 -> 1
        assert validator.features[0].end == 1000
        assert validator.features[1].start == 1501  # 1500 -> 1501
        assert validator.features[1].end == 2500

        # Verify GFF output
        output_file = output_dir / "test.gff"
        assert output_file.exists()

        with open(output_file, "r") as f:
            lines = f.readlines()
            assert lines[0].strip() == "##gff-version 3"
            # Check first feature has 1-based coordinates
            assert "\t1\t1000\t" in lines[1]

    def test_compressed_gff_workflow(self, temp_dir, output_dir):
        """
        Test compressed file handling:
        - Read gzip-compressed GFF
        - Process features
        - Write bzip2-compressed output
        """
        import gzip
        import bz2
        from validation_pkg.config_manager import FeatureConfig
        from validation_pkg.validators.feature_validator import FeatureValidator
        from validation_pkg.utils.formats import FeatureFormat, CodingType

        # Create compressed input
        gff_content = b"##gff-version 3\nchr1\ttest\tgene\t100\t500\t.\t+\t.\tID=gene1\n"
        gff_file = temp_dir / "test.gff.gz"
        with gzip.open(gff_file, "wb") as f:
            f.write(gff_content)

        feature_config = FeatureConfig(
            filename="test.gff.gz",
            filepath=gff_file,
            coding_type=CodingType.GZIP,
            detected_format=FeatureFormat.GFF
        )

        settings = FeatureValidator.Settings(coding_type='bz2')
        validator = FeatureValidator(feature_config, output_dir, settings)
        validator.validate()

        # Verify features parsed
        assert len(validator.features) == 1

        # Verify bz2 output
        output_file = output_dir / "test.gff.bz2"
        assert output_file.exists()

        with bz2.open(output_file, "rt") as f:
            content = f.read()
            assert "##gff-version 3" in content
            assert "ID=gene1" in content

    def test_replace_id_with_workflow(self, temp_dir, output_dir):
        """
        Test seqname replacement feature:
        - Replace all seqnames with custom value
        - Store original seqname in attributes
        """
        from validation_pkg.config_manager import FeatureConfig
        from validation_pkg.validators.feature_validator import FeatureValidator
        from validation_pkg.utils.formats import FeatureFormat, CodingType

        # Create GFF with multiple seqnames
        gff_file = temp_dir / "multi_chr.gff"
        with open(gff_file, "w") as f:
            f.write("##gff-version 3\n")
            f.write("scaffold_1\ttest\tgene\t100\t200\t.\t+\t.\tID=gene1\n")
            f.write("scaffold_2\ttest\tgene\t300\t400\t.\t+\t.\tID=gene2\n")
            f.write("contig_A\ttest\tgene\t500\t600\t.\t+\t.\tID=gene3\n")

        feature_config = FeatureConfig(
            filename="multi_chr.gff",
            filepath=gff_file,
            coding_type=CodingType.NONE,
            detected_format=FeatureFormat.GFF
        )

        settings = FeatureValidator.Settings(replace_id_with='chr1')
        validator = FeatureValidator(feature_config, output_dir, settings)
        validator.validate()

        # Verify all seqnames replaced
        for feature in validator.features:
            assert feature.seqname == 'chr1'

        # Check output
        output_file = output_dir / "multi_chr.gff"
        with open(output_file, "r") as f:
            content = f.read()
            # All lines should have chr1 as seqname
            for line in content.split('\n')[1:]:  # Skip header
                if line and not line.startswith('#'):
                    assert line.startswith('chr1\t')
            # Original seqnames should be in Note
            assert "Original_seqname:scaffold_1" in content
            assert "Original_seqname:scaffold_2" in content
            assert "Original_seqname:contig_A" in content

    def test_subdirectory_and_suffix_output(self, temp_dir, output_dir):
        """
        Test output customization:
        - Custom output subdirectory
        - Custom filename suffix
        """
        from validation_pkg.config_manager import FeatureConfig
        from validation_pkg.validators.feature_validator import FeatureValidator
        from validation_pkg.utils.formats import FeatureFormat, CodingType

        gff_file = temp_dir / "input.gff"
        with open(gff_file, "w") as f:
            f.write("##gff-version 3\n")
            f.write("chr1\ttest\tgene\t100\t200\t.\t+\t.\tID=gene1\n")

        feature_config = FeatureConfig(
            filename="input.gff",
            filepath=gff_file,
            coding_type=CodingType.NONE,
            detected_format=FeatureFormat.GFF
        )

        settings = FeatureValidator.Settings(
            output_subdir_name="features",
            output_filename_suffix="validated"
        )

        validator = FeatureValidator(feature_config, output_dir, settings)
        validator.validate()

        # Check output in subdirectory with suffix
        output_file = output_dir / "features" / "input_validated.gff"
        assert output_file.exists()

    def test_feature_length_property(self):
        """Test the Feature.length property."""
        from validation_pkg.validators.feature_validator import Feature

        feature = Feature(
            seqname="chr1",
            start=100,
            end=200,
            feature_type="gene"
        )

        assert feature.length == 100

        feature2 = Feature(
            seqname="chr1",
            start=500,
            end=1500,
            feature_type="CDS"
        )

        assert feature2.length == 1000


if __name__ == "__main__":
    # Run this file directly with pytest
    pytest.main([__file__, "-v", "--tb=short"])
