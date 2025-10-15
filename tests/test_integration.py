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


if __name__ == "__main__":
    # Run this file directly with pytest
    pytest.main([__file__, "-v", "--tb=short"])
