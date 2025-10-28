"""
Example 12: Parallel Processing for Multiple Files
===================================================

Demonstrates how to use file-level parallelization to process multiple
files concurrently, significantly reducing processing time.

This example shows:
1. Basic parallel processing with validate_reads()
2. Using validate_genomes() for multiple genome files
3. Error handling in parallel mode
4. Combining parallelization with other optimizations

Performance comparison:
- Sequential: 8 files @ 30s each = 240 seconds
- Parallel (4 workers): 8 files = ~65 seconds (3.7x speedup)
"""

from validation_pkg import (
    ConfigManager,
    validate_reads,
    validate_genomes,
    validate_features_list,
    ReadValidator,
    GenomeValidator,
    FeatureValidator,
    setup_logging
)

def example_1_basic_parallel_reads():
    """
    Example 1: Basic parallel processing of read files.

    Process multiple FASTQ files concurrently using max_workers.
    """
    print("=" * 80)
    print("Example 1: Basic Parallel Read Processing")
    print("=" * 80)

    # Load configuration
    config = ConfigManager.load("test_configs/parallel_reads.json")

    # Configure parallel processing
    settings = ReadValidator.Settings()
    settings = settings.update(
        max_workers=4,  # Use 4 parallel workers
        validation_level='trust',  # Fast validation
        coding_type='gz'  # Compress output
    )

    # Process all read files in parallel
    results = validate_reads(config.reads, config.output_dir, settings)

    # Check results
    print("\nResults:")
    successes = [r for r in results if r['success']]
    failures = [r for r in results if not r['success']]

    print(f"✓ Successful: {len(successes)}/{len(results)} files")
    for result in successes:
        print(f"  - {result['filename']}")

    if failures:
        print(f"✗ Failed: {len(failures)} files")
        for result in failures:
            print(f"  - {result['filename']}: {result['error']}")


def example_2_parallel_genomes():
    """
    Example 2: Process reference and modified genomes in parallel.

    Useful when you have multiple genome files to validate.
    """
    print("\n" + "=" * 80)
    print("Example 2: Parallel Genome Processing")
    print("=" * 80)

    # Load configuration
    config = ConfigManager.load("test_configs/genomes.json")

    # Collect all genome configs
    genome_list = []
    if config.ref_genome:
        genome_list.append(config.ref_genome)
    if config.mod_genome:
        genome_list.append(config.mod_genome)

    if not genome_list:
        print("No genome files found in config")
        return

    # Configure parallel processing
    settings = GenomeValidator.Settings()
    settings = settings.update(
        max_workers=2,  # 2 workers for 2 genomes
        coding_type='gz',
        plasmid_split=True
    )

    # Process in parallel
    results = validate_genomes(genome_list, config.output_dir, settings)

    # Display results
    print(f"\nProcessed {len(results)} genome files:")
    for result in results:
        status = "✓" if result['success'] else "✗"
        print(f"{status} {result['filename']}")
        if not result['success']:
            print(f"  Error: {result['error']}")


def example_3_optimal_performance():
    """
    Example 3: Combine all optimizations for maximum performance.

    Demonstrates combining:
    - File-level parallelization (max_workers)
    - Trust mode validation (validation_level='trust')
    - Parallel compression tools (pigz/pbzip2)
    - User-specified compression threads
    """
    print("\n" + "=" * 80)
    print("Example 3: Maximum Performance Configuration")
    print("=" * 80)

    config = ConfigManager.load("test_configs/large_dataset.json")

    # Get compression threads from config if available
    config_threads = config.get_threads() if hasattr(config, 'get_threads') else None

    # Configure for maximum performance
    settings = ReadValidator.Settings()
    settings = settings.update(
        max_workers=4,  # File-level parallelization
        validation_level='trust',  # Fast validation (10-15x speedup)
        coding_type='gz',  # Use gzip compression
        compression_threads=8  # Use 8 threads for pigz
    )

    print("Performance optimizations enabled:")
    print(f"  - File-level parallelization: {settings.max_workers} workers")
    print(f"  - Validation level: {settings.validation_level} (10-15x faster)")
    print(f"  - Compression threads: {settings.compression_threads} (uses pigz if available)")
    print("\nExpected speedup for 8 large FASTQ files:")
    print("  - Baseline (sequential, strict): ~240 seconds")
    print("  - With optimizations: ~15-20 seconds (12-16x faster)")

    results = validate_reads(config.reads, config.output_dir, settings, config_threads)

    print(f"\n✓ Processed {len(results)} files")


def example_4_error_handling():
    """
    Example 4: Handling errors in parallel mode.

    Shows how individual file failures don't stop other files.
    """
    print("\n" + "=" * 80)
    print("Example 4: Error Handling in Parallel Mode")
    print("=" * 80)

    config = ConfigManager.load("test_configs/mixed_quality.json")

    settings = ReadValidator.Settings()
    settings = settings.update(
        max_workers=4,
        check_invalid_chars=True  # Strict validation
    )

    results = validate_reads(config.reads, config.output_dir, settings)

    # Analyze results
    successes = [r for r in results if r['success']]
    failures = [r for r in results if not r['success']]

    print("\nProcessing Summary:")
    print(f"Total files: {len(results)}")
    print(f"Successful: {len(successes)}")
    print(f"Failed: {len(failures)}")

    if failures:
        print("\nFailed files:")
        for result in failures:
            print(f"  - {result['filename']}")
            print(f"    Error: {result['error']}")

    print("\nNote: In parallel mode, individual failures don't stop other files.")
    print("All files are processed, and errors are collected at the end.")


def example_5_sequential_fallback():
    """
    Example 5: Automatic fallback to sequential processing.

    Shows when parallel processing is not used:
    - Only 1 file
    - max_workers not set
    - max_workers = 1
    """
    print("\n" + "=" * 80)
    print("Example 5: Sequential Processing (Automatic Fallback)")
    print("=" * 80)

    config = ConfigManager.load("test_configs/single_file.json")

    # Case 1: Only one file (automatic fallback)
    settings = ReadValidator.Settings()
    settings = settings.update(max_workers=4)  # Request parallel

    print("Case 1: Single file with max_workers=4")
    print("  → Automatically falls back to sequential processing")

    results = validate_reads(config.reads, config.output_dir, settings)
    print(f"  ✓ Processed 1 file")

    # Case 2: max_workers not set (default behavior)
    print("\nCase 2: max_workers not set (None)")
    print("  → Uses sequential processing")

    settings_seq = ReadValidator.Settings()
    # max_workers is None by default
    results = validate_reads(config.reads, config.output_dir, settings_seq)
    print(f"  ✓ Processed 1 file")


def main():
    """Run all examples."""
    # Setup logging
    setup_logging(console_level="INFO")

    print("\n")
    print("╔" + "=" * 78 + "╗")
    print("║" + " " * 20 + "PARALLEL PROCESSING EXAMPLES" + " " * 30 + "║")
    print("╚" + "=" * 78 + "╝")
    print()

    try:
        example_1_basic_parallel_reads()
    except Exception as e:
        print(f"Example 1 skipped: {e}")

    try:
        example_2_parallel_genomes()
    except Exception as e:
        print(f"Example 2 skipped: {e}")

    try:
        example_3_optimal_performance()
    except Exception as e:
        print(f"Example 3 skipped: {e}")

    try:
        example_4_error_handling()
    except Exception as e:
        print(f"Example 4 skipped: {e}")

    try:
        example_5_sequential_fallback()
    except Exception as e:
        print(f"Example 5 skipped: {e}")

    print("\n" + "=" * 80)
    print("Key Takeaways:")
    print("=" * 80)
    print("1. Set max_workers=N to process N files concurrently")
    print("2. Automatic fallback to sequential for single files")
    print("3. Individual file errors don't stop other files")
    print("4. Combine with validation_level='trust' for maximum speed")
    print("5. Thread-safe logging works correctly in parallel mode")
    print("=" * 80)


if __name__ == "__main__":
    main()
