"""
Example 13: Parallel Processing for Multiple Files
===================================================

Demonstrates how to use file-level parallelization to process multiple
files concurrently, significantly reducing processing time.

This example shows:
1. Unified threads parameter (RECOMMENDED) - automatic thread splitting
2. Basic parallel processing with validate_reads()
3. Using validate_genomes() for multiple genome files
4. Manual control with max_workers and compression_threads (power users)
5. Error handling in parallel mode
6. Combining parallelization with other optimizations

Performance comparison:
- Sequential: 8 files @ 30s each = 240 seconds
- Parallel (4 workers): 8 files = ~65 seconds (3.7x speedup)
- With unified threads=8: Automatically optimizes based on file count
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

def example_1_unified_threads_parameter():
    """
    Example 1: Using the unified threads parameter (RECOMMENDED).

    The threads parameter automatically splits between:
    - File-level parallelization (max_workers)
    - Compression-level parallelization (compression_threads)

    This is the simplest and recommended approach for most users.
    """
    print("=" * 80)
    print("Example 1: Unified Threads Parameter (RECOMMENDED)")
    print("=" * 80)

    # Load configuration
    config = ConfigManager.load("test_configs/parallel_reads.json")

    # Simple: Just set threads, automatic splitting
    settings = ReadValidator.Settings()
    settings = settings.update(
        threads=8,  # Automatically splits between file and compression parallelization
        validation_level='trust',  # Fast validation
        coding_type='gz'  # Compress output
    )

    print("\nUsing threads=8:")
    print("  - Automatically splits between file and compression parallelization")
    print("  - Smart strategy based on number of files")
    print("  - For 4 files: 4 workers × 2 compression threads each")
    print("  - For 1 file: 1 worker × 8 compression threads")
    print("  - For 8+ files: 8 workers × 1 compression thread each")

    # Process all read files
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


def example_2_basic_parallel_reads():
    """
    Example 2: Manual control with max_workers (traditional approach).

    Process multiple FASTQ files concurrently using max_workers.
    Power users can still set max_workers and compression_threads separately.
    """
    print("\n" + "=" * 80)
    print("Example 2: Manual Control with max_workers")
    print("=" * 80)

    # Load configuration
    config = ConfigManager.load("test_configs/parallel_reads.json")

    # Configure parallel processing manually
    settings = ReadValidator.Settings()
    settings = settings.update(
        max_workers=4,  # Use 4 parallel workers
        compression_threads=2,  # Use 2 threads per worker for compression
        validation_level='trust',  # Fast validation
        coding_type='gz'  # Compress output
    )

    print("\nUsing manual control:")
    print("  - max_workers=4: Process 4 files concurrently")
    print("  - compression_threads=2: Each worker uses 2 threads for compression")
    print("  - Total threads: 4 × 2 = 8 threads")

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


def example_3_parallel_genomes():
    """
    Example 3: Process reference and modified genomes in parallel.

    Useful when you have multiple genome files to validate.
    Uses the unified threads parameter for simplicity.
    """
    print("\n" + "=" * 80)
    print("Example 3: Parallel Genome Processing with Unified Threads")
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

    # Configure parallel processing with unified threads
    settings = GenomeValidator.Settings()
    settings = settings.update(
        threads=4,  # Automatically splits: 2 genomes → 2 workers × 2 threads each
        coding_type='gz',
        plasmid_split=True
    )

    print(f"\nProcessing {len(genome_list)} genomes with threads=4")
    print("  - Smart split: 2 workers × 2 compression threads each")

    # Process in parallel
    results = validate_genomes(genome_list, config.output_dir, settings)

    # Display results
    print(f"\nProcessed {len(results)} genome files:")
    for result in results:
        status = "✓" if result['success'] else "✗"
        print(f"{status} {result['filename']}")
        if not result['success']:
            print(f"  Error: {result['error']}")


def example_4_optimal_performance():
    """
    Example 4: Combine all optimizations for maximum performance.

    Demonstrates combining:
    - Unified threads parameter (automatic splitting)
    - Trust mode validation (validation_level='trust')
    - Parallel compression tools (pigz/pbzip2 auto-detected)
    - Config-level threads specification
    """
    print("\n" + "=" * 80)
    print("Example 4: Maximum Performance Configuration")
    print("=" * 80)

    config = ConfigManager.load("test_configs/large_dataset.json")

    # Get compression threads from config if available
    config_threads = config.get_threads() if hasattr(config, 'get_threads') else None

    # Configure for maximum performance using unified threads
    settings = ReadValidator.Settings()
    settings = settings.update(
        threads=8,  # Automatic split: 4 files → 4 workers × 2 threads each
        validation_level='trust',  # Fast validation (10-15x speedup)
        coding_type='gz'  # Use gzip compression
    )

    print("Performance optimizations enabled:")
    print(f"  - Unified threads: {settings.threads} (automatically split)")
    print(f"  - Validation level: {settings.validation_level} (10-15x faster)")
    print(f"  - Parallel compression: pigz/pbzip2 auto-detected")
    print("\nExpected speedup for 8 large FASTQ files:")
    print("  - Baseline (sequential, strict): ~240 seconds")
    print("  - With optimizations: ~15-20 seconds (12-16x faster)")

    results = validate_reads(config.reads, config.output_dir, settings, config_threads)

    print(f"\n✓ Processed {len(results)} files")


def example_5_error_handling():
    """
    Example 5: Handling errors in parallel mode.

    Shows how individual file failures don't stop other files.
    """
    print("\n" + "=" * 80)
    print("Example 5: Error Handling in Parallel Mode")
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


def example_6_sequential_fallback():
    """
    Example 6: Automatic fallback to sequential processing.

    Shows when parallel processing is not used:
    - Only 1 file
    - max_workers not set
    - max_workers = 1
    """
    print("\n" + "=" * 80)
    print("Example 6: Sequential Processing (Automatic Fallback)")
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
        example_1_unified_threads_parameter()
    except Exception as e:
        print(f"Example 1 skipped: {e}")

    try:
        example_2_basic_parallel_reads()
    except Exception as e:
        print(f"Example 2 skipped: {e}")

    try:
        example_3_parallel_genomes()
    except Exception as e:
        print(f"Example 3 skipped: {e}")

    try:
        example_4_optimal_performance()
    except Exception as e:
        print(f"Example 4 skipped: {e}")

    try:
        example_5_error_handling()
    except Exception as e:
        print(f"Example 5 skipped: {e}")

    try:
        example_6_sequential_fallback()
    except Exception as e:
        print(f"Example 6 skipped: {e}")

    print("\n" + "=" * 80)
    print("Key Takeaways:")
    print("=" * 80)
    print("1. Use threads=N for automatic splitting (RECOMMENDED)")
    print("2. threads parameter adapts to number of files automatically")
    print("3. Power users can still use max_workers + compression_threads manually")
    print("4. Automatic fallback to sequential for single files")
    print("5. Individual file errors don't stop other files")
    print("6. Combine with validation_level='trust' for maximum speed")
    print("7. Thread-safe logging works correctly in parallel mode")
    print("=" * 80)


if __name__ == "__main__":
    main()
