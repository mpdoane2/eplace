#!/usr/bin/env python
"""
ePLACE: environmental Phylogenetic Localisation and Clade Estimation

Main command-line interface for ePLACE toolkit.
Provides unified access to database download, BLAST analysis, and grouped workflows.
"""

import sys
import argparse
import logging
from pathlib import Path
from collections import defaultdict

from .ncbi_download import setup_ncbi_database
from .blast_analysis import run_blast_search, FastaReader
from .taxonomy import (
    process_blast_results_for_taxonomy,
    rewrite_blast_hits,
    generate_classification_summary
)
from .alignment import (
    process_query_alignment_and_tree_parallel,
    process_grouped_alignment_and_tree_parallel,
    check_alignment_consistency,
    group_hits_by_group_rank,
    create_grouped_fasta_with_queries,
    IQTreeBuilder,
    concatenate_all_groups_and_build_tree
)

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)


def download_command(args):
    """Handle the download subcommand."""
    logger.info("=" * 60)
    logger.info("ePLACE Database Download")
    logger.info("=" * 60)
    
    # Download the database
    success, message = setup_ncbi_database(force_download=args.force)
    
    if success:
        logger.info(f"✓ {message}")
        return 0
    else:
        logger.error(f"✗ {message}")
        return 1


def blast_command(args):
    """Handle the blast subcommand - individual workflow."""
    # Validate input file
    if not args.query_fasta.exists():
        logger.error(f"Query FASTA file not found: {args.query_fasta}")
        return 1
    
    # Create output directory
    args.output_dir.mkdir(parents=True, exist_ok=True)
    
    # Set default output classification file if not provided
    if args.output_classification is None:
        base_name = args.query_fasta.stem
        for ext in ['.fasta', '.fa', '.fna', '.ffn', '.faa', '.frn']:
            if base_name.endswith(ext):
                base_name = base_name[:-len(ext)]
                break
        args.output_classification = args.output_dir / f"{base_name}_classification.tsv"
    if not args.output_classification.is_absolute():
        args.output_classification = args.output_dir / args.output_classification
    
    skip_existing = not args.overwrite_existing_blast

    logger.info("=" * 60)
    logger.info("ePLACE BLAST Workflow")
    logger.info("=" * 60)
    logger.info(f"Query FASTA: {args.query_fasta}")
    logger.info(f"Output directory: {args.output_dir}")
    logger.info(f"Taxonomic rank: {args.rank}")
    logger.info(f"Taxonomic rank for tree labeling: {args.tree_label_rank}")
    logger.info(f"Classification output file: {args.output_classification}")
    logger.info(f"Min identity: {args.min_identity}%")
    logger.info(f"Min coverage: {args.min_coverage}%")
    logger.info(f"Database: {args.database}")
    logger.info(f"Threads: {args.num_threads}")
    logger.info("=" * 60)
    
    # Step 1: Read query sequences
    logger.info("\n[Step 1/5] Reading query sequences...")
    try:
        sequences = FastaReader.read_fasta(args.query_fasta)
        logger.info(f"Found {len(sequences)} query sequences")
        for seq_id in sequences.keys():
            logger.info(f"  - {seq_id} ({len(sequences[seq_id])} bp)")
    except Exception:
        logger.exception(f"Error reading FASTA file")
        return 1
    
    # Step 2: Run BLAST search
    logger.info("\n[Step 2/5] Running BLAST search...")
    blast_output = args.output_dir / "blast_results.txt"
    
    try:
        success, filtered_hits = run_blast_search(
            query_fasta=args.query_fasta,
            output_file=blast_output,
            min_identity=args.min_identity,
            min_coverage=args.min_coverage,
            database=args.database,
            blastdb_path=args.blastdb_path,
            num_threads=args.num_threads,
            skip_existing=skip_existing
        )
        
        if not success:
            logger.error("BLAST search failed")
            return 1
        
        logger.info(f"BLAST search completed successfully")
        logger.info(f"Found {len(filtered_hits)} hits after filtering")
        
    except Exception:
        logger.exception(f"Error during BLAST search")
        return 1
    
    # Step 3: Group hits by query and display summary
    logger.info("\n[Step 3/5] Analyzing BLAST results...")
    hits_by_query = defaultdict(int)
    for hit in filtered_hits:
        hits_by_query[hit.query_id] += 1
    
    for query_id, count in hits_by_query.items():
        logger.info(f"  {query_id}: {count} hits")
    
    # Step 4: Extract representative sequences
    logger.info(f"\n[Step 4/5] Extracting representative sequences (rank: {args.rank})...")
    
    try:
        results = process_blast_results_for_taxonomy(
            blast_hits=filtered_hits,
            output_dir=args.output_dir,
            rank=args.rank,
            database=args.database,
            blastdb_path=args.blastdb_path
        )
        
        logger.info("\nResults:")
        for query_id, output_fasta in results.items():
            if output_fasta:
                logger.info(f"  {query_id}: {output_fasta}")
            else:
                logger.warning(f"  {query_id}: Failed to extract sequences")
        
    except Exception:
        logger.exception(f"Error extracting sequences")
        return 1

    logger.info("Rewriting the blast output file with the new annotations")
    try:
        rewrite_blast_hits(
            blast_hits=filtered_hits,
            output_file=args.output_dir / "blast_results_annotated.txt",
            header=True
        )
    except Exception:
        logger.exception(f"Error rewriting the blast hits")
        return 1

    # Step 5: Align sequences and build trees (if not skipped)
    if not args.skip_alignment:
        logger.info("\n[Step 5/5] Aligning sequences and building phylogenetic trees...")
        
        # Group hits by query for processing
        hits_by_query_map = defaultdict(list)
        for hit in filtered_hits:
            hits_by_query_map[hit.query_id].append(hit)
        
        # Process all queries to do trimming and alignment
        # and start tree building in background
        tree_jobs = []
        query_job_info = {}
        
        for query_id, query_hits in hits_by_query_map.items():
            logger.info(f"\nProcessing alignment for query: {query_id}")
            
            # Get the query directory
            safe_query_id = query_id.replace('|', '_').replace('/', '_')
            query_dir = args.output_dir / safe_query_id
            
            if not query_dir.exists():
                logger.warning(f"Query directory not found: {query_dir}")
                continue
            
            # Process alignment and tree
            try:
                result = process_query_alignment_and_tree_parallel(
                    query_id=query_id,
                    query_dir=query_dir,
                    blast_hits=query_hits,
                    taxonomic_rank=args.rank,
                    query_fasta=args.query_fasta,
                    num_threads=args.num_threads,
                    background_tree=True
                )
                
                # Log progress
                if result['alignment']:
                    logger.info(f"  ✓ Alignment: {result['alignment']}")
                if result['trimmed_fasta']:
                    logger.info(f"  ✓ Trimmed sequences: {result['trimmed_fasta']}")
                
                # Collect tree jobs for later waiting
                if result['tree_job']:
                    tree_jobs.append(result['tree_job'])
                    query_job_info[str(result['tree_file'])] = {
                        'query_id': query_id,
                        'tree_file': result['tree_file'],
                        'labeled_tree_path': result['labeled_tree_path'],
                        'blast_hits': result['blast_hits'],
                        'taxonomic_rank': result['taxonomic_rank']
                    }
                    logger.info(f"  ✓ Tree building started in background")
                    
            except Exception as e:
                logger.error(f"Error processing {query_id}: {e}")
                continue
        
        # Wait for all tree building jobs to complete
        if tree_jobs:
            logger.info(f"\n{'='*60}")
            logger.info(f"Waiting for {len(tree_jobs)} tree building jobs to complete...")
            logger.info(f"This may take a while depending on the size of the alignments.")
            logger.info(f"{'='*60}\n")
            
            tree_results = IQTreeBuilder.wait_for_tree_jobs(tree_jobs)
            
            # Relabel trees that completed successfully
            logger.info("\nRelabeling trees with taxonomic names...")
            for tree_path, success in tree_results.items():
                if success and tree_path in query_job_info:
                    job_info = query_job_info[tree_path]
                    query_id = job_info['query_id']
                    tree_file = job_info['tree_file']
                    labeled_tree = job_info['labeled_tree_path']
                    blast_hits = job_info['blast_hits']
                    taxonomic_rank = job_info['taxonomic_rank']
                    
                    logger.info(f"  Processing {query_id}...")
                    if IQTreeBuilder.relabel_tree_with_taxonomy(
                        tree_file=tree_file,
                        blast_hits=blast_hits,
                        output_tree=labeled_tree,
                        taxonomic_rank=taxonomic_rank
                    ):
                        logger.info(f"    ✓ Labeled tree: {labeled_tree}")
                    else:
                        logger.warning(f"    ! Failed to relabel tree, but unlabeled tree is available: {tree_file}")
        
        logger.info(f"\nAlignment and tree building completed for {len(hits_by_query_map)} queries")
    else:
        logger.info("\n[Step 5/5] Skipping alignment and tree building (--skip-alignment)")

    # Step 6: Generate classification summary TSV
    logger.info("\nGenerating classification summary TSV file...")
    try:
        success = generate_classification_summary(
            sequences=sequences,
            blast_hits=filtered_hits,
            output_file=args.output_classification,
            rank=args.rank,
            group_rank=args.rank,  # For individual workflow, group_rank same as rank
            tree_label_rank=args.tree_label_rank
        )
        
        if success:
            logger.info(f"✓ Classification summary: {args.output_classification}")
        else:
            logger.warning("Failed to generate classification summary")
    except Exception as e:
        logger.error(f"Error generating classification summary: {e}")

    logger.info("\n" + "=" * 60)
    logger.info("Workflow completed successfully!")
    logger.info("=" * 60)
    
    return 0


def grouped_command(args):
    """Handle the grouped subcommand - grouped workflow."""
    # Validate input file
    if not args.query_fasta.exists():
        logger.error(f"Query FASTA file not found: {args.query_fasta}")
        return 1
    
    # Create output directory
    args.output_dir.mkdir(parents=True, exist_ok=True)
    
    # Set default output classification file if not provided
    if args.output_classification is None:
        base_name = args.query_fasta.stem
        for ext in ['.fasta', '.fa', '.fna', '.ffn', '.faa', '.frn']:
            if base_name.endswith(ext):
                base_name = base_name[:-len(ext)]
                break
        args.output_classification = args.output_dir / f"{base_name}_classification.tsv"
    if not args.output_classification.is_absolute():
        args.output_classification = args.output_dir / args.output_classification
    
    skip_existing = not args.overwrite_existing_blast

    logger.info("=" * 60)
    logger.info("ePLACE Grouped BLAST Workflow")
    logger.info("=" * 60)
    logger.info(f"Query FASTA: {args.query_fasta}")
    logger.info(f"Output directory: {args.output_dir}")
    logger.info(f"Representative rank: {args.rank}")
    logger.info(f"Grouping rank: {args.group_rank}")
    logger.info(f"Tree labeling rank: {args.tree_label_rank}")
    logger.info(f"Classification output file: {args.output_classification}")
    logger.info(f"Min identity: {args.min_identity}%")
    logger.info(f"Min coverage: {args.min_coverage}%")
    logger.info(f"Database: {args.database}")
    logger.info(f"Threads: {args.num_threads}")
    logger.info("=" * 60)
    
    # Step 1: Read query sequences
    logger.info("\n[Step 1/9] Reading query sequences...")
    try:
        sequences = FastaReader.read_fasta(args.query_fasta)
        logger.info(f"Found {len(sequences)} query sequences")
        for seq_id in sequences.keys():
            logger.info(f"  - {seq_id} ({len(sequences[seq_id])} bp)")
    except Exception:
        logger.exception(f"Error reading FASTA file")
        return 1
    
    # Step 2: Run BLAST search
    logger.info("\n[Step 2/9] Running BLAST search...")
    blast_output = args.output_dir / "blast_results.txt"
    
    try:
        success, filtered_hits = run_blast_search(
            query_fasta=args.query_fasta,
            output_file=blast_output,
            min_identity=args.min_identity,
            min_coverage=args.min_coverage,
            database=args.database,
            blastdb_path=args.blastdb_path,
            num_threads=args.num_threads,
            skip_existing=skip_existing
        )
        
        if not success:
            logger.error("BLAST search failed")
            return 1
        
        logger.info(f"BLAST search completed successfully")
        logger.info(f"Found {len(filtered_hits)} hits after filtering")
        
    except Exception:
        logger.exception(f"Error during BLAST search")
        return 1
    
    # Step 3: Process taxonomy information
    logger.info(f"\n[Step 3/9] Processing taxonomy information (rank: {args.rank})...")
    
    try:
        results = process_blast_results_for_taxonomy(
            blast_hits=filtered_hits,
            output_dir=args.output_dir,
            rank=args.rank,
            database=args.database,
            blastdb_path=args.blastdb_path
        )
        
        logger.info("\nExtracted representative sequences per query:")
        for query_id, output_fasta in results.items():
            if output_fasta:
                logger.info(f"  {query_id}: {output_fasta}")
            else:
                logger.warning(f"  {query_id}: Failed to extract sequences")
        
    except Exception:
        logger.exception(f"Error extracting sequences")
        return 1

    logger.info("Rewriting the BLAST output file with the new annotations")
    try:
        rewrite_blast_hits(
            blast_hits=filtered_hits,
            output_file=args.output_dir / "blast_results_annotated.txt",
            header=True
        )
    except Exception as e:
        logger.exception(f"Error rewriting the blast hits: {e}")
        return 1

    # Step 4: Check alignment consistency
    logger.info("\n[Step 4/9] Checking alignment consistency...")
    consistency = check_alignment_consistency(
        blast_hits=filtered_hits,
        tolerance=args.alignment_tolerance
    )
    
    inconsistent_count = sum(1 for is_consistent in consistency.values() if not is_consistent)
    if inconsistent_count > 0:
        logger.warning(
            f"Found {inconsistent_count} reference sequences with inconsistent alignments "
            f"(tolerance: {args.alignment_tolerance} bp)"
        )
    else:
        logger.info("All reference sequences have consistent alignments across queries")

    # Step 5: Group hits by group_rank
    logger.info(f"\n[Step 5/9] Grouping hits by {args.group_rank}...")
    grouped_hits = group_hits_by_group_rank(filtered_hits, args.group_rank)
    
    if not grouped_hits:
        logger.error("No groups found after grouping by rank")
        return 1

    # Step 6: Create grouped FASTA files
    logger.info(f"\n[Step 6/9] Creating grouped FASTA files...")
    
    group_results = {}
    for group_tid, query_hits_map in grouped_hits.items():
        group_name = group_tid
        
        if not group_name:
            logger.warning(f"No group name found for group {group_tid}, skipping")
            continue
        
        logger.info(f"\nProcessing group: {group_name} ({group_tid})")
        logger.info(f"  Queries in this group: {', '.join(query_hits_map.keys())}")
        
        # Create group directory
        safe_group_name = group_name.replace(' ', '_').replace('/', '_').replace('|', '_')
        group_dir = args.output_dir / safe_group_name
        group_dir.mkdir(parents=True, exist_ok=True)
        
        # Create combined FASTA file
        combined_fasta = group_dir / f"{safe_group_name}_combined.fasta"
        success = create_grouped_fasta_with_queries(
            group_tid=group_tid,
            group_name=group_name,
            query_hits_map=query_hits_map,
            labeling_rank=args.rank,
            query_fasta=args.query_fasta,
            output_fasta=combined_fasta,
            database=args.database,
            blastdb_path=args.blastdb_path
        )
        
        if success:
            logger.info(f"  ✓ Created: {combined_fasta}")
            group_results[group_tid] = {
                'name': group_name,
                'dir': group_dir,
                'fasta': combined_fasta,
                'query_ids': list(query_hits_map.keys()),
                'hits': [hit for hits in query_hits_map.values() for hit in hits]
            }
        else:
            logger.error(f"  ✗ Failed to create grouped FASTA for {group_name}")

    # Step 7: Build alignments and trees for each group
    if not args.skip_alignment:
        logger.info("\n[Step 7/9] Building alignments and phylogenetic trees for each group...")
        
        # Process all groups to do trimming and alignment
        # and start tree building in background
        tree_jobs = []
        group_job_info = {}
        
        for group_tid, group_info in group_results.items():
            group_name = group_info['name']
            group_dir = group_info['dir']
            query_ids = group_info['query_ids']
            blast_hits = group_info['hits']
            
            logger.info(f"\nProcessing group: {group_name} ({group_tid})")
            
            try:
                result = process_grouped_alignment_and_tree_parallel(
                    group_name=group_name,
                    group_dir=group_dir,
                    taxonomic_rank=args.tree_label_rank,
                    blast_hits=blast_hits,
                    query_ids=query_ids,
                    num_threads=args.num_threads,
                    background_tree=True
                )
                
                # Log progress
                if result['alignment']:
                    logger.info(f"  ✓ Alignment: {result['alignment']}")
                if result['trimmed_fasta']:
                    logger.info(f"  ✓ Trimmed sequences: {result['trimmed_fasta']}")
                
                # Collect tree jobs for later waiting
                if result['tree_job']:
                    tree_jobs.append(result['tree_job'])
                    group_job_info[str(result['tree_file'])] = {
                        'group_name': group_name,
                        'tree_file': result['tree_file'],
                        'labeled_tree_path': result['labeled_tree_path'],
                        'blast_hits': result['blast_hits'],
                        'taxonomic_rank': result['taxonomic_rank']
                    }
                    logger.info(f"  ✓ Tree building started in background")
                    
            except Exception as e:
                logger.error(f"Error processing {group_name}: {e}")
                continue
        
        # Wait for all tree building jobs to complete
        if tree_jobs:
            logger.info(f"\n{'='*60}")
            logger.info(f"Waiting for {len(tree_jobs)} tree building jobs to complete...")
            logger.info(f"This may take a while depending on the size of the alignments.")
            logger.info(f"{'='*60}\n")
            
            tree_results = IQTreeBuilder.wait_for_tree_jobs(tree_jobs)
            
            # Relabel trees that completed successfully
            logger.info("\nRelabeling trees with taxonomic names...")
            for tree_path, success in tree_results.items():
                if success and tree_path in group_job_info:
                    job_info = group_job_info[tree_path]
                    group_name = job_info['group_name']
                    tree_file = job_info['tree_file']
                    labeled_tree = job_info['labeled_tree_path']
                    blast_hits = job_info['blast_hits']
                    taxonomic_rank = job_info['taxonomic_rank']
                    
                    logger.info(f"  Processing {group_name}...")
                    if IQTreeBuilder.relabel_tree_with_taxonomy(
                        tree_file=tree_file,
                        blast_hits=blast_hits,
                        output_tree=labeled_tree,
                        taxonomic_rank=taxonomic_rank
                    ):
                        logger.info(f"    ✓ Labeled tree: {labeled_tree}")
                    else:
                        logger.warning(f"    ! Failed to relabel tree, but unlabeled tree is available: {tree_file}")
        
        logger.info(f"\nAlignment and tree building completed for {len(group_results)} groups")
    else:
        logger.info("\n[Step 7/9] Skipping alignment and tree building (--skip-alignment)")

    # Step 8: Generate classification summary TSV
    logger.info("\n[Step 8/9] Generating classification summary TSV file...")
    try:
        success = generate_classification_summary(
            sequences=sequences,
            blast_hits=filtered_hits,
            output_file=args.output_classification,
            rank=args.rank,
            group_rank=args.group_rank,
            tree_label_rank=args.tree_label_rank
        )
        
        if success:
            logger.info(f"✓ Classification summary: {args.output_classification}")
        else:
            logger.warning("Failed to generate classification summary")
    except Exception as e:
        logger.error(f"Error generating classification summary: {e}")

    # Step 9: Build combined tree from all groups
    if not args.skip_alignment and group_results:
        logger.info("\n[Step 9/9] Building combined tree from all groups...")
        try:
            combined_results = concatenate_all_groups_and_build_tree(
                output_dir=args.output_dir,
                query_fasta=args.query_fasta,
                classification_file=args.output_classification,
                blast_hits=filtered_hits,
                combined_tree_label_rank=args.combined_tree_label_rank,
                num_threads=args.num_threads
            )
            
            if combined_results['alignment']:
                logger.info(f"✓ Combined alignment: {combined_results['alignment']}")
            if combined_results['tree']:
                logger.info(f"✓ Combined tree: {combined_results['tree']}")
            if combined_results['labeled_tree']:
                logger.info(f"✓ Combined labeled tree: {combined_results['labeled_tree']}")
        except Exception as e:
            logger.error(f"Error building combined tree: {e}")

    logger.info("\n" + "=" * 60)
    logger.info("Grouped workflow completed successfully!")
    logger.info(f"Processed {len(group_results)} taxonomic groups")
    logger.info("=" * 60)
    
    return 0


def main():
    """Main entry point for the ePLACE CLI."""
    parser = argparse.ArgumentParser(
        prog='eplace',
        description='ePLACE: environmental Phylogenetic Localisation and Clade Estimation',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Download NCBI database
  eplace download
  
  # Run individual BLAST workflow
  eplace blast query.fasta output_dir
  
  # Run grouped BLAST workflow
  eplace grouped query.fasta output_dir --group-rank order
  
For detailed help on each subcommand:
  eplace download --help
  eplace blast --help
  eplace grouped --help

Documentation: https://github.com/linsalrob/eplace
        """
    )
    
    parser.add_argument(
        '--version',
        action='version',
        version='%(prog)s 0.1.0'
    )
    
    subparsers = parser.add_subparsers(
        dest='command',
        title='Available commands',
        description='Use these commands to access different ePLACE functionalities',
        help='Command to run'
    )
    
    # Download subcommand
    download_parser = subparsers.add_parser(
        'download',
        help='Download NCBI BLAST database',
        description='Download and setup the NCBI core_nt BLAST database',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Download database to default location ($BLASTDB or ~/blastdb)
  eplace download
  
  # Force redownload even if database exists
  eplace download --force

Notes:
  - The download is large (several GB) and may take time
  - Database will be stored in $BLASTDB if set, otherwise ~/blastdb
  - MD5 checksums are verified automatically
        """
    )
    download_parser.add_argument(
        '--force',
        action='store_true',
        help='Force redownload even if database exists'
    )
    
    # BLAST subcommand (individual workflow)
    blast_parser = subparsers.add_parser(
        'blast',
        help='Run BLAST search with individual taxonomy analysis',
        description='Run BLAST search and extract representative sequences per query',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Basic usage with default parameters
  eplace blast query.fasta output_dir
  
  # Specify taxonomic rank and custom thresholds
  eplace blast query.fasta output_dir --rank genus --min-identity 95
  
  # Skip alignment and tree building (BLAST only)
  eplace blast query.fasta output_dir --skip-alignment
  
  # Use custom BLAST database location
  eplace blast query.fasta output_dir --blastdb-path /path/to/blastdb

Notes:
  - Creates one phylogenetic tree per query sequence
  - Default filtering: 90% identity, 80% query coverage
  - Requires BLAST+ tools and optionally MAFFT and IQTree
        """
    )
    blast_parser.add_argument(
        'query_fasta',
        type=Path,
        help='Path to query FASTA file'
    )
    blast_parser.add_argument(
        'output_dir',
        type=Path,
        help='Output directory for results'
    )
    blast_parser.add_argument(
        '--rank',
        type=str,
        default='genus',
        choices=['phylum', 'class', 'order', 'family', 'genus', 'species'],
        help='Taxonomic rank for representative selection (default: genus)'
    )
    blast_parser.add_argument(
        '--tree-label-rank',
        type=str,
        default='genus',
        choices=['phylum', 'class', 'order', 'family', 'genus', 'species'],
        help='Taxonomic rank for tree labeling (default: genus)'
    )
    blast_parser.add_argument(
        '--min-identity',
        type=float,
        default=90.0,
        help='Minimum percent identity for BLAST hits (default: 90.0)'
    )
    blast_parser.add_argument(
        '--min-coverage',
        type=float,
        default=80.0,
        help='Minimum query coverage percentage (default: 80.0)'
    )
    blast_parser.add_argument(
        '--database',
        type=str,
        default='core_nt',
        help='BLAST database name (default: core_nt)'
    )
    blast_parser.add_argument(
        '--blastdb-path',
        type=Path,
        default=None,
        help='Path to BLAST database directory'
    )
    blast_parser.add_argument(
        '--num-threads',
        type=int,
        default=1,
        help='Number of threads for BLAST and alignment (default: 1)'
    )
    blast_parser.add_argument(
        '--overwrite-existing-blast',
        action='store_true',
        help='Overwrite existing BLAST results'
    )
    blast_parser.add_argument(
        '--skip-alignment',
        action='store_true',
        help='Skip alignment and tree building steps'
    )
    blast_parser.add_argument(
        '--output-classification',
        type=Path,
        default=None,
        help='Path to output classification TSV file'
    )
    
    # Grouped subcommand
    grouped_parser = subparsers.add_parser(
        'grouped',
        help='Run BLAST search with grouped taxonomy analysis',
        description='Run BLAST search and group sequences by taxonomic rank',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Basic usage (group by class, default)
  eplace grouped query.fasta output_dir
  
  # Group by different taxonomic rank
  eplace grouped query.fasta output_dir --group-rank order
  
  # Specify both representative and grouping ranks
  eplace grouped query.fasta output_dir --rank genus --group-rank family
  
  # Skip alignment and tree building
  eplace grouped query.fasta output_dir --skip-alignment

Notes:
  - Creates one phylogenetic tree per taxonomic group
  - Groups queries that match the same taxonomic rank together
  - Useful for analyzing multiple related queries in one phylogenetic context
  - Default filtering: 90% identity, 80% query coverage
        """
    )
    grouped_parser.add_argument(
        'query_fasta',
        type=Path,
        help='Path to query FASTA file'
    )
    grouped_parser.add_argument(
        'output_dir',
        type=Path,
        help='Output directory for results'
    )
    grouped_parser.add_argument(
        '--rank',
        type=str,
        default='genus',
        choices=['phylum', 'class', 'order', 'family', 'genus', 'species'],
        help='Taxonomic rank for representative selection (default: genus)'
    )
    grouped_parser.add_argument(
        '--group-rank',
        type=str,
        default='class',
        choices=['phylum', 'class', 'order', 'family', 'genus', 'species'],
        help='Taxonomic rank for grouping sequences (default: class)'
    )
    grouped_parser.add_argument(
        '--tree-label-rank',
        type=str,
        default='genus',
        choices=['phylum', 'class', 'order', 'family', 'genus', 'species'],
        help='Taxonomic rank for tree labeling (default: genus)'
    )
    grouped_parser.add_argument(
        '--combined-tree-label-rank',
        type=str,
        default='genus',
        choices=['phylum', 'class', 'order', 'family', 'genus', 'species'],
        help='Taxonomic rank for tree labeling for the combined tree (default: genus)'
    )
    grouped_parser.add_argument(
        '--min-identity',
        type=float,
        default=90.0,
        help='Minimum percent identity for BLAST hits (default: 90.0)'
    )
    grouped_parser.add_argument(
        '--min-coverage',
        type=float,
        default=80.0,
        help='Minimum query coverage percentage (default: 80.0)'
    )
    grouped_parser.add_argument(
        '--database',
        type=str,
        default='core_nt',
        help='BLAST database name (default: core_nt)'
    )
    grouped_parser.add_argument(
        '--blastdb-path',
        type=Path,
        default=None,
        help='Path to BLAST database directory'
    )
    grouped_parser.add_argument(
        '--num-threads',
        type=int,
        default=1,
        help='Number of threads for BLAST and alignment (default: 1)'
    )
    grouped_parser.add_argument(
        '--overwrite-existing-blast',
        action='store_true',
        help='Overwrite existing BLAST results'
    )
    grouped_parser.add_argument(
        '--skip-alignment',
        action='store_true',
        help='Skip alignment and tree building steps'
    )
    grouped_parser.add_argument(
        '--alignment-tolerance',
        type=int,
        default=50,
        help='Maximum coordinate difference for alignment consistency (default: 50)'
    )
    grouped_parser.add_argument(
        '--output-classification',
        type=Path,
        default=None,
        help='Path to output classification TSV file'
    )
    
    # Parse arguments
    args = parser.parse_args()
    
    # If no command provided, show help
    if not args.command:
        parser.print_help()
        return 1
    
    # Route to appropriate command handler
    if args.command == 'download':
        return download_command(args)
    elif args.command == 'blast':
        return blast_command(args)
    elif args.command == 'grouped':
        return grouped_command(args)
    else:
        parser.print_help()
        return 1


if __name__ == '__main__':
    sys.exit(main())
