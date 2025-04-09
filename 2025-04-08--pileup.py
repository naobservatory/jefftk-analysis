#!/usr/bin/env python3

# Claude wrote most of this, with my tweaks and suggestions.

import os
import json
import subprocess
import numpy as np
from PIL import Image
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from matplotlib.patches import Rectangle
from Bio.SeqIO.FastaIO import SimpleFastaParser
from collections import defaultdict

IN_1 = "reads_ff_1.fastq"
IN_2 = "reads_ff_2.fastq"
contigs_fname = "contigs.fasta"
seeds_fname = "seeds.fasta"
db_name = "db-contigs"
aligned_fname = "aligned.sam"
annotations_fname = "annotations.json"
attested_regions_fname = "attested_regions.fasta"  # New output file for attested regions

# Duplicate removal parameters
DUPLICATE_TOLERANCE = 3  # Allow up to this many bp difference for start/end positions

# Minimum read coverage for attested regions
MIN_COVERAGE_THRESHOLD = 3  # Require at least this many deduplicated seed-containing reads

# Minimum length for attested regions to output (to avoid tiny fragments)
MIN_REGION_LENGTH = 25  # Minimum length of attested region to output

# Color constants (RGBA)
BLUE = (0, 0, 255, 255)      # Matching regions
RED = (255, 0, 0, 255)       # Mismatches
BLACK = (0, 0, 0, 255)       # Insertions and deletions
LIGHT_GREY = (200, 200, 200, 255)  # Gaps between paired reads
GREEN = (0, 255, 0, 255)     # Seed regions
PURPLE = (128, 0, 128, 255)  # Non-seed portion of reads that contain the seed
WHITE = (255, 255, 255, 255)  # Background

# Define colors for annotation regions
ANNOTATION_COLORS = [
    '#8dd3c7', '#ffffb3', '#bebada', '#fb8072', '#80b1d3',
    '#fdb462', '#b3de69', '#fccde5', '#d9d9d9', '#bc80bd'
]

# Color for seed annotations
SEED_COLOR = '#00cc00'  # Bright green

def run_alignment():
    """Run bowtie2 alignment if needed"""
    contigs = {}
    with open(contigs_fname) as inf:
        for title, seq in SimpleFastaParser(inf):
            contigs[title] = seq

    if not os.path.exists("%s.1.bt2" % db_name):
        print("Building bowtie2 index...")
        subprocess.check_call([
            "bowtie2-build", contigs_fname, db_name])

    if not os.path.exists("aligned.sam"):
        print("Aligning reads...")
        subprocess.check_call([
            "bowtie2",
            "-x", db_name,
            "-1", IN_1,
            "-2", IN_2,
            "--local",
            "--threads", "8",
            "--no-sq",
            "--no-head",
            "--no-unal",
            "-S", aligned_fname,
        ])

    return contigs

def load_seeds():
    """Load seed sequences from FASTA file"""
    seeds = {}
    if not os.path.exists(seeds_fname):
        print(f"Warning: Seeds file {seeds_fname} not found")
        return seeds

    with open(seeds_fname) as inf:
        for title, seq in SimpleFastaParser(inf):
            seeds[title] = seq

    print(f"Loaded {len(seeds)} seed sequences")
    return seeds

def find_seed_positions(genome_name, genome_seq, seed_seq):
    """Find positions of seed sequence in genome"""
    if not seed_seq:
        return []

    positions = []
    seed_len = len(seed_seq)

    # Search for exact matches of seed in genome
    for i in range(len(genome_seq) - seed_len + 1):
        if genome_seq[i:i+seed_len] == seed_seq:
            positions.append((i, i + seed_len - 1))

    if not positions:
        print(f"Warning: Seed for {genome_name} not found in the genome sequence")
    else:
        print(f"Found seed for {genome_name} at positions: {positions}")

    return positions

def load_annotations():
    """Load genome annotations from JSON file"""
    if not os.path.exists(annotations_fname):
        print(f"Warning: Annotations file {annotations_fname} not found")
        return {}

    with open(annotations_fname, 'r') as f:
        annotations = json.load(f)

    # Convert 1-indexed to 0-indexed
    for genome_id, regions in annotations.items():
        for i, region in enumerate(regions):
            # Convert 1-indexed to 0-indexed
            region[0] -= 1
            annotations[genome_id][i] = region

    return annotations

def parse_cigar(cigar, pos, genome_name):
    """Parse CIGAR string to get alignment details"""
    aligned_positions = []
    aligned_colors = []
    insertions = []
    deletions = []

    genome_pos = pos
    read_pos = 0

    i = 0
    while i < len(cigar):
        # Find the number
        num_start = i
        while i < len(cigar) and cigar[i].isdigit():
            i += 1
        if i == len(cigar):
            break

        num = int(cigar[num_start:i])
        op = cigar[i]
        i += 1

        if op == 'M' or op == '=' or op == 'X':  # Match or mismatch
            for j in range(num):
                aligned_positions.append(genome_pos)
                if op == 'X':  # Mismatch
                    aligned_colors.append(RED)
                else:  # Match
                    aligned_colors.append(BLUE)
                genome_pos += 1
                read_pos += 1
        elif op == 'I':  # Insertion
            insertions.append((genome_pos, num))
            read_pos += num
        elif op == 'D':  # Deletion
            deletions.append((genome_pos, num))
            genome_pos += num
        elif op == 'S' or op == 'H':  # Soft or hard clip
            if op == 'S':
                read_pos += num
        elif op == 'N':  # Skipped region
            genome_pos += num

    return {
        'positions': aligned_positions,
        'colors': aligned_colors,
        'insertions': insertions,
        'deletions': deletions,
        'start_pos': pos,
        'end_pos': genome_pos
    }

def collect_read_pairs(genome_name):
    """Collect all read pairs aligned to the specified genome"""
    read_pairs = {}

    with open(aligned_fname) as inf:
        for line in inf:
            if line.startswith('@'):  # Skip header lines
                continue

            fields = line.strip().split('\t')
            qname = fields[0]
            flag = int(fields[1])
            rname = fields[2]
            pos = int(fields[3]) - 1  # Convert to 0-based
            cigar = fields[5]

            # Skip if not aligned to current genome
            if rname != genome_name:
                continue

            # Parse which read in the pair (1 or 2)
            is_first = flag & 0x40
            is_second = flag & 0x80
            is_reverse = flag & 0x10

            # Store alignments by read pair name
            if qname not in read_pairs:
                read_pairs[qname] = {'read1': None, 'read2': None}

            # Parse CIGAR string
            alignment_data = parse_cigar(cigar, pos, genome_name)
            alignment_data['is_reverse'] = is_reverse

            if is_first:
                read_pairs[qname]['read1'] = alignment_data
            else:
                read_pairs[qname]['read2'] = alignment_data

    # Filter incomplete pairs and sort by start position
    sorted_pairs = []
    for qname, pair_data in read_pairs.items():
        read1 = pair_data.get('read1')
        read2 = pair_data.get('read2')

        # Skip incomplete pairs
        if read1 is None or read2 is None:
            continue

        read1_start = min(read1['positions']) if read1['positions'] else read1['start_pos']
        read2_start = min(read2['positions']) if read2['positions'] else read2['start_pos']
        pair_start = min(read1_start, read2_start)

        # Add fragment boundaries information
        read1_end = max(read1['positions']) if read1['positions'] else read1['end_pos']
        read2_end = max(read2['positions']) if read2['positions'] else read2['end_pos']
        pair_end = max(read1_end, read2_end)

        pair_data['fragment_start'] = pair_start
        pair_data['fragment_end'] = pair_end

        sorted_pairs.append((qname, pair_data, pair_start))

    # Sort by start position
    sorted_pairs.sort(key=lambda x: x[2])

    return sorted_pairs

def remove_duplicates(sorted_pairs):
    """
    Remove sequencing duplicates by identifying read pairs where fragments
    align to the same positions (within DUPLICATE_TOLERANCE bp)
    """
    if not sorted_pairs:
        return []

    print(f"Removing duplicates from {len(sorted_pairs)} read pairs...")

    # Group read pairs by their fragment coordinates (with tolerance)
    fragments = defaultdict(list)

    for qname, pair_data, pair_start in sorted_pairs:
        # Round the fragment start/end to bins of size DUPLICATE_TOLERANCE
        # This creates "equivalence classes" for fragment coordinates
        binned_start = pair_data['fragment_start'] // DUPLICATE_TOLERANCE
        binned_end = pair_data['fragment_end'] // DUPLICATE_TOLERANCE

        # Use the binned coordinates as a key
        fragments[(binned_start, binned_end)].append((qname, pair_data, pair_start))

    # Select one representative from each group (the first one in each bin)
    deduplicated_pairs = []

    for fragment_group in fragments.values():
        # Take the first read pair from each group
        deduplicated_pairs.append(fragment_group[0])

    print(f"Removed {len(sorted_pairs) - len(deduplicated_pairs)} duplicates, {len(deduplicated_pairs)} unique fragments remain")

    # Re-sort by start position
    deduplicated_pairs.sort(key=lambda x: x[2])

    return deduplicated_pairs

def filter_seed_containing_reads(sorted_pairs, seed_positions):
    """Filter read pairs to only include those that overlap with seed regions"""
    if not seed_positions:
        return []

    # Create a set of seed positions for quick lookup
    seed_pos_set = set()
    for start, end in seed_positions:
        for pos in range(start, end + 1):
            seed_pos_set.add(pos)

    # Filter read pairs that overlap with seed regions
    seed_containing_pairs = []

    for qname, pair_data, pair_start in sorted_pairs:
        read1 = pair_data['read1']
        read2 = pair_data['read2']

        # Check if either read overlaps with seed regions
        contains_seed = False

        if read1 and read1['positions']:
            for pos in read1['positions']:
                if pos in seed_pos_set:
                    contains_seed = True
                    break

        if not contains_seed and read2 and read2['positions']:
            for pos in read2['positions']:
                if pos in seed_pos_set:
                    contains_seed = True
                    break

        if contains_seed:
            seed_containing_pairs.append((qname, pair_data, pair_start))

    print(f"Found {len(seed_containing_pairs)} read pairs containing seed regions")
    return seed_containing_pairs

def create_pileup_image(genome_name, genome_length, sorted_pairs, seed_positions, filter_by_seed=False, remove_dups=True):
    """Create a pixel-based pileup image with seed regions highlighted"""
    # Apply duplicate removal if requested
    if remove_dups:
        sorted_pairs = remove_duplicates(sorted_pairs)
        dedup_suffix = "_dedup"
    else:
        dedup_suffix = ""

    # Filter pairs if needed
    if filter_by_seed:
        pairs_to_use = filter_seed_containing_reads(sorted_pairs, seed_positions)
        output_suffix = f"_seed_containing{dedup_suffix}"
    else:
        pairs_to_use = sorted_pairs
        output_suffix = dedup_suffix

    # Create a blank white image
    height = len(pairs_to_use)
    width = genome_length

    if height == 0:
        print(f"No suitable read pairs found for genome {genome_name}{output_suffix}")
        return None

    print(f"Creating pileup image for {genome_name}{output_suffix} with {height} read pairs")

    # Use numpy for efficiency - one row per read pair
    img_array = np.full((height, width, 4), WHITE, dtype=np.uint8)

    # Create seed positions mask for quick lookup
    seed_mask = np.zeros(width, dtype=bool)
    for start, end in seed_positions:
        seed_mask[start:end+1] = True

    # First pass: Identify reads containing seed regions
    seed_containing_read_pairs = set()
    for y, (qname, pair_data, _) in enumerate(pairs_to_use):
        read1 = pair_data['read1']
        read2 = pair_data['read2']

        # Check if either read overlaps with seed regions
        if read1 and read1['positions']:
            for pos in read1['positions']:
                if 0 <= pos < width and seed_mask[pos]:
                    seed_containing_read_pairs.add(qname)
                    break

        if read2 and read2['positions'] and qname not in seed_containing_read_pairs:
            for pos in read2['positions']:
                if 0 <= pos < width and seed_mask[pos]:
                    seed_containing_read_pairs.add(qname)
                    break

    # Second pass: Color read pairs based on seed content
    for y, (qname, pair_data, _) in enumerate(pairs_to_use):
        read1 = pair_data['read1']
        read2 = pair_data['read2']

        # Check if this read pair contains seed
        read_pair_contains_seed = qname in seed_containing_read_pairs

        # Plot read1
        if read1 and read1['positions']:
            for pos, color in zip(read1['positions'], read1['colors']):
                if 0 <= pos < width:
                    if seed_mask[pos]:
                        # Seed regions in green
                        img_array[y, pos] = GREEN
                    elif read_pair_contains_seed:
                        # Non-seed portion of seed-containing read in purple
                        img_array[y, pos] = PURPLE
                    else:
                        # Normal alignment colors
                        img_array[y, pos] = color

        # Plot read2
        if read2 and read2['positions']:
            for pos, color in zip(read2['positions'], read2['colors']):
                if 0 <= pos < width:
                    if seed_mask[pos]:
                        # Seed regions in green
                        img_array[y, pos] = GREEN
                    elif read_pair_contains_seed:
                        # Non-seed portion of seed-containing read in purple
                        img_array[y, pos] = PURPLE
                    else:
                        # Normal alignment colors
                        img_array[y, pos] = color

        # Plot gap between reads if they don't overlap
        if read1 and read2:
            read1_end = max(read1['positions']) if read1['positions'] else read1['end_pos']
            read2_start = min(read2['positions']) if read2['positions'] else read2['start_pos']
            read2_end = max(read2['positions']) if read2['positions'] else read2['end_pos']
            read1_start = min(read1['positions']) if read1['positions'] else read1['start_pos']

            if read1_end < read2_start:
                for pos in range(read1_end, read2_start):
                    if 0 <= pos < width:
                        # Use purple for gaps in seed-containing read pairs, otherwise light grey
                        if read_pair_contains_seed and not filter_by_seed:
                            img_array[y, pos] = PURPLE  # Make gaps purple too for seed-containing reads
                        else:
                            img_array[y, pos] = LIGHT_GREY
            elif read2_end < read1_start:
                for pos in range(read2_end, read1_start):
                    if 0 <= pos < width:
                        # Use purple for gaps in seed-containing read pairs, otherwise light grey
                        if read_pair_contains_seed and not filter_by_seed:
                            img_array[y, pos] = PURPLE  # Make gaps purple too for seed-containing reads
                        else:
                            img_array[y, pos] = LIGHT_GREY

        # Plot insertions and deletions
        for pos, length in read1['insertions'] if read1 else []:
            if 0 <= pos < width:
                # Keep black for insertions, even in seed regions
                img_array[y, pos] = BLACK

        for pos, length in read2['insertions'] if read2 else []:
            if 0 <= pos < width:
                # Keep black for insertions, even in seed regions
                img_array[y, pos] = BLACK

        for pos, length in read1['deletions'] if read1 else []:
            for p in range(pos, min(pos + length, width)):
                if p >= 0:
                    # Keep black for deletions, even in seed regions
                    img_array[y, p] = BLACK

        for pos, length in read2['deletions'] if read2 else []:
            for p in range(pos, min(pos + length, width)):
                if p >= 0:
                    # Keep black for deletions, even in seed regions
                    img_array[y, p] = BLACK

    # Convert numpy array to PIL Image and save
    img = Image.fromarray(img_array)

    # If the image is very wide, we'll need to save it as a large image
    # PIL has dimension limits, so we'll check and resize if needed
    max_dimension = 65500  # PIL's maximum image dimension

    if width > max_dimension:
        print(f"Warning: Genome {genome_name} is too wide ({width} px). Rescaling width to {max_dimension}.")
        new_width = max_dimension
        new_height = int(height * (new_width / width))
        img = img.resize((new_width, new_height), Image.NEAREST)

    output_file = f"{genome_name}{output_suffix}_pileup.png"
    img.save(output_file)
    print(f"Saved pileup image to {output_file}")

    return pairs_to_use

def context_positions(positions, min_context):
    if not min_context:
        return positions

    return positions[min_context:-min_context]

def calculate_coverage(genome_length, sorted_pairs, min_context=0):
    """
    Calculate the coverage at each position along the genome.

    Args:
        genome_length: Length of the genome
        sorted_pairs: List of read pairs
        min_context: Minimum number of bases of context required

    Returns:
        Numpy array of coverage counts at each position
    """
    coverage = np.zeros(genome_length, dtype=np.int32)

    for _, pair_data, _ in sorted_pairs:
        read1 = pair_data['read1']
        read2 = pair_data['read2']

        # Get the observed positions from both reads
        observed_positions = set()

        if read1 and read2 and read1['positions'] and read2['positions'] and \
           min(read2['positions']) < max(read1['positions']):
            observed_positions = context_positions(
                sorted(set(read1['positions'] + read2['positions'])),
                min_context)
        else:
            if read1 and read1['positions']:
                observed_positions.update(context_positions(
                    read1['positions'], min_context))

            if read2 and read2['positions']:
                observed_positions.update(context_positions(
                    read2['positions'], min_context))

        # Increment coverage for each observed position
        for pos in observed_positions:
            assert 0 <= pos < genome_length
            coverage[pos] += 1

    return coverage

def add_seed_annotations(annotations, genome_name, seed_positions):
    """Add seed regions as annotations"""
    if not seed_positions:
        return

    if genome_name not in annotations:
        annotations[genome_name] = []

    # Add each seed position as a new annotation
    for i, (start, end) in enumerate(seed_positions):
        annotations[genome_name].append([start, end, f"Seed {i+1}"])

def create_coverage_plots(genome_name, genome_length, sorted_pairs, annotations, seed_positions, filter_by_seed=False, remove_dups=True):
    """Create coverage plots for the genome with annotations and seed regions"""
    # Apply duplicate removal if requested
    if remove_dups:
        sorted_pairs = remove_duplicates(sorted_pairs)
        dedup_suffix = "_dedup"
    else:
        dedup_suffix = ""

    # Filter pairs if needed
    if filter_by_seed:
        pairs_to_use = filter_seed_containing_reads(sorted_pairs, seed_positions)
        output_suffix = f"_seed_containing{dedup_suffix}"
        plot_title_suffix = " (Seed-Containing Reads Only, Duplicates Removed)" if remove_dups else " (Seed-Containing Reads Only)"
    else:
        pairs_to_use = sorted_pairs
        output_suffix = dedup_suffix
        plot_title_suffix = " (Duplicates Removed)" if remove_dups else ""

    if not pairs_to_use:
        print(f"No suitable read pairs found for genome {genome_name}{output_suffix}, skipping coverage plots")
        return

    print(f"Creating coverage plots for {genome_name}{output_suffix}")

    # Calculate coverage with different context requirements
    coverage_0bp = calculate_coverage(genome_length, pairs_to_use, min_context=0)
    coverage_25bp = calculate_coverage(genome_length, pairs_to_use, min_context=25)
    coverage_50bp = calculate_coverage(genome_length, pairs_to_use, min_context=50)

    # Create figure with three subplots
    plt.figure(figsize=(15, 10))
    gs = GridSpec(3, 1, figure=plt.gcf(), hspace=0.4)

    # Get genome annotations
    genome_annotations = annotations.get(genome_name, [])

    # Create a mapping of unique annotation labels to colors
    unique_annotations = {}
    for region in genome_annotations:
        annotation = region[2]
        # Check if it's a seed annotation
        if annotation.startswith("Seed"):
            unique_annotations[annotation] = SEED_COLOR
        elif annotation not in unique_annotations:
            unique_annotations[annotation] = ANNOTATION_COLORS[len(unique_annotations) % len(ANNOTATION_COLORS)]

    # Plot 1: Basic coverage (0bp context)
    ax1 = plt.subplot(gs[0])
    ax1.plot(range(genome_length), coverage_0bp, 'b-', linewidth=1)
    ax1.set_title(f"Coverage for {genome_name}{plot_title_suffix} (0bp context)")
    ax1.set_ylabel("Read pairs")
    ax1.set_xlim(0, genome_length)
    ax1.grid(True, linestyle='--', alpha=0.7)

    # Add annotations to plot 1
    add_annotations_to_plot(ax1, genome_annotations, unique_annotations, coverage_0bp)

    # Plot 2: Coverage with 25bp context
    ax2 = plt.subplot(gs[1])
    ax2.plot(range(genome_length), coverage_25bp, 'g-', linewidth=1)
    ax2.set_title(f"Coverage for {genome_name}{plot_title_suffix} (25bp context)")
    ax2.set_ylabel("Read pairs")
    ax2.set_xlim(0, genome_length)
    ax2.grid(True, linestyle='--', alpha=0.7)

    # Add annotations to plot 2
    add_annotations_to_plot(ax2, genome_annotations, unique_annotations, coverage_25bp)

    # Plot 3: Coverage with 50bp context
    ax3 = plt.subplot(gs[2])
    ax3.plot(range(genome_length), coverage_50bp, 'r-', linewidth=1)
    ax3.set_title(f"Coverage for {genome_name}{plot_title_suffix} (50bp context)")
    ax3.set_xlabel("Genome Position")
    ax3.set_ylabel("Read pairs")
    ax3.set_xlim(0, genome_length)
    ax3.grid(True, linestyle='--', alpha=0.7)

    # Add annotations to plot 3
    add_annotations_to_plot(ax3, genome_annotations, unique_annotations, coverage_50bp)

    # Add a legend for annotations if there are any
    if unique_annotations:
        handles = []
        for annotation, color in unique_annotations.items():
            handle = plt.Rectangle((0, 0), 1, 1, color=color, alpha=0.3, label=annotation)
            handles.append(handle)

        # Place legend at the bottom of the figure
        plt.figlegend(handles=handles, loc='lower center', ncol=min(len(handles), 3),
                      bbox_to_anchor=(0.5, 0.02), frameon=True)

        # Add some extra space at the bottom for the legend
        plt.subplots_adjust(bottom=0.15)

    # Adjust layout and save
    plt.tight_layout()
    output_file = f"{genome_name}{output_suffix}_coverage.png"
    plt.savefig(output_file, dpi=180)
    plt.close()

    print(f"Saved coverage plots to {output_file}")

def add_annotations_to_plot(ax, genome_annotations, unique_annotations, coverage):
    """Add annotation regions to the plot"""
    max_height = np.max(coverage) if np.any(coverage) else 1

    for region in genome_annotations:
        start_pos, end_pos, annotation = region
        color = unique_annotations.get(annotation, 'gray')

        # Add colored rectangle for the annotation region
        rect = Rectangle((start_pos, 0), end_pos - start_pos + 1, max_height,
                          facecolor=color, alpha=0.3, edgecolor='none')
        ax.add_patch(rect)

        # Add text label in the middle of the region
        mid_pos = (start_pos + end_pos) / 2

        # For seed regions, place the label at the top
        if annotation.startswith("Seed"):
            ax.text(mid_pos, max_height * 0.95, annotation,
                   horizontalalignment='center', verticalalignment='center',
                   fontsize=8, rotation=0, bbox=dict(facecolor='white', alpha=0.7, boxstyle='round,pad=0.2'))
        else:
            # For regular annotations, place in the middle
            ax.text(mid_pos, max_height * 0.5, annotation,
                   horizontalalignment='center', verticalalignment='center',
                   fontsize=8, rotation=0, bbox=dict(facecolor='white', alpha=0.7, boxstyle='round,pad=0.2'))

def find_attested_regions(genome_name, genome_seq, sorted_pairs, seed_positions):
    """
    Find regions of the genome that are attested by at least MIN_COVERAGE_THRESHOLD
    deduplicated seed-containing reads

    Returns:
        List of (start, end, coverage) tuples of attested regions
    """
    # We need deduplicated seed-containing reads
    deduped_pairs = remove_duplicates(sorted_pairs)
    seed_pairs = filter_seed_containing_reads(deduped_pairs, seed_positions)

    if not seed_pairs:
        print(f"No seed-containing read pairs found for {genome_name}")
        return []

    # Calculate coverage from the filtered read pairs
    coverage = calculate_coverage(len(genome_seq), seed_pairs)

    # Find regions with coverage above threshold
    attested_regions = []
    in_region = False
    start = 0

    for pos in range(len(genome_seq)):
        if coverage[pos] >= MIN_COVERAGE_THRESHOLD:
            if not in_region:
                # Start a new region
                in_region = True
                start = pos
        else:
            if in_region:
                # End current region
                end = pos - 1
                if end - start + 1 >= MIN_REGION_LENGTH:
                    # Only include regions that meet minimum length
                    avg_coverage = np.mean(coverage[start:end+1])
                    attested_regions.append((start, end, avg_coverage))
                in_region = False

    # Check if we were in a region at the end
    if in_region:
        end = len(genome_seq) - 1
        if end - start + 1 >= MIN_REGION_LENGTH:
            avg_coverage = np.mean(coverage[start:end+1])
            attested_regions.append((start, end, avg_coverage))

    print(f"Found {len(attested_regions)} attested regions for {genome_name}")
    return attested_regions

def write_attested_regions_fasta(contigs, output_fasta=attested_regions_fname):
    """
    Write attested regions from all genomes to a FASTA file
    """
    with open(output_fasta, 'w') as outf:
        for genome_name, genome_seq in sorted(contigs.items()):
            # Find seed positions for this genome
            seed_seq = seeds.get(genome_name)
            seed_positions = find_seed_positions(genome_name, genome_seq, seed_seq)

            # Collect read pairs and find attested regions
            sorted_pairs = collect_read_pairs(genome_name)
            attested_regions = find_attested_regions(genome_name, genome_seq, sorted_pairs, seed_positions)

            # Write each attested region to the FASTA file
            for i, (start, end, coverage) in enumerate(attested_regions):
                region_id = f"{genome_name}_region_{i+1}_{start+1}-{end+1}_cov{coverage:.1f}"
                region_seq = genome_seq[start:end+1]
                outf.write(f">{region_id}\n")

                # Write sequence with 60 characters per line
                for j in range(0, len(region_seq), 60):
                    outf.write(f"{region_seq[j:j+60]}\n")

            if not attested_regions:
                print(f"No attested regions found for {genome_name} with coverage >= {MIN_COVERAGE_THRESHOLD}")

    print(f"Wrote attested regions to {output_fasta}")

def main():
    # Run alignment if needed and get contigs
    contigs = run_alignment()

    # Load seeds
    global seeds
    seeds = load_seeds()

    # Load annotations
    annotations = load_annotations()
    print(f"Loaded annotations for {len(annotations)} genomes")

    # Process each genome
    for genome_name, genome_seq in sorted(contigs.items()):
        print(f"Processing genome: {genome_name}")
        genome_length = len(genome_seq)

        # Find seed positions in this genome
        seed_seq = seeds.get(genome_name)
        seed_positions = find_seed_positions(genome_name, genome_seq, seed_seq)

        # Add seed regions as annotations
        add_seed_annotations(annotations, genome_name, seed_positions)

        # Collect read pairs aligned to this genome
        sorted_pairs = collect_read_pairs(genome_name)

        if not sorted_pairs:
            print(f"No read pairs found for genome {genome_name}")
            continue

        # Create all the pileup and coverage visualizations
        # 1. Main pileup with duplicate removal
        create_pileup_image(genome_name, genome_length, sorted_pairs, seed_positions,
                           filter_by_seed=False, remove_dups=True)

        # 2. Main coverage plots with duplicate removal
        create_coverage_plots(genome_name, genome_length, sorted_pairs, annotations, seed_positions,
                             filter_by_seed=False, remove_dups=True)

        # 3. Seed-containing reads pileup with duplicate removal
        create_pileup_image(genome_name, genome_length, sorted_pairs, seed_positions,
                           filter_by_seed=True, remove_dups=True)

        # 4. Seed-containing reads coverage plots with duplicate removal
        create_coverage_plots(genome_name, genome_length, sorted_pairs, annotations, seed_positions,
                             filter_by_seed=True, remove_dups=True)

    # Generate the attested regions FASTA file
    write_attested_regions_fasta(contigs)

if __name__ == "__main__":
    main()
