#!/bin/bash

# AME Motif Enrichment Analysis Script
# This script performs motif enrichment analysis on CNE regions using AME (Analysis of Motif Enrichment)
# from the MEME suite, running inside a Docker container
# 
# Usage: bash meme_motif_enrichment.sh <vertebrata_cne.bed> <teleostei_cne.bed> <genome.fa> [output_dir] [motif_db] [docker_image]
#
# Dependencies: Docker, docker image with MEME suite and bedtools
# Example docker image: memesuite/memesuite:latest
# Example motif database: /usr/local/share/meme/motif_databases/hocomoco/hocomoco11_core_HUMAN_mono_homer_format.motif

set -e

# Color codes for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

# Function to print colored output
print_status() {
    echo -e "${GREEN}[$(date '+%Y-%m-%d %H:%M:%S')]${NC} $1"
}

print_error() {
    echo -e "${RED}[ERROR]${NC} $1" >&2
}

print_warning() {
    echo -e "${YELLOW}[WARNING]${NC} $1"
}

# Check inputs
if [ $# -lt 3 ]; then
    print_error "Usage: bash meme_motif_enrichment.sh <vertebrata_cne.bed> <teleostei_cne.bed> <genome.fa> [output_dir] [motif_db] [docker_image]"
    exit 1
fi

VERTEBRATA_BED=$1
TELEOSTEI_BED=$2
GENOME_FA=$3
OUTPUT_DIR=${4:-.}
MOTIF_DB=${5:-}
DOCKER_IMAGE=${6:-memesuite/memesuite:latest}

# Get absolute paths
VERTEBRATA_BED=$(cd "$(dirname "$VERTEBRATA_BED")" && pwd)/$(basename "$VERTEBRATA_BED")
TELEOSTEI_BED=$(cd "$(dirname "$TELEOSTEI_BED")" && pwd)/$(basename "$TELEOSTEI_BED")
GENOME_FA=$(cd "$(dirname "$GENOME_FA")" && pwd)/$(basename "$GENOME_FA")
OUTPUT_DIR=$(mkdir -p "$OUTPUT_DIR" && cd "$OUTPUT_DIR" && pwd)

# Convert motif database path to absolute if it's a relative path
if [ ! -z "$MOTIF_DB" ] && [ ! -d "$(dirname "$MOTIF_DB")" ]; then
    # If directory doesn't exist, it might be a relative path
    if [[ "$MOTIF_DB" != /* ]]; then
        MOTIF_DB=$(cd "$(dirname "$MOTIF_DB")" && pwd)/$(basename "$MOTIF_DB")
    fi
fi

# Validate input files
for file in "$VERTEBRATA_BED" "$TELEOSTEI_BED" "$GENOME_FA"; do
    if [ ! -f "$file" ]; then
        print_error "Input file not found: $file"
        exit 1
    fi
done

# Check for Docker
if ! command -v docker &> /dev/null; then
    print_error "Docker is not installed or not in PATH"
    exit 1
fi

# Check if Docker image exists, pull if necessary
print_status "Checking Docker image: $DOCKER_IMAGE"
if ! docker image inspect "$DOCKER_IMAGE" &> /dev/null; then
    print_status "Pulling Docker image: $DOCKER_IMAGE"
    docker pull "$DOCKER_IMAGE"
fi

# --- Resolve motif DB: host file vs container-internal path ---

MOTIF_DB_HOST=""          # absolute host path if motif DB is on host
MOTIF_DB_CONT="$MOTIF_DB" # path we will pass to ame inside container

if [ -n "$MOTIF_DB" ] && [ -f "$MOTIF_DB" ]; then
    # motif DB is a host file -> mount it into an absolute container path
    MOTIF_DB_HOST="$(realpath "$MOTIF_DB")"
    MOTIF_DB_CONT="/motifs/$(basename "$MOTIF_DB_HOST")"
fi


# Create output directory
mkdir -p "$OUTPUT_DIR"
print_status "Output directory: $OUTPUT_DIR"

# Create subdirectories for organization
VERTEBRATA_DIR="$OUTPUT_DIR/vertebrata_motifs"
TELEOSTEI_DIR="$OUTPUT_DIR/teleostei_motifs"
BACKGROUND_DIR="$OUTPUT_DIR/background"
COMPARISON_DIR="$OUTPUT_DIR/motif_comparison"

mkdir -p "$VERTEBRATA_DIR" "$TELEOSTEI_DIR" "$BACKGROUND_DIR" "$COMPARISON_DIR"

# Extract sequences from BED files
print_status "Extracting sequences from vertebrata CNEs..."
docker run --rm \
    -v "$GENOME_FA":"$GENOME_FA":ro \
    -v "$GENOME_FA.fai":"$GENOME_FA.fai":ro \
    -v "$VERTEBRATA_BED":"$VERTEBRATA_BED":ro \
    -v "$VERTEBRATA_DIR":"$VERTEBRATA_DIR" \
    "$DOCKER_IMAGE" \
    bed2fasta -o "$VERTEBRATA_DIR/vertebrata_cne.fa" "$VERTEBRATA_BED" "$GENOME_FA"

print_status "Extracting sequences from teleostei CNEs..."
docker run --rm \
    -v "$GENOME_FA":"$GENOME_FA":ro \
    -v "$GENOME_FA.fai":"$GENOME_FA.fai":ro \
    -v "$TELEOSTEI_BED":"$TELEOSTEI_BED":ro \
    -v "$TELEOSTEI_DIR":"$TELEOSTEI_DIR" \
    "$DOCKER_IMAGE" \
    bed2fasta -o "$TELEOSTEI_DIR/teleostei_cne.fa" "$TELEOSTEI_BED" "$GENOME_FA"

# Validate extracted sequences
if [ ! -s "$VERTEBRATA_DIR/vertebrata_cne.fa" ]; then
    print_error "No sequences extracted for vertebrata CNEs. Check BED file coordinates."
    exit 1
fi

if [ ! -s "$TELEOSTEI_DIR/teleostei_cne.fa" ]; then
    print_error "No sequences extracted for teleostei CNEs. Check BED file coordinates."
    exit 1
fi

print_status "Vertebrata sequences: $(grep -c '^>' "$VERTEBRATA_DIR/vertebrata_cne.fa") sequences"
print_status "Teleostei sequences: $(grep -c '^>' "$TELEOSTEI_DIR/teleostei_cne.fa") sequences"

# Run AME (Analysis of Motif Enrichment) for vertebrata CNEs
print_status "Running AME on vertebrata CNEs..."
docker run --rm \
    -v "$VERTEBRATA_DIR":"$VERTEBRATA_DIR" \
    ${MOTIF_DB_HOST:+-v "$MOTIF_DB_HOST":"$MOTIF_DB_CONT":ro} \
    "$DOCKER_IMAGE" \
    ame \
    --method ranksum \
    --oc "$VERTEBRATA_DIR/ame_output" \
    "$VERTEBRATA_DIR/vertebrata_cne.fa" \
    "$MOTIF_DB_CONT"

print_status "✓ Vertebrata motif enrichment analysis complete"

# Run AME (Analysis of Motif Enrichment) for teleostei CNEs
print_status "Running AME on teleostei CNEs..."
docker run --rm \
    -v "$TELEOSTEI_DIR":"$TELEOSTEI_DIR" \
    ${MOTIF_DB_HOST:+-v "$MOTIF_DB_HOST":"$MOTIF_DB_CONT":ro} \
    "$DOCKER_IMAGE" \
    ame \
    --method ranksum \
    --oc "$TELEOSTEI_DIR/ame_output" \
    "$TELEOSTEI_DIR/teleostei_cne.fa" \
    "$MOTIF_DB_CONT"

print_status "✓ Teleostei motif enrichment analysis complete"

# Compare enrichment results between datasets
if [ -f "$VERTEBRATA_DIR/ame_output/ame.tsv" ] && [ -f "$TELEOSTEI_DIR/ame_output/ame.tsv" ]; then
    print_status "Comparing enrichment results between datasets..."
    
    # Create a comparison report
    cat > "$COMPARISON_DIR/enrichment_comparison.txt" << 'EOFCOMP'
Motif Enrichment Comparison: Vertebrata vs Teleostei CNEs
==========================================================
EOFCOMP
    
    echo "" >> "$COMPARISON_DIR/enrichment_comparison.txt"
    echo "Top 10 enriched motifs in Vertebrata CNEs:" >> "$COMPARISON_DIR/enrichment_comparison.txt"
    head -11 "$VERTEBRATA_DIR/ame_output/ame.tsv" >> "$COMPARISON_DIR/enrichment_comparison.txt"
    
    echo "" >> "$COMPARISON_DIR/enrichment_comparison.txt"
    echo "Top 10 enriched motifs in Teleostei CNEs:" >> "$COMPARISON_DIR/enrichment_comparison.txt"
    head -11 "$TELEOSTEI_DIR/ame_output/ame.tsv" >> "$COMPARISON_DIR/enrichment_comparison.txt"
    
    print_status "✓ Enrichment comparison complete"
fi

# Generate summary report
print_status "Generating summary report..."
REPORT="$OUTPUT_DIR/motif_enrichment_summary.txt"
cat > "$REPORT" << EOF
MEME Motif Enrichment Analysis Summary
======================================
Generated: $(date)

Input Files:
  Vertebrata CNEs: $VERTEBRATA_BED
  Teleostei CNEs: $TELEOSTEI_BED
  Reference Genome: $GENOME_FA

Analysis Details:
  - Analysis Method: AME (Analysis of Motif Enrichment)
  - Motif Database: $MOTIF_DB
  - Docker image: $DOCKER_IMAGE

Output Files:
  Vertebrata Enrichment Results: $VERTEBRATA_DIR/ame_output/
  Teleostei Enrichment Results: $TELEOSTEI_DIR/ame_output/
  Enrichment Comparison: $COMPARISON_DIR/enrichment_comparison.txt
  Extracted Sequences: 
    - $VERTEBRATA_DIR/vertebrata_cne.fa
    - $TELEOSTEI_DIR/teleostei_cne.fa

Next Steps:
1. Review ame.html in each ame_output/ directory for interactive results
2. Check ame.tsv for detailed enrichment p-values and fold changes
3. Compare enrichment_comparison.txt for top motifs in each dataset
4. Investigate specific enriched motifs for functional roles in fin development

For detailed results, see:
  - ame.html files in each ame_output directory
  - ame.tsv for complete enrichment statistics (p-value, E-value, fold change)
  - enrichment_comparison.txt for top 10 motifs comparison
EOF

print_status "Summary report saved to: $REPORT"

# Fix permissions so user can modify output directories without sudo
print_status "Fixing output directory permissions..."
chmod -R u+rwX,g+rX,o+rX "$OUTPUT_DIR"
find "$OUTPUT_DIR" -type d -exec chmod u+rwx,g+rx,o+rx {} \;
find "$OUTPUT_DIR" -type f -exec chmod u+rw,g+r,o+r {} \;

print_status "✓ All analyses complete!"
print_status "Results directory: $OUTPUT_DIR"
