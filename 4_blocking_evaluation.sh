#!/bin/bash
# Unified Blocking Primer Evaluation with Taxonomic Information
# Exit on error
set -e

# Gene-specific configurations
declare -A GENE_CONFIGS

# 12S rRNA Configuration
GENE_CONFIGS["12S_blocking_primer"]="CTATGCTCAGCCCTAAACATAGATAGCTTACATAACAAAACTATCTGC"
GENE_CONFIGS["12S_blocking_name"]="12SV5_blocker"
GENE_CONFIGS["12S_f_primer"]="TTAGATACCCCACTATGC"
GENE_CONFIGS["12S_r_primer"]="TAGAACAGGCTCCTCTAG"
GENE_CONFIGS["12S_base_dir"]="/nas/collabs04/refseq/mockup_sequence/12S"
GENE_CONFIGS["12S_primer_path"]="primer_evaluation/12S_12SV5F_12SV5R"

# 16S rRNA Configuration  
GENE_CONFIGS["16S_blocking_primer"]="XXX"
GENE_CONFIGS["16S_blocking_name"]="Otter_Vert16S_blocker"
GENE_CONFIGS["16S_f_primer"]="AGACGAGAAGACCCYdTGGAGCTT"
GENE_CONFIGS["16S_r_primer"]="GATCCAACATCGAGGTCGTAA"
GENE_CONFIGS["16S_base_dir"]="/nas/collabs04/refseq/mockup_sequence/16S"
GENE_CONFIGS["16S_primer_path"]="primer_evaluation/16S_Vert-16S-eDNA-F1_Vert-16S-eDNA-R1"

# COI Configuration
GENE_CONFIGS["COI_blocking_primer"]="XXX"
GENE_CONFIGS["COI_blocking_name"]="coi175f345r_blocker"
GENE_CONFIGS["COI_f_primer"]="GGAGGCTTTGGMAAYTGRYT"
GENE_CONFIGS["COI_r_primer"]="TAGAGGRGGGTARACWGTYCA"
GENE_CONFIGS["COI_base_dir"]="/nas/collabs04/refseq/mockup_sequence/COI"
GENE_CONFIGS["COI_primer_path"]="primer_evaluation/COI_coi.175f_coi.345r"

# Available classes
CLASSES=(
    "Actinopterygii"
    "Amphibia"
    "Arachnida"
    "Aves"
    "Bivalvia"
    "Clitellata"
    "Gastropoda"
    "Insecta"
    "Malacostraca"
    "Mammalia"
    "Reptile"
)

# Function to show usage
show_usage() {
    echo "Usage: $0 <gene> [options]"
    echo ""
    echo "Required argument:"
    echo "  gene                Gene type: 12S, 16S, or COI"
    echo ""
    echo "Options:"
    echo "  --blocking-primer   Custom blocking primer sequence"
    echo "  --blocking-name     Custom blocking primer name"  
    echo "  --classes           Comma-separated list of classes (default: all)"
    echo "  -h, --help          Show this help message"
    echo ""
    echo "Available genes and configurations:"
    echo "  12S: 12SV5 blocker for 12S rRNA"
    echo "  16S: Otter Vert16S blocker for 16S rRNA"  
    echo "  COI: COI 175f/345r blocker for COI"
    echo ""
    echo "Available classes:"
    echo "  ${CLASSES[*]}"
    echo ""
    echo "Examples:"
    echo "  $0 12S"
    echo "  $0 16S --classes Actinopterygii,Mammalia"
    echo "  $0 COI --blocking-primer CUSTOM_SEQUENCE --blocking-name Custom_blocker"
}

# Function to log messages
log_message() {
    local gene=$1
    local class=$2
    local message=$3
    local blocking_eval_dir="${GENE_CONFIGS[${gene}_base_dir]}/primer_blocking_evaluation"
    local blocking_name="${GENE_CONFIGS[${gene}_blocking_name]}"
    local logfile="$blocking_eval_dir/${blocking_name}_${class}/${blocking_name}_${class}_blocking_evaluation.log"
    mkdir -p "$(dirname "$logfile")"
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] $message" | tee -a "$logfile"
}

# Function to analyze blocking efficiency with taxonomy
analyze_blocking_efficiency_with_taxonomy() {
    local gene=$1
    local class=$2
    local fasta_file=$3
    local taxonomy_file=$4
    local output_dir=$5

    local blocking_primer="${GENE_CONFIGS[${gene}_blocking_primer]}"
    local f_primer="${GENE_CONFIGS[${gene}_f_primer]}"
    local r_primer="${GENE_CONFIGS[${gene}_r_primer]}"
    local blocking_name="${GENE_CONFIGS[${gene}_blocking_name]}"

    log_message "$gene" "$class" "Analyzing blocking efficiency with taxonomic mapping..."

    # Analyze mismatches and map to taxonomy using Python script
    python3 << EOF
from Bio import SeqIO
import pandas as pd
import numpy as np
from collections import defaultdict
import re

def count_mismatches(primer, target_seq):
    """Count mismatches between primer and target sequence (forward direction only)"""
    primer = primer.upper()
    target_seq = target_seq.upper()

    # IUPAC code mapping
    iupac_codes = {
        'A': ['A'], 'T': ['T'], 'U': ['T'], 'G': ['G'], 'C': ['C'],
        'R': ['A', 'G'], 'Y': ['C', 'T'], 'S': ['G', 'C'], 'W': ['A', 'T'],
        'K': ['G', 'T'], 'M': ['A', 'C'], 'B': ['C', 'G', 'T'],
        'D': ['A', 'G', 'T'], 'H': ['A', 'C', 'T'], 'V': ['A', 'C', 'G'],
        'N': ['A', 'T', 'G', 'C']
    }

    def bases_match(primer_base, target_base):
        primer_bases = iupac_codes.get(primer_base, [primer_base])
        target_bases = iupac_codes.get(target_base, [target_base])
        return any(pb in target_bases or tb in primer_bases
                  for pb in primer_bases for tb in target_bases)

    # Find best match position in forward direction only
    min_mismatches = float('inf')
    best_position = -1

    for i in range(len(target_seq) - len(primer) + 1):
        target_window = target_seq[i:i + len(primer)]
        mismatches = sum(1 for p, t in zip(primer, target_window)
                        if not bases_match(p, t))

        if mismatches < min_mismatches:
            min_mismatches = mismatches
            best_position = i

    if len(target_seq) < len(primer):
        mismatches = sum(1 for p, t in zip(primer, target_seq)
                        if not bases_match(p, t))
        mismatches += len(primer) - len(target_seq)
        if mismatches < min_mismatches:
            min_mismatches = mismatches
            best_position = 0

    return min_mismatches, best_position

def parse_taxonomy_string(tax_string):
    """Parse QIIME2 taxonomy string into components"""
    if pd.isna(tax_string) or tax_string == '':
        return {'Kingdom': 'Unknown', 'Phylum': 'Unknown', 'Class': 'Unknown',
                'Order': 'Unknown', 'Family': 'Unknown', 'Genus': 'Unknown', 'Species': 'Unknown'}

    # Remove confidence scores (numbers in parentheses)
    clean_tax = re.sub(r'\([^)]*\)', '', tax_string)

    # Split by semicolon and clean up
    levels = [level.strip() for level in clean_tax.split(';')]

    # Map to taxonomic levels
    taxonomy = {'Kingdom': 'Unknown', 'Phylum': 'Unknown', 'Class': 'Unknown',
                'Order': 'Unknown', 'Family': 'Unknown', 'Genus': 'Unknown', 'Species': 'Unknown'}

    level_names = ['Kingdom', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species']

    for i, level in enumerate(levels[:len(level_names)]):
        if level and level != '' and not level.startswith('__'):
            # Remove prefixes like 'k__', 'p__', etc.
            clean_level = re.sub(r'^[kpcofgs]__', '', level)
            if clean_level and clean_level != '':
                taxonomy[level_names[i]] = clean_level

    return taxonomy

# Read sequences
sequences = []
seq_ids = []
for record in SeqIO.parse("$fasta_file", "fasta"):
    sequences.append(str(record.seq))
    seq_ids.append(record.id)

# Read taxonomy
taxonomy_dict = {}
try:
    taxonomy_df = pd.read_csv("$taxonomy_file", sep='\t', header=0)
    # Assume columns are: Feature ID, Taxon, Confidence (or similar)
    if 'Feature ID' in taxonomy_df.columns and 'Taxon' in taxonomy_df.columns:
        for _, row in taxonomy_df.iterrows():
            taxonomy_dict[row['Feature ID']] = row['Taxon']
    elif len(taxonomy_df.columns) >= 2:
        # Use first two columns if headers are different
        for _, row in taxonomy_df.iterrows():
            taxonomy_dict[row.iloc[0]] = row.iloc[1]
    print(f"Loaded taxonomy for {len(taxonomy_dict)} sequences")
except Exception as e:
    print(f"Warning: Could not load taxonomy file: {e}")
    print("Proceeding without taxonomic information")

blocking_primer = "$blocking_primer"

print(f"\\n=== BLOCKING PRIMER ANALYSIS WITH TAXONOMY ===")
print(f"Gene: $gene")
print(f"Class: $class")
print(f"PCR Primers used for extraction: $f_primer / $r_primer")
print(f"Blocking primer for analysis: {blocking_primer}")
print(f"Blocking primer length: {len(blocking_primer)} bp")
print(f"Sequences to analyze: {len(sequences)}")
print(f"Analysis: Testing blocking primer binding on [primer_binding_region--amplicon--primer_binding_region] regions")

# Analyze each sequence
results = []
for i, seq in enumerate(sequences):
    mismatches, position = count_mismatches(blocking_primer, seq)

    # Get taxonomy for this sequence
    tax_string = taxonomy_dict.get(seq_ids[i], '')
    taxonomy = parse_taxonomy_string(tax_string)

    result = {
        'sequence_id': seq_ids[i],
        'sequence_length': len(seq),
        'mismatches': mismatches,
        'position': position,
        'blocking_category': 'excellent' if mismatches <= 1 else 'moderate' if 2 <= mismatches <= 3 else 'poor',
        'kingdom': taxonomy['Kingdom'],
        'phylum': taxonomy['Phylum'],
        'class': taxonomy['Class'],
        'order': taxonomy['Order'],
        'family': taxonomy['Family'],
        'genus': taxonomy['Genus'],
        'species': taxonomy['Species'],
        'full_taxonomy': tax_string
    }
    results.append(result)

# Calculate overall statistics
total = len(results)
if total == 0:
    print("No sequences to analyze!")
    exit(1)

excellent = sum(1 for r in results if r['mismatches'] <= 1)
moderate = sum(1 for r in results if 2 <= r['mismatches'] <= 3)
poor = sum(1 for r in results if r['mismatches'] >= 4)

print(f"\\n=== OVERALL RESULTS FOR {str('$class').upper()} ===")
print(f"Total sequences: {total}")
print(f"Excellent blocking (0-1 mismatches): {excellent} ({excellent/total*100:.1f}%)")
print(f"Moderate blocking (2-3 mismatches): {moderate} ({moderate/total*100:.1f}%)")
print(f"Poor blocking (4+ mismatches): {poor} ({poor/total*100:.1f}%)")

# Save detailed results with taxonomy
df = pd.DataFrame(results)
df.to_csv("$output_dir/${blocking_name}_${class}_detailed_results_with_taxonomy.csv", index=False)

# Create taxonomic summaries
print(f"\\n=== BLOCKING EFFICIENCY BY TAXONOMY ===")

# Summary by Genus
if 'genus' in df.columns:
    genus_summary = df.groupby(['genus', 'blocking_category']).size().unstack(fill_value=0)
    genus_summary['total'] = genus_summary.sum(axis=1)
    genus_summary['excellent_pct'] = (genus_summary.get('excellent', 0) / genus_summary['total'] * 100).round(1)
    genus_summary['moderate_pct'] = (genus_summary.get('moderate', 0) / genus_summary['total'] * 100).round(1)
    genus_summary['poor_pct'] = (genus_summary.get('poor', 0) / genus_summary['total'] * 100).round(1)

    # Sort by excellent percentage (descending)
    genus_summary = genus_summary.sort_values('excellent_pct', ascending=False)
    genus_summary.to_csv("$output_dir/${blocking_name}_${class}_genus_summary.csv")

    print("\\nTop 10 Genera by Excellent Blocking:")
    print(genus_summary[['total', 'excellent_pct', 'moderate_pct', 'poor_pct']].head(10).to_string())

# Summary by Family
if 'family' in df.columns:
    family_summary = df.groupby(['family', 'blocking_category']).size().unstack(fill_value=0)
    family_summary['total'] = family_summary.sum(axis=1)
    family_summary['excellent_pct'] = (family_summary.get('excellent', 0) / family_summary['total'] * 100).round(1)
    family_summary['moderate_pct'] = (family_summary.get('moderate', 0) / family_summary['total'] * 100).round(1)
    family_summary['poor_pct'] = (family_summary.get('poor', 0) / family_summary['total'] * 100).round(1)

    family_summary = family_summary.sort_values('excellent_pct', ascending=False)
    family_summary.to_csv("$output_dir/${blocking_name}_${class}_family_summary.csv")

    print("\\nTop 10 Families by Excellent Blocking:")
    print(family_summary[['total', 'excellent_pct', 'moderate_pct', 'poor_pct']].head(10).to_string())

# Show specific examples of excellent and poor blocking
print(f"\\n=== SPECIFIC EXAMPLES ===")

excellent_examples = df[df['blocking_category'] == 'excellent'].head(5)
if len(excellent_examples) > 0:
    print("\\nExcellent Blocking Examples:")
    for _, row in excellent_examples.iterrows():
        print(f"  {row['genus']} {row['species']}: {row['mismatches']} mismatches")

moderate_examples = df[df['blocking_category'] == 'moderate'].head(5)
if len(moderate_examples) > 0:
    print("\\nModerate Blocking Examples:")
    for _, row in moderate_examples.iterrows():
        print(f"  {row['genus']} {row['species']}: {row['mismatches']} mismatches")

poor_examples = df[df['blocking_category'] == 'poor'].head(5)
if len(poor_examples) > 0:
    print("\\nPoor Blocking Examples:")
    for _, row in poor_examples.iterrows():
        print(f"  {row['genus']} {row['species']}: {row['mismatches']} mismatches")

# Save summary for overall report
summary = {
    'gene': '$gene',
    'class': '$class',
    'blocking_primer': blocking_primer,
    'primer_length': len(blocking_primer),
    'total_sequences': total,
    'excellent_blocking': excellent,
    'moderate_blocking': moderate,
    'poor_blocking': poor,
    'excellent_percentage': round(excellent/total*100, 2) if total > 0 else 0,
    'moderate_percentage': round(moderate/total*100, 2) if total > 0 else 0,
    'poor_percentage': round(poor/total*100, 2) if total > 0 else 0,
    'mean_mismatches': round(np.mean([r['mismatches'] for r in results]), 2) if results else 0,
    'median_mismatches': round(np.median([r['mismatches'] for r in results]), 2) if results else 0
}

summary_df = pd.DataFrame([summary])
summary_df.to_csv("$output_dir/${blocking_name}_${class}_summary.csv", index=False)

print(f"\\nFiles saved:")
print(f"- Detailed results with taxonomy: ${blocking_name}_${class}_detailed_results_with_taxonomy.csv")
print(f"- Summary for overall report: ${blocking_name}_${class}_summary.csv")
print(f"- Genus summary: ${blocking_name}_${class}_genus_summary.csv")
print(f"- Family summary: ${blocking_name}_${class}_family_summary.csv")

EOF

    log_message "$gene" "$class" "Blocking efficiency analysis with taxonomy completed"
}

# Main blocking evaluation function
evaluate_blocking_primer_for_class() {
    local gene=$1
    local class=$2

    local blocking_primer="${GENE_CONFIGS[${gene}_blocking_primer]}"
    local blocking_name="${GENE_CONFIGS[${gene}_blocking_name]}"
    local f_primer="${GENE_CONFIGS[${gene}_f_primer]}"
    local r_primer="${GENE_CONFIGS[${gene}_r_primer]}"
    local base_dir="${GENE_CONFIGS[${gene}_base_dir]}"
    local primer_path="${GENE_CONFIGS[${gene}_primer_path]}"

    echo "Processing blocking primer for gene: $gene, class: $class"
    echo "Blocking primer: $blocking_primer"

    # Set output directory
    local blocking_eval_dir="$base_dir/primer_blocking_evaluation"
    local output_dir="$blocking_eval_dir/${blocking_name}_${class}"
    mkdir -p "$output_dir"

    # Input files
    local derep3_sequence="$primer_path}_${class}_Thailand/${gene}_*_${class}_Thailand_sequence_derep3.qza"
    local derep1_sequence="$primer_path}_${class}_Thailand/${gene}_*_${class}_Thailand_sequence_derep1.qza"
    local taxonomy="$primer_path}_${class}_Thailand/${gene}_*_${class}_Thailand_taxonomy_derep1.qza"

    # Find actual files (handle wildcard patterns)
    derep3_sequence=$(find "$base_dir" -name "*_${class}_Thailand_sequence_derep3.qza" | head -1)
    derep1_sequence=$(find "$base_dir" -name "*_${class}_Thailand_sequence_derep1.qza" | head -1)
    taxonomy=$(find "$base_dir" -name "*_${class}_Thailand_taxonomy_derep1.qza" | head -1)

    # Check if input files exist
    if [[ ! -f "$derep3_sequence" ]]; then
        echo "✗ derep3 file not found for $class"
        return 1
    fi

    if [[ ! -f "$derep1_sequence" ]]; then
        echo "✗ derep1 file not found for $class"
        return 1
    fi

    if [[ ! -f "$taxonomy" ]]; then
        echo "⚠ taxonomy file not found for $class (proceeding without taxonomy)"
        taxonomy=""
    fi

    echo "✓ Input files found"
    log_message "$gene" "$class" "Starting blocking evaluation with primer binding region extraction and taxonomy"

    # Step 1: Export derep3 to get sequence IDs
    log_message "$gene" "$class" "Extracting sequence IDs from derep3"
    qiime tools export \
        --input-path "$derep3_sequence" \
        --output-path "$output_dir/derep3_export"

    if [[ ! -f "$output_dir/derep3_export/dna-sequences.fasta" ]]; then
        echo "✗ Failed to export derep3 sequences"
        return 1
    fi

    grep "^>" "$output_dir/derep3_export/dna-sequences.fasta" | sed 's/^>//' > "$output_dir/amplifiable_seq_ids.txt"
    local amplifiable_count=$(wc -l < "$output_dir/amplifiable_seq_ids.txt")
    echo "✓ Found $amplifiable_count sequences that passed primer filtering"

    if [ $amplifiable_count -eq 0 ]; then
        echo "✗ No sequences found in derep3"
        return 1
    fi

    # Step 2: Export derep1 sequences
    log_message "$gene" "$class" "Extracting full-length sequences from derep1"
    qiime tools export \
        --input-path "$derep1_sequence" \
        --output-path "$output_dir/derep1_export"

    if [[ ! -f "$output_dir/derep1_export/dna-sequences.fasta" ]]; then
        echo "✗ Failed to export derep1 sequences"
        return 1
    fi

    # Step 3: Export taxonomy if available
    local taxonomy_file=""
    if [[ -n "$taxonomy" ]]; then
        log_message "$gene" "$class" "Extracting taxonomy information"
        qiime tools export \
            --input-path "$taxonomy" \
            --output-path "$output_dir/taxonomy_export"

        if [[ -f "$output_dir/taxonomy_export/taxonomy.tsv" ]]; then
            taxonomy_file="$output_dir/taxonomy_export/taxonomy.tsv"
            echo "✓ Taxonomy information exported"
        else
            echo "⚠ Failed to export taxonomy, proceeding without it"
            taxonomy_file=""
        fi
    fi

    # Step 4: Extract primer-binding regions
    log_message "$gene" "$class" "Extracting [primer_binding_region--amplicon--primer_binding_region] regions"
    python3 << EOF
from Bio import SeqIO

def extract_with_primer_binding_regions(full_seq, trimmed_seq, f_primer_length, r_primer_length):
    """Extract regions with primer binding contexts"""
    full_seq_str = str(full_seq).upper()
    trimmed_seq_str = str(trimmed_seq).upper()

    # Find where trimmed sequence aligns in full sequence
    align_pos = full_seq_str.find(trimmed_seq_str)

    if align_pos >= 0:
        # Extract with primer binding regions:
        start_pos = max(0, align_pos - f_primer_length)
        end_pos = min(len(full_seq_str), align_pos + len(trimmed_seq_str) + r_primer_length)

        extracted_region = full_seq_str[start_pos:end_pos]
        
        actual_f_region_length = align_pos - start_pos
        actual_r_region_length = end_pos - (align_pos + len(trimmed_seq_str))

        return extracted_region, {
            'start_pos': start_pos,
            'end_pos': end_pos,
            'align_pos': align_pos,
            'actual_f_region_length': actual_f_region_length,
            'actual_r_region_length': actual_r_region_length,
            'amplicon_length': len(trimmed_seq_str),
            'total_length': end_pos - start_pos
        }

    return None, None

# Read amplifiable sequence IDs
amplifiable_ids = set()
with open("$output_dir/amplifiable_seq_ids.txt", 'r') as f:
    for line in f:
        amplifiable_ids.add(line.strip())

print(f"Looking for {len(amplifiable_ids)} amplifiable sequence IDs...")

# Primer lengths for extraction
f_primer_length = len("$f_primer")
r_primer_length = len("$r_primer")

print(f"Forward primer: $f_primer (length: {f_primer_length} bp)")
print(f"Reverse primer: $r_primer (length: {r_primer_length} bp)")
print(f"Extraction strategy: [real {f_primer_length}bp + amplicon + real {r_primer_length}bp]")

# Read sequences
trimmed_sequences = {}
for record in SeqIO.parse("$output_dir/derep3_export/dna-sequences.fasta", "fasta"):
    if record.id in amplifiable_ids:
        trimmed_sequences[record.id] = record

full_sequences = {}
for record in SeqIO.parse("$output_dir/derep1_export/dna-sequences.fasta", "fasta"):
    if record.id in amplifiable_ids:
        full_sequences[record.id] = record

print(f"Found {len(trimmed_sequences)} trimmed sequences from derep3")
print(f"Found {len(full_sequences)} full sequences from derep1")

# Extract primer binding regions + amplicon
extracted_sequences = []
successful_extractions = 0
extraction_details = []

for seq_id in amplifiable_ids:
    if seq_id in trimmed_sequences and seq_id in full_sequences:
        trimmed_seq = trimmed_sequences[seq_id]
        full_seq = full_sequences[seq_id]

        extracted_region, positions = extract_with_primer_binding_regions(
            full_seq.seq, trimmed_seq.seq, f_primer_length, r_primer_length
        )

        if extracted_region and positions:
            new_record = full_seq[:]
            new_record.seq = type(full_seq.seq)(extracted_region)
            new_record.description = f"{full_seq.description} | [{positions['actual_f_region_length']}bp + {positions['amplicon_length']}bp + {positions['actual_r_region_length']}bp] region {positions['start_pos']}-{positions['end_pos']} | total_length {positions['total_length']}"

            extracted_sequences.append(new_record)
            extraction_details.append({
                'seq_id': seq_id,
                'start_pos': positions['start_pos'],
                'end_pos': positions['end_pos'],
                'align_pos': positions['align_pos'],
                'actual_f_region_length': positions['actual_f_region_length'],
                'actual_r_region_length': positions['actual_r_region_length'],
                'amplicon_length': positions['amplicon_length'],
                'total_length': positions['total_length'],
                'target_f_length': f_primer_length,
                'target_r_length': r_primer_length
            })
            successful_extractions += 1

print(f"\\nExtraction results:")
print(f"✓ Successfully extracted: {successful_extractions} sequences")

if successful_extractions > 0:
    # Save extracted sequences
    with open("$output_dir/extracted_primer_amplicon_primer.fasta", 'w') as output:
        SeqIO.write(extracted_sequences, output, "fasta")

    # Save extraction details
    import pandas as pd
    details_df = pd.DataFrame(extraction_details)
    details_df.to_csv("$output_dir/extraction_details.csv", index=False)

    print(f"\\nSaved [primer_binding_region--amplicon--primer_binding_region] sequences")
else:
    print("\\nERROR: No sequences were successfully extracted!")
    exit(1)

EOF

    # Check if extraction was successful
    if [[ ! -f "$output_dir/extracted_primer_amplicon_primer.fasta" ]]; then
        echo "✗ Failed to extract regions"
        return 1
    fi

    local extracted_count=$(grep -c "^>" "$output_dir/extracted_primer_amplicon_primer.fasta")
    echo "✓ Successfully extracted $extracted_count regions"

    # Step 5: Analyze blocking efficiency
    log_message "$gene" "$class" "Running blocking primer analysis"
    analyze_blocking_efficiency_with_taxonomy "$gene" "$class" "$output_dir/extracted_primer_amplicon_primer.fasta" "$taxonomy_file" "$output_dir"

    echo "✓ Blocking evaluation completed for $gene - $class"
    echo "  → Analyzed $extracted_count regions with taxonomic information"
}

# Function to generate comprehensive summary report
generate_summary_report() {
    local gene=$1
    local classes_array=("${@:2}")
    
    local base_dir="${GENE_CONFIGS[${gene}_base_dir]}"
    local blocking_name="${GENE_CONFIGS[${gene}_blocking_name]}"
    local blocking_eval_dir="$base_dir/primer_blocking_evaluation"

    echo "Generating comprehensive blocking efficiency summary report for $gene..."

    local summary_dir="$blocking_eval_dir/summary_reports"
    mkdir -p "$summary_dir"

    # Create summary CSV file
    local summary_csv="$summary_dir/blocking_efficiency_summary.csv"
    echo "gene,class,blocking_primer,primer_length,total_sequences,excellent_blocking,moderate_blocking,poor_blocking,excellent_percentage,moderate_percentage,poor_percentage,mean_mismatches,median_mismatches" > "$summary_csv"

    # Collect all individual summary files
    for class in "${classes_array[@]}"; do
        local individual_summary="$blocking_eval_dir/${blocking_name}_${class}/${blocking_name}_${class}_summary.csv"
        if [ -f "$individual_summary" ]; then
            tail -n +2 "$individual_summary" >> "$summary_csv"
        fi
    done

    # Generate visualization and final report
    python3 << EOF
import pandas as pd
import numpy as np
import sys

# Check matplotlib availability
try:
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    print("✓ Matplotlib imported successfully")
except ImportError as e:
    print(f"✗ Error importing matplotlib: {e}")
    sys.exit(1)

# Read the summary data
try:
    df = pd.read_csv("$summary_csv")
    print(f"✓ Successfully loaded data for {len(df)} taxonomic classes")
except Exception as e:
    print(f"✗ Error reading summary file: {e}")
    sys.exit(1)

if len(df) == 0:
    print("✗ No data found in summary file")
    sys.exit(1)

# Print overall statistics
print(f"\\n=== BLOCKING PRIMER EVALUATION SUMMARY ===")
print(f"Gene: $gene")
print(f"Blocking primer: {df['blocking_primer'].iloc[0]}")
print(f"Blocking primer length: {df['primer_length'].iloc[0]} bp")
print(f"Taxonomic classes evaluated: {len(df)}")
print(f"Total sequences analyzed: {df['total_sequences'].sum():,}")

# Overall statistics
total_seqs = df['total_sequences'].sum()
total_excellent = df['excellent_blocking'].sum()
total_moderate = df['moderate_blocking'].sum()
total_poor = df['poor_blocking'].sum()

print(f"\\n=== OVERALL BLOCKING EFFICIENCY ===")
print(f"Excellent blocking (0-1 mismatches): {total_excellent:,} ({(total_excellent/total_seqs)*100:.1f}%)")
print(f"Moderate blocking (2-3 mismatches): {total_moderate:,} ({(total_moderate/total_seqs)*100:.1f}%)")
print(f"Poor blocking (4+ mismatches): {total_poor:,} ({(total_poor/total_seqs)*100:.1f}%)")

# Create visualization
df_sorted = df.sort_values('class', ascending=True)
fig, axes = plt.subplots(2, 2, figsize=(16, 12))

blocking_primer = df['blocking_primer'].iloc[0]
primer_length = df['primer_length'].iloc[0]

fig.suptitle(f'$gene Blocking Efficiency Analysis\\nPrimer: {blocking_primer}\\nLength: {primer_length} bp',
             fontsize=14, fontweight='bold', y=0.95)

# Plot 1: Blocking efficiency by class
ax1 = axes[0, 0]
y_pos = np.arange(len(df_sorted))

bars1 = ax1.barh(y_pos, df_sorted['excellent_percentage'],
                 label='Excellent (0-1 mm)', color='#2E8B57', alpha=0.8)
bars2 = ax1.barh(y_pos, df_sorted['moderate_percentage'],
                 left=df_sorted['excellent_percentage'],
                 label='Moderate (2-3 mm)', color='#FF8C00', alpha=0.8)
bars3 = ax1.barh(y_pos, df_sorted['poor_percentage'],
                 left=df_sorted['excellent_percentage'] + df_sorted['moderate_percentage'],
                 label='Poor (4+ mm)', color='#DC143C', alpha=0.8)

ax1.set_xlabel('Percentage of Sequences')
ax1.set_ylabel('Taxonomic Class')
ax1.set_title('Blocking Efficiency by Taxonomic Class')
ax1.set_yticks(y_pos)
ax1.set_yticklabels(df_sorted['class'])
ax1.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
ax1.set_xlim(0, 100)

# Plot 2: Number of sequences by class
ax2 = axes[0, 1]
bars = ax2.bar(range(len(df_sorted)), df_sorted['total_sequences'],
               color='skyblue', alpha=0.8, edgecolor='navy', linewidth=0.5)
ax2.set_xlabel('Taxonomic Class')
ax2.set_ylabel('Number of Sequences')
ax2.set_title('Number of Sequences Analyzed by Class')
ax2.set_xticks(range(len(df_sorted)))
ax2.set_xticklabels(df_sorted['class'], rotation=45, ha='right')
ax2.grid(axis='y', alpha=0.3)

for i, bar in enumerate(bars):
    height = bar.get_height()
    ax2.text(bar.get_x() + bar.get_width()/2., height + height*0.01,
             f'{int(height):,}', ha='center', va='bottom', fontsize=9, fontweight='bold')

# Plot 3: Mean and Median mismatches by class
ax3 = axes[1, 0]
x_pos = np.arange(len(df_sorted))
width = 0.35

bars1 = ax3.bar(x_pos - width/2, df_sorted['mean_mismatches'], width,
                label='Mean mismatches', color='#8B4B8B', alpha=0.8)
bars2 = ax3.bar(x_pos + width/2, df_sorted['median_mismatches'], width,
                label='Median mismatches', color='#B8860B', alpha=0.8)

ax3.set_xlabel('Taxonomic Class')
ax3.set_ylabel('Number of Mismatches')
ax3.set_title('Mean and Median Mismatches by Class')
ax3.set_xticks(x_pos)
ax3.set_xticklabels(df_sorted['class'], rotation=45, ha='right')
ax3.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
ax3.grid(axis='y', alpha=0.3)

# Plot 4: Overall pie chart
ax4 = axes[1, 1]
sizes = [total_excellent, total_moderate, total_poor]
colors = ['#2E8B57', '#FF8C00', '#DC143C']
labels = ['Excellent (0-1 mm)', 'Moderate (2-3 mm)', 'Poor (4+ mm)']

non_zero_sizes = []
non_zero_colors = []
non_zero_labels = []
for size, color, label in zip(sizes, colors, labels):
    if size > 0:
        non_zero_sizes.append(size)
        non_zero_colors.append(color)
        non_zero_labels.append(label)

if non_zero_sizes:
    wedges, texts, autotexts = ax4.pie(non_zero_sizes, colors=non_zero_colors,
                                       autopct='%1.1f%%', startangle=90,
                                       textprops={'fontsize': 10, 'fontweight': 'bold'})

    legend_labels = []
    for size, label in zip(non_zero_sizes, non_zero_labels):
        legend_labels.append(f'{label}: {size:,} sequences')

    ax4.legend(wedges, legend_labels, title="Blocking Categories",
               loc="center left", bbox_to_anchor=(1, 0, 0.5, 1), fontsize=9)

    ax4.set_title(f'Overall Blocking Efficiency\\n({total_seqs:,} total sequences)',
                  fontsize=12, fontweight='bold')

plt.tight_layout()
plt.subplots_adjust(top=0.88, right=0.85)

# Save visualization
output_file = "$summary_dir/blocking_efficiency_analysis.png"
try:
    plt.savefig(output_file, dpi=300, bbox_inches='tight',
                facecolor='white', edgecolor='none')
    print(f"✓ Visualization saved to: {output_file}")
except Exception as e:
    print(f"✗ Error saving visualization: {e}")

plt.close()

# Print detailed results
print(f"\\n=== DETAILED RESULTS BY CLASS ===")
display_df = df_sorted[['class', 'total_sequences', 'excellent_percentage', 'moderate_percentage', 'poor_percentage', 'mean_mismatches', 'median_mismatches']]
print(display_df.to_string(index=False, float_format='%.1f'))

# Save sorted results
df_sorted.to_csv("$summary_dir/blocking_efficiency_detailed.csv", index=False)

print(f"\\n✓ Analysis complete! Results saved to: $summary_dir")

EOF

    echo "Comprehensive summary report generated in $summary_dir"
}

# Main execution function
main() {
    # Check for help or no arguments
    if [ $# -eq 0 ] || [ "$1" = "-h" ] || [ "$1" = "--help" ]; then
        show_usage
        exit 0
    fi

    # Get gene type (required first argument)
    local gene=$(echo "$1" | tr '[:lower:]' '[:upper:]')
    if [[ ! -n "${GENE_CONFIGS[${gene}_blocking_primer]}" ]]; then
        echo "Error: Invalid gene type '$1'. Must be one of: 12S, 16S, COI"
        show_usage
        exit 1
    fi

    # Parse additional arguments
    local classes_to_process=("${CLASSES[@]}")
    local custom_blocking_primer=""
    local custom_blocking_name=""

    shift # Remove gene argument
    while [ $# -gt 0 ]; do
        case $1 in
            --blocking-primer)
                if [ $# -lt 2 ]; then
                    echo "Error: --blocking-primer requires a sequence"
                    exit 1
                fi
                custom_blocking_primer="$2"
                shift 2
                ;;
            --blocking-name)
                if [ $# -lt 2 ]; then
                    echo "Error: --blocking-name requires a name"
                    exit 1
                fi
                custom_blocking_name="$2"
                shift 2
                ;;
            --classes)
                if [ $# -lt 2 ]; then
                    echo "Error: --classes requires a comma-separated list"
                    exit 1
                fi
                IFS=',' read -ra classes_to_process <<< "$2"
                shift 2
                ;;
            *)
                echo "Error: Unknown option $1"
                show_usage
                exit 1
                ;;
        esac
    done

    # Apply custom blocking primer if provided
    if [ -n "$custom_blocking_primer" ]; then
        GENE_CONFIGS["${gene}_blocking_primer"]="$custom_blocking_primer"
        echo "Using custom blocking primer: $custom_blocking_primer"
    fi

    if [ -n "$custom_blocking_name" ]; then
        GENE_CONFIGS["${gene}_blocking_name"]="$custom_blocking_name"
        echo "Using custom blocking name: $custom_blocking_name"
    fi

    echo "=== UNIFIED BLOCKING PRIMER EVALUATION ==="
    echo "Gene: $gene"
    echo "Blocking primer: ${GENE_CONFIGS[${gene}_blocking_primer]}"
    echo "Blocking name: ${GENE_CONFIGS[${gene}_blocking_name]}"
    echo "PCR primers: ${GENE_CONFIGS[${gene}_f_primer]} / ${GENE_CONFIGS[${gene}_r_primer]}"
    echo "Base directory: ${GENE_CONFIGS[${gene}_base_dir]}"
    echo "Classes to process: ${classes_to_process[*]}"
    echo ""

    # Validate base directory
    local base_dir="${GENE_CONFIGS[${gene}_base_dir]}"
    if [ ! -d "$base_dir" ]; then
        echo "Error: Base directory not found: $base_dir"
        exit 1
    fi

    # Check dependencies
    if ! command -v qiime &> /dev/null; then
        echo "Error: QIIME2 not found"
        exit 1
    fi

    if ! python3 -c "import pandas, numpy, Bio" 2>/dev/null; then
        echo "Error: Required Python packages not found (pandas, numpy, biopython)"
        exit 1
    fi

    # Process each class
    local successful_classes=()
    for class in "${classes_to_process[@]}"; do
        echo "Processing: $class"
        echo "------------------------"

        if evaluate_blocking_primer_for_class "$gene" "$class"; then
            echo "✓ Completed: $class"
            successful_classes+=("$class")
        else
            echo "✗ Failed: $class"
        fi
        echo ""
    done

    # Generate comprehensive summary report if any classes succeeded
    if [ ${#successful_classes[@]} -gt 0 ]; then
        echo ""
        echo "Generating comprehensive summary report..."
        generate_summary_report "$gene" "${successful_classes[@]}"

        echo ""
        echo "=== SUMMARY ==="
        echo "Gene: $gene"
        echo "Successfully processed classes: ${#successful_classes[@]}"
        echo "Results saved in: ${GENE_CONFIGS[${gene}_base_dir]}/primer_blocking_evaluation"
        echo ""
        echo "Key output files:"
        echo "- blocking_efficiency_analysis.png (main visualization)"
        echo "- blocking_efficiency_summary.csv (overall data by class)"
        echo "- blocking_efficiency_detailed.csv (classes sorted by performance)"
        echo ""
        echo "Individual class directories contain:"
        echo "- extracted_primer_amplicon_primer.fasta (extracted regions for analysis)"
        echo "- *_detailed_results_with_taxonomy.csv (detailed blocking results)"
        echo "- *_genus_summary.csv (blocking efficiency by genus)"
        echo "- *_family_summary.csv (blocking efficiency by family)"
    else
        echo "✗ No classes were successfully processed"
        exit 1
    fi
}

# Run the script
main "$@"