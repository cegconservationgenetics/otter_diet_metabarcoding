#!/usr/bin/env python3

import matplotlib
matplotlib.use('Agg')  # Use non-interactive backend
import matplotlib.pyplot as plt
from Bio import SeqIO
import os
import subprocess
import glob
import traceback
import sys
import math

# Gene-specific configurations
GENE_CONFIGS = {
    "12s": {
        "input_dir": "/mnt/d/BI/crabs/rawdata/12s",
        "fasta_output_dir": "/mnt/d/BI/crabs/12s",
        "plots_output_dir": "/mnt/d/BI/crabs/12s/crabs_12s_output",
        "primer_pairs": {
            "12S_tele02": ("AAACTCGTGCCAGCCACC", "GGGTATCTAATCCCAGTTTG"),
            "12S_teleo": ("ACACCGCCCGTCACTCT", "CTTCCGGTACACTTACCRTG"),
            "12S_MiFish-U": ("GTCGGTAAAACTCGTGCCAGC", "CATAGTGGGGTATCTAATCCCAGTTTG"),
            "12S_12SV5": ("TTAGATACCCCACTATGC", "TAGAACAGGCTCCTCTAG")
        }
    },
    "16s": {
        "input_dir": "/mnt/d/BI/crabs/rawdata/16s",
        "fasta_output_dir": "/mnt/d/BI/crabs/16s",
        "plots_output_dir": "/mnt/d/BI/crabs/16s/crabs_16s_output",
        "primer_pairs": {
            "16S_16SMAV": ("CCAACATCGAGGTCRYAA", "ARTTACYNTAGGGATAACAG"),
            "16S_Vert16SeDNA": ("AGACGAGAAGACCCYDTGGAGCTT", "GATCCAACATCGAGGTCGTAA"),
            "16S_MarVer3": ("AGACGAGAAGACCCTRTG", "GGATTGCGCTGTTATCCC")
        }
    },
    "coi": {
        "input_dir": "/mnt/d/BI/crabs/rawdata/coi",
        "fasta_output_dir": "/mnt/d/BI/crabs/coi",
        "plots_output_dir": "/mnt/d/BI/crabs/coi/crabs_coi_output",
        "primer_pairs": {
            "COI_VF2FishR1": ("TCAACCAACCACAAAGACATTGGCAC", "TAGACTTCTGGGTGGCCAAAGAATCA"),
            "COI_coi175f345r": ("GGAGGCTTTGGMAAYTGRYT", "TAGAGGRGGGTARACWGTYCA"),
            "COI_L2513H2714": ("GCCTGTTTACCAAAAACATCA", "CTCCATAGGGTCTTCTCGTCTT")
        }
    }
}

# Default mismatch percentage
mismatch_percentage = 10.0

def show_usage():
    """Display usage information"""
    print("Usage: python3 crabs_insilico_pcr.py <gene> [options]")
    print("")
    print("Required argument:")
    print("  gene              Gene type: 12s, 16s, or coi")
    print("")
    print("Options:")
    print("  --primer <name>   Use specific primer (default: all primers)")
    print("  --mismatch <n>    Mismatch percentage (default: 10.0)")
    print("  --pcr-only        Run only in silico PCR (skip plotting)")
    print("  --plot-only       Generate plots only (skip PCR)")
    print("  --no-clean        Keep empty FASTA files")
    print("  -h, --help        Show this help message")
    print("")
    print("Available genes and primers:")
    for gene, config in GENE_CONFIGS.items():
        print(f"  {gene.upper()}: {', '.join(config['primer_pairs'].keys())}")
    print("")
    print("Examples:")
    print("  python3 crabs_insilico_pcr.py 12s")
    print("  python3 crabs_insilico_pcr.py 16s --primer 16S_MarVer3")
    print("  python3 crabs_insilico_pcr.py coi --mismatch 15.0 --plot-only")

def calculate_mismatches(forward, reverse, percentage):
    """Calculate allowed mismatches based on primer length and percentage"""
    # Calculate average primer length
    avg_length = (len(forward) + len(reverse)) / 2
    
    # Calculate mismatches as percentage of average length
    mismatches = math.ceil(avg_length * percentage / 100)
    
    return mismatches

def check_and_remove_empty_fastas(fasta_output_dir):
    """Check for empty FASTA files and remove them"""
    print("\n=== Checking for Empty FASTA Files ===")
    
    # Find all FASTA files in the output directory
    fasta_files = glob.glob(os.path.join(fasta_output_dir, "*_mismatch*.fasta"))
    
    if not fasta_files:
        print("No FASTA files found to check.")
        return
    
    print(f"Found {len(fasta_files)} FASTA files to check.")
    
    empty_files = []
    for fasta_file in fasta_files:
        try:
            # Try to count sequences
            sequence_count = sum(1 for _ in SeqIO.parse(fasta_file, "fasta"))
            
            if sequence_count == 0:
                print(f"Empty FASTA file found: {os.path.basename(fasta_file)}")
                empty_files.append(fasta_file)
        except Exception as e:
            print(f"Error checking {fasta_file}: {e}")
    
    # Remove empty files
    if empty_files:
        print(f"\nRemoving {len(empty_files)} empty FASTA files:")
        for empty_file in empty_files:
            try:
                os.remove(empty_file)
                print(f"Removed: {os.path.basename(empty_file)}")
            except Exception as e:
                print(f"Error removing {empty_file}: {e}")
    else:
        print("No empty FASTA files found.")
    
    print(f"FASTA file check complete. Removed {len(empty_files)} empty files.")

def insilico_pcr_individual_files(gene, primer_name=None):
    """Run in silico PCR on individual class files"""
    print(f"\n=== Starting In Silico PCR for {gene.upper()} ===")
    
    config = GENE_CONFIGS[gene]
    input_dir = config["input_dir"]
    fasta_output_dir = config["fasta_output_dir"]
    
    # Create output directory if it doesn't exist
    os.makedirs(fasta_output_dir, exist_ok=True)
    
    # Determine which primers to use
    primers_to_use = {}
    if primer_name and primer_name in config["primer_pairs"]:
        primers_to_use[primer_name] = config["primer_pairs"][primer_name]
    else:
        primers_to_use = config["primer_pairs"]
    
    if not primers_to_use:
        print(f"Error: Primer '{primer_name}' not found in {gene.upper()} primer pairs.")
        return []
        
    print(f"Using {len(primers_to_use)} primer pairs: {', '.join(primers_to_use.keys())}")
    print(f"Mismatch tolerance: {mismatch_percentage}%")
    
    # Get all FASTA files in the input directory
    fasta_files = glob.glob(os.path.join(input_dir, "*_sequences.fasta"))
    processed_files = []
    
    if not fasta_files:
        print(f"No FASTA files found in {input_dir}")
        return []
        
    print(f"Found {len(fasta_files)} FASTA files to process")
    
    # Process each file individually with each primer pair
    for input_fasta in fasta_files:
        class_name = os.path.basename(input_fasta).replace("_sequences.fasta", "")
        
        for primer_name, (forward, reverse) in primers_to_use.items():
            # Calculate mismatches based on primer length and percentage
            mismatches = calculate_mismatches(forward, reverse, mismatch_percentage)
            
            # New file naming pattern: primer_class_mismatchN
            output_fasta = os.path.join(fasta_output_dir, f"{primer_name}_{class_name}_mismatch{mismatches}.fasta")
            
            print(f"\nProcessing {class_name} with primer {primer_name}...")
            print(f"Forward primer: {forward} ({len(forward)} bp)")
            print(f"Reverse primer: {reverse} ({len(reverse)} bp)")
            print(f"Allowing {mismatches} mismatches ({mismatch_percentage}% of average primer length)")
            
            command = [
                "crabs", "insilico_pcr",
                "-i", input_fasta,
                "-o", output_fasta,
                "-f", forward,
                "-r", reverse,
                "-e", str(mismatches)
            ]
            
            try:
                subprocess.run(command, check=True)
                
                # Verify the output file exists
                if os.path.exists(output_fasta):
                    # Check if the output FASTA has any sequences
                    sequence_count = sum(1 for _ in SeqIO.parse(output_fasta, "fasta"))
                    
                    if sequence_count > 0:
                        processed_files.append((class_name, primer_name, output_fasta, mismatches))
                        print(f"Successfully processed {class_name} with primer {primer_name} - Found {sequence_count} sequences")
                    else:
                        print(f"No sequences found for {class_name} with primer {primer_name}")
                else:
                    print(f"Warning: Output file {output_fasta} was not created")
                    
            except subprocess.CalledProcessError as e:
                print(f"Error processing {class_name} with primer {primer_name}: {e}")
            except Exception as e:
                print(f"Unexpected error processing {class_name} with primer {primer_name}: {e}")
                print(traceback.format_exc())
    
    print(f"\nIn silico PCR completed. Successfully processed {len(processed_files)} file-primer combinations.")
    return processed_files

def plot_sequences(gene, processed_files=None, primer_name=None):
    """Generate plots for sequence lengths"""
    print(f"\n=== Starting Plot Generation for {gene.upper()} ===")
    
    config = GENE_CONFIGS[gene]
    fasta_output_dir = config["fasta_output_dir"]
    plots_output_dir = config["plots_output_dir"]
    
    # Create output directory if it doesn't exist
    os.makedirs(plots_output_dir, exist_ok=True)
    
    # If no processed files were passed, find files in the FASTA output directory
    if not processed_files:
        print("No processed files provided, searching FASTA output directory...")
        
        # Find all files matching the pattern [primer]_[class]_mismatch*.fasta
        if primer_name:
            fasta_pattern = f"{primer_name}_*_mismatch*.fasta"
        else:
            fasta_pattern = "*_*_mismatch*.fasta"
        
        fasta_files = glob.glob(os.path.join(fasta_output_dir, fasta_pattern))
        processed_files = []
        
        for file_path in fasta_files:
            basename = os.path.basename(file_path)
            name_parts = basename.replace(".fasta", "").split("_mismatch")
            if len(name_parts) == 2:
                name_parts2 = name_parts[0].split("_", 1)  # Split at first underscore
                if len(name_parts2) == 2:
                    primer = name_parts2[0]
                    class_name = name_parts2[1]
                    mismatches = int(name_parts[1])
                    
                    # Only process files with matching primer if specified
                    if primer_name and primer != primer_name:
                        continue
                    
                    # Check if the file has sequences
                    sequence_count = sum(1 for _ in SeqIO.parse(file_path, "fasta"))
                    if sequence_count > 0:
                        processed_files.append((class_name, primer, file_path, mismatches))
    
    if not processed_files:
        print(f"No valid FASTA files found in {fasta_output_dir}")
        return
        
    print(f"Generating plots for {len(processed_files)} FASTA files")
    
    plots_created = 0
    
    # Group files by primer for consistent plots
    by_primer = {}
    for class_name, primer, file_path, mismatches in processed_files:
        if primer not in by_primer:
            by_primer[primer] = []
        by_primer[primer].append((class_name, file_path, mismatches))
    
    # Process each primer group
    for primer, files in by_primer.items():
        print(f"\nGenerating plots for primer {primer}...")
        
        for class_name, file_path, mismatches in files:
            try:
                print(f"Processing plot for {class_name} with primer {primer}...")
                
                # Check if the file exists and is readable
                if not os.path.exists(file_path):
                    print(f"Error: File {file_path} does not exist. Skipping.")
                    continue
                    
                # Extract sequence lengths for this class
                lengths = []
                
                for record in SeqIO.parse(file_path, 'fasta'):
                    lengths.append(len(record.seq))
                
                print(f"Found {len(lengths)} sequences for {class_name}")
                
                # Skip if no sequences
                if not lengths:
                    print(f"No sequences found for {class_name}, skipping plot.")
                    continue
                
                # Create plot
                plt.figure(figsize=(15, 8))
                seq_nums = range(1, len(lengths) + 1)
                lengths_1000 = [(l if l < 1000 else 1000) for l in lengths]
                colors = ['blue' if l < 1000 else 'red' for l in lengths]
                
                plt.scatter(lengths_1000, seq_nums, c=colors, alpha=0.6)
                plt.title(f'{primer}: Dot Plot of Sequence Lengths - {class_name} (Mismatch: {mismatches})')
                plt.xlabel('Sequence Length (bp)')
                plt.ylabel('Sequence Number')
                plt.xticks(range(0, 1001, 50))
                plt.xlim(0, 1000)
                plt.grid(True, linestyle='--', alpha=0.7)
                
                # Add statistics
                stats = (f'Total sequences: {len(lengths)}\n'
                         f'≤1000bp: {sum(1 for l in lengths if l <= 1000)}\n'
                         f'>1000bp: {sum(1 for l in lengths if l > 1000)}\n'
                         f'Minimum length: {min(lengths)}\n'
                         f'Maximum length: {max(lengths)}\n'
                         f'Mean length: {sum(lengths)/len(lengths):.1f}')
                
                plt.text(1.02, 0.95, stats, transform=plt.gca().transAxes,
                         bbox=dict(facecolor='white', edgecolor='black'))
                plt.scatter([], [], c='blue', label='≤1000bp')
                plt.scatter([], [], c='red', label='>1000bp')
                plt.legend()
                plt.tight_layout()
                
                # Save the plot - using the same naming pattern as FASTA files
                plot_path = os.path.join(plots_output_dir, f'{primer}_{class_name}_mismatch{mismatches}_dot_plot.png')
                plt.savefig(plot_path, dpi=300)
                plt.close()
                
                # Verify the plot was created
                if os.path.exists(plot_path):
                    print(f"Plot saved to {plot_path}")
                    plots_created += 1
                else:
                    print(f"Warning: Plot file {plot_path} was not created")
                
            except Exception as e:
                print(f"Error generating plot for {class_name}:")
                print(traceback.format_exc())
    
    print(f"\nPlot generation completed. Created {plots_created} plots.")

def main():
    """Main execution function"""
    global mismatch_percentage
    
    try:
        # Check for help or no arguments
        if len(sys.argv) < 2 or sys.argv[1] in ['-h', '--help']:
            show_usage()
            return 0
        
        # Get gene type (required first argument)
        gene = sys.argv[1].lower()
        if gene not in GENE_CONFIGS:
            print(f"Error: Invalid gene type '{gene}'. Must be one of: {', '.join(GENE_CONFIGS.keys())}")
            show_usage()
            return 1
        
        # Parse command-line arguments
        run_pcr = True
        run_plot = True
        primer_name = None
        clean_empty = True
        
        # Check for specific flags
        i = 2
        while i < len(sys.argv):
            if sys.argv[i] == '--plot-only':
                run_pcr = False
                print("Running in plot-only mode")
            elif sys.argv[i] == '--pcr-only':
                run_plot = False
                print("Running in PCR-only mode")
            elif sys.argv[i] == '--primer':
                if i + 1 < len(sys.argv):
                    primer_name = sys.argv[i + 1]
                    print(f"Using primer: {primer_name}")
                    i += 1
                else:
                    print("Error: --primer flag requires a primer name")
                    return 1
            elif sys.argv[i] == '--mismatch':
                if i + 1 < len(sys.argv):
                    try:
                        mismatch_percentage = float(sys.argv[i + 1])
                        print(f"Using mismatch percentage: {mismatch_percentage}%")
                        i += 1
                    except ValueError:
                        print(f"Error: Invalid mismatch percentage: {sys.argv[i + 1]}")
                        return 1
                else:
                    print("Error: --mismatch flag requires a percentage value")
                    return 1
            elif sys.argv[i] == '--no-clean':
                clean_empty = False
                print("Keeping empty FASTA files")
            i += 1
        
        processed_files = []
        
        # Run in silico PCR if requested
        if run_pcr:
            processed_files = insilico_pcr_individual_files(gene, primer_name)
        
        # Clean empty FASTA files if requested
        if clean_empty:
            check_and_remove_empty_fastas(GENE_CONFIGS[gene]["fasta_output_dir"])
        
        # Run plot generation if requested
        if run_plot:
            plot_sequences(gene, processed_files if processed_files else None, primer_name)
        
        print(f"\n=== Processing complete for {gene.upper()}! ===")
        print(f"FASTA files have been saved to: {GENE_CONFIGS[gene]['fasta_output_dir']}")
        print(f"Plot files have been saved to: {GENE_CONFIGS[gene]['plots_output_dir']}")
        
    except Exception as e:
        print("\n=== An error occurred during execution: ===")
        print(traceback.format_exc())
        return 1
    
    return 0

if __name__ == "__main__":
    sys.exit(main())