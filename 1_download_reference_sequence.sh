#!/bin/bash

# Define gene regions search terms
declare -A gene_queries
gene_queries["12S"]='(12S[Title] OR 12s[Title] OR "s-rRNA"[Title] OR rrn12[Title] OR rrnS[Title] OR rrns[Title] OR srRNA[Title] OR "small subunit ribosomal RNA"[Title] OR "small ribosomal RNA subunit"[Title] OR mtDNA[Title] OR mitochondrial[Title] OR mitochondrion[Title] OR "mt-rRNA"[Title] OR "mt-rrnS"[Title])'
gene_queries["16S"]='(16S[Title] OR 16s[Title] OR rrn16[Title] OR rrnA[Title] OR rrnB[Title] OR lrRNA[Title] OR "l-rRNA"[Title] OR mtDNA[Title] OR mitochondrial[Title] OR mitochondrion[Title] OR "mt-rRNA"[Title] OR "mt-rrnL"[Title])'
gene_queries["COI"]='(COI[Title] OR CO1[Title] OR "cytochrome oxidase I"[Title] OR "cytochrome c oxidase subunit I"[Title] OR "cytochrome c oxidase I"[Title] OR cox1[Title] OR coxi[Title] OR "COX subunit 1"[Title] OR mtDNA[Title] OR mitochondrial[Title] OR mitochondrion[Title] OR "mt-co1"[Title] OR "mt-coI"[Title] OR "mt-rRNA"[Title])'

# Define taxonomic groups
groups=("Actinopterygii" "Amphibia" "Mammalia" "Aves" "Crocodylia" "Lepidosauria" "Testudines" "Insecta" "Bivalvia" "Gastropoda" "Clitellata" "Malacostraca" "Insecta")

# Process each gene region in sequence
for gene in "12S" "16S" "COI"; do
    gene_query="${gene_queries[$gene]}"

    # Process each taxonomic group
    for group in "${groups[@]}"; do
        echo "Processing $gene for $group..."

        qiime rescript get-ncbi-data \
            --p-query "$group[Organism] AND 80[SLEN]:30000[SLEN] AND ($gene_query) NOT \"environmental sample\"[Title] NOT \"environmental samples\"[Title] NOT environmental[Title] NOT uncultured[Title] NOT unclassified[Title] NOT unidentified[Title] NOT unverified[Title] NOT PREDICTED[Title]" \
            --p-ranks kingdom phylum class order family genus species \
            --p-rank-propagation \
            --p-n-jobs 5 \
            --o-sequences "${gene}_12smetazoa_sequence.qza" \
            --o-taxonomy "${gene}_12smetazoa_taxonomy.qza" \
            --verbose > "${gene}_12smetazoa_qiime_download.log" 2>&1

        # Add a small delay to avoid overwhelming the NCBI servers
        sleep 10
    done
done
