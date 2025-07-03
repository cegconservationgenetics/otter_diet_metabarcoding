# otter_diet_metabarcoding
DNA Metabarcoding-based diet diversity and niche differences of otters in Thailand

🔹 1. Retrieve Reference Sequences from NCBI Using QIIME 2

This step downloads curated reference sequences for common mitochondrial genes (e.g., 12S, 16S, COI) across major taxonomic groups relevant to metabarcoding studies.
🧬 Target Genes

The gene-specific query strings below target commonly used mitochondrial markers, using multiple aliases to ensure broad capture of relevant sequences from NCBI:

gene_queries["12S"] = '(12S[Title] OR 12s[Title] OR "s-rRNA"[Title] OR rrn12[Title] OR rrnS[Title] OR rrns[Title] OR srRNA[Title] OR "small subunit ribosomal RNA"[Title] OR "small ribosomal RNA subunit"[Title] OR mtDNA[Title] OR mitochondrial[Title] OR mitochondrion[Title] OR "mt-rRNA"[Title] OR "mt-rrnS"[Title])'

gene_queries["16S"] = '(16S[Title] OR 16s[Title] OR rrn16[Title] OR rrnA[Title] OR rrnB[Title] OR lrRNA[Title] OR "l-rRNA"[Title] OR mtDNA[Title] OR mitochondrial[Title] OR mitochondrion[Title] OR "mt-rRNA"[Title] OR "mt-rrnL"[Title])'

gene_queries["COI"] = '(COI[Title] OR CO1[Title] OR "cytochrome oxidase I"[Title] OR "cytochrome c oxidase subunit I"[Title] OR "cytochrome c oxidase I"[Title] OR cox1[Title] OR coxi[Title] OR "COX subunit 1"[Title] OR mtDNA[Title] OR mitochondrial[Title] OR mitochondrion[Title] OR "mt-co1"[Title] OR "mt-coI"[Title] OR "mt-rRNA"[Title])'

🧫 Target Taxonomic Groups

We focus on diverse eukaryotic groups of ecological interest:
