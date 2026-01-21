# Evolutionary Analysis of Cancer-Associated Loci in Mammals

This repository contains the computational pipeline, statistical scripts, and visualization tools used to analyze the evolutionary conservation of cancer-associated mutation sites across mammalian species.

---

## 1. Data Acquisition

### Cancer Mutation Data
* COSMIC Database: Cancer mutation data (v99 GRCh38) were downloaded using the following scripted API:
  - Download Genome Screens Mutant TSV: Use curl with the COSMIC V99 endpoint.
    curl -H "Authorization: Basic XXX" "https://cancer.sanger.ac.uk/api/mono/products/v1/downloads/scripted?path=grch38/cosmic/v99/Cosmic_GenomeScreensMutant_Tsv_v99_GRCh38.tar&bucket=downloads"
  - Download Sample TSV: Use curl with the COSMIC V99 Sample endpoint.
    curl -H "Authorization: Basic XXX" "https://cancer.sanger.ac.uk/api/mono/products/v1/downloads/scripted?path=grch38/cosmic/v99/Cosmic_Sample_Tsv_v99_GRCh38.tar&bucket=downloads"

* cBioPortal: Additional cancer site data were retrieved from cBioPortal.
  - get_mutation_from_cbioportal.R: Script to automate data extraction from the cBioPortal database.

### Genomic Resources
* Multiple Sequence Alignment (MSA): 62 mammalian species alignments were obtained from the UCSC Genome Browser (hg38 100-way).
* Reference Genome: Homo sapiens (GRCh38/hg38).
  http://hgdownload.soe.ucsc.edu/goldenPath/hg38/multiz100way/

---

## 2. Computational Pipeline & Quality Control

### Step 1: Sequence Selection and Alignment Processing
To ensure high-quality cross-species comparison, we implemented a human-centric filtering strategy:
* Gap Filtering: Sites are excluded if the human residue is missing or if >70% of the other species exhibit gaps at that position.
* seq_select.py: Extracts genomic positions and amino acid sequences for specific loci.
  Command: python seq_select.py -i Homo_sapiens.GRCh38.dna_sm.toplevel.fa -l coding.mutation.aa.sort.tsv.pylist -o gene.fasta --detail

### Step 2: Evolutionary Rate Estimation
We employed both phylogeny-aware and model-independent methods to calculate site-specific conservation:
* Rate4Site: Normalizes evolutionary rates by considering phylogenetic branch lengths.
  Command: rate4site -s GeneX.aa.aln -t hg38.100way.nh > out.res
* RAxML: Used for maximum likelihood analysis to calculate site-specific selection pressure.
  Command: raxmlHPC -fe -t all58.nwk -m PROTGAMMAJTTX -s GeneX.raxml58.fa -n RAxML58
* free-model: 
  A model-independent approach used to quantify site-specific conservation without assuming a predefined amino acid substitution matrix. This script calculates conservation scores based on the distribution of amino acids and physicochemical properties directly from the Multiple Sequence Alignment (MSA).
  Command: python conservation_score_free-model.py

### Step 3: Functional Enrichment
* polysel.R: Performs PolySel enrichment analysis to detect polygenic selection signals in specific pathways.
  
### Step 4: Branch Model Analysis (PAML)
To compare selection pressures between long-lived and short-lived mammalian species:
* Code: `paml codeml` using the branch model (model = 2, NSsites = 0).
* Objective: Identify genes where cancer pathogenic sites show significantly different selection pressures (dN/dS) between longevity groups.
* Statistical Criteria: Significance is assessed via 95% Confidence Intervals (Point Estimate ± 1.96 * SE). If the interval excludes zero, the null hypothesis is rejected.

### Step 5: Site-Specific Evolutionary Rates
To evaluate the functional importance of human cancer-associated sites, two complementary likelihood-based approaches were used:
1. LEISR (HyPhy): A site-rate estimator to calculate purifying selection.
   Command: `hyphy leisr --alignment GeneX.msa --tree species.nwk`
2. Site Models (PAML): To test for selection constraints on specific loci.
   Command: `paml codeml` (using M1a/M2a or M7/M8 models).

### Step 6: Sensitivity & Robustness Checks
* Gap Filtering: Sites with excessive gaps (>70%) or low alignment confidence were excluded.
* Validation: Mapped positions were cross-referenced with cBioPortal data.
* Result: 86.8% overlap with cBioPortal and 69.5% retention from COSMIC after rigorous filtering.

### Step 7: Intersection Analysis of Epigenetic Modifications and Cancer Sites

This section describes the integration of DNA 5mC and RNA m6A modification data with human-centric cancer mutation sites across 21 mammalian species.

#### 7.1 Homology Mapping of Cancer Sites
* **Objective**: To identify homologous positions of human cancer-associated sites in 20 other mammalian species (21 species total).
* **Tool**: **Clustal Omega** was utilized to perform multi-species sequence alignment and precise coordinate mapping.

#### 7.2 Coordinate Transformation (BED Mapping)
To intersect multi-species modification data (5mC/m6A) with human cancer sites, non-human genomic coordinates were converted into human-equivalent coordinates (GRCh38/hg38).
* **Script**: `map_species_to_human_bed.py` (formerly `A_to_human_bed.py`)
* **Function**: This script transforms species-specific BED files into human-reference BED format based on orthologous mapping, facilitating cross-species genomic comparisons.
* **Usage**:
  ```bash
  python map_species_to_human_bed.py --input species_specific.bed --map homology_table.txt --output human_coords.bed
---

## 3. Statistical Analysis and Visualization

The following list maps the scripts to the corresponding figures in the manuscript:

* Figure 1: plot_tree.r - Visualizes the mammalian phylogenetic tree.
* Figure 1: plot_age_of_cancer_sample.r - Shows the age distribution of cancer samples.
* Figure 1: plot_Transition_transversion_ratio.r - Calculates and plots the Ti/Tv ratio.
* Figure 1: plot_Rratio.py - Distribution of R-ratio scores.
* Figure all: box.r - Boxplots comparing conservation scores across groups.
* Figure 3: plot_polysel.r - Visualizes PolySel enrichment results.
* Figure 1-4: The source data for the boxplots in Figures 1–4 are all stored in the data folder.
* Statistics: PGLS_cancerRisk.R - Phylogenetic Generalized Least Squares (PGLS) analysis for cancer risk vs. lifespan.
* Supplementary: tissue_site_distribution.r - Pie chart for tissue-specific cancer site distribution.

---

## 4. Environment Requirements

* R: >= 4.0 (Dependencies: ggplot2, ape, caper)
* Python: >= 3.7 (Dependencies: biopython, pandas)
* Bioinformatics Tools: Rate4Site, RAxML-HPC

---

## Contact
For any queries regarding the scripts or data, please contact Xiaoxiao Zhang at zhangxiaoxiao@ioz.ac.cn.

