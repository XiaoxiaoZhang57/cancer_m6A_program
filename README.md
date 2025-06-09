# cancer_m6A_program

#1. Download cancer sites data

curl -H "Authorization: Basic XXX" "https://cancer.sanger.ac.uk/api/mono/products/v1/downloads/scripted?path=grch38/cosmic/v99/Cosmic_GenomeScreensMutant_Tsv_v99_GRCh38.tar&bucket=downloads"

curl -H "Authorization: Basic XXX" "https://cancer.sanger.ac.uk/api/mono/products/v1/downloads/scripted?path=grch38/cosmic/v99/Cosmic_Sample_Tsv_v99_GRCh38.tar&bucket=downloads"


#2. Download of 62 species genome alignment

http://hgdownload.soe.ucsc.edu/goldenPath/hg38/multiz100way/


#3. Calculate the value of rate4site

rate4site -s GeneX.aa.aln -t hg38.100way.nh >out.res

#4. Check the positions of cancer mutation sites in the genome

python seq_select.py -i Homo_sapiens.GRCh38.dna_sm.toplevel.fa -l coding.mutation.aa.sort.tsv.pylist -o gene.fasta --detail


#5. Correspond the rate4site sites in humans with cancer sites

