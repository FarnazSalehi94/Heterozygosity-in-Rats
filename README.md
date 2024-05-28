# Heterozygosity-in-Rats


#First thing is writing script to aligning primary, alternate assembly and reference using pggb

```
#!/bin/bash
set -x
set -v

samtools=/lizardfs/salehi/micromamba/bin/samtools

# Input filename
input1=$1
input2=$2
reference=$3

# output filename
base1=$(basename "$input1" .fa.gz)

output1=$4
output2=$5
outputfinal=$6


# First command

zcat "$input1" | awk '/^>/{sub(">", ">'$base1'#"1"#")}1' > "$output1"


# Second command

zcat "$input2" | awk '/^>/{sub(">", ">'$base1'#"2"#")}1' > "$output2"


#Third command combine them together

(cat "$output1" "$output2" "$reference") | bgzip > "$outputfinal"


#Forth command samtools

$samtools faidx "$outputfinal"
sbatch -p workers -c 48 --wrap 'pggb -D /scratch/PGGB -i "$outputfinal" -t 8 -o outputvcf2.2 -V GRCr8'
```

#Next step is using wfmash to match each ID in gff format of grcr8 to name of the genes
#First it is needed to download *.gff and *.fna from this link ```https://ftp.ncbi.nlm.nih.gov/genomes/all/annotation_releases/10116/GCF_036323735.1-RS_2024_02/```
#The other file we need is Grcr8.fa.gz

```
#!/bin/bash


wfmash=/gnu/store/8kma8dr0g1h6kz04dkcs336c9zn2q7gm-wfmash-0.12.5-1+0222f7c/bin/wfmash



# Input files
target_fa=$1
query_fa=$2
output=$3


# Run wfmash with the specified parameters
$wfmash -t 48  $1 $2 > $3

```
```bash script.sh GRCr8.fa.gz GCF_036323735.1_GRCr8_genomic.fna.gz outputwfmash.paf```

#Then we use this command to update the name in our gff file: 

```awk -F'\t' 'BEGIN {OFS="\t"} NR==FNR{a[$1]=$6;next} $1 in a{$1=a[$1]}1' outputwfmash.paf GCF_036323735.1_GRCr8_genomic.gff > updated_gff_file.gff```

Next script : extracting Heterozigosity, and related genes.

```
#!/bin/bash
set -x
set -v



input1=$1

base1=$(basename "$input1"combinedfinalGRCr8.fa.gz.bf3285f.eb0f3d3.11fe66b.smooth.final.grcr8.vcf)
base2=$(basename "$input1"combinedfinalGRCr8vcf.2)
output1="$input1"heterovariant.vcf
output2="$input1"heterovariants.bed
output3="$input1"bedtoolsintersectwao.bed
output4="$input1"gene_list.txt
output5="$input1"gene_list_without_prefix.txt
output6="$input1"unique_gene_list.txt

#extracting heterozigocity from vcf file
bcftools view  -i 'GT="0|1" || GT="1|0" || GT="0|2" || GT="2|0"' /lizardfs/salehi/Ass_Ref/$input1/pggb/$base2/$base1  > /lizardfs/salehi/Ass_Ref/$input1/pggb/enrich/$output1                                  

cut -f 1,2,4,10 /lizardfs/salehi/Ass_Ref/$input1/pggb/enrich/$output1 > /lizardfs/salehi/Ass_Ref/$input1/pggb/enrich/$output2                                                                                  

bedtools intersect -wao -a $output1 -b /lizardfs/salehi/Ass_Ref/BN_Lx_Cub/pggb/BN_Lx_CubcombinedfinalGRCr8vcf.2/test/updated_gff_file4.gff > /lizardfs/salehi/Ass_Ref/$input1/pggb/enrich/$output3             


grep -o 'gene-[^;]\+' /lizardfs/salehi/Ass_Ref/$input1/pggb/enrich/$output3 > /lizardfs/salehi/Ass_Ref/$input1/pggb/enrich/$output4                                                                            

awk '{sub(/^gene-/, "", $0); print}' /lizardfs/salehi/Ass_Ref/$input1/pggb/enrich/$output4 > /lizardfs/salehi/Ass_Ref/$input1/pggb/enrich/$output5                                                             


sort /lizardfs/salehi/Ass_Ref/$input1/pggb/enrich/$output5 | uniq > /lizardfs/salehi/Ass_Ref/$input1/pggb/enrich/$output6                                                                                      


echo This is something $(date)
```


Script for copy files: 
```
#!/bin/bash
set -x
set -v



input1=$1
base1=$(basename "$input1"unique_gene_list.txt)
base2=$(basename "$input1"gene_list_without_prefix.txt)

scp salehi@octopus01:/lizardfs/salehi/Ass_Ref/$input1/pggb/enrich/$base1 salehi@octopus01:/lizardfs/salehi/Ass_Ref/$input1/pggb/enrich/$base2  /Users/erikgarrison/Desktop/GO/Enrichment\ Analysis/$input1
```





#Intersect the het region calls with genes from grcr8

```bedtools intersect -wao -a heterovarinats.vcf -b updated_gff_file.gff > bedtoolsintersectwao.bed```




#Given gene list for Gene Ontology

Enrichment Analysis

we use four different application gProfiler, DAVID, EnrichR, and ToppGene.

for EnrichR, and ToppGene, they do not accept rat gene list; therefore we need to translate our rat gene list to human. 


For this, we use following R script:


```
######homologene#####
install.packages("homologene")
library(homologene)

setwd("/Users/.../")



# Install and load the homologene package if not already installed
if (!requireNamespace("homologene", quietly = TRUE)) {
  install.packages("homologene")
}
library(homologene)

# Specify the path to the gene list file
gene_list_file = "unique_gene_list.txt" # Update with the correct path if needed

# Check if the file exists before attempting to read it
if (!file.exists(gene_list_file)) {
  stop(paste("File not found:", gene_list_file))
}

# Read the list of rat genes from the input file
rat_genes = readLines(gene_list_file)

# Function to retrieve human orthologs
get_human_orthologs <- function(rat_genes) {
  orthologs <- homologene(rat_genes, inTax = 10116, outTax = 9606)
  orthologs_df <- as.data.frame(orthologs)
  return(orthologs_df)
}

# Retrieve human orthologs
human_orthologs_df = get_human_orthologs(rat_genes)

# View the results
print(human_orthologs_df)


# Extract the first column of the data frame
first_column <- human_orthologs_df[, 1]

# Save the first column to a file
output_file <- "human_orthologs.txt"
write.table(first_column, file = output_file, quote = FALSE, row.names = FALSE, col.names = FALSE)

```







