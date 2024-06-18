Heterozygosity-in-Rats


First Step:

We need assemblies (alternate and primary) and reference (latest version of the reference in your species).
The file types are FASTA.
Align all these using PGGB.

```
#!/bin/bash
set -x
set -v

samtools=/lizardfs/salehi/micromamba/bin/samtools

# Input filename
input1=$1 #for primary file
input2=$2 #for alternate file
reference=$3 #the file you select as reference for your specious 

# output filename
base1=$(basename "$input1" .fa.gz) # if your file is not zipped you can remove .gz

#editing the name of contigs like #name_of_assembly#type of strain(1 or 2)#name of the contig
output1=$4 #primary file after editing name of the contigs 
output2=$5 #alternate file after editing name of the contigs

outputfinal=$6 #aligned file for all three parametes together


# First command

zcat "$input1" | awk '/^>/{sub(">", ">'$base1'#"1"#")}1' > "$output1"


# Second command

zcat "$input2" | awk '/^>/{sub(">", ">'$base1'#"2"#")}1' > "$output2"


#Third command combine them together

(cat "$output1" "$output2" "$reference") | bgzip > "$outputfinal"


#Forth command samtools

$samtools faidx "$outputfinal" #need to index file before running pggb

sbatch -p workers -c 48 --wrap 'pggb -D /scratch/PGGB -i "$outputfinal" -t 48 -o outputvcf2.2 -V grcr8'
```

Second Step:

Editing column name;


Downloading GFF and FNA format files from these links:

NCBI Link 1 : https://www.ncbi.nlm.nih.gov/datasets/genome/GCA_036323735.1/
NCBI Link 2 : https://ftp.ncbi.nlm.nih.gov/genomes/all/annotation_releases/10116/GCF_036323735.1-RS_2024_02/


In addition, we need our reference FASTA file.


```
#!/bin/bash


wfmash=/gnu/store/8kma8dr0g1h6kz04dkcs336c9zn2q7gm-wfmash-0.12.5-1+0222f7c/bin/wfmash

bcftools=/home/salehi/.guix-profile/bin/bcftools

# Input files
target_fa=$4 #GRCr8.fa.gz
query_fa=$5 #GCF_036323735.1_GRCr8_genomic.fna.gz

$bcftools index -c $query 

# Run wfmash with the specified parameters
$wfmash -t 48  $target $query > wfmashoutput.paf

# Update the name in gff file:
awk -F'\t' 'BEGIN {OFS="\t"} NR==FNR{a[$1]=$6;next} $1 in a{$1=a[$1]}1' wfmashoutput.paf GCF_036323735.1_GRCr8_genomic.gff > updated_gff_file4.gff
```
There are different types of gene biotypes in the GFF file. To extract the gene biotype "protein coding," you need to look at your file and select the desired gene biotype.

Extracting protein-coding genes.

```grep 'gene_biotype=protein_coding' updated_gff_file4.gff > gff_proteincoding.gff```

Third Step: 


Intersecting the first and second steps to get a gene list related to heterozygosity.

Extracting heterozygosity regions to get a unique gene list.

Now, we have two separate files: one from the first step, which is a VCF file including heterozygosity sites, and the other as a reference genome from second step which we extract protein-coding genes.

Now, we need to intersect these two files to identify the genes related to and causing these heterozygosity regions.


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
bcftools view  -i 'GT="0|1" || GT="1|0" || GT="0|2" || GT="2|0"' /lizardfs/salehi/Ass_Ref/$input1/pggb/$base2/$base1  > $output1                                  
cut -f 1,2,4,10 $output1 > $output2                                                                                  

bedtools intersect -wao -a $output1 -b gff_proteincoding.gff > $output3             

grep -o 'gene-[^;]\+' $output3 > $output4                                                                            

awk '{sub(/^gene-/, "", $0); print}' $output4 > $output5                                                             

sort $output5 | uniq > $output6                                                                                      

echo This is something $(date)
```



Enrichment Analysis:

gProfiler

```
cat > automatedEAgprofiler.R << 'EOF'
# Load the necessary libraries
library(gprofiler2)
library(gridExtra)
library(gprofiler2)
library(ggplot2)
library(dplyr)
library(cowplot)
library(gridExtra)
library(tibble)
library(ggplotify)




# Read the gene list from the file
gene_list <- scan("BN_Lx_Cubunique_gene_list.txt", what = "", sep = "\n")

# Perform the enrichment analysis using g:Profiler
gprofiler_results <- gost(query = gene_list, organism = "rnorvegicus")

# Generate the main plot using gostplot
main_plot <- gostplot(gprofiler_results, capped = FALSE, interactive = FALSE)



# Extract the relevant data for the table
plot_data <- gprofiler_results$result %>%
  as_tibble() %>%
  mutate(log_p_value = -log10(p_value)) %>%
  arrange(desc(log_p_value))

# Extract top results for the table
top_results <- plot_data %>%
  group_by(source) %>%
  top_n(3, log_p_value) %>%
  ungroup() %>%
  arrange(source, desc(log_p_value)) %>%
  dplyr::select(source, term_id, term_name, log_p_value)

# Create the table plot with adjusted size and white background
# Create the table plot with adjusted size and white background
table_theme <- ttheme_default(
  core = list(bg_params = list(fill = "white", col = "black")),
  colhead = list(bg_params = list(fill = "white", col = "black")),
  rowhead = list(bg_params = list(fill = "white", col = "black"))
)

table_plot <- tableGrob(top_results, rows = NULL, theme = table_theme)



# Combine the main plot and the table plot using ggdraw and draw_plot
combined_plot <- ggdraw() +
  draw_plot(main_plot, x = 0, y = 0.35, width = 1, height = 0.65) +
  draw_plot(table_plot, x = 0, y = 0, width = 1, height = 0.25)



# Save the combined plot as a PNG file
ggsave("gpr_results.png", plot = combined_plot, width = 16, height = 12)

# Print a message indicating that the process is complete
cat("Enrichment analysis complete and result saved as gprofiler_combined_results.png\n")

EOF
```

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










