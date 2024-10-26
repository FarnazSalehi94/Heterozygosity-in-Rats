Heterozygosity-in-Rats

# Mapping 
## First Step:

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

## Second Step:

### Editing column name;


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

## Third Step: 


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
bcftools view  -i 'GT="0|1" || GT="1|0" || GT="0|2" || GT="2|0"' outputpggb.vcf  > $output1                                  
cut -f 1,2,4,10 $output1 > $output2                                                                                  

bedtools intersect -wao -a $output1 -b gff_proteincoding.gff > $output3             

grep -o 'gene-[^;]\+' $output3 > $output4                                                                            

awk '{sub(/^gene-/, "", $0); print}' $output4 > $output5                                                             

sort $output5 | uniq > $output6                                                                                      

echo This is something $(date)
```

# Heatmap


# Coverage

# Segmental Duplications


# Het Rate

In this section, for calculation the per-base pair rate of heterozygosity for each gene, the following steps are needed:

__1.__ Counting the overlaps for each gene.

__2.__ Dividing the count by the length of each gene.


```#!/bin/bash

# Check if the assembly name, path, and destination are provided
if [ "$#" -ne 3 ]; then
            echo "Usage: $0 <assembly_name> <path> <destination>"
                exit 1
fi

# Set the assembly name, path, and destination
assembly=$1
path=$2
destination=$3

# Create the 'hetrate' directory in the provided path and move to it
mkdir -p "$path/hetrate"
cd "$path/hetrate" || { echo "Failed to change directory to $path/hetrate"; exit 1; }

# first step : Het Rate counts
awk '{count[$1]++} END {for (gene in count) print gene, count[gene]}' "$path/${assembly}gene_list_without_prefix.txt" > "${assembly}gene_counts.txt"

# second step: Gene Length
awk 'BEGIN {OFS="\t"} $3 == "gene" {split($9, a, ";"); for (i in a) {if (a[i] ~ /^Name=/) {split(a[i], b, "="); print b[2], $1, $4, $5}}}' /lizardfs/salehi/Ass_Ref/BN_Lx_Cub/pggb/BN_Lx_CubcombinedfinalGRCr8vcf.2/test/gff_proteincoding.gff > "${assembly}gene_positions.txt"

# third : merging the files
join -1 1 -2 1 <(sort "${assembly}gene_counts.txt") <(sort "${assembly}gene_positions.txt") > "${assembly}merged_output.txt"

# fourth step: measuring het rate
awk 'BEGIN { OFS="\t" } 
{
          if (NF >= 5) {  # Check if there are at least 5 fields
                      total_variants = $2;
                          region_length = $5 - $4;
                              het_rate = total_variants / region_length;
                                  printf "%s\t%s\t%.8f\n", $0, region_length, het_rate; 
                                    } else {
                                        print "Error: Not enough fields in line:", $0;  # Error handling
                                          }
                          }' "${assembly}merged_output.txt" > "${assembly}merged_with_hetrate.txt"

                  # fifth step: sorting based on the Het Rate
                  sort -k 7,7nr "${assembly}merged_with_hetrate.txt" -o "${assembly}merged_with_hetrate_sorted.txt"

                  # header
                  awk 'BEGIN { print "Contig\tStart\tEnd\tGeneName\tCount\tLength\tHetRate" } { print $3, $4, $5, $1, $2, $6, $7 }' OFS="\t" "${assembly}merged_with_hetrate_sorted.txt" > "${assembly}.hetrate.bed"

                  # sixth step : top 100
                  awk 'NR==1 { print; next } { print $0 | "sort -k7,7nr" }' "${assembly}.hetrate.bed" | head -n 100 > "${assembly}top100.bed"

                  # Copy the top100 file to the specified destination
                  cp "${assembly}top100.bed" "$destination"

                  echo "Processing complete. Output files prefixed with '${assembly}_', and '${assembly}top100.bed' copied to '$destination'."

```


# Enrichment Analysis

## gProfiler:

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
gene_list <- scan("input.txt", what = "", sep = "\n")

# Perform the enrichment analysis using g:Profiler
gprofiler_results <- gost(query = gene_list, organism = "rnorvegicus") #replace this with your organism

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

### Homologene

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



## EnrichR:


```vi automatedEnrichR.R```

```
#!/usr/bin/env Rscript

# Load necessary libraries
library(enrichR)
library(ggplot2)
library(dplyr)

# Define function to perform enrichment analysis and save results
perform_enrichment_analysis <- function(gene_list_file) {
  # Read gene list
  genes <- readLines(gene_list_file)
  
  # Define databases for enrichment analysis
  databases <- c("KEGG_2021_Human", "Reactome_2022", "GO_Biological_Process_2023", 
                 "GO_Cellular_Component_2023", "GO_Molecular_Function_2023")

                 
  
  # Perform enrichment analysis
  enrich_results <- tryCatch({
    enrichr(genes, databases)
  }, error = function(e) {
    message("Error in retrieving results: ", e)
    return(NULL)
  })
  
  if (is.null(enrich_results)) {
    return(NULL)
  }
  

  
  # Save results to CSV files
  for (db in names(enrich_results)) {
    write.csv(enrich_results[[db]], paste0(db, "_enrichment_results.csv"), row.names = FALSE)
  }
  
  # Function to plot top 10 enrichment results with gradient color
  plot_top_10 <- function(enrichment_data, title) {
    top_10 <- enrichment_data %>%
      arrange(P.value) %>%
      head(10)
    
    top_10$Term <- factor(top_10$Term, levels = rev(top_10$Term))
    
    ggplot(top_10, aes(x = Term, y = -log10(P.value), fill = -log10(P.value))) +
      geom_bar(stat = "identity") +
      coord_flip() +
      scale_fill_gradient(low = "#c6dbef", high = "#08306b") +
      labs(title = title, x = "Pathway/Term", y = "-log10(P.value)") +
      theme_minimal(base_size = 14) +
      theme(plot.background = element_rect(fill = "white"),
            panel.background = element_rect(fill = "white"),
            panel.grid.major = element_line(color = "grey90"),
            panel.grid.minor = element_line(color = "grey90"),
            plot.title = element_text(hjust = 0.5, face = "bold"),
            axis.title.x = element_text(face = "bold"),
            axis.title.y = element_text(face = "bold"),
            legend.position = "none")
  }
  
  # Plot and save top 10 results for each database
  for (db in names(enrich_results)) {
    if (nrow(enrich_results[[db]]) > 0) {
      p <- plot_top_10(enrich_results[[db]], paste0("Top 10 Enriched Terms in ", db))
      ggsave(paste0(db, "_top10.png"), plot = p, width = 15, height = 5)
    } else {
      message(paste("No results for", db))
    }
  }
}

# Check if the correct number of arguments is provided
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1) {
  stop("Usage: Rscript enrichment_analysis.R <gene_list_file>")
}

# Perform enrichment analysis
perform_enrichment_analysis(args[[1]])
```

```Rscript automatedEnrichR.R inputfile.txt```










