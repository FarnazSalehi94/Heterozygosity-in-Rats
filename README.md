# Heterozygosity-in-Rats


#First thing is writing script to aligning primary, alternate assembly and reference using pggb

```#!/bin/bash
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

#The next step should be extracting heterozigocity from vcf file

```bcftools view -i 'GT="1|0" || GT="0|1"' combinedfinal.fa.gz.bf3285f.eb0f3d3.867196c.smooth.final.grcr8.vcf > output.vcf```

#convert the file to the bed format

```awk '!/^#/ && $10 ~ /0\/1/ {print $1"\t"$2 - 1"\t"$2"\t"$4"\t"$10}' heterovarinats.vcf > variants.bed```


#Next step is using wfmash to match each ID in gff format of grcr8 to name of the genes

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


#Intersect the het region calls with genes from grcr8

```bedtools intersect -wao -a heterovarinats.vcf -b updated_gff_file.gff > bedtoolsintersectwao.bed```

#Given gene list for Gene Ontology




