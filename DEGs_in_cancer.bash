#!/bin/bash

### Step 1: Load required modules
module load anaconda3/2023.09-0

### Step 2: Download GTEX and TCGA GEM files for bladder cancer
echo "Downloading data files..."
wget -O BLCA_GTEX.gz https://ndownloader.figshare.com/files/9150181
wget -O BLCA_TCGA.gz https://ndownloader.figshare.com/files/9150184

### Step 3: Decompress the downloaded files
echo "Decompressing files..."
gunzip -f BLCA_GTEX.gz
gunzip -f BLCA_TCGA.gz

# Verify files exist and are not empty
if [ ! -s BLCA_GTEX ] || [ ! -s BLCA_TCGA ]; then
    echo "ERROR: Downloaded files are empty or failed to decompress"
    exit 1
fi

### Step 4: Clone the mergegem repository
echo "Cloning mergegem repository..."
if [ ! -d "mergegem" ]; then
    git clone https://github.com/feltus/mergegem.git
fi

### Step 5: Merge GEM files
echo "Merging GEM files..."
python ./mergegem/mergegem.py BLCA_GTEX BLCA_TCGA bladder-gtex-tcga.txt

# Verify merge was successful
if [ ! -s bladder-gtex-tcga.txt ]; then
    echo "ERROR: Merge failed - output file is empty"
    echo "Checking if input files are valid..."
    head -n 5 BLCA_GTEX
    head -n 5 BLCA_TCGA
    exit 1
else
    echo "Merge successful!"
fi

### Step 6: Preprocess the GEM file
echo "Preprocessing the GEM file..."

# Convert RSEM decimal values into integers
sed -e 's/\.[0-9]*//g' -e 's/ *$//' bladder-gtex-tcga.txt > bladder-gtex-tcga-integer.txt
echo "Converted decimal values to integers"

# Remove duplicate gene rows
awk '!a[$1]++' bladder-gtex-tcga-integer.txt > bladder-gtex-tcga-integer-unique.txt
echo "Removed duplicate gene rows"

# Replace dashes with underscores to comply with DESeq2 requirements
sed 's/-/_/g' bladder-gtex-tcga-integer-unique.txt > bladder-gtex-tcga-clean.txt
echo "Replaced dashes with underscores"

# Convert to CSV format
sed 's/\s/,/g' bladder-gtex-tcga-clean.txt > bladder-gtex-tcga-clean.csv
echo "Converted to CSV format"

# Verify output file exists and is not empty
if [ ! -s bladder-gtex-tcga-clean.csv ]; then
    echo "ERROR: Preprocessing failed - CSV output file is empty"
    exit 1
else
    echo "Preprocessing successful!"
    echo "Lines in final CSV file: $(wc -l < bladder-gtex-tcga-clean.csv)"
fi

### Step 7: Prepare DESeq2 group definition file
echo "Preparing DESeq2 group definition file..."

# Check if labels file was created during merge
if [ ! -s bladder-gtex-tcga.labels.txt ]; then
    echo "ERROR: Labels file not found. Check if mergegem created it properly."
    exit 1
fi

# Convert dashes to underscores in the labels file
sed 's/-/_/g' bladder-gtex-tcga.labels.txt > bladder-gtex-tcga-dash.labels.txt

# Add comparison label to each row
sed 's/\s/,/g' bladder-gtex-tcga-dash.labels.txt | sed 's/$/,GTEX_BLCA_TCGA/' > bladder-gtex-tcga.comparison.tmp

# Add header row
echo "Sample,Group,Comparison" > bladder-gtex-tcga.comparison.csv
grep -v "^Sample" bladder-gtex-tcga.comparison.tmp >> bladder-gtex-tcga.comparison.csv

echo "DESeq2 group definition file created successfully"

### Step 8: Generate the R script for DESeq2 analysis
cat > run_deseq2.R << 'EOF'
# Load required libraries
library(DESeq2)
library(ggplot2)
library(pheatmap)
library(RColorBrewer)
library(EnhancedVolcano)

# Read in the count matrix
countData <- read.csv("bladder-gtex-tcga-clean.csv", header=TRUE, row.names=1)

# Read in the sample information
colData <- read.csv("bladder-gtex-tcga.comparison.csv", header=TRUE, row.names=1)

# Ensure the row names in colData match the column names in countData
if (!all(colnames(countData) %in% rownames(colData))) {
  stop("Column names in count matrix do not match row names in sample info")
}
# Reorder colData to match countData
colData <- colData[colnames(countData),]

# Create DESeq2 dataset
dds <- DESeqDataSetFromMatrix(countData = countData,
                              colData = colData,
                              design = ~ Group)

# Run DESeq2
dds <- DESeq(dds)

# Extract results
res <- results(dds, contrast=c("Group", "TCGA", "GTEX"))

# Order results by adjusted p-value
resOrdered <- res[order(res$padj),]

# Write results to CSV
write.csv(as.data.frame(resOrdered), 
          file="bladder-gtex-tcga-deseq2-results.csv")

# Create MA plot
png("bladder-gtex-tcga-ma-plot.png", width=800, height=600)
plotMA(res, main="MA Plot", ylim=c(-5,5))
dev.off()

# Create Volcano plot
png("bladder-gtex-tcga-volcano-plot.png", width=1000, height=800)
EnhancedVolcano(res,
                lab = rownames(res),
                x = 'log2FoldChange',
                y = 'padj',
                title = 'TCGA vs GTEX',
                pCutoff = 0.05,
                FCcutoff = 1,
                pointSize = 2.0,
                labSize = 3.0)
dev.off()

# Count significant genes
sig_genes <- sum(res$padj < 0.05, na.rm=TRUE)
cat("Number of significantly differentially expressed genes:", sig_genes, "\n")

# Count upregulated and downregulated genes
up_reg <- sum(res$padj < 0.05 & res$log2FoldChange > 0, na.rm=TRUE)
down_reg <- sum(res$padj < 0.05 & res$log2FoldChange < 0, na.rm=TRUE)
cat("Up-regulated genes:", up_reg, "\n")
cat("Down-regulated genes:", down_reg, "\n")

# Extract top 50 differentially expressed genes for heatmap
res_df <- as.data.frame(resOrdered)
top_genes <- head(rownames(res_df), 50)
top_gene_counts <- countData[top_genes, ]

# Generate heatmap of top genes
png("bladder-gtex-tcga-heatmap.png", width=1200, height=800)
pheatmap(log2(top_gene_counts + 1), 
         scale="row", 
         annotation_col=colData[,c("Group"), drop=FALSE],
         main="Top 50 DEGs",
         fontsize_row=8,
         fontsize_col=6)
dev.off()

print("DESeq2 analysis complete. Results saved to bladder-gtex-tcga-deseq2-results.csv")
EOF

echo "R script for DESeq2 analysis created as 'run_deseq2.R'"
echo "To run the analysis, execute: Rscript run_deseq2.R"

echo "Script completed successfully!"
