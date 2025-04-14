#install DESeq2
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("DESeq2")

#set working directory
setwd("/scratch/kemeeha/DEG_analysis_finding_DEGs_in_cancer")

#1. extract a sub-matrix of counts for each group
subgem <- function(gem, anot, group ){
  datalist = list()
  subanot = subset(anot, Comparison == group)
  for (id in subanot$Sample) {
    ind = which(colnames(gem) == id)
    genes = gem[0]
    exp = gem[,ind, drop=FALSE]
    datalist[[id]] <- exp
  }
  subcounts = cbind(genes, datalist)
  return(subcounts)
}

#2. extract a subset of the sample annotation matrix for each group
subanot <- function(anot, group){
  datalist = list()
  print(str(group))
  subanot = subset(anot, Comparison == group)
  rownames(subanot) <- subanot$Sample  # Ensures rownames match sample IDs
  print(str(subanot))
  return(subanot)
}

#3. run DESeq2
run_deseq <- function(counts, annotation){
  # Match column names in counts with sample IDs in annotation
  print("Counts column names:")
  print(colnames(counts))
  print("Annotation row names:")
  print(rownames(annotation))
  
  common_samples <- intersect(colnames(counts), rownames(annotation))
  if (length(common_samples) == 0) {
    stop("No matching sample names between counts and annotation.")
  }

  counts <- counts[, common_samples]
  annotation <- annotation[common_samples, ]

  dds <- DESeqDataSetFromMatrix(countData = counts,
                                colData = annotation,
                                design = ~ Group)

  #Filter and normalize genes with low total counts across all samples (Edit as needed): 
  dds <- dds[rowSums(counts(dds)) >= 50,]
  dds <- DESeq(dds)
  norm = fpm(dds)

  #Sort the columns in the FPM data frame (Edit as needed): 
  conditionA = which(annotation$Group == "BLCA_GTEX")  #EDIT
  conditionB = which(annotation$Group == "BLCA_TCGA")  #EDIT
  norm = subset(norm, select=c(conditionA, conditionB))
  print(str(norm))

  #Retrieve the results (Edit groupIDs as needed): 
  res <- results(dds, contrast=c("Group", "BLCA_GTEX", "BLCA_TCGA"))
  print(summary(res))
  res <- cbind(res, norm) # Add FPM values to results for easy visualization
  resultsNames(dds)
  return(res)
}

#4. Print results to a file
main <- function(countfile, anotfile, outfile){
  outname = outfile
  counts = read.delim(countfile, sep=',', header=TRUE, row.names='Hugo_Symbol') #EDIT
  samples = read.delim(anotfile, sep=',', row.names = NULL, check.names=FALSE)
  groups = unique(samples$Comparison)
  for (t in groups){
    subcounts = subgem(counts, samples, t)
    subannotation = subanot(samples, t)
    results = run_deseq(subcounts, subannotation)
    #Filter and sort results table
    f_results = subset(results, padj < 0.05)
    o_results = f_results[order(f_results$padj),]
    write.csv(o_results, outname, row.names = TRUE)
  }
  return(results)
}

#5. run the analysis with the input file names and make an output file name
main('bladder-gtex-tcga-clean.csv', 'bladder-gtex-tcga.comparison.csv', 'bladder-gtex-tcga-degs.csv')
