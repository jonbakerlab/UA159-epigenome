## DESeq2 
## Updated for multiple sample (pairwise) comparisons RNA-seq 6-3-26

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(version = "3.22")

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("DESeq2")
BiocManager::install("genefilter")

# Set working directory
setwd("~/path/to/working/dir/")
getwd()

# Import & pre-process ----------------------------------------------------

# Import data from featureCounts

countdata <- read.table("featureCount.txt", header=TRUE, row.names=1)

# Convert to matrix
#countdata = as.matrix(UAB_deSeq_input)
countdata <- as.matrix(countdata)
head(countdata)

## Parse GTF annotation to map gene IDs to display names
gtf_file <- "{ref.gtf}"
display_map <- NULL
if (file.exists(gtf_file)) {
  gtf <- read.table(gtf_file, sep="\t", header=FALSE, quote="", stringsAsFactors=FALSE, comment.char="#")
  if (ncol(gtf) >= 9) {
    attrs <- gtf[,9]
    parse_attrs <- function(s) {
      parts <- unlist(strsplit(s, ";"))
      parts <- trimws(parts)
      res <- list()
      for (p in parts) {
        if (nchar(p) == 0) next
        kv <- unlist(strsplit(p, " ", perl=TRUE))
        key <- kv[1]
        val <- paste(kv[-1], collapse=" ")
        val <- gsub('^"|"$', "", trimws(val))
        res[[key]] <- val
      }
      return(res)
    }
    parsed <- lapply(attrs, parse_attrs)
    gene_ids <- vapply(parsed, function(x) if (!is.null(x$gene_id)) x$gene_id else if (!is.null(x$ID)) x$ID else NA_character_, character(1))
    gene_field <- vapply(parsed, function(x) {
      if (!is.null(x$gene)) return(x$gene)
      if (!is.null(x$gene_name)) return(x$gene_name)
      if (!is.null(x$Name)) return(x$Name)
      return(NA_character_)
    }, character(1))
    old_locus <- vapply(parsed, function(x) {
      if (!is.null(x$old_locus_tag)) return(x$old_locus_tag)
      if (!is.null(x$old_locus)) return(x$old_locus)
      if (!is.null(x$locus_tag)) return(x$locus_tag)
      return(NA_character_)
    }, character(1))
    anno_df <- data.frame(gene_id=gene_ids, gene=gene_field, old_locus=old_locus, stringsAsFactors=FALSE)
    ## keep first occurrence per gene_id
    anno_df <- anno_df[!is.na(anno_df$gene_id), ]
    anno_df <- anno_df[!duplicated(anno_df$gene_id), ]
    display_map <- setNames(ifelse(!is.na(anno_df$gene) & nzchar(anno_df$gene), anno_df$gene, ifelse(!is.na(anno_df$old_locus) & nzchar(anno_df$old_locus), anno_df$old_locus, anno_df$gene_id)), anno_df$gene_id)
  }
} else {
  message("GTF file not found; will use gene IDs as display names")
}

# Assign samples
(condition <- factor(c(
  rep("{sample1}", {#replicates}), 
  rep("{sample2}", {#replicates}),
  rep("{sample3}", {#replicates}),
  rep("{sample4}", {#replicates}),
  rep("{sampleN}", {#replicates})
)))

# DESeq analysis----------------------------------------------------

library(DESeq2)


# Create a coldata frame and instantiate the DESeqDataSet 
(coldata <- data.frame(row.names=colnames(countdata), condition))
dds <- DESeqDataSetFromMatrix(countData=countdata, colData=coldata, design=~condition)
dds

#set a reference level
dds$condition <- relevel(dds$condition, "{reference_condition}")

# Run the DESeq pipeline
dds <- DESeq(dds)

# Plot dispersions
png("qc-dispersions.png", 1000, 1000, pointsize=20)
plotDispEsts(dds, main="Dispersion plot")
dev.off()

# Regularized log transformation for clustering/heatmaps, etc
rld <- rlogTransformation(dds)
head(assay(rld))
hist(assay(rld))

# Colors for plots 
## Use RColorBrewer
library(RColorBrewer)
(mycols <- brewer.pal(8, "Dark2")[1:length(unique(condition))])

# Sample distance heatmap
sampleDists <- as.matrix(dist(t(assay(rld))))
library(gplots)
png("qc-heatmap-samples.png", w=1000, h=1000, pointsize=20)
heatmap.2(as.matrix(sampleDists), key=F, trace="none",
          col=colorpanel(100, "black", "white"),
          ColSideColors=mycols[condition], RowSideColors=mycols[condition],
          margin=c(11, 11), main="Sample Distance Matrix")
dev.off()

# Principal components analysis
## Could do with built-in DESeq2 function:
## DESeq2::plotPCA(rld, intgroup="condition") + theme_bw()
rld_pca <- function (rld, intgroup = "condition", ntop = 500, colors=NULL, legendpos="bottomleft", main="PCA Biplot", textcx=1, ...) {
  require(genefilter)
  require(calibrate)
  require(RColorBrewer)
  rv = rowVars(assay(rld))
  select = order(rv, decreasing = TRUE)[seq_len(min(ntop, length(rv)))]
  pca = prcomp(t(assay(rld)[select, ]))
  fac = factor(apply(as.data.frame(colData(rld)[, intgroup, drop = FALSE]), 1, paste, collapse = " : "))
  if (is.null(colors)) {
    if (nlevels(fac) >= 3) {
      colors = brewer.pal(nlevels(fac), "Paired")
    }   else {
      colors = c("black", "red")
    }
  }
  pc1var <- round(summary(pca)$importance[2,1]*100, digits=1)
  pc2var <- round(summary(pca)$importance[2,2]*100, digits=1)
  pc1lab <- paste0("PC1 (",as.character(pc1var),"%)")
  pc2lab <- paste0("PC1 (",as.character(pc2var),"%)")
  plot(PC2~PC1, data=as.data.frame(pca$x), bg=colors[fac], pch=21, xlab=pc1lab, ylab=pc2lab, main=main, ...)
  with(as.data.frame(pca$x), textxy(PC1, PC2, labs=rownames(as.data.frame(pca$x)), cex=textcx))
  legend(legendpos, legend=levels(fac), col=colors, pch=20)
  #     rldyplot(PC2 ~ PC1, groups = fac, data = as.data.frame(pca$rld),
  #            pch = 16, cerld = 2, aspect = "iso", col = colours, main = draw.key(key = list(rect = list(col = colours),
  #                                                                                         terldt = list(levels(fac)), rep = FALSE)))
}
png("qc-pca.png", 1000, 1000, pointsize=20)
rld_pca(rld, colors=mycols, intgroup="condition", xlim=c(-50, 50))
dev.off()


## DE results: pairwise comparisons between all conditions
## Function to draw MA plot for a comparison
maplot <- function (res, name="comparison", thresh=0.05, labelsig=TRUE, textcx=1, ...) {
  with(res, plot(baseMean, log2FoldChange, pch=20, cex=.5, log="x", ...))
  with(subset(res, padj<thresh), points(baseMean, log2FoldChange, col="red", pch=20, cex=1.5))
  if (labelsig) {
    require(calibrate)
    with(subset(res, padj<thresh), textxy(baseMean, log2FoldChange, labs=name, cex=textcx, col=2))
  }
}

## Volcano plot with labeled significant genes
volcanoplot <- function (res, name="comparison", lfcthresh=2, sigthresh=0.01, main="Volcano Plot", legendpos="bottomright", labelsig=TRUE, textcx=1, xlim=NULL, ...) {
  # determine x-axis limits from data if not provided
  lfc <- res$log2FoldChange
  if (is.null(xlim)) {
    finite_lfc <- lfc[is.finite(lfc)]
    if (length(finite_lfc) == 0) {
      xlim <- c(-2, 2)
    } else {
      xmin <- min(finite_lfc)
      xmax <- max(finite_lfc)
      rng <- xmax - xmin
      pad <- ifelse(rng == 0, 1, rng * 0.05)
      xlim <- c(xmin - pad, xmax + pad)
    }
  }
  with(res, plot(log2FoldChange, -log10(pvalue), pch=20, main=main, xlim=xlim, ...))
  with(subset(res, padj<sigthresh ), points(log2FoldChange, -log10(pvalue), pch=20, col="red", ...))
  with(subset(res, abs(log2FoldChange)>lfcthresh), points(log2FoldChange, -log10(pvalue), pch=20, col="orange", ...))
  with(subset(res, padj<sigthresh & abs(log2FoldChange)>lfcthresh), points(log2FoldChange, -log10(pvalue), pch=20, col="green", ...))
  if (labelsig) {
    require(calibrate)
    with(subset(res, padj<sigthresh & abs(log2FoldChange)>lfcthresh), textxy(log2FoldChange, -log10(pvalue), labs = name, cex=textcx, ...))
  }
  legend(legendpos, xjust=1, yjust=1, legend=c(paste("FDR<",sigthresh,sep=""), paste("|LogFC|>",lfcthresh,sep=""), "both"), pch=20, col=c("red","orange","green"))
}

## Significance cutoffs for collecting differential genes
sig_padj <- 0.05
sig_lfc <- 1

## Run pairwise comparisons for all condition pairs and save results/plots
conds <- levels(coldata$condition)
pairs <- t(combn(conds, 2))
sig_genes_all_ids <- character()
sig_by_comp_ids <- list()
for (i in seq_len(nrow(pairs))) {
  a <- pairs[i, 1]
  b <- pairs[i, 2]
  comp_name <- paste0(b, "_vs_", a)
  message("Running comparison: ", comp_name)
  res <- results(dds, contrast=c("condition", a, b))

  ## Basic summaries
  print(table(res$padj < 0.05))

  ## Order by adjusted p-value
  res <- res[order(res$padj), ]

  ## Merge with normalized counts
  resdata <- merge(as.data.frame(res), as.data.frame(counts(dds, normalized=TRUE)), by="row.names", sort=FALSE)
  names(resdata)[1] <- "gene"

  ## Collect significant genes for combined clustermap and per-comparison lists
  sig_idx <- which(!is.na(res$padj) & res$padj < sig_padj & abs(res$log2FoldChange) > sig_lfc)
  sig_genes_cur <- character()
  if (length(sig_idx) > 0) {
    sig_genes_cur <- rownames(res)[sig_idx]
    sig_genes_all_ids <- c(sig_genes_all_ids, sig_genes_cur)
  }
  sig_by_comp_ids[[comp_name]] <- sig_genes_cur

  ## Create subdirectory for this comparison
  outdir <- file.path(getwd(), comp_name)
  if (!dir.exists(outdir)) dir.create(outdir, recursive=TRUE, showWarnings=FALSE)

  ## Write list of significant genes for this comparison (display names)
  if (!is.null(display_map) && length(sig_genes_cur) > 0) {
    sig_genes_cur_disp <- ifelse(!is.na(display_map[sig_genes_cur]), display_map[sig_genes_cur], sig_genes_cur)
  } else {
    sig_genes_cur_disp <- sig_genes_cur
  }
  write.table(sig_genes_cur_disp, file=file.path(outdir, paste0("sig_genes_", comp_name, ".txt")), row.names=FALSE, col.names=FALSE, quote=FALSE)

  ## Add gene_id and display_name to results, and write results for this pair into its subdirectory
  resdata$gene_id <- resdata$gene
  if (!is.null(display_map)) {
    resdata$display_name <- ifelse(!is.na(display_map[resdata$gene_id]), display_map[resdata$gene_id], resdata$gene_id)
  } else {
    resdata$display_name <- resdata$gene_id
  }
  ## replace first column with display name for readability but keep gene_id
  resdata$gene <- resdata$display_name
  write.csv(resdata, file=file.path(outdir, paste0("deseq2-results_", comp_name, ".csv")), row.names=FALSE)

  ## p-value histogram
  png(filename=file.path(outdir, paste0("pvalue_hist_", comp_name, ".png")), width=800, height=600)
  hist(res$pvalue, breaks=50, col="grey", main=paste("P-value histogram:", comp_name))
  dev.off()

  ## MA plot
  png(filename=file.path(outdir, paste0("diffexpr-maplot_", comp_name, ".png")), 1500, 1000, pointsize=20)
  maplot(resdata, name=comp_name, main=paste("MA Plot", comp_name))
  dev.off()

  ## Volcano plot (PDF)
  pdf(file=file.path(outdir, paste0("diffexpr-volcanoplot_", comp_name, ".pdf")), width=12, height=10, pointsize=20)
  volcanoplot(resdata, name=comp_name, lfcthresh=1, sigthresh=0.01, textcx=.8, main=paste("Volcano", comp_name))
  dev.off()

  ## Per-comparison clustermap for genes significant in this comparison
  if (length(sig_genes_cur) > 0) {
    
    tryCatch({
      
      heat_mat_cur <- assay(rld)[sig_genes_cur, , drop=FALSE]
      heat_mat_cur_z <- t(scale(t(heat_mat_cur)))
      ## PDF output and red-black-green color palette (portrait/taller for readability)
      colormap <- colorRampPalette(c("red","black","green"))(100)
      ## Choose per-cell sizes in points (pt). These can be adjusted if needed.
      cellwidth_pt <- 12
      cellheight_pt <- 12
      ncols <- max(ncol(heat_mat_cur_z), 1)
      nrows <- max(nrow(heat_mat_cur_z), 1)
      ## Estimate rowname label width (pts) based on max characters and fontsize
      max_label_chars <- max(nchar(rownames(heat_mat_cur_z)))
      fontsize_row <- ceiling(cellheight_pt)
      label_width_pt <- max(60, min(2000, fontsize_row * max_label_chars * 0.6))
      ## Add margins (pts) for dendrograms, legend, etc.
      extra_width_pt <- 200
      extra_height_pt <- 200
      total_width_pt <- label_width_pt + (ncols * cellwidth_pt) + extra_width_pt
      total_height_pt <- (nrows * cellheight_pt) + extra_height_pt
      pdf_w <- max(4, total_width_pt / 72)
      pdf_h <- max(4, total_height_pt / 72)
      pdf(file.path(outdir, paste0("clustermap_sig_genes_", comp_name, ".pdf")),
          width=pdf_w, height=pdf_h)
      ## Map gene_ids to display names for row labels; if duplicates, append gene_id in parentheses
      gene_ids_cur <- rownames(heat_mat_cur_z)
      if (!is.null(display_map)) {
        disp_names_cur <- ifelse(!is.na(display_map[gene_ids_cur]),
                                 display_map[gene_ids_cur],
                                 gene_ids_cur)
      } else {
        disp_names_cur <- gene_ids_cur
      }
      dup <- duplicated(disp_names_cur) | duplicated(disp_names_cur, fromLast=TRUE)
      if (any(dup)) {
        disp_names_cur[dup] <- paste0(disp_names_cur[dup], " (", gene_ids_cur[dup], ")")
      }
      rownames(heat_mat_cur_z) <- disp_names_cur
      pheatmap::pheatmap(
        heat_mat_cur_z,
        annotation_col=as.data.frame(colData(rld)[, "condition", drop=FALSE]),
        show_rownames=TRUE,
        main=paste0("Clustermap: ", comp_name),
        color=colormap,
        cellwidth=cellwidth_pt,
        cellheight=cellheight_pt,
        fontsize_row=fontsize_row
      )
      dev.off()
      
    }, error = function(e) {
      message("Skipping heatmap for ", comp_name, " due to error: ", e$message)
    })
    
  } else {
    message("Skipping heatmap for ", comp_name, " (<2 significant genes)")
  }
}

sig_genes_all_ids <- unique(sig_genes_all_ids)
message("Total unique significant genes across comparisons: ", length(sig_genes_all_ids))

if (length(sig_genes_all_ids) == 0) {
  message("No genes passed the significance cutoffs (padj<", sig_padj, " & |LFC|>", sig_lfc, ")")
  
} else {
  
  heat_mat <- assay(rld)[sig_genes_all_ids, , drop=FALSE]
  
  ## ✅ Row scaling
  heat_mat_z <- t(scale(t(heat_mat)))
  
  ## ✅ REMOVE BAD ROWS (CRITICAL FIX)
  bad_rows <- !complete.cases(heat_mat_z)
  if (any(bad_rows)) {
    message("Removing ", sum(bad_rows), " rows with NA/NaN values before clustering")
    heat_mat_z <- heat_mat_z[!bad_rows, , drop=FALSE]
  }
  
  ## ✅ FINAL SAFETY CHECK
  if (nrow(heat_mat_z) < 2) {
    message("Skipping global clustermap: fewer than 2 valid genes after filtering")
    
  } else {
    
    gene_ids_all <- rownames(heat_mat_z)
    
    if (!is.null(display_map)) {
      disp_names_all <- ifelse(!is.na(display_map[gene_ids_all]),
                               display_map[gene_ids_all],
                               gene_ids_all)
    } else {
      disp_names_all <- gene_ids_all
    }
    
    dup_all <- duplicated(disp_names_all) | duplicated(disp_names_all, fromLast=TRUE)
    if (any(dup_all)) {
      disp_names_all[dup_all] <- paste0(disp_names_all[dup_all], " (", gene_ids_all[dup_all], ")")
    }
    
    annot <- as.data.frame(colData(rld)[, "condition", drop=FALSE])
    
    colormap <- colorRampPalette(c("red","black","green"))(100)
    
    cellwidth_pt <- 12
    cellheight_pt <- 12
    
    ncols_all <- max(ncol(heat_mat_z), 1)
    nrows_all <- max(nrow(heat_mat_z), 1)
    
    max_label_chars_all <- max(nchar(disp_names_all))
    fontsize_row_all <- ceiling(cellheight_pt)
    
    label_width_pt_all <- max(60, min(4000, fontsize_row_all * max_label_chars_all * 0.6))
    
    extra_width_pt <- 200
    extra_height_pt <- 200
    
    total_width_pt <- label_width_pt_all + (ncols_all * cellwidth_pt) + extra_width_pt
    total_height_pt <- (nrows_all * cellheight_pt) + extra_height_pt
    
    pdf_w <- max(4, total_width_pt / 72)
    pdf_h <- max(4, total_height_pt / 72)
    
    pdf_file <- "clustermap_sig_genes.pdf"
    
    ## ✅ SAFE PLOTTING
    tryCatch({
      
      pdf(pdf_file, width=pdf_w, height=pdf_h)
      
      rownames(heat_mat_z) <- disp_names_all
      
      pheatmap::pheatmap(
        heat_mat_z,
        annotation_col=annot,
        show_rownames=TRUE,
        main=paste0("Clustermap: genes (padj<", sig_padj, ", |LFC|>", sig_lfc, ")"),
        color=colormap,
        cellwidth=cellwidth_pt,
        cellheight=cellheight_pt,
        fontsize_row=fontsize_row_all
      )
      
      message("Clustermap successfully generated: ", pdf_file)
      
    }, error = function(e) {
      message("Failed to generate global clustermap: ", e$message)
      
    }, finally = {
      dev.off()
    })
  }
}
## Write master summary: which comparisons each gene (by gene_id) was significant in
if (length(sig_genes_all_ids) > 0) {
  comps <- names(sig_by_comp_ids)
  all_gene_ids <- sort(unique(unlist(sig_by_comp_ids)))
  summary_df <- data.frame(gene_id=all_gene_ids, stringsAsFactors=FALSE)
  if (!is.null(display_map)) {
    summary_df$display_name <- ifelse(!is.na(display_map[summary_df$gene_id]), display_map[summary_df$gene_id], summary_df$gene_id)
  } else {
    summary_df$display_name <- summary_df$gene_id
  }
  for (comp in comps) {
    summary_df[[comp]] <- summary_df$gene_id %in% sig_by_comp_ids[[comp]]
  }
  ## also add semicolon-separated list of comparisons per gene (using comp names)
  summary_df$comparisons <- apply(summary_df[ , comps, drop=FALSE], 1, function(row) paste(comps[which(row)], collapse=";"))
  write.csv(summary_df, file="genes_by_comparison_summary.csv", row.names=FALSE)
}


#Heat map for the most abundant BGCs

dds <- DESeqDataSetFromMatrix(countData=countdata, colData=coldata, design=~1)
dds
dds <- estimateSizeFactors(dds)
normalizeddata <- counts(dds, normalized=TRUE)
write.csv(normalizeddata, file="normalizeddata.csv")
library("pheatmap")
select <- order(rowMeans(counts(dds,normalized=TRUE)),
                decreasing=TRUE)[1:20]
nt <- normTransform(dds) # defaults to log2(x+1)
log2.norm.counts <- assay(nt)[select,]
df <- as.data.frame(colData(dds)[,c("condition")]) 
pheatmap(log2.norm.counts, cluster_rows=FALSE, show_rownames=TRUE, cluster_cols=FALSE)

