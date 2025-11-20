suppressPackageStartupMessages({
    library("edgeR")
    library("ggplot2")
    library("dplyr")
    library("ggrepel")
    library("pheatmap")
    library("RColorBrewer")
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) %% 2 != 0) {
    stop("Arguments must be provided as --flag value pairs.")
}

parse_arg <- function(flag, default = NULL) {
    if (flag %in% args) {
        idx <- match(flag, args)
        if (idx == length(args)) stop(paste("Missing value for", flag))
        return(args[idx + 1])
    }
    default
}

counts_file <- parse_arg("--counts", "results/counts/gene_counts.txt")
metadata_file <- parse_arg("--metadata", "metadata/samples.csv")
output_dir <- parse_arg("--outdir", "results/edger")
lfc_cutoff <- as.numeric(parse_arg("--lfc", "1"))
pvalue_cutoff <- as.numeric(parse_arg("--pvalue", "0.05"))

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

metadata <- read.csv(metadata_file, stringsAsFactors = FALSE)
required_cols <- c("sample_id", "condition")
if (!all(required_cols %in% colnames(metadata))) {
    stop("Metadata file must contain 'sample_id' and 'condition' columns.")
}

counts_raw <- read.table(counts_file, header = TRUE, sep = "\t", check.names = FALSE, stringsAsFactors = FALSE)
meta_cols <- c("Chr", "Start", "End", "Strand", "Length")
sample_cols <- setdiff(colnames(counts_raw), c("Geneid", meta_cols))

clean_names <- function(x) {
    x <- gsub("^.*/", "", x)
    x <- gsub("_sorted\\.bam$", "", x)
    x <- gsub("\\.bam$", "", x)
    x <- gsub("\\.fastq$", "", x)
    x
}

sample_cols_clean <- clean_names(sample_cols)
counts_matrix <- counts_raw[, sample_cols, drop = FALSE]
colnames(counts_matrix) <- sample_cols_clean
rownames(counts_matrix) <- counts_raw$Geneid

metadata$sample_id <- as.character(metadata$sample_id)
metadata$condition <- as.character(metadata$condition)

missing_samples <- setdiff(metadata$sample_id, colnames(counts_matrix))
if (length(missing_samples) > 0) {
    stop(paste("Samples missing from count matrix:", paste(missing_samples, collapse = ", ")))
}

counts_matrix <- counts_matrix[, metadata$sample_id, drop = FALSE]

if (ncol(counts_matrix) < 2) {
    stop("Counts file must contain at least two samples.")
}

group <- factor(metadata$condition)

dge <- DGEList(counts = counts_matrix, group = group)
keep <- filterByExpr(dge, group = group)
dge <- dge[keep, , keep.lib.sizes = FALSE]
dge <- calcNormFactors(dge)
dge <- estimateDisp(dge)
logcpm <- cpm(dge, log = TRUE, prior.count = 1)

contrasts <- levels(group)
if (length(contrasts) != 2) {
    stop("Exactly two experimental conditions are required for this script.")
}

et <- exactTest(dge, pair = contrasts)
res <- topTags(et, n = Inf)$table %>% arrange(PValue)

res <- res %>%
    mutate(significant = case_when(
        PValue < pvalue_cutoff & logFC > lfc_cutoff ~ "Upregulated",
        PValue < pvalue_cutoff & logFC < -lfc_cutoff ~ "Downregulated",
        TRUE ~ "Not Significant"
    ))

write.csv(res, file = file.path(output_dir, "edger_results.csv"))

top_label_genes <- res %>%
    filter(significant != "Not Significant") %>%
    head(15)

volcano <- ggplot(res, aes(x = logFC, y = -log10(PValue), color = significant)) +
    geom_point(alpha = 0.6, size = 1.5) +
    geom_text_repel(
        data = top_label_genes,
        aes(label = rownames(top_label_genes)),
        size = 3.5,
        max.overlaps = 20
    ) +
    scale_color_manual(values = c("Not Significant" = "gray70",
                                  "Upregulated" = "#d73027",
                                  "Downregulated" = "#4575b4")) +
    geom_vline(xintercept = c(-lfc_cutoff, lfc_cutoff), col = "black", linetype = "dashed") +
    geom_hline(yintercept = -log10(pvalue_cutoff), col = "black", linetype = "dotted") +
    labs(
        title = sprintf("Differential Expression (%s vs %s)", contrasts[2], contrasts[1]),
        x = "Log2 Fold Change",
        y = "-Log10 P-Value",
        color = "Significance"
    ) +
    theme_minimal(base_size = 14)

ggsave(file.path(output_dir, "volcano_plot.png"), plot = volcano, width = 8, height = 6, dpi = 300)

# Principal Component Analysis (top 500 most variable genes)
row_vars <- apply(logcpm, 1, var)
top_var_genes <- names(sort(row_vars, decreasing = TRUE))[seq_len(min(length(row_vars), 500))]
pca_data <- prcomp(t(logcpm[top_var_genes, , drop = FALSE]), scale. = FALSE)
variance <- (pca_data$sdev^2) / sum(pca_data$sdev^2)

pca_df <- data.frame(
    sample = colnames(logcpm),
    PC1 = pca_data$x[, 1],
    PC2 = pca_data$x[, 2]
) %>%
    left_join(metadata, by = c("sample" = "sample_id"))

pca_plot <- ggplot(pca_df, aes(x = PC1, y = PC2, color = condition, label = sample)) +
    geom_point(size = 4) +
    geom_text_repel(size = 3.5, max.overlaps = 50) +
    labs(
        title = "PCA of logCPM expression (top 500 variable genes)",
        x = sprintf("PC1 (%.1f%%)", variance[1] * 100),
        y = sprintf("PC2 (%.1f%%)", variance[2] * 100),
        color = "Condition"
    ) +
    theme_minimal(base_size = 14)

ggsave(file.path(output_dir, "pca_plot.png"), plot = pca_plot, width = 7, height = 5, dpi = 300)

# Heatmap for top differentially expressed genes
heatmap_file <- file.path(output_dir, "top50_heatmap.png")
top_heatmap_genes <- head(rownames(res), 50)

if (length(top_heatmap_genes) >= 2) {
    heatmap_data <- logcpm[top_heatmap_genes, , drop = FALSE]
    annotation_col <- data.frame(
        condition = group,
        row.names = colnames(heatmap_data)
    )
    palette_length <- 50
    my_colors <- colorRampPalette(rev(brewer.pal(n = 9, name = "RdBu")))(palette_length)
    pheatmap(
        heatmap_data,
        scale = "row",
        show_rownames = TRUE,
        fontsize_row = 6,
        clustering_distance_rows = "correlation",
        clustering_distance_cols = "euclidean",
        annotation_col = annotation_col,
        color = my_colors,
        filename = heatmap_file,
        width = 6,
        height = 8
    )
} else {
    file.create(heatmap_file)
}

deg_summary <- res %>%
    group_by(significant) %>%
    summarise(count = n(), .groups = "drop")
write.csv(deg_summary, file = file.path(output_dir, "deg_summary.csv"), row.names = FALSE)

significant_genes <- res %>%
    filter(significant != "Not Significant") %>%
    rownames()
writeLines(significant_genes, file.path(output_dir, "significant_genes.txt"))

writeLines(capture.output(sessionInfo()), file.path(output_dir, "session_info.txt"))

message("edgeR analysis completed successfully.")
message("Results saved to: ", file.path(output_dir, "edger_results.csv"))
