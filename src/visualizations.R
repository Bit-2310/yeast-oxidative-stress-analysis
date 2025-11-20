suppressPackageStartupMessages({
    library("ggplot2")
    library("dplyr")
    library("readr")
})

args <- commandArgs(trailingOnly = TRUE)

parse_arg <- function(flag, default = NULL) {
    if (flag %in% args) {
        idx <- match(flag, args)
        if (idx == length(args)) stop(paste("Missing value for", flag))
        return(args[idx + 1])
    }
    default
}

results_file <- parse_arg("--results", "results/edger/edger_results.csv")
output_dir <- parse_arg("--outdir", "results/edger/visualizations")

if (!file.exists(results_file)) {
    stop("Differential expression table not found: ", results_file)
}

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

de <- read_csv(results_file, show_col_types = FALSE)

if (!all(c("logFC", "logCPM", "PValue", "significant") %in% colnames(de))) {
    stop("The results file must contain logFC, logCPM, PValue, and significant columns.")
}

ma_plot <- ggplot(de, aes(x = logCPM, y = logFC, color = significant)) +
    geom_point(alpha = 0.6, size = 1.4) +
    scale_color_manual(values = c(
        "Upregulated" = "#d73027",
        "Downregulated" = "#4575b4",
        "Not Significant" = "gray70"
    )) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "black", size = 0.3) +
    labs(
        title = "MA Plot",
        x = "logCPM",
        y = "log2 Fold Change",
        color = "Category"
    ) +
    theme_minimal(base_size = 14)

ggsave(file.path(output_dir, "ma_plot.png"), ma_plot, width = 7, height = 5, dpi = 300)

top_up <- de %>%
    filter(significant == "Upregulated", PValue < 0.05) %>%
    arrange(PValue) %>%
    slice_head(n = 15)
top_down <- de %>%
    filter(significant == "Downregulated", PValue < 0.05) %>%
    arrange(PValue) %>%
    slice_head(n = 15)

top_genes <- bind_rows(top_up, top_down) %>%
    mutate(direction = ifelse(logFC > 0, "Up", "Down")) %>%
    mutate(gene = factor(row_number(), labels = .$`...1`))

if (nrow(top_genes) > 0) {
    top_genes <- top_genes %>% mutate(gene_label = coalesce(`...1`, gene))
    ordered <- top_genes %>% arrange(direction, ifelse(direction == "Up", -logFC, logFC))
    ordered$gene_label <- factor(ordered$gene_label, levels = ordered$gene_label)

    bar_plot <- ggplot(ordered, aes(x = gene_label, y = logFC, fill = direction)) +
        geom_col() +
        coord_flip() +
        scale_fill_manual(values = c("Up" = "#d73027", "Down" = "#4575b4")) +
        labs(
            title = "Top Differentially Expressed Genes",
            x = NULL,
            y = "log2 Fold Change",
            fill = "Direction"
        ) +
        theme_minimal(base_size = 13)

    ggsave(file.path(output_dir, "top_de_bar.png"), bar_plot, width = 7, height = 6, dpi = 300)
} else {
    warning("No genes met the p-value threshold for the bar plot.")
}

message("Additional visualizations saved to ", output_dir)
