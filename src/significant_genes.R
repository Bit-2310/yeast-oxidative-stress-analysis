suppressPackageStartupMessages(library("dplyr"))

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

results_file <- parse_arg("--input", "results/edger/edger_results.csv")
output_file <- parse_arg("--output", "results/edger/significant_genes.txt")
pvalue_cutoff <- as.numeric(parse_arg("--pvalue", "0.05"))
lfc_cutoff <- as.numeric(parse_arg("--logfc", "1"))

if (!file.exists(results_file)) {
    stop("Results file not found: ", results_file)
}

res <- read.csv(results_file, row.names = 1, check.names = FALSE)

significant_genes <- res %>%
    filter(PValue < pvalue_cutoff, abs(logFC) > lfc_cutoff) %>%
    rownames()

if (length(significant_genes) == 0) {
    warning("No genes meet the supplied thresholds. Creating an empty file.")
    file.create(output_file)
} else {
    writeLines(significant_genes, output_file)
}

cat("Significant gene list saved to:", output_file, "\n")
