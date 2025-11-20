suppressPackageStartupMessages(library("biomaRt"))

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

input_file <- parse_arg("--input", "results/edger/significant_genes.txt")
output_file <- parse_arg("--output", "results/edger/converted_genes.txt")

if (!file.exists(input_file)) {
    stop("Input gene list not found: ", input_file)
}

gene_list <- readLines(input_file)
gene_list <- unique(na.omit(gene_list))
if (!length(gene_list)) {
    stop("Input gene list is empty.")
}

mart <- useMart("fungi_mart", dataset = "scerevisiae_eg_gene")

converted <- getBM(
    attributes = c("external_gene_name", "systematic_name"),
    filters = "external_gene_name",
    values = gene_list,
    mart = mart
)

if (nrow(converted) == 0) {
    warning("No matches were returned by BioMart.")
} else {
    write.table(
        converted,
        file = output_file,
        row.names = FALSE,
        col.names = TRUE,
        quote = FALSE,
        sep = "\t"
    )
    cat("Converted gene list saved to:", output_file, "\n")
}
