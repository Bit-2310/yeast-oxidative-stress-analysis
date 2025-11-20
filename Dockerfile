FROM rocker/tidyverse:4.3.2

# Install dependencies
RUN apt-get update && apt-get install -y \
    fastqc \
    hisat2 \
    samtools \
    subread \
    wget \
    curl \
    && rm -rf /var/lib/apt/lists/*

# Install R libraries
RUN R -e "BiocManager::install(c('edgeR', 'ggplot2', 'dplyr', 'ggrepel', 'pheatmap', 'RColorBrewer'))"

# Set working directory
WORKDIR /project

CMD ["/bin/bash"]
