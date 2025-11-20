SHELL := /bin/bash
PROJECT_ROOT := $(abspath $(dir $(lastword $(MAKEFILE_LIST))))
DATA_DIR := $(PROJECT_ROOT)/data
RESULTS_DIR := $(PROJECT_ROOT)/results
SAMPLESHEET := $(PROJECT_ROOT)/metadata/samples.csv

.PHONY: all help download-data prepare-reference docker-image run clean-results

all: help

help:
	@echo "Available targets:"
	@echo "  download-data     Download PRJNA184040 FASTQs listed in metadata/samples.csv"
	@echo "  prepare-reference Download sacCer3 genome + GTF and build HISAT2 index"
	@echo "  docker-image      Build the R/edgeR image used by edgeRAnalysis"
	@echo "  run               Execute the complete Nextflow pipeline"
	@echo "  clean-results     Remove generated results (keeps raw data)"

download-data:
	bash scripts/download_fastqs.sh

prepare-reference:
	bash scripts/prepare_reference.sh

docker-image:
	docker build -t yeast-oxidative-stress-analysis .

run:
	nextflow run src/main.nf -profile $${PROFILE:-conda} -resume \
		--data_dir $(DATA_DIR) \
		--results_dir $(RESULTS_DIR) \
		--samplesheet $(SAMPLESHEET)

clean-results:
	rm -rf $(RESULTS_DIR)/qc \
		$(RESULTS_DIR)/trimmed \
		$(RESULTS_DIR)/alignment \
		$(RESULTS_DIR)/counts \
		$(RESULTS_DIR)/edger
