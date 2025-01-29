# Connor Murray
# Started 12.2.2024
# analyzing TOPchef eQTLs and perform colocalization
# module load gcc/11.4.0 openmpi/4.1.4 R/4.3.1; R

# Libraries
library(data.table)
library(coloc)
library(tidyverse)
library(locuscomparer)
library(foreach)

# Working directory 
setwd("/standard/vol185/cphg_Manichaikul/users/csm6hg/nextflow_dna/")

### Datasets & Setup ###

# Coloc files
coloc <- list.files(path = "output/coloc/", pattern = ".txt", full.names = T)

# Coloc run - levin et al 2022
dt1 <- rbindlist(lapply(coloc, function(t) {fread(t, header = T)}))
