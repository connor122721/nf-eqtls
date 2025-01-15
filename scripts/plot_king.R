#!/usr/bin/env Rscript
# Connor Murray
# Started 11.10.2024

# Libraries
library(tidyverse)
library(data.table)

# Command line arguments
args <- commandArgs(trailingOnly = TRUE)
king_input <- args[1]
metadata_file <- args[2]

# Metadata
meta <- data.table(fread(metadata_file, header = TRUE))

# Kinship
kin <- data.table(fread(king_input, header = T))

# Add relatedness boundaries
kin <- data.table(kin %>% 
            mutate(rel=case_when(Kinship >= 0.3536 ~ "Full-sibling",
                             Kinship >= 0.1768 & Kinship < 0.3536 ~ "1st-degree",
                             Kinship >= 0.0884 & Kinship < 0.1768 ~ "2nd-degree",
                             Kinship >= 0.0442 & Kinship < 0.0884 ~ "3rd-degree",
                             Kinship >= 0.0221 & Kinship < 0.0442 ~ "4th-degree",
                             TRUE ~ "Unrelated")))

# Reorder
kin$rel <- factor(kin$rel, levels=c("Full-sibling", 
                                    "1st-degree", 
                                    "2nd-degree", 
                                    "3rd-degree", 
                                    "4th-degree",
                                    "Unrelated"))
table(kin$rel)

# Common theme element
themei <- {
  theme_bw() + 
    theme(axis.title.x = element_text(face = "bold", size = 16),
          axis.text.x = element_text(face = "bold", size = 14),
          axis.title.y = element_text(face = "bold", size = 16),
          axis.text.y = element_text(face = "bold", size = 14),
          legend.text = element_text(face = "bold", size = 14),
          legend.title = element_text(face = "bold", size = 16)) 
}

# Plot kinship and IBS0
kinPlot <- {
  kin %>% 
    ggplot(., 
            aes(x=IBS0,
                y=Kinship,
                color=rel)) +
    geom_point(alpha=0.6) +
    labs(x="Identity-by-state", y="Kinship", color="") +
    themei +
    theme(aspect.ratio = 1)
}

# Output
ggsave(plot = kinPlot, filename = "kinship.pdf")

# Remove related individuals
fullSibs <- kin[rel=="Full-sibling"]$ID1
firstDegree <- kin[rel=="1st-degree"]$ID1
table(kin[!ID1 %in% c(fullSibs, firstDegree)]$rel)

# Related Individuals
write.table(x = unique(fullSibs, firstDegree), 
            file = "relatedIndividuals.txt", 
            quote = F, 
            sep = "\t", 
            row.names = F, 
            col.names = F)
