#!/usr/bin/env Rscript
# permanova_validation.R
# Independent R cross-check of the primary 4-24 h PERMANOVA reported in the paper,
# using vegan::adonis2 with restricted (within-animal) permutations.
#
# Mirrors reproduce_all.py: family-level relative abundances of the 18 depth-passed
# 4-24 h lymph samples, after union (control-prevalence + genus-catalogue) decontamination.
# Permutation p-values are Monte Carlo and may differ in the last digit from the Python run.
#
# Usage: Rscript permanova_validation.R [path/to/otu_taxon.xls]
#   (defaults to data/otu_taxon.xls next to this script)
# Requires: vegan

suppressMessages(library(vegan))

args <- commandArgs(trailingOnly = TRUE)
here <- tryCatch(dirname(sub("^--file=", "",
        grep("^--file=", commandArgs(FALSE), value = TRUE))), error = function(e) ".")
otu_path <- if (length(args) >= 1) args[1] else file.path(here, "data", "otu_taxon.xls")
if (!file.exists(otu_path)) stop(sprintf("OTU table not found at: %s", otu_path))

df  <- read.delim(otu_path, row.names = 1, check.names = FALSE, stringsAsFactors = FALSE)
tax <- df[["taxonomy"]]
C   <- as.matrix(df[, setdiff(colnames(df), "taxonomy")])
storage.mode(C) <- "double"

rank_of <- function(t, r) {
  m <- regmatches(t, regexpr(paste0(r, "__([^;]+)"), t))
  if (length(m) == 0) "unclassified" else sub(paste0(r, "__"), "", m)
}
family <- vapply(tax, rank_of, "", r = "f")
genus  <- tolower(vapply(tax, rank_of, "", r = "g"))

kitome <- c('acetobacter','ralstonia','pseudomonas','acinetobacter','phyllobacterium','bradyrhizobium',
 'sphingomonas','sphingobium','methylobacterium','stenotrophomonas','delftia','herbaspirillum','burkholderia',
 'burkholderia-caballeronia-paraburkholderia','comamonas','cutibacterium','propionibacterium','rhodococcus',
 'microbacterium','brevundimonas','chryseobacterium','pseudoxanthomonas','massilia','rhizobium','agrobacterium',
 'variovorax','novosphingobium','afipia','caulobacter','pelomonas','methylophilus','oxalobacter','solirubrobacter','devosia')

cols   <- colnames(C); depth <- colSums(C)
ctrl   <- cols[grepl("^con_", cols) | cols == "water"]
fecal  <- cols[grepl("_FB$", cols)]; balf <- cols[grepl("_FP$", cols)]
later  <- cols[grepl("^(1|3|5|6|9)_(4|8|12|24)hLB$", cols)]
later_qc <- later[depth[later] >= 5000]
bio    <- c(later_qc, fecal, balf)
animal <- sub("_.*", "", later_qc)
timeph <- sub("^[0-9]+_([0-9]+)h.*", "\\1", later_qc)

# union decontamination: control-prevalence (one-sided Fisher p<0.05) OR genus catalogue
prev <- character(0)
for (o in rownames(C)) {
  cp <- sum(C[o, ctrl]  > 0); bp <- sum(C[o, bio] > 0)
  m  <- matrix(c(cp, length(ctrl) - cp, bp, length(bio) - bp), nrow = 2, byrow = TRUE)
  if (fisher.test(m, alternative = "greater")$p.value < 0.05) prev <- c(prev, o)
}
kit   <- rownames(C)[genus %in% kitome]
flag  <- union(prev, kit)
clean <- C[setdiff(rownames(C), flag), , drop = FALSE]
famcl <- tolower(family[rownames(clean)]); famcl <- family[rownames(clean)]

# family-level relative abundance for the 4-24 h lymph samples
fam_counts <- rowsum(clean[, later_qc, drop = FALSE], group = famcl)
rel <- sweep(fam_counts, 2, colSums(fam_counts), "/")
D   <- vegdist(t(rel), method = "bray")

meta <- data.frame(animal = factor(animal), timepoint = factor(timeph))
set.seed(1)
h   <- how(blocks = meta$animal, nperm = 1999)     # restrict permutations within animal
fit <- adonis2(D ~ timepoint, data = meta, permutations = h)
print(fit)
cat(sprintf("\nPERMANOVA (R/vegan) 4-24 h: pseudo-F=%.2f  R2=%.2f  p=%.3f\n",
            fit$F[1], fit$R2[1], fit$`Pr(>F)`[1]))
cat("Compare with reproduce_all.py: pseudo-F=1.21, R2=0.21, p=0.136 (Python).\n")
