#!/usr/bin/env Rscript
argv <- commandArgs(trailingOnly = TRUE)

if(length(argv) < 3) {
    q()
}

ld.file <- argv[1]   # e.g., 'LD/fourier_ls-all.bed'
gwas.file <- argv[2] # e.g., 'BBJ_SUMMARY/BBJ.AG.autosome.txt.gz'
out.dir <- argv[3]   # e.g., 'gwas_data/'

options(stringsAsFactors = FALSE)
source('util.R')
library(dplyr)
library(readr)

dir.create(out.dir, recursive = TRUE, showWarnings = FALSE)

ld.info.tab <- read.ld.info(ld.idx = NULL, ld.file = ld.file)

rm.str <- function(x, pat) gsub(x, pattern = pat, replacement = '')

gwas.name <- function(x) {
    basename(x) %>%
        rm.str(pat = 'autosome[.]txt[.]gz') %>%
            rm.str(pat = 'BBJ[.]') %>%
                rm.str(pat = '[.]')
}

gwas.tab <- read_tsv(gwas.file) %>%
    select(SNP, CHR, POS, REF, ALT, Frq, BETA, SE, LOG10P)

ref.cols <- c('CHR', 'SNP', 'miss', 'POS', 'REF.1KG', 'ALT.1KG')
ref.types <- 'iciicc'
ref.var.tab <- '1KG_EAS/chr' %&&% 1:22 %&&% '.bim' %>%
    lapply(FUN = read_tsv, col_names = ref.cols, col_types = ref.types) %>%
        bind_rows()

valid.stat.tab <- ref.var.tab %>%
    left_join(gwas.tab %>% select(-SNP)) %>%
        na.omit() %>%
            mutate(BETA = if_else(REF == REF.1KG, BETA, -BETA)) %>%
                select(-REF, -ALT) %>%
                    rename(REF = REF.1KG, ALT = ALT.1KG) %>%
                        filter(Frq > 1e-2) %>%
                            select(CHR, SNP, POS, REF, ALT, BETA, SE, LOG10P)

rm(gwas.tab)
gc()

for(rr in 1:nrow(ld.info.tab)) {

    ld.idx <- ld.info.tab$ld.idx[rr]
    chr <- ld.info.tab$chr.input[rr]
    lb <- ld.info.tab$lb.input[rr]
    ub <- ld.info.tab$ub.input[rr]

    out.file <- out.dir %&&% '/' %&&% ld.idx %&&% '_' %&&% gwas.name(gwas.file) %&&% '.gz'

    if(!file.exists(out.file)) {
        valid.stat.tab %>%
            filter(CHR == chr, POS >= lb, POS <= ub) %>%
                write_tsv(path = out.file)
        gc()
        log.msg('Wrote %s file\n', out.file)
    }
}

