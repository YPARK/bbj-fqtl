#!/usr/bin/env Rscript

argv <- commandArgs(trailingOnly = TRUE)

if(length(argv) < 4) {
    q()
}

result.dir <- argv[1]               # e.g., result.dir = 'result/20180606/fgwas/50/'
ld.file <- argv[2]                  # e.g., 'LD/fourier_ls-all.bed'
pip.cutoff <- as.numeric(argv[3])   # e.g., pip.cutoff = 0.9
out.file <- argv[4]                 # e.g., out.file = 'temp'

source('util.R')
library(dplyr)
library(readr)
library(tidyr)

dir.create(dirname(out.file), recursive = TRUE, showWarnings = FALSE)

read.summary <- function(ld.idx, result.dir, lodds.cutoff) {

    .collapse <- function(...) paste(..., collapse = '|')

    snp.file <- result.dir %&&% '/' %&&% ld.idx %&&% '.snp-factor.gz'
    trait.file <- result.dir %&&% '/' %&&% ld.idx %&&% '.trait-factor.gz'
    z.file <- result.dir %&&% '/' %&&% ld.idx %&&% '.zscore.gz'
    var.file <- result.dir %&&% '/' %&&% ld.idx %&&% '.var.gz'
    snp.tab <- suppressMessages(read_tsv(snp.file))

    valid.factors <- snp.tab %>% group_by(factor) %>% slice(which.max(lodds)) %>%
        filter(lodds > lodds.cutoff) %>%
            select(factor) %>% unlist()

    if(length(valid.factors) > 0) {

        trait.tab <- suppressMessages(read_tsv(trait.file))
        z.tab <- suppressMessages(read_tsv(z.file))
        var.tab <- suppressMessages(read_tsv(var.file, col_types = 'iiiiiiccd'))

        var.tot <- var.tab %>% filter(factor == 'total')
        var.trait <- var.tab %>% filter(factor != 'total') %>%
            mutate(factor = as.integer(factor)) %>%
                filter(factor %in% valid.factors)

        var.summary <-
            var.trait %>% select(factor, trait, var)

        trait.summary <- trait.tab %>% filter(factor %in% valid.factors) %>%
            left_join(var.summary) %>%
                group_by(chr, ld.idx, LB, UB, factor) %>%
                    summarize(trait = .collapse(trait),
                              trait.theta = .collapse(theta),
                              trait.theta.se = .collapse(theta.se),
                              trait.lodds = .collapse(lodds),
                              var = .collapse(var))

        traits <- trait.tab %>% select(trait) %>% unique() %>% unlist()

        z.summary <- snp.tab %>%
            filter(lodds > lodds.cutoff) %>%
                select(snp.loc, factor) %>%
                    left_join(z.tab %>% select(-chr, -rs, -a1, -a2, -plink.pos)) %>%
                        gather_(key_col= 'snp.best.trait',
                                value_col = 'snp.best.z',
                                gather_col = traits)

        z.summary <- z.summary %>%
            group_by(snp.loc, factor) %>%
                slice(which.max(abs(snp.best.z))) %>%
                    select(starts_with('snp'), factor)

        snp.summary <- snp.tab %>% filter(lodds > lodds.cutoff) %>%
            left_join(z.summary) %>%
                group_by(chr, ld.idx, LB, UB, factor) %>%
                    summarize(snp = .collapse(rs),
                              snp.loc = .collapse(snp.loc),
                              snp.theta = .collapse(theta),
                              snp.theta.se = .collapse(theta.se),
                              snp.lodds = .collapse(lodds),
                              snp.best.trait = .collapse(snp.best.trait),
                              snp.best.z = .collapse(snp.best.z))

        tot.summary <- snp.summary %>% left_join(trait.summary)
        return(tot.summary %>% as.data.frame())
    } else {
        return(NULL)
    }
}

ld.info <- read.ld.info(ld.idx = NULL, ld.file)

out.tab <-
    ld.info$ld.idx %>%
        lapply(FUN = read.summary,
               result.dir = result.dir,
               lodds.cutoff = log(pip.cutoff) - log(1 - pip.cutoff)) %>%
                   bind_rows()

write_tsv(out.tab, out.file)
