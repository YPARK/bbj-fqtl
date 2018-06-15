#!/usr/bin/env Rscript
argv <- commandArgs(trailingOnly = TRUE)

if(length(argv) < 5) {
    q()
}

ld.idx <- as.integer(argv[1])
ld.file <- argv[2]                   # e.g., 'LD/fourier_ls-all.bed'
result.dir <- argv[3]                # e.g., 'temp'
lodds.cutoff <- as.numeric(argv[4])  # e.g., 0
out.file <- argv[5]                  # e.g., 'temp.txt.gz'

if(file.exists(out.file)){
    q()
}

options(stringsAsFactors = FALSE)
source('util.R')

result.hdr <- result.dir %&&% '/' %&&% ld.idx
snp.factor.file <- result.hdr %&&% '.snp-factor.gz'
trait.factor.file <- result.hdr %&&% '.trait-factor.gz'
zscore.file <- result.hdr %&&% '.zscore.gz'
resid.file <- result.hdr %&&% '.resid.gz'
var.file <- result.hdr %&&% '.var.gz'

.files <- c(snp.factor.file, trait.factor.file, zscore.file, resid.file, var.file)

if(!all(sapply(.files, file.exists))) {
    log.msg('Missing files exists:\n%s\n', paste(.files, collapse = '\n'))
    q()
}

################################################################
source('util.fgwas.R')
library(dplyr)
library(readr)
library(zqtl)

ld.info <- read.ld.info(ld.idx, ld.file)

temp.dir <- system('mkdir -p /broad/hptmp/ypp/bbj-fgwas/' %&&% result.hdr %&&%
                   '; mktemp -d /broad/hptmp/ypp/bbj-fgwas/' %&&% result.hdr %&&%
                   '/temp.XXXXXXXX',
                   intern = TRUE,
                   ignore.stderr = TRUE)
plink.hdr <- '1KG_EAS/chr' %&&% ld.info$chr.input
plink <- subset.plink(plink.hdr = plink.hdr,
                      chr = ld.info$chr.input,
                      plink.lb = ld.info$lb.input,
                      plink.ub = ld.info$ub.input,
                      temp.dir)
if(dir.exists(temp.dir)) system('rm -r ' %&&% temp.dir)

if(is.null(plink)) {
    log.msg('Failed to read plink\n')
    q()
}

cor.factor.gwas <- function(kk, snp.tab, gwas.tab, traits, DVt) {

    theta.k <- snp.tab %>%
        filter(factor == kk) %>%
            select(snp.loc, theta, theta.se, lodds) %>%
                mutate(theta.se = theta.se + 1e-4)

    gwas.tab.k <- gwas.tab %>%
        left_join(theta.k, by = 'snp.loc')

    Z.gwas <- gwas.tab.k %c% traits %>%
        as.matrix() %>%
            scale(scale = FALSE, center = TRUE)

    Z.factor <- gwas.tab.k %>%
        mutate(z.factor = theta / theta.se) %>%
            select(z.factor) %>%
                as.matrix()

    eta.factor <- DVt %*% Z.factor

    num.vec <- apply(sweep(Z.gwas, 1, Z.factor, `*`), 2, sum) %>%
        round(digits = 2)

    denom <- sqrt(sum(eta.factor^2))

    zz <- (num.vec / denom) %>% round(digits = 2)

    ret <- data.frame(trait = traits, dot = num.vec, z = zz, factor = kk)                      
}

snp.tab <- suppressMessages(read_tsv(snp.factor.file))

factors <- snp.tab %>%
    filter(lodds > lodds.cutoff) %>%
        select(factor) %>% unique() %>%
            unlist()

if(length(factors) < 1) {
    write_tsv(data.frame(), out.file)
    q()
}

trait.tab <- suppressMessages(read_tsv(trait.factor.file))
resid.tab <- suppressMessages(read_tsv(resid.file))
z.tab <- suppressMessages(read_tsv(zscore.file))
var.tab <- suppressMessages(read_tsv(var.file, col_types = 'iiiiccd'))

opt <- list(eigen.tol = 1e-2, do.stdize = TRUE)
x.svd.out <- plink$BED %c% z.tab$plink.pos %>%
    zqtl::take.ld.svd(options = opt)

Vt <- x.svd.out$V.t
D <- x.svd.out$D
DVt <- sweep(Vt, 1, D, `*`)

traits <- trait.tab$trait %>% unique()

out.z.tab <- lapply(factors, cor.factor.gwas,
                    snp.tab = snp.tab, gwas.tab = z.tab,
                    traits = traits, DVt = DVt) %>%
                        bind_rows()

out.resid.tab <- lapply(factors, cor.factor.gwas,
                        snp.tab = snp.tab, gwas.tab = resid.tab,
                        traits = traits, DVt = DVt) %>%
                            bind_rows()

.tab1 <- var.tab %>% filter(factor == 'total') %>% select(trait, var) %>% rename(var.g = var)
.tab2 <- var.tab %>% filter(factor == 'conf') %>% select(trait, var) %>% rename(var.c = var)
var.tot.tab <- left_join(.tab1, .tab2, by = 'trait')

temp.tab <- bind_rows(out.z.tab %>% mutate(gwas = 'obs'),
                      out.resid.tab %>% mutate(gwas = 'resid'))

temp.trait <- trait.tab %>%
    select(trait, factor, theta, theta.se, lodds) %>%
        mutate(theta = round(theta, 2),
               theta.se = round(theta.se, 2),
               lodds = round(lodds, 2))

out.tab <- var.tab %>% filter(factor %in% factors) %>%
    mutate(factor = as.integer(factor)) %>%
        right_join(temp.tab) %>%
            left_join(var.tot.tab) %>%
                left_join(temp.trait)

write_tsv(out.tab, out.file)
