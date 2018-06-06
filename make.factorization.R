#!/usr/bin/env Rscript
argv <- commandArgs(trailingOnly = TRUE)

if(length(argv) < 6) {
    q()
}

ld.idx <- as.integer(argv[1])
ld.file <- argv[2]            # e.g., 'LD/fourier_ls-all.bed'
data.dir <- argv[3]           # e.g., 'gwas_data'
trait.names <- argv[4]        # e.g., 'AG:ALP:ALT:APTT'
K <- as.integer(argv[5])      # e.g., 10
out.hdr <- argv[6]            # e.g., 'temp'

## trait.names <- 'AG:ALP:ALT:APTT:AST:Alb:BS:BUN:Baso:CK:CRP:Ca:Cl:DBP:EA:EF:Eosino:FS:Fbg:GGT:HDL-C:Hb:HbA1c:Ht:IVS:K:LDH:LDL-C:LVDd:LVDs:LVM:LVMI:Lym:MAP:MCH:MCHC:MCV:Mono:NAP:Na:Neutro:P:PP:PT:PW:Plt:RBC:RWT:SBP:TBil:TC:TG:TP:UA:WBC:ZTT:eGFR:sCr'

options(stringsAsFactors = FALSE)
source('util.R')

dir.create(dirname(out.hdr), recursive = TRUE, showWarnings = FALSE)
conf.ind.file <- out.hdr %&&% '.conf-ind.gz'
conf.loading.file <- out.hdr %&&% '.conf-trait.gz'
conf.lodds.file <- out.hdr %&&% '.conf-lodds.gz'

.files <- c(conf.ind.file, conf.loading.file, conf.lodds.file)

if(all(sapply(.files, file.exists))) {
    log.msg('Files exists:\n%s\n', paste(.files, collapse = '\n'))
    q()
}

################################################################
source('util.fgwas.R')
library(dplyr)
library(readr)

ld.info <- read.ld.info(ld.idx, ld.file)

temp.dir <- system('mkdir -p /broad/hptmp/ypp/ukbb-fgwas/' %&&% out.hdr %&&%
                   '; mktemp -d /broad/hptmp/ypp/ukbb-fgwas/' %&&% out.hdr %&&%
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

traits <- strsplit(trait.names, split = ':') %>%
    (function(x) x[[1]])

gwas.files <- traits %>%
    (function(x) data.dir %&&% '/' %&&% ld.idx %&&% '_' %&&% x %&&% '.gz')

if(!all(file.exists(gwas.files))) {
    missing.files <- gwas.files %>%
        sapply(FUN = function(x) ifelse(file.exists(x), '', x)) %>%
            paste(collapse = '\n')
    log.msg('Missing GWAS files :\n%s\n', missing.files)
    q()
}

gwas.stat.tab <- gwas.files %>% read.gwas.files(plink.obj = plink)
if(is.null(gwas.stat.tab)) {
    log.msg('Failed to read GWAS files\n')
    q()
}

gwas.data <- construct.data.matrix(gwas.stat.tab, plink, traits)
if(is.null(gwas.data)) {
    log.msg('Failed to construct GWAS data\n')
    q()
}

################################################################
n <- nrow(gwas.data$X)
K <- max(min(c(length(traits) - 1, ncol(gwas.data$X) - 1, K, n)), 1)

## Just run factorization to check if there were any confounders
vb.opt <- list(pi.ub = -1, pi.lb = -3, tau = -5, do.hyper = TRUE,
               right.nn = FALSE, do.stdize = TRUE, do.rescale = TRUE,
               eigen.tol = 1e-2, gammax = 1e2, vbiter = 3500,
               svd.init = TRUE, jitter = 0.1,
               tol = 1e-8, rate = 1e-2, k = K)

z.out <- fit.zqtl.factorize(effect = gwas.data$beta, effect.se = gwas.data$se,
                            X = gwas.data$X, options = vb.opt)

log.msg('Finished factorization estimation\n\n')

colnames(z.out$param.indiv$theta) <- 'factor.' %&&% 1:K
Z.conf.ind.tab <- data.frame(iid = plink$FAM$iid, signif(z.out$param.indiv$theta, 4))

colnames(z.out$param.trait$theta) <- 'factor.' %&&% 1:K
Z.conf.loading.tab <- data.frame(traits, signif(z.out$param.trait$theta, 4))

conf.lodds <- signif(as.numeric(z.out$param.trait$lodds), 4)
Z.conf.lodds.tab <- data.frame(conf.factor = 1:K, conf.lodds)

write_tsv(Z.conf.loading.tab, path = conf.loading.file)
write_tsv(Z.conf.ind.tab, path = conf.ind.file)
write_tsv(Z.conf.lodds.tab, path = conf.lodds.file)

log.msg('Successfully finished fQTL\n\n')
