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


TODO: generate null data and save




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

TODO: generate the null and write them down

## use this: make.zqtl.null(X, beta.mat, se.mat, eig.tol, is.indep = FALSE, stdize = TRUE)
