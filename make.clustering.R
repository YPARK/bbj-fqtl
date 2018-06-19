#!/usr/bin/env Rscript

argv <- commandArgs(trailingOnly = TRUE)

if(length(argv) < 3) {
    q()
}

stat.file <- argv[1] # e.g., stat.file = 'result/20180613/fgwas-50-05.txt.gz'
K.max <- as.integer(argv[2]) # e.g., K.max = 100
out.hdr <- argv[3]

library(readr)
library(dplyr)
library(tidyr)
source('util.R')
library(ClusterR)

clust.file <- out.hdr %&&% '.clust.txt.gz'
snp.file <- out.hdr %&&% '.snp.txt.gz'
trait.file <- out.hdr %&&% '.trait.txt.gz'

.files <- c(snp.file, trait.file, clust.file)
if(all(sapply(.files, file.exists))) {
    log.msg('Files exist: %s\n', paste(.files, sep = ', '))
    q()
}

stat.tab <- read_tsv(stat.file)

.split <- function(...) {
    unlist(..., use.names = FALSE) %>%
        strsplit(split = '[|]') %>% (function(x) x[[1]])
}
.split.d <- function(...) .split(...) %>% as.numeric()
.split.i <- function(...) .split(...) %>% as.integer()

.split.trait <- function(tab) {
    data.frame(trait = .split(tab$trait),
               trait.theta = .split.d(tab$trait.theta),
               trait.se = .split.d(tab$trait.theta.se),
               trait.lodds = .split.d(tab$trait.lodds),
               gwas.6 = .split.d(tab$gwas.6),
               resid.6 = .split.d(tab$resid.6))
}

.split.snp <- function(tab) {
    data.frame(snp = .split(tab$snp),
               snp.theta = .split.d(tab$snp.theta),
               snp.theta.se = .split.d(tab$snp.theta.se),
               snp.best.z = .split.d(tab$snp.best.z),
               snp.best.trait = .split(tab$snp.best.trait),
               snp.loc = .split.i(tab$snp.loc))
}

trait.tab <- stat.tab %>%
    group_by(ld.idx, factor) %>%
    do(.split.trait(.)) %>%
    select(ld.idx, factor, trait, starts_with('trait'),
           starts_with('gwas'), starts_with('resid'))

trait.mat <- trait.tab %>%
    select(ld.idx, factor, trait, trait.lodds) %>%
    group_by(ld.idx, factor) %>%
    spread(key = trait, value = trait.lodds) %>%
    as.data.frame()

snp.tab <- stat.tab %>%
    select(ld.idx, factor, starts_with('snp')) %>%
    group_by(ld.idx, factor) %>%
    do(.split.snp(.))

M <- trait.mat[, -(1:2)] %>% as.matrix()
lodds.trunc <- log(0.9999) - log(0.0001)
M <- pmin(pmax(M, -lodds.trunc), lodds.trunc)

################################################################
## run clustering
K.max <- min(max(nrow(M) - 2, 2), K.max)
opt <- Optimal_Clusters_KMeans(M, max_clusters = K.max,
                               num_init = 3,
                               plot_clusters = FALSE,
                               verbose = TRUE,
                               criterion = 'BIC',
                               seed = 13)

K <- which.min(opt)

clust.out <- KMeans_rcpp(M, clusters = K,
                         num_init = 3,
                         max_iters = 200,
                         verbose = TRUE,
                         seed = 13)

################################################################
## reorder centroids
C <- clust.out$centroid
ko <- row.order(C > 0)
to.2 <- row.order(t(C) > 0)    # trait order for pair-wise map
to <- apply(C[ko, ] > 0, 2, function(x) median(which(x))) %>%
    order(decreasing = TRUE)

clust.order.tab <- 
    tibble(k = as.integer(ko)) %>%
    mutate(cluster = 1:n())

trait.order.tab <-
    tibble(trait = colnames(M),
           trait.o = colnames(M)[to],
           trait.o2 = colnames(M)[to.2])

colnames(clust.out$centroid) <- trait.order.tab$trait

factor2clust.tab <- trait.mat[, 1:2] %>%
    mutate(k = clust.out$clusters) %>%
    left_join(clust.order.tab)

centroid.tab <- cbind(k = 1:K, clust.out$centroid) %>%
    as.data.frame() %>%
    gather(key = trait, value = lodds, trait.order.tab$trait) %>%
    mutate(pip = 1/(1+exp(-lodds))) %>%
    mutate(trait = factor(trait, trait.order.tab$trait.o)) %>%
    left_join(clust.order.tab)    

################################################################
## reordered data
snp.membership <- snp.tab %>% 
    left_join(factor2clust.tab) %>%
    arrange(cluster) %>%
    select(-k)

trait.membership <- trait.tab %>%
    mutate(trait = factor(trait, trait.order.tab$trait.o)) %>%
    mutate(pip = 1/(1+exp(-trait.lodds))) %>%
    left_join(factor2clust.tab) %>%
    arrange(cluster) %>%
    select(-k) %>%
    mutate(col = ld.idx %&&% '_' %&&% factor) %>%
    as.data.frame()

factor.o <- trait.membership %>%
    rename(row = trait) %>%
    mutate(weight = sign(trait.theta) * pip) %>%
    select(row, col, weight) %>%
    order.pair()

trait.membership <-
    trait.membership %>%
    mutate(col = factor(col, factor.o$cols))

clust.stat <- 
    snp.membership %>% group_by(cluster) %>%
    summarize(n.snp = n()) %>%
    left_join(factor2clust.tab %>% group_by(cluster) %>%
              summarize(n.factor = n()))

trait.stat <- trait.tab %>%
    filter(trait.lodds > 0) %>%
    left_join(factor2clust.tab) %>%
    group_by(trait) %>%
    summarize(n.cluster = n())

################################################################
## Save clustering results

write_tsv(trait.membership, path = trait.file)
write_tsv(snp.membership, path = snp.file)
write_tsv(clust.stat, path = clust.file)

