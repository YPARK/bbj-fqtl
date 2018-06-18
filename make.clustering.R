#!/usr/bin/env Rscript
## Identify clusters of tissue patterns

library(readr)
library(dplyr)
library(tidyr)
source('util.R')
library(ClusterR)

stat.file <- 'result/20180613/fgwas-50-09.txt.gz'
K.max <- 50

################################################################

stat.tab <- read_tsv(stat.file)

.split <- function(...) unlist(..., use.names = FALSE) %>%
    strsplit(split = '[|]') %>% (function(x) x[[1]])
.split.d <- function(...) .split(...) %>% as.numeric()
.split.i <- function(...) .split(...) %>% as.integer()

.split.trait <- function(tab) {
    data.frame(trait = .split(tab$trait),
               trait.theta = .split.d(tab$trait.theta),
               trait.se = .split.d(tab$trait.theta.se),
               trait.lodds = .split.d(tab$trait.lodds),
               gwas.6 = .split.d(tab$gwas.6))
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
    select(ld.idx, factor, trait, starts_with('trait'), starts_with('gwas'))

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

opt <- Optimal_Clusters_KMeans(M, max_clusters = K.max,
                               num_init = 3,
                               plot_clusters = FALSE,
                               verbose = TRUE,
                               criterion = 'BIC',
                               seed = 13)

K <- which.min(opt)

clust.out <- KMeans_rcpp(M, clusters = K, num_init = 3, max_iters = 200,
                         verbose = TRUE, seed = 13)

################################################################
## reorder centroids
C <- clust.out$centroid
ko <- row.order(C > 0)         # cluster order
to.2 <- row.order(t(C) > 0)    # trait order for pair-wise map
to <- order(apply(C[ko, ] > 0, 2, function(x) max(which(x))), decreasing = TRUE)

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

centroid.tab <- cbind(k = 1:K, clust.out$centroid) %>% as.data.frame() %>%
    gather(key = trait, value = lodds, trait.order.tab$trait) %>%
    mutate(pip = 1/(1+exp(-lodds))) %>%
    mutate(trait = factor(trait, trait.order.tab$trait.o)) %>%
    left_join(clust.order.tab)    

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
    snp.membership %>% group_by(cluster) %>% summarize(n.snp = n()) %>%
    left_join(factor2clust.tab %>% group_by(cluster) %>% summarize(n.factor = n())) %>%
    left_join(trait.membership %>% group_by(cluster) %>% summarize(n.trait = n()))

trait.stat <- trait.tab %>%
    filter(trait.lodds > 0) %>%
    left_join(factor2clust.tab) %>%
    group_by(trait) %>%
    summarize(n.cluster = n())

################################################################

.aes <- aes(x = col, y = trait, fill = pip)

plt <- gg.plot(trait.membership, .aes) +
    theme(axis.text.x = element_blank(), panel.spacing = unit(0.1, 'lines')) +
    geom_tile(color = 'gray') +
    facet_grid(.~cluster, space = 'free', scales = 'free') +
    scale_fill_gradientn('PIP', colors = c('#FFFFFF', '#FF0000'),
                         breaks = c(0, .05, .25, .5, .75, .95, 1))

plt <- plt + xlab('Factors') + ylab('Traits')

plt

.aes <- aes(x = col, y = trait, fill = gwas.6)

plt <- gg.plot(trait.membership, .aes) +
    theme(axis.text.x = element_blank(), panel.spacing = unit(0.1, 'lines')) +
    geom_tile(color = 'gray') +
    facet_grid(.~cluster, space = 'free', scales = 'free') +
    scale_fill_gradientn('# p < 1e-6',
                         colors = c('#FFFFFF', '#FFAA00', '#FF0000'),
                         breaks = c(1, 10, 100, 1000),
                         trans = 'log10', na.value = '#FFFFFF')

plt <- plt + xlab('Factors') + ylab('Traits')

plt

################################################################

TODO



################################################################
## save trait order
## write_tsv(data.frame(trait = traits.ordered) %>% mutate(trait.order = 1:n()),
##           path = gzfile('result/ukbb-fqtl-traits-order.txt.gz'))
##
## write_tsv(centroid.tab, path = gzfile(centroid.out.file))

################################################################
plt <-
    gg.plot(centroid.tab, aes(x = cluster, y = trait, fill = pip)) +
    geom_tile(color = 'gray') +
    scale_fill_gradientn('PIP', colors = c('#FFFFFF', '#FF0000'),
                         breaks = c(0, .05, .25, .5, .75, .95, 1))

plt <- plt +
    xlab('Clustered genomic regions') + ylab(nrow(trait.order.tab) %&&% ' traits') +
    theme(legend.background = element_rect(fill = 'white', color = 'white'),
          axis.text.x = element_text(size = 4))

plt

log.msg('Done\n')
