#!/usr/bin/env Rscript

rm(list = ls(all = T))

library(qvalue)
library(data.table)
library(ggplot2)

### read in the results of the v-eQTL mapping
results.veqtl <- data.frame(fread(file = "veqtl.results.best", header = T))

### calculate P values adjusted for multiple testing across the genes, controlling the FDR and print results from paper (proportion of genes estimated to have a v-eQTL and number of significant v-eQTL)
qvalues.veqtl <- qvalue(results.veqtl$BETA)
print(c(1 - qvalues.veqtl$pi0, sum(qvalues.veqtl$qvalues < 0.05)))
results.veqtl <- data.frame(results.veqtl, qvalues = qvalue(results.veqtl$BETA)$qvalues)

### Print the most signficant association
head(results.veqtl[order(results.veqtl$BETA), ], 1)

### Extract gene and SNP information to plot an example v-eQTL later
expression.veqtl <- as.vector(unlist(fread("grep 'ENSG00000075234.12' geuvadis.linc.protein.cov.50PC.bed | cut -f5-")))
names(expression.veqtl) <- as.character(unlist(fread("head -1 geuvadis.linc.protein.cov.50PC.bed | cut -f5-", header = F)))

snp.veqtl <- data.frame(fread("tabix -h snps.vcf.gz 22:46616732-46616732 | grep -v '##' | cut -f10-"))

expression.veqtl <- expression.veqtl[names(expression.veqtl)%in%colnames(snp.veqtl)]
stopifnot(sum(names(expression.veqtl)!=colnames(snp.veqtl)) == 0)
snp.veqtl <- as.vector(unlist(snp.veqtl))

### code heterozygotes consistently
snp.veqtl <- ifelse(snp.veqtl == "1|0", "0|1", snp.veqtl)

snp.veqtl <- as.numeric(factor(snp.veqtl, levels = c("0|0", "0|1", "1|1")))

### This repeats previous steps but with variance effects suggesting parent of origin rather than additive
results.poeqtl <- data.frame(fread(file = "poveqtl.results.best", header = T))

qvalues.poeqtl <- qvalue(results.poeqtl$BETA)
print(c(1 - qvalues.poeqtl$pi0, sum(qvalues.poeqtl$qvalues < 0.05)))
results.poeqtl <- data.frame(results.poeqtl, qvalues = qvalue(results.poeqtl$BETA)$qvalues)

head(results.poeqtl[order(results.poeqtl$BETA), ])

expression.poeqtl <- as.vector(unlist(fread("grep 'ENSG00000196126.6' geuvadis.linc.protein.cov.50PC.bed | cut -f5-")))
names(expression.poeqtl) <- as.character(unlist(fread("head -1 geuvadis.linc.protein.cov.50PC.bed | cut -f5-", header = F)))

snp.poeqtl <- data.frame(fread("tabix -h snps.vcf.gz 6:32454727-32454727 | grep -v '##' | cut -f10-"))

expression.poeqtl <- expression.poeqtl[names(expression.poeqtl)%in%colnames(snp.poeqtl)]
stopifnot(sum(names(expression.poeqtl)!=colnames(snp.poeqtl)) == 0)
snp.poeqtl <- as.vector(unlist(snp.poeqtl))

snp.poeqtl <- ifelse(snp.poeqtl == "1|0", "0|1", snp.poeqtl)
### this ensures poeqtl genotypes have different values to veqtl genotypes, except T/T which exists for both SNPs
snp.poeqtl <- as.numeric(factor(snp.poeqtl, levels = c("0|0", "0|1", "1|1"))) - 2

Figure.1 <- data.frame(Expression = c(expression.veqtl, expression.poeqtl),
                       Genotypes = factor(c(snp.veqtl, snp.poeqtl), levels = -1:3, labels = c("TG/TG", "TG/T", "T/T", "T/C", "C/C")),
                       Type = factor(rep(c("v-eQTL", "Parent of Origin"), each = 445)), labels = c("v-eQTL", "Parent of Origin"))

ggplot(Figure.1, aes(x = Genotypes, y = Expression)) + geom_boxplot() + facet_wrap(~Type, scales = 'free') + theme_bw() + theme(axis.text = element_text(size = rel(2)), axis.title = element_text(size = rel(2)), strip.text = element_text(size = rel(2)))

### Store beta adjusted P values to produce qq plot
Supp.1 <- data.frame(Observed = -log10(c(sort(results.veqtl$BETA), sort(results.poeqtl$BETA))),
                     Expected = -log10(rep(seq(from = 1 / 14201, to = 14200 / 14201, by = 1 / 14201), 2)),
                     Type = factor(rep(c("v-eQTL", "Parent of Origin"), each = 14200)), labels = c("v-eQTL", "Parent of Origin"))

X11()

ggplot(Supp.1, aes(x = Expected, y = Observed)) + geom_point() + facet_wrap(~ Type) + geom_abline(intercept = 0, slope = 1) + theme_bw() + theme(axis.text = element_text(size = rel(2)), axis.title = element_text(size = rel(2)), strip.text = element_text(size = rel(2))) + labs(x = expression(-log[10] ~ "Expected P value"), y = expression(-log[10] ~ "Observed P value"))

Supp.2 <- data.frame(Genotypes = factor(snp.veqtl, labels = c("T/T", "T/C", "C/C")),
                     Expression = expression.veqtl,
                     Residuals = residuals(lm(expression.veqtl ~ factor(snp.veqtl))),
                     Distance = residuals(lm(expression.veqtl ~ factor(snp.veqtl)))^2)

panel.1 <- ggplot(Supp.2, aes(x = Genotypes, y = Expression)) + geom_boxplot() + theme_bw() + theme(axis.text = element_text(size = rel(2)), axis.title = element_text(size = rel(2))) + labs(title = "Expression vs Genotype")

panel.2 <- ggplot(Supp.2, aes(x = Genotypes, y = Residuals)) + geom_boxplot() + theme_bw() + theme(axis.text = element_text(size = rel(2)), axis.title = element_text(size = rel(2))) + labs(title = "Expression (removing eQTL effect) vs Genotype")

panel.3 <- ggplot(Supp.2, aes(x = Genotypes, y = Distance)) + geom_boxplot() + geom_smooth(aes(x = as.numeric(Genotypes)), method = 'lm', formula = y ~ x + I(x^2), se = FALSE) + theme_bw() + theme(axis.text = element_text(size = rel(2)), axis.title = element_text(size = rel(2))) + labs(title = "Distance vs Genotype")

## to have multiple panels in the same figure, use the multiplot function defined here: http://www.cookbook-r.com/Graphs/Multiple_graphs_on_one_page_(ggplot2)/

multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)

  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)

  numPlots = length(plots)

  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                    ncol = cols, nrow = ceiling(numPlots/cols))
  }

 if (numPlots==1) {
    print(plots[[1]])

  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))

    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))

      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}

X11()
multiplot(panel.1, panel.2, panel.3, cols = 3)
