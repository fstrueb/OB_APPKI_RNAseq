library(tidyverse)
library(rtracklayer)
library(edgeR)
mywd  <- '/home/fstruebi/projects/OB_RNAseq/' # substitute your folder with STAR counts
source('~/projects/ggplot_theme_FLS.R') # substitute your favorite ggplot2 theme

metad <- readRDS(paste0(mywd, 'resources/OB_RNAseq_metad.rds'))
fdrcutoff <-  0.1 # minimum cutoff for adjusted p-values to be considered significant
logfccutoff <- 0.5 # minimum cutoff for log2-fold changes to be considered significant

#### get mapping stats from FastQC ####
fqstats <- readRDS(paste0(mywd, 'resources/fqstats.rds'))
names(fqstats)[6:length(names(fqstats))] <- paste0('FastQC_', names(fqstats)[6:length(names(fqstats))])
fqstats <- fqstats %>% pivot_wider(names_from = 'read', values_from = starts_with('FastQC')) %>% 
    dplyr::filter(FLS_ID != 'Undetermined') %>% 
    mutate(FLS_id = as.numeric(FLS_ID)) %>% 
    dplyr::select(-where(is.character))

#### get counts from STAR output ####
allcounts <- lapply(seq(1:16), function(x) {
    read_delim(paste0(mywd, 'aligned/', x, '_ReadsPerGene.out.tab'), delim = '\t', col_names = FALSE) %>% 
        mutate(FLS_id = x)
}) %>% data.table::rbindlist(.) %>% 
    pivot_wider(id_cols = FLS_id, names_from = 'X1', values_from = 'X4') %>% 
    left_join(., fqstats, by = 'FLS_id') %>% 
    pivot_longer(-FLS_id) %>% 
    pivot_wider(id_cols = name, names_from = 'FLS_id', values_from = 'value')

#### get GFF for gene ID to symbol mapping ####
gtf <- import.gff('/earth/public_data/gtf/gencode.vM32.primary_assembly.annotation.gtf') %>% # substitute your gtf file
    plyranges::filter(type == 'gene') 
mycounts <- allcounts %>% 
    dplyr::filter(grepl('ENSMUSG', name))
## if you want to see the duplicated symbols run:
# setdiff(mycounts$name, gtf$gene_id) 
# table(duplicated(gtf$gene_name)) # number of duplicated gene names
# gtf[duplicated(gtf$gene_name)]
dupe_gene_symbols <- gtf[duplicated(gtf$gene_name)]$gene_name

#### replace non-unique/duplicated gene symbol names by appending ENSG ID ####
myanno <- as.data.frame(gtf) %>% 
    group_by(gene_name) %>% 
    mutate(dupe_id = row_number()) %>% 
    mutate(gene_symbol = ifelse(dupe_id > 1, paste0(gene_name, '_', gene_id), gene_name)) %>%
    dplyr::ungroup() %>% 
    dplyr::select(seqnames, start, end, width, strand, gene_id, gene_symbol, gene_type)
mycounts$name <- myanno[match(myanno$gene_id, mycounts$name),]$gene_symbol

#### make count matrix ####
count_mat <- matrix(as.matrix(mycounts[,-1]), dimnames = list(mycounts$name, names(mycounts[,-1])), ncol = 16)

#### make design matrix ####
mygroups <- factor(paste(metad$group, metad$sex, sep = '.'))
design <- model.matrix(~0+mygroups)
colnames(design) <- levels(mygroups)

#### run edgeR ####
# make DGEList:
dgelist <- DGEList(counts = count_mat, group = metad$group)
# filter genes with low expression:
genestokeep <- filterByExpr(dgelist, min.count = 5, min.prop = 0.5)
dgelist <- dgelist[genestokeep, , keep.lib.sizes = FALSE] %>% 
    calcNormFactors(., method = 'TMM') %>% 
    estimateDisp(., design = design)

# Apply GLM fit
dgefit <- glmQLFit(dgelist, design)
# make contrast
mycontrasts <- makeContrasts(
    KI.vs.WT = (KI.male + KI.female) - (WT.male + WT.female), levels = design)
# edgeR tables
groupdiffs <- glmQLFTest(dgefit, contrast = mycontrasts[,1])
res_groupdiffs <- topTags(groupdiffs, n = Inf) %>% as.data.frame() %>%
    rownames_to_column() %>% 
    dplyr::select(gene_symbol = rowname, logFC, FDR, PValue, logCPM) %>% 
    left_join(., myanno) %>% 
    mutate(comparison = colnames(mycontrasts)[1]) %>% 
    mutate(significant = ifelse(FDR < fdrcutoff & abs(logFC) > logfccutoff, 'yes', 'no'))

#### Volcano Plot ####
plt <- res_groupdiffs %>% 
    ggplot(aes(x = logFC, y = -log10(PValue), color = significant, text = gene_symbol)) +
    geom_point(alpha = 0.7) + 
    labs(title = 'DE: KI vs. WT') +
    theme_FLS()
plt

#### GO annotation ####
library(enrichR)
mydbs <- listEnrichrDbs() %>% 
    dplyr::filter(grepl('GO.*2023', libraryName))
edger_all_anno <- enrichr(res_groupdiffs, %>% dplyr::pull(gene_symbol), mydbs = mydbs$libraryName)
plotEnrich(edger_all_anno, y = 'Ratio', orderBy = 'Combined.Score', numChar = 60, title = 'GO annotations for DE genes (edgeR)')


