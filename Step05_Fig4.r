## BrainSpan
# ------------
library(clusterProfiler)
library(ClusterGVis)
library(ComplexHeatmap)
library(org.Hs.eg.db)
library(org.Mm.eg.db)
library(colorspace)
library(cols4all)
library(MASS)
library(class)
library(cluster)
library(impute)
library(Hmisc)
library(dplyr)
library(RColorBrewer)
library(ggnewscale)
library(cowplot)
options(stringsAsFactors = F)
dat = read.csv('/n02dat01/users/ypwang/Gradient/GeneticGradient/Data/GCIsig/BS_cere_891*34.csv',
      row.names=1,header=TRUE, check.names=FALSE)
colnames(dat)
dat_num = aggregate(dat, by=list(type=dat$age),mean)
dat_num = dat_num[order(dat_num$log_age_weeks),]
df = t(dplyr::select(dat_num, -c( "type","log_age_weeks","age","age_in_weeks" )))
colnames(df) = dat_num$log_age_weeks
rownames(df)
df = df[which(rowSums(df) > 0),]

RdBu = read.csv('/n02dat01/users/ypwang/Gradient/GeneticGradient/Fig/Color/RdBu_Hex_200.csv', row.names = 1, header=TRUE, check.names=FALSE)
output_dir = "/n02dat01/users/ypwang/Gradient/GeneticGradient/Data/GCIsig/Temporal/"
setwd(output_dir)
library(ggplot2)
mytheme = theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(), panel.border = element_blank(),
                panel.background = element_blank(), axis.line = element_line(colour = "black"),
                axis.title.x = element_text(size = 18,  face ="plain"),
                axis.title.y = element_text(size = 18,  face ="plain"),
                axis.text.x  = element_text(size = 15,  face ="plain", angle = 45, hjust = 1),
                axis.text.y  = element_text(size = 15,  face ="plain"),
                plot.margin = unit(c(0.2, 0.2, 0.2, 0.2), "inches"),
                plot.title = element_text(size = 20,  face ="bold", hjust = 0.5),
                axis.line.x  = element_line(size = 0.8),
                axis.line.y  = element_line(size = 0.8),
                legend.position = 'right', legend.text = element_text(size = 15,  face ="plain"),
                legend.title = element_blank()) 
pdf('Determinative_cluster.pdf', height = 9, width = 20, onefile = F)
# print(getClusters(exp = t(scale(t(df),scale=T)))) # apply(df,1, scale)
print(factoextra::fviz_nbclust(x = t(scale(t(df),scale=T)), stats::kmeans, method = 'wss', print.summary=False)
# +factoextra::fviz_nbclust(x = t(scale(t(df),scale=T)), stats::kmeans, method = "gap_stat",print.summary=False ) 
+mytheme
+factoextra::fviz_nbclust(x = t(scale(t(df),scale=T)), stats::kmeans, method = 'silhouette',print.summary=False)
+mytheme) 
# library(cluster)
# print(factoextra::fviz_nbclust(x = t(scale(t(df),scale=T)),  
# cluster::clara, method = 'silhouette', print.summary=T))
# a = t(exp)
# b = t(scale(t(df),scale=T))
# c = b[!is.na(b)]
dev.off()

# 2 cluster
j = 'mfuzz'
ck <- clusterData(exp = df, cluster.method = j,cluster.num = 2,seed = 123)
enrich.data = ck$wide.res

visCluster(object = ck,plot.type = "line", add.mline = TRUE)
seed =123
fromType = "SYMBOL"
toType = "ENTREZID"
OrgDb = org.Hs.eg.db
organism = "hsa"
pvalueCutoff = 0.05
topn = 10
readable = TRUE

x = 1
tmp <- enrich.data %>% dplyr::filter(cluster == unique(enrich.data$cluster)[x])
gene.ent <- clusterProfiler::bitr(tmp$gene,
                                    fromType = "SYMBOL",
                                    toType = c("ENTREZID"),
                                    OrgDb = org.Hs.eg.db)
tartget.gene = unlist(gene.ent[, "ENTREZID"])
ego <- clusterProfiler::enrichGO(gene       = tartget.gene,
                              keyType       = toType,
                              OrgDb         = OrgDb,
                              ont           = "ALL",
                              pAdjustMethod = "BH",
                              pvalueCutoff  = 1,
                              qvalueCutoff  = 0.2,
                              readable      = readable)
ekegg <- clusterProfiler::enrichKEGG(gene   = tartget.gene,
                              keyType       = "kegg",
                              organism      = "hsa",
                              universe      = NULL,
                              pvalueCutoff  = 0.2,
                              pAdjustMethod = "BH",
                              qvalueCutoff  = 1)
df_go <- data.frame(clusterProfiler::setReadable(ego,OrgDb = OrgDb,keyType = toType)) %>%
dplyr::filter(pvalue < pvalueCutoff) %>%
dplyr::mutate(group = paste("C",unique(enrich.data$cluster)[x],sep = '')) %>%
dplyr::arrange(pvalue)
df_kegg <- data.frame(clusterProfiler::setReadable(ekegg,OrgDb = OrgDb,keyType = toType)) %>%
dplyr::filter(pvalue < pvalueCutoff) %>%
dplyr::mutate(group = paste("C",unique(enrich.data$cluster)[x],sep = '')) %>%
dplyr::arrange(pvalue) 
df_1 = rbind(df_go[,colnames(df_kegg)], df_kegg)%>%dplyr::arrange(pvalue)
df_1_top = df_1 %>% dplyr::slice_head(n = topn)
df_1_top = df_1_top[,c('group','Description','pvalue','GeneRatio')] 

x = 2
tmp <- enrich.data %>% dplyr::filter(cluster == unique(enrich.data$cluster)[x])
gene.ent <- clusterProfiler::bitr(tmp$gene,
                                    fromType = "SYMBOL",
                                    toType = c("ENTREZID"),
                                    OrgDb = org.Hs.eg.db)
tartget.gene = unlist(gene.ent[, "ENTREZID"])
ego <- clusterProfiler::enrichGO(gene       = tartget.gene,
                              keyType       = toType,
                              OrgDb         = OrgDb,
                              ont           = "ALL",
                              pAdjustMethod = "BH",
                              pvalueCutoff  = 1,
                              qvalueCutoff  = 0.2,
                              readable      = readable)
ekegg <- clusterProfiler::enrichKEGG(gene   = tartget.gene,
                              keyType       = "kegg",
                              organism      = "hsa",
                              universe      = NULL,
                              pvalueCutoff  = 0.2,
                              pAdjustMethod = "BH",
                              qvalueCutoff  = 1)

df_go <- data.frame(clusterProfiler::setReadable(ego,OrgDb = OrgDb,keyType = toType)) %>%
dplyr::filter(pvalue < pvalueCutoff) %>%
dplyr::mutate(group = paste("C",unique(enrich.data$cluster)[x],sep = '')) %>%
dplyr::arrange(pvalue)
df_kegg <- data.frame(clusterProfiler::setReadable(ekegg,OrgDb = OrgDb,keyType = toType)) %>%
dplyr::filter(pvalue < pvalueCutoff) %>%
dplyr::mutate(group = paste("C",unique(enrich.data$cluster)[x],sep = '')) %>%
dplyr::arrange(pvalue) 
      
df_2 = rbind(df_go[,colnames(df_kegg)], df_kegg)%>%dplyr::arrange(pvalue)
df_2_top = df_2 %>% dplyr::slice_head(n = topn)
df_2_top = df_2_top[,c('group','Description','pvalue','GeneRatio')] 

enrich.all = rbind(df_1_top,df_2_top)
write.csv(enrich.all, file = "enrich_cluster2_20231129.csv")

annt_top = HeatmapAnnotation(age = anno_text(dat_num$type),
                             log_age_weeks = anno_barplot(dat_num$log_age_weeks)
                              )
mycol = c('#02274a','#832404' ) # ,'#8ecae6''#eb9172', 

pdf(paste0('Heatmap_cluster2.pdf'), height = 9, width = 12, onefile = F)
visCluster_ping(object = ck, plot.type = "both",  # line.side = "left", 
            ht.col = "RdBu", mline.col = '#808080',
            row_dend_gp = gpar(col = "#636363"),
            HeatmapAnnotation = annt_top, # column_names_rot = 60,  markGenes = mark, markGenes.side = "left", genes.gp = c('italic',fontsize = 10), 
            column_names_gp = gpar(fontsize = 0), 
            show_heatmap_legend = F,
            sample.group = rep(c("C1","C2","C3","C4","C5"),c(10,3,4,5,6)), # c("fetal","infant","child","adolescent","adult")
            column_title = c("Fetal","Infant","Child","Adolescent","Adult"),
            show_row_dend = T, annoTerm.data = enrich.all, go.size = "pval" , # add.bar = T,
            ctAnno.col = mycol, # line.col = 'grey',# ggsci::pal_d3()(5),line.col = mycol,
            go.col = c(rep(mycol[1],10),rep(mycol[2],10))
            )   go.size = "pval"
dev.off()   

j = 'mfuzz'
expression = read.csv(paste0(data_dir,'expression_4mm_317_cere2extra.csv'), header=TRUE, check.names=FALSE)
report = read.csv(paste0(data_dir,'report_4mm_317_cere2extra.csv'), header=TRUE, check.names=FALSE)

ind = 'pos'
for (ind in c('pos','neg')){
      print(ind)
      GCI_cur = get(paste0('GCI_',ind))
      expr_cur = expression[,GCI_cur$gene] # "PROS1" 
      R_cur = lapply(colnames(expr_cur), function(x){
            cor.test(report$gradient1, expr_cur[,x])$estimate})
      assign(paste0('R_', ind), R_cur)
      print(length(R_cur[R_cur>0]))
}
length(R_pos[R_pos>0])/length(R_pos)  # 0.9180672
length(R_neg[R_neg<0])/length(R_neg)  # 0.9288321

ck <- clusterData(exp = df, cluster.method = j,cluster.num = 2, seed = 123)
enrich.data = ck$wide.res
tmp <- enrich.data %>% dplyr::filter(cluster == unique(enrich.data$cluster)[1])
print(length(intersect(tmp$gene, GCI_neg$gene))/length(tmp$gene))
# 238:243  238/481
tmp <- enrich.data %>% dplyr::filter(cluster == unique(enrich.data$cluster)[2]) # 0.4948025
print(length(intersect(tmp$gene, GCI_pos$gene))/length(tmp$gene)) # 0.4327628
# 232:177 177/409
length(intersect(rownames(df), GCI_pos$gene)) # 420
length(intersect(rownames(df), GCI_neg$gene)) # 470
for (c_num in c(2,5)){
      ck <- clusterData(exp = df, cluster.method = j,cluster.num = c_num,seed = 123)
      enrich.data = ck$wide.res
      for (x in 1:c_num){
            print(x)
            tmp <- enrich.data %>% dplyr::filter(cluster == unique(enrich.data$cluster)[x])
            print(length(intersect(tmp$gene, GCI_pos$gene))/length(tmp$gene))
      }
}

# 8 cluster
j = 'mfuzz'
ck <- clusterData(exp = df, cluster.method = j,cluster.num = 8,seed = 123)
enrich.data = ck$wide.res

visCluster(object = ck,plot.type = "line", add.mline = TRUE)
seed =123
fromType = "SYMBOL"
toType = "ENTREZID"
OrgDb = org.Hs.eg.db
organism = "hsa"
pvalueCutoff = 0.05
topn = 5
readable = TRUE

x = 1
tmp <- enrich.data %>% dplyr::filter(cluster == unique(enrich.data$cluster)[x])
gene.ent <- clusterProfiler::bitr(tmp$gene,
                                    fromType = "SYMBOL",
                                    toType = c("ENTREZID"),
                                    OrgDb = org.Hs.eg.db)
tartget.gene = unlist(gene.ent[, "ENTREZID"])
ego <- clusterProfiler::enrichGO(gene       = tartget.gene,
                              keyType       = toType,
                              OrgDb         = OrgDb,
                              ont           = "ALL",
                              pAdjustMethod = "BH",
                              pvalueCutoff  = 1,
                              qvalueCutoff  = 0.2,
                              readable      = readable)
ekegg <- clusterProfiler::enrichKEGG(gene   = tartget.gene,
                              keyType       = "kegg",
                              organism      = "hsa",
                              universe      = NULL,
                              pvalueCutoff  = 0.2,
                              pAdjustMethod = "BH",
                              qvalueCutoff  = 1)
df_go <- data.frame(clusterProfiler::setReadable(ego,OrgDb = OrgDb,keyType = toType)) %>%
dplyr::filter(pvalue < pvalueCutoff) %>%
dplyr::mutate(group = paste("C",unique(enrich.data$cluster)[x],sep = '')) %>%
dplyr::arrange(pvalue)
df_kegg <- data.frame(clusterProfiler::setReadable(ekegg,OrgDb = OrgDb,keyType = toType)) %>%
dplyr::filter(pvalue < pvalueCutoff) %>%
dplyr::mutate(group = paste("C",unique(enrich.data$cluster)[x],sep = '')) %>%
dplyr::arrange(pvalue) 
df_1 = rbind(df_go[,colnames(df_kegg)], df_kegg)%>%dplyr::arrange(pvalue)
# df_1_top = df_1 %>% dplyr::slice_head(n = topn)
df_1_top = df_1[1:4,c('group','Description','pvalue','GeneRatio')] 

x = 2
tmp <- enrich.data %>% dplyr::filter(cluster == unique(enrich.data$cluster)[x])
gene.ent <- clusterProfiler::bitr(tmp$gene,
                                    fromType = "SYMBOL",
                                    toType = c("ENTREZID"),
                                    OrgDb = org.Hs.eg.db)
tartget.gene = unlist(gene.ent[, "ENTREZID"])
ego <- clusterProfiler::enrichGO(gene       = tartget.gene,
                              keyType       = toType,
                              OrgDb         = OrgDb,
                              ont           = "ALL",
                              pAdjustMethod = "BH",
                              pvalueCutoff  = 1,
                              qvalueCutoff  = 0.2,
                              readable      = readable)
ekegg <- clusterProfiler::enrichKEGG(gene   = tartget.gene,
                              keyType       = "kegg",
                              organism      = "hsa",
                              universe      = NULL,
                              pvalueCutoff  = 0.2,
                              pAdjustMethod = "BH",
                              qvalueCutoff  = 1)

df_go <- data.frame(clusterProfiler::setReadable(ego,OrgDb = OrgDb,keyType = toType)) %>%
dplyr::filter(pvalue < pvalueCutoff) %>%
dplyr::mutate(group = paste("C",unique(enrich.data$cluster)[x],sep = '')) %>%
dplyr::arrange(pvalue)
df_kegg <- data.frame(clusterProfiler::setReadable(ekegg,OrgDb = OrgDb,keyType = toType)) %>%
dplyr::filter(pvalue < pvalueCutoff) %>%
dplyr::mutate(group = paste("C",unique(enrich.data$cluster)[x],sep = '')) %>%
dplyr::arrange(pvalue) 
      
df_2 = rbind(df_go[,colnames(df_kegg)], df_kegg)%>%dplyr::arrange(pvalue)
df_2_top = df_2[1:4,c('group','Description','pvalue','GeneRatio')] 

x = 3
tmp <- enrich.data %>% dplyr::filter(cluster == unique(enrich.data$cluster)[x])
gene.ent <- clusterProfiler::bitr(tmp$gene,
                                    fromType = "SYMBOL",
                                    toType = c("ENTREZID"),
                                    OrgDb = org.Hs.eg.db)
tartget.gene = unlist(gene.ent[, "ENTREZID"])
ego <- clusterProfiler::enrichGO(gene       = tartget.gene,
                              keyType       = toType,
                              OrgDb         = OrgDb,
                              ont           = "ALL",
                              pAdjustMethod = "BH",
                              pvalueCutoff  = 1,
                              qvalueCutoff  = 0.2,
                              readable      = readable)
ekegg <- clusterProfiler::enrichKEGG(gene   = tartget.gene,
                              keyType       = "kegg",
                              organism      = "hsa",
                              universe      = NULL,
                              pvalueCutoff  = 0.2,
                              pAdjustMethod = "BH",
                              qvalueCutoff  = 1)

df_go <- data.frame(clusterProfiler::setReadable(ego,OrgDb = OrgDb,keyType = toType)) %>%
dplyr::filter(pvalue < pvalueCutoff) %>%
dplyr::mutate(group = paste("C",unique(enrich.data$cluster)[x],sep = '')) %>%
dplyr::arrange(pvalue)
df_kegg <- data.frame(clusterProfiler::setReadable(ekegg,OrgDb = OrgDb,keyType = toType)) %>%
dplyr::filter(pvalue < pvalueCutoff) %>%
dplyr::mutate(group = paste("C",unique(enrich.data$cluster)[x],sep = '')) %>%
dplyr::arrange(pvalue) 
      
df_3 = rbind(df_go[,colnames(df_kegg)], df_kegg)%>%dplyr::arrange(pvalue)
df_3_top = df_3[1:4,c('group','Description','pvalue','GeneRatio')] 

x = 4
tmp <- enrich.data %>% dplyr::filter(cluster == unique(enrich.data$cluster)[x])
gene.ent <- clusterProfiler::bitr(tmp$gene,
                                    fromType = "SYMBOL",
                                    toType = c("ENTREZID"),
                                    OrgDb = org.Hs.eg.db)
tartget.gene = unlist(gene.ent[, "ENTREZID"])
ego <- clusterProfiler::enrichGO(gene       = tartget.gene,
                              keyType       = toType,
                              OrgDb         = OrgDb,
                              ont           = "ALL",
                              pAdjustMethod = "BH",
                              pvalueCutoff  = 1,
                              qvalueCutoff  = 0.2,
                              readable      = readable)
ekegg <- clusterProfiler::enrichKEGG(gene   = tartget.gene,
                              keyType       = "kegg",
                              organism      = "hsa",
                              universe      = NULL,
                              pvalueCutoff  = 0.2,
                              pAdjustMethod = "BH",
                              qvalueCutoff  = 1)

df_go <- data.frame(clusterProfiler::setReadable(ego,OrgDb = OrgDb,keyType = toType)) %>%
dplyr::filter(pvalue < pvalueCutoff) %>%
dplyr::mutate(group = paste("C",unique(enrich.data$cluster)[x],sep = '')) %>%
dplyr::arrange(pvalue)
df_kegg <- data.frame(clusterProfiler::setReadable(ekegg,OrgDb = OrgDb,keyType = toType)) %>%
dplyr::filter(pvalue < pvalueCutoff) %>%
dplyr::mutate(group = paste("C",unique(enrich.data$cluster)[x],sep = '')) %>%
dplyr::arrange(pvalue) 
      
df_4 = rbind(df_go[,colnames(df_kegg)], df_kegg)%>%dplyr::arrange(pvalue)
df_4_top = df_4[1:4,c('group','Description','pvalue','GeneRatio')] 

x = 5
tmp <- enrich.data %>% dplyr::filter(cluster == unique(enrich.data$cluster)[x])
gene.ent <- clusterProfiler::bitr(tmp$gene,
                                    fromType = "SYMBOL",
                                    toType = c("ENTREZID"),
                                    OrgDb = org.Hs.eg.db)
tartget.gene = unlist(gene.ent[, "ENTREZID"])
ego <- clusterProfiler::enrichGO(gene       = tartget.gene,
                              keyType       = toType,
                              OrgDb         = OrgDb,
                              ont           = "ALL",
                              pAdjustMethod = "BH",
                              pvalueCutoff  = 1,
                              qvalueCutoff  = 1,
                              readable      = readable)
ekegg <- clusterProfiler::enrichKEGG(gene   = tartget.gene,
                              keyType       = "kegg",
                              organism      = "hsa",
                              universe      = NULL,
                              pvalueCutoff  = 1,
                              pAdjustMethod = "BH",
                              qvalueCutoff  = 1)

df_go <- data.frame(clusterProfiler::setReadable(ego,OrgDb = OrgDb,keyType = toType)) %>%
dplyr::filter(pvalue < pvalueCutoff) %>%
dplyr::mutate(group = paste("C",unique(enrich.data$cluster)[x],sep = '')) %>%
dplyr::arrange(pvalue)
df_kegg <- data.frame(clusterProfiler::setReadable(ekegg,OrgDb = OrgDb,keyType = toType)) %>%
dplyr::filter(pvalue < pvalueCutoff) %>%
dplyr::mutate(group = paste("C",unique(enrich.data$cluster)[x],sep = '')) %>%
dplyr::arrange(pvalue) 
      
df_5 = rbind(df_go[,colnames(df_kegg)], df_kegg)%>%dplyr::arrange(pvalue)
df_5_top = df_5[c(1,2,4,6),c('group','Description','pvalue','GeneRatio')] 

x = 6
tmp <- enrich.data %>% dplyr::filter(cluster == unique(enrich.data$cluster)[x])
gene.ent <- clusterProfiler::bitr(tmp$gene,
                                    fromType = "SYMBOL",
                                    toType = c("ENTREZID"),
                                    OrgDb = org.Hs.eg.db)
tartget.gene = unlist(gene.ent[, "ENTREZID"])
ego <- clusterProfiler::enrichGO(gene       = tartget.gene,
                              keyType       = toType,
                              OrgDb         = OrgDb,
                              ont           = "ALL",
                              pAdjustMethod = "BH",
                              pvalueCutoff  = 1,
                              qvalueCutoff  = 0.2,
                              readable      = readable)
ekegg <- clusterProfiler::enrichKEGG(gene   = tartget.gene,
                              keyType       = "kegg",
                              organism      = "hsa",
                              universe      = NULL,
                              pvalueCutoff  = 0.2,
                              pAdjustMethod = "BH",
                              qvalueCutoff  = 1)

df_go <- data.frame(clusterProfiler::setReadable(ego,OrgDb = OrgDb,keyType = toType)) %>%
dplyr::filter(pvalue < pvalueCutoff) %>%
dplyr::mutate(group = paste("C",unique(enrich.data$cluster)[x],sep = '')) %>%
dplyr::arrange(pvalue)
df_kegg <- data.frame(clusterProfiler::setReadable(ekegg,OrgDb = OrgDb,keyType = toType)) %>%
dplyr::filter(pvalue < pvalueCutoff) %>%
dplyr::mutate(group = paste("C",unique(enrich.data$cluster)[x],sep = '')) %>%
dplyr::arrange(pvalue)  
df_6 = rbind(df_go[,colnames(df_kegg)], df_kegg)%>%dplyr::arrange(pvalue)
df_6_top = df_6[1:4,c('group','Description','pvalue','GeneRatio')] 

x = 7
tmp <- enrich.data %>% dplyr::filter(cluster == unique(enrich.data$cluster)[x])
gene.ent <- clusterProfiler::bitr(tmp$gene,
                                    fromType = "SYMBOL",
                                    toType = c("ENTREZID"),
                                    OrgDb = org.Hs.eg.db)
tartget.gene = unlist(gene.ent[, "ENTREZID"])
ego <- clusterProfiler::enrichGO(gene       = tartget.gene,
                              keyType       = toType,
                              OrgDb         = OrgDb,
                              ont           = "ALL",
                              pAdjustMethod = "BH",
                              pvalueCutoff  = 1,
                              qvalueCutoff  = 0.2,
                              readable      = readable)
ekegg <- clusterProfiler::enrichKEGG(gene   = tartget.gene,
                              keyType       = "kegg",
                              organism      = "hsa",
                              universe      = NULL,
                              pvalueCutoff  = 0.2,
                              pAdjustMethod = "BH",
                              qvalueCutoff  = 1)
df_go <- data.frame(clusterProfiler::setReadable(ego,OrgDb = OrgDb,keyType = toType)) %>%
dplyr::filter(pvalue < pvalueCutoff) %>%
dplyr::mutate(group = paste("C",unique(enrich.data$cluster)[x],sep = '')) %>%
dplyr::arrange(pvalue)
df_kegg <- data.frame(clusterProfiler::setReadable(ekegg,OrgDb = OrgDb,keyType = toType)) %>%
dplyr::filter(pvalue < pvalueCutoff) %>%
dplyr::mutate(group = paste("C",unique(enrich.data$cluster)[x],sep = '')) %>%
dplyr::arrange(pvalue) 
df_7 = rbind(df_go[,colnames(df_kegg)], df_kegg)%>%dplyr::arrange(pvalue)
df_7_top = df_7[1:4,c('group','Description','pvalue','GeneRatio')] 

x = 8
tmp <- enrich.data %>% dplyr::filter(cluster == unique(enrich.data$cluster)[x])
gene.ent <- clusterProfiler::bitr(tmp$gene,
                                    fromType = "SYMBOL",
                                    toType = c("ENTREZID"),
                                    OrgDb = org.Hs.eg.db)
tartget.gene = unlist(gene.ent[, "ENTREZID"])
ego <- clusterProfiler::enrichGO(gene       = tartget.gene,
                              keyType       = toType,
                              OrgDb         = OrgDb,
                              ont           = "ALL",
                              pAdjustMethod = "BH",
                              pvalueCutoff  = 1,
                              qvalueCutoff  = 0.2,
                              readable      = readable)
ekegg <- clusterProfiler::enrichKEGG(gene   = tartget.gene,
                              keyType       = "kegg",
                              organism      = "hsa",
                              universe      = NULL,
                              pvalueCutoff  = 0.2,
                              pAdjustMethod = "BH",
                              qvalueCutoff  = 1)
df_go <- data.frame(clusterProfiler::setReadable(ego,OrgDb = OrgDb,keyType = toType)) %>%
dplyr::filter(pvalue < pvalueCutoff) %>%
dplyr::mutate(group = paste("C",unique(enrich.data$cluster)[x],sep = '')) %>%
dplyr::arrange(pvalue)
df_kegg <- data.frame(clusterProfiler::setReadable(ekegg,OrgDb = OrgDb,keyType = toType)) %>%
dplyr::filter(pvalue < pvalueCutoff) %>%
dplyr::mutate(group = paste("C",unique(enrich.data$cluster)[x],sep = '')) %>%
dplyr::arrange(pvalue)   
df_8 = rbind(df_go[,colnames(df_kegg)], df_kegg)%>%dplyr::arrange(pvalue)
df_8_top = df_8[1:4,c('group','Description','pvalue','GeneRatio')] 


enrich.all = rbind(df_1_top,df_2_top,df_3_top,df_4_top,df_5_top,df_6_top,df_7_top,df_8_top)
write.csv(enrich.all, file = "enrich_cluster8_20231130.csv")

annt_top = HeatmapAnnotation(age = anno_text(dat_num$type),
                             log_age_weeks = anno_barplot(dat_num$log_age_weeks)
                              )
mycol = c('#4361EE', '#7209B7', '#3A0CA3','#f725c6' , "#4CC9F0")
mycol = c('#dabfff', '#907ad6', '#4f518c','#2c2a4a' , "#7fdeff")
mycol = c('#360568', '#5b2a86', '#7785ac','#9ac6c5' , "#a5e6ba")
mycol = c('#97dffc', '#858ae3', '#613dc1','#4e148c' , "#2c0735")

# mycol = c(  "#e76f51" ,'#f4a261', '#e9c46a','#2a9d8f','#264653')
mycol = c('#5188ba', '#eb9172', '#8ecae6','#832404', '#e86a59','#02274a', '#2c0735', '#f725c6') # '#8ecae6',
# "#1F77B4FF" "#FF7F0EFF" "#2CA02CFF" "#D62728FF" "#9467BDFF"
# mark <- rownames(df)[sample(1:nrow(df),28, replace = F)]
# png(paste0('cluster_heatmap_nosplit_',j,'.png'), width = 3000, height = 3000, units = 'px', res = 300)

source("/n02dat01/users/ypwang/Gradient/Need_Packages/ClusterGVis-main/R/visCluster_ping.R")
enrich.all$pvalue = enrich.all$pvalue*0.6
pdf(paste0('Heatmap_cluster8.pdf'), height = 9, width = 12, onefile = F)
visCluster_ping(object = ck, plot.type = "both",  # line.side = "left", 
            ht.col = "RdBu", mline.col = '#808080',
            row_dend_gp = gpar(col = "#636363"),
            HeatmapAnnotation = annt_top, # column_names_rot = 60,  markGenes = mark, markGenes.side = "left", genes.gp = c('italic',fontsize = 10), 
            column_names_gp = gpar(fontsize = 0), 
            show_heatmap_legend = F,
            sample.group = rep(c("C1","C2","C3","C4","C5"),c(10,3,4,5,6)), # c("fetal","infant","child","adolescent","adult")
            column_title = c("Fetal","Infant","Child","Adolescent","Adult"),
            show_row_dend = T, annoTerm.data = enrich.all, go.size = "pval" , # add.bar = T,
            ctAnno.col = mycol, # line.col = 'grey',# ggsci::pal_d3()(5),line.col = mycol,
            go.col = c(rep(mycol[1],4),rep(mycol[2],4), 
                      rep(mycol[3],4),rep(mycol[4],4),
                      rep(mycol[5],4),rep(mycol[6],4),
                      rep(mycol[7],4),rep(mycol[8],4))
            ) # go.size = "pval"
dev.off()   

# ------------

## dN/dS
# ------------
library(dplyr)
library(biomaRt)
library(reshape2)
library(ggplot2)
library(cowplot)
library(ggpubr)
library(ggrepel)
library(dplyr)
library(tidyr)
library(tibble)
library(ggplot2)
library(ggplotify)
library(cowplot)
library(readxl)
library(RColorBrewer)
library(ggpubr)
library(Hmisc)
library(ggrepel)
library(reshape2)
library(gglayer)
require(gglayer)
library(PupillometryR)
## set up directories
out_dir = '/n02dat01/users/ypwang/Gradient/GeneticGradient/Data/GCI2dNdS/'
RdBu = read.csv('/n02dat01/users/ypwang/Gradient/GeneticGradient/Fig/Color/RdBu_Hex_200.csv', row.names = 1, header=TRUE, check.names=FALSE)
setwd(out_dir)
data_dir <- "/n02dat01/users/ypwang/Gradient/GeneticGradient/Data/"
output_dir = "/n02dat01/users/ypwang/Gradient/GeneticGradient/Data/GCI2DEG/"
setwd(output_dir)
gene_FG1=read.csv(paste0(data_dir,'PLS_Betas_cere2extra_15624_permutation.csv'), header=TRUE, check.names=FALSE)
colnames(gene_FG1) = c('hgnc_symbol', "GCI","p_perm","p_perm_fdr0.05","type","color")  
hsapiens_GFList <- gene_FG1$hgnc_symbol

# Biomart
# For the human data, get the homology gene dN/ds
# need to change the clash into global mode
ensemblhsapiens = useMart(host = 'https://jan2020.archive.ensembl.org', 
                        biomart = 'ENSEMBL_MART_ENSEMBL', 
                        dataset = "hsapiens_gene_ensembl")

GF_info = getBM(attributes = c('hgnc_symbol','ensembl_gene_id' ), 
                    filters = 'hgnc_symbol', 
                    values = hsapiens_GFList, 
                    mart = ensemblhsapiens)
GFList = merge(x=GF_info, y=gene_FG1, by='hgnc_symbol', all.y = T) # dplyr::left_join(gene_FG1,GF_info)

filters = listFilters(ensemblhsapiens)
attributes = listAttributes(ensemblhsapiens)
a = listAttributes(ensemblhsapiens, page="feature_page")
species = listDatasets(useMart("ensembl"))
for (species in c('mmusculus','ptroglodytes')){ # celegans, 'fish',
    print(species)
    hsp <- getBM(attributes = c('ensembl_gene_id', paste0(species,'_homolog_ensembl_gene'),  
                              paste0(species,'_homolog_dn'), paste0(species,'_homolog_ds'), 
                              paste0(species,'_homolog_orthology_type')), 
                filters = 'hgnc_symbol', 
                values = hsapiens_GFList, 
                mart = ensemblhsapiens)

    hsp_1 <- subset(hsp, hsp[paste0(species,'_homolog_orthology_type')] == 'ortholog_one2one')
    hsp_1$dnds_ratios <- hsp_1[,paste0(species,'_homolog_dn')]/hsp_1[,paste0(species,'_homolog_ds')]
    # colnames(hsp_1)[6] = paste0(species,'_dnds_ratios')
    hsp_1 = hsp_1[is.finite(hsp_1$dnds_ratios),] 

    GF_dnds = merge(x=GFList, y=hsp_1, by='ensembl_gene_id', all.x = T)
    GF_dnds = GF_dnds[is.finite(GF_dnds$dnds_ratios),] 
    print(cor.test(GF_dnds[,'GCI'],GF_dnds[,'dnds_ratios'],na.action = "na.exclude" ))
    print(cor.test(GF_dnds[GF_dnds$p_perm_fdr0.05=='True','GCI'],
                GF_dnds[GF_dnds$p_perm_fdr0.05=='True','dnds_ratios'],na.action = "na.exclude" ))
    write.csv(hsp_1, file = paste0(species,"_dnds.csv"))
    assign(paste0('hsp_',species), hsp_1)
    assign(paste0('dnds_',species),GF_dnds)
}

## Is dN/dS higher in rats compared to GCI distribution in monkeys and rats?
colnames(hsp_mmusculus)[6] = 'mmusculus_dnds_ratios'
colnames(hsp_ptroglodytes)[6] = 'ptroglodytes_dnds_ratios'
GF_dnds = merge(x=GFList, y=hsp_mmusculus, by='ensembl_gene_id', all.x = T)
GF_dnds = merge(x=GF_dnds, y=hsp_ptroglodytes, by='ensembl_gene_id',all.x = T)
GF_dnds_FG1 = GF_dnds[GF_dnds$p_perm_fdr0.05=='True',]
plotData <- rbind(data.frame(Gene=GF_dnds_FG1$hgnc_symbol,
                             condition=rep('Mouse', length(GF_dnds_FG1$hgnc_symbol)), 
                             value = as.numeric(GF_dnds_FG1$mmusculus_dnds_ratios)),
                  data.frame(Gene=GF_dnds_FG1$hgnc_symbol,
                             condition=rep('Chimpanzee', length(GF_dnds_FG1$hgnc_symbol)), 
                             value = as.numeric(GF_dnds_FG1$ptroglodytes_dnds_ratios)))
names(plotData) <- c("Gene","condition","value")

library(ggnewscale)
top10 <- filter(plotData, value > 1) %>%
           distinct(Gene, .keep_all = T) %>%
           top_n(10, abs(value))
anno_har = plotData[plotData$Gene%in%c('RWDD3', 'HNMT','TBL1XR1'),]
anno = rbind(top10 %>% top_n(5, abs(value)), anno_har[anno_har$condition=='Chimpanzee',])

dNdS = ggplot(plotData, aes(x = condition, y = value)) + 
    geom_point(data=(plotData %>% dplyr::filter(value > 1)), position = position_jitter(0.03), 
                      size = 4, alpha = 0.9, aes(x = condition, y = value, color = value)) +
    scale_color_gradientn(colours = RdBu[120:199,1]) +
    geom_point(data=(plotData %>% dplyr::filter(value <= 1)), position = position_jitter(0.15), 
                      size = 4, alpha = 0.45, aes(x = condition, y = value), color="gray") +
    geom_hline(yintercept = 1, linetype = "dashed", size = 0.2, alpha = 1)  + 
    scale_y_continuous("dN/dS ratio",expand = c(0,0), limits = c(0,3.4)) +
    annotate("rect", xmin = 0.4, xmax = Inf, ymin = 0, ymax = 1, alpha = .2, fill="gray80")  + 
    annotate("text", x = .47, y = 0.5, label="\U2190 More conserved", color="gray20", size=3,fontface = 'italic')+ 
    annotate("text", x = .47, y = 1.5, label="Less conserved \U2192", color="gray20", size=3,fontface = 'italic')+
    new_scale_color() +
    geom_flat_violin(aes(x = condition, y = value, fill = condition, color = condition), 
                     position = position_nudge(x = 0.25, y = 0), adjust = 2, 
                     alpha = 0.9, trim = TRUE, scale = "width") +
    scale_color_manual(values = rev(c("#5188ba", "#eb9172"))) +
    scale_fill_manual(values = rev(c("#5188ba", "#eb9172"))) +
    stat_compare_means(aes(label = ..p.signif..), method = "wilcox.test",  
                       comparisons = list(c("Mouse","Chimpanzee"))) + 
    geom_text_repel(data=anno, aes(label=Gene), 
                    nudge_x=-.2, color="gray50", fill="white",size=2,seed = 1) + 
    xlab("") +
    coord_flip()+
    theme_cowplot() + 
    scale_shape_identity() +
    theme(legend.position = "none", plot.title = element_text(size = 20), 
          axis.title = element_text(size = 21), axis.text = element_text(size = 15), 
          axis.text.x = element_text(angle = 0, hjust = 0, vjust = 0)) 

png(file = 'dNdS.png', width = 3000, height = 1800, units = 'px', res = 300)
    # CairoPDF(figure_out, height=2, width=3)
print(dNdS)
dev.off()
write.csv(file = 'GCI_dNdS.csv', as.data.frame(plotData))
# ----------------

## GCIsig2DEGpatient
# ----------------
library(tidyverse)
library(ggcorrplot)
library(ggpubr)
library(cowplot)
library(showtext)
library(Cairo)
library(psych)
# library(customLayout)
library(RColorBrewer)

## set up directories
base_dir = '/n02dat01/users/ypwang/Gradient/GeneticGradient/Data/GCI2DEG'
RdBu = read.csv('/n02dat01/users/ypwang/Gradient/GeneticGradient/Fig/Color/RdBu_Hex_200.csv', row.names = 1, header=TRUE, check.names=FALSE)

# read gene-wise spatial correlations to MDD effect maps
data_dir <- "/n02dat01/users/ypwang/Gradient/GeneticGradient/Data/"
output_dir = "/n02dat01/users/ypwang/Gradient/GeneticGradient/Data/GCI2DEG/"
setwd(output_dir)
gene_FG1=read.csv(paste0(data_dir,'PLS_Betas_cere2extra_15624_permutation.csv'), header=TRUE, check.names=FALSE)
colnames(gene_FG1) = c('gene', "GCI","p_perm","p_perm_fdr0.05","type","color")  
gene_FG1_sig = gene_FG1[gene_FG1[['p_perm_fdr0.05']]=='True',] # 0.05 1024
colnames(gene_FG1_sig)
gene_FG1     = gene_FG1[order(gene_FG1$GCI),]
gene_FG1_sig = gene_FG1_sig[order(gene_FG1_sig$GCI),]

# load diff expr 
table_dir  = '/n02dat01/users/ypwang/Gradient/Need_Packages/Shared-molecular-neuropathology-across-major-psychiatric-disorders-parallels-polygenic-overlap-master/results/tables'
for (dis in c('AAD', 'ASD', 'BD', 'MDD', 'SCZ')){
    data = read_csv(paste0(table_dir, '/Microarray_',dis,'_metaanalysis_092017.csv'))
    colnames(data)[1]  =  "ensembl_gene_id"
    assign(paste0('gandal_',dis),data)
}
gandal_meta = read_csv(paste0(table_dir, '/Manuscript/TableS1 - Microarray Meta-Analyses.csv'))
gandal_AAD = merge(x=gandal_AAD, y=gandal_meta[c("ensembl_gene_id",'hgnc_symbol')], 
              by="ensembl_gene_id", all.x= TRUE)
colnames(gandal_AAD) = colnames(gandal_ASD)

output_dir = "/n02dat01/users/ypwang/Gradient/GeneticGradient/Data/GCI2DEG/"
setwd(output_dir)

for (gci in c('FG1', 'FG1_sig')){
    
    print(gci)
    for (dis in c('AAD', 'ASD', 'BD', 'MDD', 'SCZ')){
        print(dis)
        data = merge(x=get(paste0('gene_', gci)), y=get(paste0('gandal_',dis)), by.x='gene', by.y='symbol')
        assign(paste0(gci,'_',dis),data)
        assign(paste0(dis,'_cor'),cor.test(data$beta, data$GCI))
    }
    # organize for plotting
    genewise_cor = rbind(data.frame(cor=MDD_cor$estimate, pheno='MDD', l95=MDD_cor$conf.int[1], u95=MDD_cor$conf.int[2], pval=MDD_cor$p.value),
                         data.frame(cor=BD_cor$estimate, pheno='BD', l95=BD_cor$conf.int[1], u95=BD_cor$conf.int[2], pval=BD_cor$p.value),
                         data.frame(cor=ASD_cor$estimate, pheno='ASD', l95=ASD_cor$conf.int[1], u95=ASD_cor$conf.int[2], pval=ASD_cor$p.value),
                         data.frame(cor=AAD_cor$estimate, pheno='AAD', l95=AAD_cor$conf.int[1], u95=AAD_cor$conf.int[2], pval=AAD_cor$p.value),
                         data.frame(cor=SCZ_cor$estimate, pheno='SCZ', l95=SCZ_cor$conf.int[1], u95=SCZ_cor$conf.int[2], pval=SCZ_cor$p.value))

    # barplot
    p = ggplot(genewise_cor, aes(x=pheno, y=cor, fill=pheno)) +
                geom_bar(stat="identity", color="black", position=position_dodge(), width=0.85) +
                geom_errorbar(aes(ymin=l95, ymax=u95), width=.4, position=position_dodge(.9)) +
                theme_classic() +
                # scale_fill_manual(values=c(RdBu[200,1],'grey69','grey69','grey69','grey69')) +
                ggtitle(paste0(gci,'_nobins')) +
                scale_y_continuous(limit=c(-0.1,0.15), breaks=seq(-0.1,0.15, 0.05), expand=c(0,0))
    figure_out = paste0(gci,'_meta_nobins.png')
    png(file = figure_out, width = 2000, height = 1800, units = 'px', res = 300)
    # CairoPDF(figure_out, height=2, width=3)
    print(p)
    dev.off()
    print(p)
}
# Figure 4
gci = 'FG1'
sp = 40
color = c(RdBu[abs(MDD_decile_cor$rho[1,2])*100,1],#  "#3c8abe"
          RdBu[BD_decile_cor$rho[1,2]*175+100,1], # "#eb9172","#f5a886","#ce4f45","#700320"
          RdBu[ASD_decile_cor$rho[1,2]*175+100,1], # "#f8bb9e"
          RdBu[SCZ_decile_cor$rho[1,2]*175+100,1], #  "#f19e7d"
          'grey69')
names(color) = c('MDD', 'BD', 'ASD', 'SCZ', 'AAD')

print(sp)
for (dis in c('AAD', 'ASD', 'BD', 'MDD', 'SCZ')){
      print(dis)
      data = get(paste0(gci,'_',dis))
      data = data[!is.na(data$beta),]
      data = data[order(data$beta),]

      data$group = as.numeric(cut_number(1:nrow(data), sp))
      by_decile         = data %>% group_by(group) %>% summarise(GCI=mean(GCI), log2FC=mean(beta))
      decile_cor        = cor.test(by_decile$GCI, by_decile$log2FC, method='spearman')
      decile_cor = cor.ci(cbind(by_decile$log2FC, by_decile$GCI), method='spearman', plot=F)
      assign(paste0(dis,'_decile_cor'), decile_cor)

      if (dis %in% c('BD', 'MDD', 'SCZ')){y_lim = 0.25}else {y_lim = 0.5}
      plot_point = ggplot(data=by_decile, aes(x=GCI, y=log2FC)) +
                        geom_point(shape=21, show.legend = FALSE, linewidth=22, color = color[dis], fill = color[dis]) +
                        geom_smooth(size=.5, fill='gray69', color='black', show.legend = FALSE, 
                                    linetype = 'dashed', fullrange = TRUE, method = lm, formula = y~x, se = FALSE) +
                        ggtitle(paste0(dis)) +
                        labs(x = 'GCI', y = base::expression(paste(log[2], 'FC')))+
                        stat_cor(method = 'spearman', label.x = -0.00009) +
                        # scale_x_continuous(limits=c(-0.0006,0.0006), expand=c(0,0), breaks=seq(-0.0006,0.0006,0.0002)) +
                        scale_x_continuous(limits=c(-0.0001,0.0001), expand=c(0,0), breaks=seq(-0.0001,0.0001,0.00005)) +
                        scale_y_continuous(limits=c(-y_lim,y_lim), expand=c(0,0), breaks=seq(-y_lim,y_lim,0.25)) +
                        scale_color_manual(values=RdBu[200,1]) +
                        theme_classic() +
                        theme(axis.text = element_text(size = 8,color='black'),
                              axis.title = element_text(size = 12),
                              legend.position = "none", 
                              plot.title = element_text(size = 15, hjust = 0.5, color = 'black'),
                              axis.ticks=element_line(color='black'))
      assign(paste0(dis,'_point'), plot_point)
}
decile_genewise_cor = rbind(data.frame(cor=MDD_decile_cor$rho[1,2], pheno='MDD', l95=MDD_decile_cor$ci$lower, u95=MDD_decile_cor$ci$upper, pval= MDD_decile_cor$p),
                        data.frame(cor=BD_decile_cor$rho[1,2], pheno='BD', l95=BD_decile_cor$ci$lower, u95=BD_decile_cor$ci$upper, pval=BD_decile_cor$p),
                        data.frame(cor=ASD_decile_cor$rho[1,2], pheno='ASD', l95=ASD_decile_cor$ci$lower, u95=ASD_decile_cor$ci$upper, pval=ASD_decile_cor$p),
                        data.frame(cor=SCZ_decile_cor$rho[1,2], pheno='SCZ', l95=SCZ_decile_cor$ci$lower, u95=SCZ_decile_cor$ci$upper, pval=SCZ_decile_cor$p),
                        data.frame(cor=AAD_decile_cor$rho[1,2], pheno='AAD', l95=AAD_decile_cor$ci$lower, u95=AAD_decile_cor$ci$upper, pval=AAD_decile_cor$p))
decile_genewise_cor$pheno = factor(decile_genewise_cor$pheno, levels=rev(decile_genewise_cor$pheno))
dis_meta = ggplot(decile_genewise_cor, aes(x=pheno, y=cor, fill=pheno)) +
            geom_bar(stat="identity", color="black", position=position_dodge(), width=0.85) +
            geom_errorbar(aes(ymin=l95, ymax=u95), width=.3, position=position_dodge(.9)) +
            # ggtitle(paste0(gci,'_',sp)) +
            labs(x = 'Disease', y = 'R')+
            coord_flip() +
            scale_fill_manual(values=rev(c(color))) +
            theme_classic() +
            scale_y_continuous(limit=c(-1,1), breaks=seq(-1,1,.5), expand=c(0,0)) +
            theme(axis.text = element_text(size = 8,color='black'),
                  axis.title = element_text(size = 12),
                  legend.position = "none", 
                  plot.title = element_text(size = 15, hjust = 0.5),
                  axis.ticks=element_line(color='black'))

dNdS = ggplot(plotData, aes(x = condition, y = value)) + 
    geom_point(data=(plotData %>% dplyr::filter(value > 1)), position = position_jitter(0.03), 
                      size = 4, alpha = 0.9, aes(x = condition, y = value, color = value)) +
    scale_color_gradientn(colours = RdBu[120:199,1]) +
    geom_point(data=(plotData %>% dplyr::filter(value <= 1)), position = position_jitter(0.15), 
                      size = 4, alpha = 0.45, aes(x = condition, y = value), color="gray") +
    geom_hline(yintercept = 1, linetype = "dashed", size = 0.2, alpha = 1)  + 
    scale_y_continuous("dN/dS ratio",expand = c(0,0), limits = c(0,3.4)) +
    annotate("rect", xmin = 0.4, xmax = Inf, ymin = 0, ymax = 1, alpha = .2, fill="gray80")  + 
    annotate("text", x = .47, y = 0.5, label="\U2190 More conserved", color="gray20", size=4,fontface = 'italic')+ 
    annotate("text", x = .47, y = 1.5, label="Less conserved \U2192", color="gray20", size=4,fontface = 'italic')+
    new_scale_color() +
    geom_flat_violin(aes(x = condition, y = value, fill = condition, color = condition), 
                     position = position_nudge(x = 0.25, y = 0), adjust = 2, 
                     alpha = 0.9, trim = TRUE, scale = "width") +
    scale_color_manual(values = rev(c("#5188ba", "#eb9172"))) +
    scale_fill_manual(values = rev(c("#5188ba", "#eb9172"))) +
    stat_compare_means(aes(label = ..p.signif..), method = "wilcox.test",  
                       comparisons = list(c("Mouse","Chimpanzee"))) + 
    geom_text_repel(data=anno, aes(label=Gene), 
                    nudge_x=-.2, color="gray50", fill="white",size=3,seed = 1) + 
    xlab("") +
    coord_flip()+
    theme_cowplot() + 
    scale_shape_identity() +
    theme(legend.position = "none", plot.title = element_text(size = 20),  
          axis.text.y = element_text(size = 15, angle = 90, hjust = 0.5),
          axis.text.x = element_text(size = 8),
          axis.title.x = element_text(size = 12)
      ) 

gg <- ggdraw() +     
draw_plot(dis_meta,  0,    1/2, 2/10, 1/2) + # 在母图上半部，占母图比例1/2  
draw_plot(MDD_point, 2/10,  1/2, 2/10, 1/2) + # 在母图左下角，占母图比例1/4  
draw_plot(BD_point,  4/10,  1/2, 2/10, 1/2) + # 在母图右下角，占母图比例1/4 
draw_plot(ASD_point, 0,    0,   2/10, 1/2) + 
draw_plot(SCZ_point, 2/10,  0,   2/10, 1/2) + 
draw_plot(AAD_point, 4/10,  0,   2/10, 1/2) +
draw_plot(dNdS, 6.1/10,0, 3.9/10,1)
print(gg) 
figure_out = paste0(gci,'_dNdS+DEG.png')
png(file = figure_out, width = 3000, height = 1000, units = 'px', res = 200)
print(gg)
dev.off()

# Supple Figure
gci = 'FG1'
color = c(RdBu[abs(MDD_decile_cor$rho[1,2])*100,1],#  "#3c8abe"
          RdBu[BD_decile_cor$rho[1,2]*175+100,1], # "#eb9172","#f5a886","#ce4f45","#700320"
          RdBu[ASD_decile_cor$rho[1,2]*175+100,1], # "#f8bb9e"
          RdBu[SCZ_decile_cor$rho[1,2]*175+100,1], #  "#f19e7d"
          'grey69')
names(color) = c('MDD', 'BD', 'ASD', 'SCZ', 'AAD')

for (sp in c(20,40,60,80,100)){
      print(sp)
      for (dis in c('AAD', 'ASD', 'BD', 'MDD', 'SCZ')){
            print(dis)
            data = get(paste0(gci,'_',dis))
            data = data[!is.na(data$beta),]
            data = data[order(data$beta),]

            data$group = as.numeric(cut_number(1:nrow(data), sp))
            by_decile         = data %>% group_by(group) %>% summarise(GCI=mean(GCI), log2FC=mean(beta))
            decile_cor        = cor.test(by_decile$GCI, by_decile$log2FC, method='spearman')
            decile_cor = cor.ci(cbind(by_decile$log2FC, by_decile$GCI), method='spearman', plot=F)
            assign(paste0(dis,'_decile_cor'), decile_cor)

            if (dis %in% c('BD', 'MDD', 'SCZ')){y_lim = 0.25}else {y_lim = 0.5}
            plot_point = ggplot(data=by_decile, aes(x=GCI, y=log2FC)) +
                              geom_point(shape=21, show.legend = FALSE, linewidth=22, color = color[dis], fill = color[dis]) +
                              geom_smooth(size=.5, fill='gray69', color='black', show.legend = FALSE, 
                                          linetype = 'dashed', fullrange = TRUE, method = lm, formula = y~x, se = FALSE) +
                              ggtitle(paste0(dis)) +
                              labs(x = 'GCI', y = base::expression(paste(log[2], 'FC')))+
                              stat_cor(method = 'spearman', label.x = -0.00009) +
                              # scale_x_continuous(limits=c(-0.0006,0.0006), expand=c(0,0), breaks=seq(-0.0006,0.0006,0.0002)) +
                              scale_x_continuous(limits=c(-0.0001,0.0001), expand=c(0,0), breaks=seq(-0.0001,0.0001,0.00005)) +
                              scale_y_continuous(limits=c(-y_lim,y_lim), expand=c(0,0), breaks=seq(-y_lim,y_lim,0.25)) +
                              scale_color_manual(values=RdBu[200,1]) +
                              theme_classic() +
                              theme(axis.text = element_text(size = 8,color='black'),
                                    axis.title = element_text(size = 12),
                                    legend.position = "none", 
                                    plot.title = element_text(size = 15, hjust = 0.5, color = 'black'),
                                    axis.ticks=element_line(color='black'))
            assign(paste0(dis,'_point'), plot_point)
      }
      decile_genewise_cor = rbind(data.frame(cor=MDD_decile_cor$rho[1,2], pheno='MDD', l95=MDD_decile_cor$ci$lower, u95=MDD_decile_cor$ci$upper, pval= MDD_decile_cor$p),
                              data.frame(cor=BD_decile_cor$rho[1,2], pheno='BD', l95=BD_decile_cor$ci$lower, u95=BD_decile_cor$ci$upper, pval=BD_decile_cor$p),
                              data.frame(cor=ASD_decile_cor$rho[1,2], pheno='ASD', l95=ASD_decile_cor$ci$lower, u95=ASD_decile_cor$ci$upper, pval=ASD_decile_cor$p),
                              data.frame(cor=SCZ_decile_cor$rho[1,2], pheno='SCZ', l95=SCZ_decile_cor$ci$lower, u95=SCZ_decile_cor$ci$upper, pval=SCZ_decile_cor$p),
                              data.frame(cor=AAD_decile_cor$rho[1,2], pheno='AAD', l95=AAD_decile_cor$ci$lower, u95=AAD_decile_cor$ci$upper, pval=AAD_decile_cor$p))
      decile_genewise_cor$pheno = factor(decile_genewise_cor$pheno, levels=rev(decile_genewise_cor$pheno))
      dis_meta = ggplot(decile_genewise_cor, aes(x=pheno, y=cor, fill=pheno)) +
                  geom_bar(stat="identity", color="black", position=position_dodge(), width=0.85) +
                  geom_errorbar(aes(ymin=l95, ymax=u95), width=.3, position=position_dodge(.9)) +
                  # ggtitle(paste0(gci,'_',sp)) +
                  labs(x = 'Disease', y = 'R')+
                  coord_flip() +
                  scale_fill_manual(values=rev(c(color))) +
                  theme_classic() +
                  scale_y_continuous(limit=c(-1,1), breaks=seq(-1,1,.5), expand=c(0,0)) +
                  theme(axis.text = element_text(size = 8,color='black'),
                        axis.title = element_text(size = 12),
                        legend.position = "none", 
                        plot.title = element_text(size = 15, hjust = 0.5),
                        axis.ticks=element_line(color='black'))
      gg <- ggdraw() +     
      draw_plot(dis_meta,  0,    1/2, 1/3.1, 1/2) + # 在母图上半部，占母图比例1/2  
      draw_plot(MDD_point, 1/3,  1/2, 1/3.1, 1/2) + # 在母图左下角，占母图比例1/4  
      draw_plot(BD_point,  2/3,  1/2, 1/3.1, 1/2) + # 在母图右下角，占母图比例1/4 
      draw_plot(ASD_point, 0,    0,   1/3.1, 1/2) + 
      draw_plot(SCZ_point, 1/3,  0,   1/3.1, 1/2) + 
      draw_plot(AAD_point, 2/3,  0,   1/3.1, 1/2) 
      print(gg) 
      figure_out = paste0(gci,'_multi_sp',sp,'.png')
      png(file = figure_out, width = 1800, height = 1000, units = 'px', res = 200)
      print(gg)
      dev.off()
}
# ----------------