
library(seurat) #4.3.0.1
library(ggplot2) #3.4.3
library(dplyr) #1.1.3
library(colorRamp2) #0.1.0
library(PupillometryR) #geom_flat_violin #0.0.5


#*******************************#
#* Preliminaries: prepare colors 
#*******************************#
pathToData="./multiSpecies_limbRegeneration/data/scRNAseq_files/" #points to GitHub folder (https://github.com/BICC-UNIL-EPFL/multiSpecies_limbRegeneration/)
filename=paste0(pathToData,"colors_scRNAseq.xlsx")
sheets <- readxl::excel_sheets(filename)
colors_scRNAseq.l <- lapply(sheets, function(X) {
  d=data.frame(readxl::read_excel(filename, sheet = X))
  setNames(d[,2],d[,1])
})
names(colors_scRNAseq.l)=sheets

pathToData="./multiSpecies_limbRegeneration/data/multiomics_files/" #points to the GitHub folder (https://github.com/BICC-UNIL-EPFL/multiSpecies_limbRegeneration/)
filename=paste0(pathToData,"colors_multiomics.xlsx")
sheets <- readxl::excel_sheets(filename)
colors_multiomics.l <- lapply(sheets, function(X) {
  d=data.frame(readxl::read_excel(filename, sheet = X))
  setNames(d[,2],d[,1])
})
names(colors_multiomics.l)=sheets

col_featDist <- c("#a6cee3", "#1f78b4", "#b2df8a",
          "#33a02c", "#fb9a99", "#e31a1c",
          "#fdbf6f", "#ff7f00", "#cab2d6",
          "#6a3d9a", "#ffff99", "#b15928")


#*******************************#
#* Preliminaries: load main multiomics seurat object 
#*******************************#
pathToGEO="files_GSE292215/" # points to the data from GEO 
pathToProcessed="./multiSpecies_limbRegeneration/data/multiomics_files/processed" #points to the GitHub folder (https://github.com/BICC-UNIL-EPFL/multiSpecies_limbRegeneration/)

expression_matrix <- ReadMtx(
  mtx = paste0(pathToGEO,"/matrix_RNA.mtx"), features = paste0(pathToGEO,"/genes_RNA.tsv"),
  cells = paste0(pathToGEO,"/barcodes_RNA.tsv"), feature.column = 1
)

metadatafile = paste0(pathToProcessed,"/metadata.csv.gz")
metadata.df=read.csv(gzfile(metadatafile), row.names=1)
metadata.df$condition_TP=factor(metadata.df$condition_TP, levels=c("WT_0h","HPA_6h","HPA_24h","HPA_72h"))

seu_multiomics <- CreateSeuratObject(counts = expression_matrix, meta.data = metadata.df)
seu_multiomics <- NormalizeData(seu_multiomics)
seu_multiomics <- FindVariableFeatures(seu_multiomics)
seu_multiomics <- ScaleData(seu_multiomics)

filename=paste0(pathToProcessed,"/multiomics_UMAP_coords.xlsx")
sheets <- readxl::excel_sheets(filename)
umap_coords.l <- lapply(sheets, function(X) {
  d=data.frame(readxl::read_excel(filename, sheet = X))
  rownames(d)=d[,1]; d=d[,-1]; return(d)
})
names(umap_coords.l)=sheets
seu_multiomics@reductions=list()
seu_multiomics@reductions[["umap.rna"]]=CreateDimReducObject(embeddings = umap_coords.l$UMAP_rna %>% as.matrix(), key = "rnaUMAP_" , assay = "RNA")
seu_multiomics@reductions[["umap.atac"]]=CreateDimReducObject(embeddings = umap_coords.l$UMAP_atac %>% as.matrix(), key = "atacUMAP_" , assay = "ATAC")
seu_multiomics@reductions[["wnn.umap"]]=CreateDimReducObject(embeddings = umap_coords.l$UMAP_wnn %>% as.matrix(), key = "wnnUMAP_" , assay = "ATAC")

# re-create ATAC assay
frags.l=readRDS(paste0(pathToProcessed,"frags_l.Rds")) 
# note: frags.l[[i]]@path should be changed accordingly to point to the fragment files available on GEO (GSE292215)
# e.g., frags.l[["6h"]]@path="6HPA_ATAC_fragments.tsv.gz" (!! also need the index .tbi files - create it with tabix - e.g. tabix -p bed 6HPA_ATAC_fragments.tsv.gz )
peaks_reduce=readRDS(paste0(pathToProcessed,"peaks_reduce.Rds")) #(1.3Mb)
featMatrix_peaksReduced  <- FeatureMatrix(fragments = frags.l, features = peaks_reduce, cells = colnames(seu_multiomics))
gtf.gr=readRDS(paste0(pathToProcessed,"ATAC_annotations_gtf_gr.rds")) 
seu_multiomics[["ATAC"]] <- CreateChromatinAssay(counts = featMatrix_peaksReduced, fragments = frags.l, genome = "mm10", annotation = gtf.gr)
#seu_multiomics <- RunTFIDF(seu_multiomics)
#seu_multiomics <- FindTopFeatures(seu_multiomics, min.cutoff = 'q0')
#seu_multiomics <- RunSVD(seu_multiomics)
#seuratObj_mutimodal <- FindMultiModalNeighbors(seu_multiomics, reduction.list = list("pca", "lsi"), dims.list = list(1:50, 2:50))




#*******************************#
#* Preliminaries: load main scRNA-seq seurat objects 
#*******************************#

pathToGEO="files_GSE292212/" # points to the data from GEO 
pathToData="./multiSpecies_limbRegeneration/data/scRNAseq_files/" #points to the GitHub folder (https://github.com/BICC-UNIL-EPFL/multiSpecies_limbRegeneration/)

loadSeuratObj <- function(curDataset){
  expression_matrix <- ReadMtx(
    mtx = paste0(pathToGEO,curDataset,"/matrix_",curDataset,".mtx"), features = paste0(pathToGEO,curDataset,"/genes_",curDataset,".tsv"),
    cells = paste0(pathToGEO,curDataset,"/barcodes_",curDataset,".tsv"), feature.column = 1
  )
  
  metadatafile = paste0(pathToData,curDataset,"/metadata.csv")
  metadata.df=read.csv(metadatafile, row.names=1)
  if(curDataset=="VehCtrl_AERfactors"){
  metadata.df$orig.ident=factor(metadata.df$orig.ident, levels=c("Mouse_VehCtrl96w","Mouse_AERfactors96w","Mouse_AERfactorsALI"))
  }

  seu <- CreateSeuratObject(counts = expression_matrix, meta.data = metadata.df)
  seu <- NormalizeData(seu)
  seu <- FindVariableFeatures(seu)
  seu <- ScaleData(seu)

  umapfile = paste0(pathToData,curDataset,"/intUMAP_coords.csv")
  umap.df = read.csv(umapfile, row.names=1)  
  seu@reductions=list()
  seu@reductions[["umap.int"]]=CreateDimReducObject(embeddings = umap.df %>% as.matrix(), key = "intUMAP_", assay = "integrated")
  
  seu@assays$integrated=seu@assays$RNA # to avoid error when call subset()
  
  return(seu)
}

seuratObj_AER=loadSeuratObj("VehCtrl_AERfactors") 
seuratObj_Mouse=loadSeuratObj("Mouse_ALI_96w") 
seuratObj_Xen=loadSeuratObj("Xen_ALI_96w") 





#*******************************#
#* Preliminaries: prepare AUCell scores
#*******************************#
pathToData="./multiSpecies_limbRegeneration/data/scRNAseq_files/" #points to the GitHub folder (https://github.com/BICC-UNIL-EPFL/multiSpecies_limbRegeneration/)
filename=paste0(pathToData,"data_AUCellScores.xlsx")
sheets <- readxl::excel_sheets(filename)
data_AUCellScores.l <- lapply(sheets, function(X) {
  data.frame(readxl::read_excel(filename, sheet = X))
})
names(data_AUCellScores.l)=sheets


#*******************************#
#* Preliminaries: prepare gene lists
#*******************************#
pathToData="./multiSpecies_limbRegeneration/data/scRNAseq_files/" #points to the GitHub folder (https://github.com/BICC-UNIL-EPFL/multiSpecies_limbRegeneration/)
filename=paste0(pathToData,"geneLists_aer_blastema.xlsx")
sheets <- readxl::excel_sheets(filename)
aer_blastema.l <- lapply(sheets, function(X) {
  data.frame(readxl::read_excel(filename, sheet = X))
})
names(aer_blastema.l)=sheets


order_CellTypeExt_Xen=c(paste("CT_",1:15,sep=""),"Basal Ectoderm","Surface Ectoderm","Differentiating Ectoderm","Muscle","Blood","Immune","Endothelial","Pericytes","Schwann","Neuron")



#*******************************#
########## Figure 3 #########
#*******************************#

## Fig 3B: umap cell-types (VehCtrl_AERfactors)
seuratObj=seuratObj_AER
p_cellTypeExt <- DimPlot(seuratObj, reduction="umap.int", group.by = "CellTypeExt", cols=colors_scRNAseq.l$cellTypeExt_Mouse,order=TRUE) + ggtitle("Xen VehCtrl_AER", subtitle="CellTypeExt") +  guides(color = guide_legend(override.aes = list(size=4), ncol=1) )

## Fig 3D: dot plot AER genes (VehCtrl_AERfactors)
seuratObj=seuratObj_AER
genesForDotplots.df=aer_blastema.l$aer_blastema_Mouse
genesetName="aer_genes" 
Idents(seuratObj)="CellTypeExt" 
curCellType = "Basal Ectoderm"
seu = subset(seuratObj, idents =curCellType)
Idents(seu)="orig.ident"
d=DotPlot(seu, features=genesForDotplots.df[[,genesetName]] )
d$data$id=factor(d$data$id, levels=rev(c("Mouse_VehCtrl96w","Mouse_AERfactors96w","Mouse_AERfactorsALI")))
d$data=d$data %>% group_by(features.plot) %>% mutate("maxExp"=max(avg.exp)) %>% mutate("avgRelativeToMax"=avg.exp/maxExp)
g <- ggplot(d$data, aes(x=features.plot,y=id, size=pct.exp, col=avgRelativeToMax)) + geom_point(aes(fill=avgRelativeToMax),color='grey75', shape=21)
g <- g + theme_light() + scale_fill_gradient(name = "avgRelativeToMax",  low = "snow2", high = "#7d0025") #, limits=c(0,10)
g <- g + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
g <- g + ylab("")
g <- g + ggtitle(genesetName, subtitle=curCellType)

## Fig 3E: individual Violin plots (VehCtrl_AERfactors)
seuratObj=seuratObj_AER
genesToPlot.l=apply(aer_blastema.l$aer_blastema_Mouse,2, function(v){intersect(unique(v), rownames(seuratObj@assays$RNA@data))})
curGeneSet="aer_genes"
genesToPlot=genesToPlot.l[[curGeneSet]]
Idents(seuratObj)="CellTypeExt"
curCellType = "Basal Ectoderm"
seu = subset(seuratObj, idents=curCellType)
d=seu@assays$RNA@data
plots.l=list()
for(curGene in genesToPlot){
  dd=merge(data.frame("cellID"=colnames(d),"geneExp"=d[curGene,]),seu@meta.data, by.x="cellID",by.y=0, all.x=TRUE)
  g <- ggplot(dd , aes(x=orig.ident,y=geneExp)) +
    geom_violin(aes(fill=orig.ident), scale = "width") + geom_boxplot(width=0.1, color="grey33", alpha=0.2) +
    scale_fill_manual(values=colors_scRNAseq.l$specie_conditions) +
    theme_light() + xlab("") + ylab("Normalized expression") + ggtitle("", subtitle="") +
    theme(axis.text.x = element_text(angle=45,hjust=1,vjust=1), plot.title = element_text(hjust = 0.5, face="bold")) + NoLegend()
  g <- g + ggtitle(curGene)
  
  plots.l[[curGene]]=g
}

## Fig 3F: dot plot blastema genes (VehCtrl_AERfactors)
seuratObj=seuratObj_AER
genesForDotplots.df=aer_blastema.l$aer_blastema_Mouse
genesetName="blastema_genes" 
Idents(seuratObj)="CellTypeExt_CTgrp" 
seu = subset(seuratObj, idents = "CT_Sox9neg")
Idents(seu)="orig.ident"
d=DotPlot(seu, features=genesForDotplots.df[,genesetName] )
d$data$id=factor(d$data$id, levels=rev(c("Mouse_VehCtrl96w","Mouse_AERfactors96w","Mouse_AERfactorsALI")))
d$data=d$data %>% group_by(features.plot) %>% mutate("maxExp"=max(avg.exp)) %>% mutate("avgRelativeToMax"=avg.exp/maxExp)
g <- ggplot(d$data, aes(x=features.plot,y=id, size=pct.exp, col=avgRelativeToMax)) + geom_point(aes(fill=avgRelativeToMax),color='grey75', shape=21)
g <- g + theme_light() + scale_fill_gradient(name = "avgRelativeToMax",  low = "snow2", high = "#7d0025") #, limits=c(0,10)
g <- g + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
g <- g + ylab("")
g <- g + ggtitle(genesetName, subtitle="Sox9 low CT cells")


## Fig 3G: individual Violin plots (VehCtrl_AERfactors)
seuratObj=seuratObj_AER
genesToPlot.l=apply(aer_blastema.l$aer_blastema_Mouse,2, function(v){intersect(unique(v), rownames(seuratObj@assays$RNA@data))})
curGeneSet="blastema_genes"
genesToPlot=genesToPlot.l[[curGeneSet]] 
Idents(seuratObj)="CellTypeExt_CTgrp"
seu = subset(seuratObj, idents="CT_Sox9neg")
d=seu@assays$RNA@data
plots.l=list()
for(curGene in genesToPlot){
  dd=merge(data.frame("cellID"=colnames(d),"geneExp"=d[curGene,]),seu@meta.data, by.x="cellID",by.y=0, all.x=TRUE)
  g <- ggplot(dd , aes(x=orig.ident,y=geneExp)) +
    geom_violin(aes(fill=orig.ident), scale = "width") + geom_boxplot(width=0.1, color="grey33", alpha=0.2) +
    scale_fill_manual(values=colors_scRNAseq.l$specie_conditions) +
    theme_light() + xlab("") + ylab("Normalized expression") + ggtitle("", subtitle="") +
    theme(axis.text.x = element_text(angle=45,hjust=1,vjust=1), plot.title = element_text(hjust = 0.5, face="bold")) + NoLegend()
  g <- g + ggtitle(curGene)
  
  plots.l[[curGene]]=g
}



#*******************************#
########## Figure S9 #########
#*******************************#

## Fig S9B: QC plots (VehCtrl_AERfactors)
pathToData="./multiSpecies_limbRegeneration/data/scRNAseq_files/" #points to the GitHub folder (https://github.com/BICC-UNIL-EPFL/multiSpecies_limbRegeneration/)
infile=paste0(pathToData,"someNumbers.txt")
toPlot.df=t(read.delim(infile)) %>% data.frame()
toPlot.df$sampleName=rownames(toPlot.df)
toPlot.df$sampleName=factor(toPlot.df$sampleName, levels=c("Mouse_96w","Mouse_ALI","Xen_96w","Xen_ALI","Mouse_VehCtrl96w","Mouse_AERfactors96w","Mouse_AERfactorsALI"))
g <- ggplot(toPlot.df, aes(x=sampleName, y=median_gene_per_cell, fill=sampleName, label = median_gene_per_cell)) + geom_bar(stat="identity")
g <- g + scale_fill_manual(values=colors_scRNAseq.l$specie_conditions)
g <- g + theme_light() + ggtitle("Median gene per cell", subtitle="all samples")
g + geom_text(aes(label = median_gene_per_cell), vjust = -0.2) + NoLegend()

infile=paste0(pathToData,"nFeature_nCount.xlsx")
toPlot.df = readxl::read_excel(infile, sheet = "Mouse_AERfactors") %>% as.data.frame # sheet names are:  "Mouse_AERfactors" "Mouse_96w_ALI"    "Xen_96w_ALI"  
#Replace y=nCount_RNA by nFeature_RNA and percent.mt accordingly
g <- ggplot(toPlot.df, aes(x=sampleName, y=nCount_RNA, fill=sampleName)) 
g <- g + geom_flat_violin(
  position = position_nudge(x = 0, y = 0),
  adjust = 2, trim = TRUE, scale = "width", color = NA, alpha = 0.7
) 
g <- g + geom_boxplot(
  position = position_nudge(x = -.2, y = 0),
  outlier.shape = NA, width = .2, lwd = .2,
) 
g <- g + theme_light() 
g <- g + scale_fill_manual(values=colors_scRNAseq.l$specie_conditions)
g + ggtitle("nCount_RNA") + NoLegend()


## Fig S9C: umap samples (VehCtrl_AERfactors)
seuratObj=seuratObj_AER
p_cond <- DimPlot(seuratObj, reduction=curRed, group.by = "orig.ident", cols=colors_scRNAseq.l$specie_conditions_transparent,label=FALSE, order=TRUE)
p_cond


## Fig S9D: umap extended cell-types (VehCtrl_AERfactors)
seuratObj=seuratObj_AER
p_cellTypeExt <- DimPlot(seuratObj, reduction=curRed, group.by = "CellTypeExt",cols=colors_scRNAseq.l$cellTypeExt_Mouse,order=TRUE) + ggtitle("Mouse ALI_96w", subtitle="CellTypeExt") +  guides(color = guide_legend(override.aes = list(size=4), ncol=1) )
p_cellTypeExt


## Fig S9E: barplot extended cell-types (VehCtrl_AERfactors)
seuratObj=seuratObj_AER
ct=table(seuratObj$CellTypeExt,seuratObj$orig.ident)
ct.df=data.frame(ct[which(rowSums(ct)>0),])
colnames(ct.df)=c("cellType","condition","Freq")
ct.df$cellType=factor(ct.df$cellType, levels=order_CellTypeExt_Mouse);
colToUse=colors_scRNAseq.l$cellTypeExt_Mouse

g_cond_cellTypeExt <- ggplot(ct.df, aes(fill=cellType, y=Freq, x=condition)) + 
  geom_bar(position="fill", stat="identity") + theme_light() + scale_fill_manual(values=colToUse) + ylab("Population %") + ggtitle("Cell-Types Ext per condition",subtitle=seuratObj@project.name)
g_cond_cellTypeExt


## Fig S9F: heatmap marker genes (VehCtrl_AERfactors)
seuratObj=seuratObj_AER
genesToPlot=c("Msx1","Prrx1","Grem1", "Fgf10", "Ptch1" ,"Acta2","Thy1","Sox5", "Sox6","Sox9","Gdf5","Runx3","Pax7","Myf5","Pax3","Myog","Des","Rbm24","Myl1","Krt5","Krt14","Krt1","Epcam","Trp63","Postn","Lox","Frem2","Eya2","Sox7","Arg1","Zeb2","Csf1r","Spi1","Cd38","Cxcr4","Gata3","Dll4","Cdh4","Pecam1")

I.l=split(rownames(seuratObj@meta.data), f=seuratObj@meta.data[,"CellTypeExt"])
I.ll=lapply(I.l, function(v){sample(v,size=min(length(v),200))}) # down-sample to 200 cells per cell-type 
seurat_subset = subset(seuratObj, cells=unlist(I.ll))
seurat_subset$CellTypeExt = droplevels(seurat_subset$CellTypeExt)
DefaultAssay(seurat_subset)="RNA";
seurat_subset=ScaleData(seurat_subset, features=genesToPlot)

colToUse=colors_scRNAseq.l$cellTypeExt_Mouse
p <- DoHeatmap(seurat_subset, assay="RNA", features = genesToPlot, group.by = "CellTypeExt", group.colors = colToUse, 
               disp.min = -1, raster = FALSE, size = 3) 
p <- p + scale_fill_gradientn(limits = c(-1,3), colors = CustomPalette(low = "lightblue", high = "orangered3", mid = "snow", k=50), 
                              values = scales::rescale(c(-1, 0, 3))  ,na.value = "white") 
p <- p + theme(axis.text.y = element_text(size = 10)) + ggtitle(seuratObj@project.name, subtitle="CellTypeExt")
print(p + guides(color = "none" ) )


## Fig S9G: barplot cell-cycle vs. sample & cell-type (VehCtrl_AERfactors) 
seuratObj=seuratObj_AER
seuratObj$condition_CellType=paste(seuratObj$orig.ident, seuratObj$CellType,sep="_")
ct=table(seuratObj$condition_CellType,seuratObj$Phase)
ct.df=data.frame(ct[which(rowSums(ct)>0),])
order_grp=c("Mouse_VehCtrl96w","Mouse_AERfactors96w","Mouse_AERfactorsALI")

colnames(ct.df)=c("condition_CellType","Phase","Freq")
ct.df$cellType=sapply(as.character(ct.df$condition_CellType), function(s){unlist(strsplit(s,"_"))[3]})
ct.df$condition=sapply(as.character(ct.df$condition_CellType), function(s){paste(unlist(strsplit(s,"_"))[1:2],collapse="_")})
ct.df$condition=factor(ct.df$condition, levels=order_grp)
ct.df$Phase=factor(ct.df$Phase, levels=rev(c("G1", "G2M", "S")))
N=length(unique(ct.df$cellType))
N_cond=length(order_grp)
ct.df$condition_CellType=factor(ct.df$condition_CellType, levels=paste(rep(order_grp,each=N),rep(unique(ct.df$cellType),N_cond),sep="_"))
g <- ggplot(ct.df %>% filter(cellType %in% names(table(ct.df$cellType))), aes(fill=Phase, y=Freq, x=condition_CellType)) + 
  geom_bar(position="fill", stat="identity") + scale_fill_manual(values=colors_scRNAseq.l$Phase) 
g <- g + theme_light() + theme(axis.text.x = element_text(color = "black", size = 12, angle = 30, hjust = 1))
g_cellType_Phase <- g +  ggtitle("Phase",subtitle=seuratObj@project.name)
g_cellType_Phase <- g_cellType_Phase + scale_x_discrete(labels = as.character(sapply(order_grp, function(s){unlist(strsplit(s,"_"))[2]})), name = "")
g_cellType_Phase <- g_cellType_Phase + facet_grid(.~cellType, scales = "free_x")



#*******************************#
########## Figure S10 #########
#*******************************#
## Fig S10A: Dot plot blastema genes (VehCtrl_AERfactors)
seuratObj=seuratObj_AER
genesForDotplots.df=aer_blastema.l$aer_blastema_Mouse
genesetName="blastema_genes" 

Idents(seuratObj)="CellType" 
seu = subset(seuratObj, idents = "CT")
Idents(seu)="orig.ident"
d=DotPlot(seu, features=genesForDotplots.df[,genesetName] )
d$data$id=factor(d$data$id, levels=rev(c("Mouse_VehCtrl96w","Mouse_AERfactors96w","Mouse_AERfactorsALI")))
d$data=d$data %>% group_by(features.plot) %>% mutate("maxExp"=max(avg.exp)) %>% mutate("avgRelativeToMax"=avg.exp/maxExp)
g <- ggplot(d$data, aes(x=features.plot,y=id, size=pct.exp, col=avgRelativeToMax)) + geom_point(aes(fill=avgRelativeToMax),color='grey75', shape=21)
g <- g + theme_light() + scale_fill_gradient(name = "avgRelativeToMax",  low = "snow2", high = "#7d0025") #, limits=c(0,10)
g <- g + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
g <- g + ylab("")
g <- g + ggtitle(genesetName, subtitle="CT")


Idents(seuratObj)="CellTypeExt_CTgrp" 
seu = subset(seuratObj, idents = "CT_Sox9pos")
Idents(seu)="orig.ident"
d=DotPlot(seu, features=genesForDotplots.df[,genesetName] )
d$data$id=factor(d$data$id, levels=rev(c("Mouse_VehCtrl96w","Mouse_AERfactors96w","Mouse_AERfactorsALI")))
d$data=d$data %>% group_by(features.plot) %>% mutate("maxExp"=max(avg.exp)) %>% mutate("avgRelativeToMax"=avg.exp/maxExp)
g <- ggplot(d$data, aes(x=features.plot,y=id, size=pct.exp, col=avgRelativeToMax)) + geom_point(aes(fill=avgRelativeToMax),color='grey75', shape=21)
g <- g + theme_light() + scale_fill_gradient(name = "avgRelativeToMax",  low = "snow2", high = "#7d0025") #, limits=c(0,10)
g <- g + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
g <- g + ylab("")
g <- g + ggtitle(genesetName, subtitle="CT_Sox9pos")


## Fig S10B: Violin plot Fgf10 and Shh genes (VehCtrl_AERfactors)
seuratObj=seuratObj_AER
Idents(seuratObj)="CellTypeExt"
seu = subset(seuratObj, idents="Basal Ectoderm")
genesToPlot=c("Fgf10","Shh")
d=seu@assays$RNA@data
plots.l=list()
for(curGene in genesToPlot){
  dd=merge(data.frame("cellID"=colnames(d),"geneExp"=d[curGene,]),seu@meta.data, by.x="cellID",by.y=0, all.x=TRUE)
  g <- ggplot(dd , aes(x=orig.ident,y=geneExp)) +
    geom_violin(aes(fill=orig.ident), scale = "width") + geom_boxplot(width=0.1, color="grey33", alpha=0.2) +
    scale_fill_manual(values=colors_scRNAseq.l$specie_conditions) +
    theme_light() + xlab("") + ggtitle("", subtitle="") +
    theme(axis.text.x = element_text(angle=45,hjust=1,vjust=1), plot.title = element_text(hjust = 0.5, face="bold")) + NoLegend()
  g <- g + ggtitle(curGene)
  
  plots.l[[curGene]]=g
}


## Fig S10C: Dot plot AER genes (VehCtrl_AERfactors)
seuratObj=seuratObj_AER
genesForDotplots.df=aer_blastema.l$aer_blastema_Mouse
genesetName="aer_genes" 

Idents(seuratObj)="CellTypeExt_CTgrp" 
seu = subset(seuratObj, idents = "CT_Sox9neg")
Idents(seu)="orig.ident"
d=DotPlot(seu, features=genesForDotplots.df[[genesetName]] )
d$data$id=factor(d$data$id, levels=rev(c("Mouse_VehCtrl96w","Mouse_AERfactors96w","Mouse_AERfactorsALI")))
g <- ggplot(d$data, aes(x=features.plot,y=id, size=pct.exp, col=avg.exp)) + geom_point(aes(fill=avg.exp),color='grey75', shape=21)
g <- g + theme_light() + scale_fill_gradient(name = "avg exp",  low = "snow2", high = "#7d0025") #, limits=c(0,10)
g <- g + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
g <- g + ylab("")
g <- g + ggtitle(genesetName, subtitle="CT_Sox9neg")



#*******************************#
########## Figure 4 #########
#*******************************#

## Fig 4A: UMAP RNA colored by cell-types (scMultiome)
# Code: see Fig S11D


## Fig 4B & 4C: mean profiles AER and random genes (scMultiome)
pathToData="./multiSpecies_limbRegeneration/data/multiomics_files/" #points to the GitHub folder (https://github.com/BICC-UNIL-EPFL/multiSpecies_limbRegeneration/)
infile=paste0(pathToData,"compareSignals_AER_random_data.xlsx")
toPlot.df=readxl::read_xlsx(infile) %>% as.data.frame()
toPlot.df$geneList=factor(toPlot.df$geneList, levels=c("AERgenes","200random"))
curLabels=seq(-1000, 1000, length=21)
g <- ggplot(toPlot.df %>% filter(Tissue=="AER_and_SurfaceEctoderm_signacNorm"), aes(x=pos_n, y=meanNormScore, col=condTP)) + geom_point(size=0.2) + geom_line() 
g <- g + scale_color_manual(values=colors_multiomics.l$condition_TP)
g <- g + geom_vline(xintercept = 21, col="cornflowerblue", linetype="dashed") + theme_light() 
g <- g + scale_x_continuous("TSS position", breaks=seq(1,41,by=2), labels=curLabels) 
g <- g + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), panel.grid.minor = element_blank(),panel.grid.major.x = element_blank())
g <- g + facet_grid(Tissue~geneList, scales="free_y")


## Fig 4D: coverage plots (scMultiome)
Idents(seu_multiomics)="MainCellTypesExt"
seu_multiomics@active.assay="ATAC"
seuratObj_subset=subset(seu_multiomics, ident="AER_and_SurfaceEctoderm")
plots.l=list()
for(curGene in c("Fgf8","Sp6","Sp8","Msx2")){
  plots.l[[curGene]] = CoveragePlot(seuratObj_subset, region = curGene, features = curGene, group.by= "condition_TP", extend.upstream = 1500, extend.downstream = 1000) & scale_fill_manual(values=colors_multiomics.l$condition_TP)
  }

cowplot::plot_grid(plotlist = plots.l)


## Fig 4H: umap and barplot cell populations (scMultiome)
# Code: see Fig S14C and S14D



#*******************************#
########## Figure S11 #########
#*******************************#

## Fig S11A: QC raw (scMultiome)
pathToData="./multiSpecies_limbRegeneration/data/multiomics_files/" #points to the GitHub folder (https://github.com/BICC-UNIL-EPFL/multiSpecies_limbRegeneration/)
infile=paste0(pathToData,"someNumbers_multiomics.txt")
toPlot.df=t(read.delim(infile)) %>% data.frame()
toPlot.df$sampleName=rownames(toPlot.df)
toPlot.df$sampleName=factor(toPlot.df$sampleName,levels=c("WT_0h","HPA_6h","HPA_24h","HPA_72h"))
g <- ggplot(toPlot.df, aes(x=sampleName, y=median_gene_per_cell, fill=sampleName, label = median_gene_per_cell)) + geom_bar(stat="identity")
g <- g + scale_fill_manual(values=colors_multiomics.l$condition_TP)
g <- g + theme_light() + ggtitle("Median gene per cell", subtitle="all samples")
g + geom_text(aes(label = median_gene_per_cell), vjust = -0.2) + NoLegend()

toPlot.df=read.csv("./QC_raw/nFeature_nCount_multiomics.csv")
toPlot.df$sampleName=factor(toPlot.df$sampleName,levels=c("WT_0h","HPA_6h","HPA_24h","HPA_72h"))
#Replace y=nCount_RNA by nFeature_RNA and percent.mt accordingly
g <- ggplot(toPlot.df, aes(x=sampleName, y=nCount_RNA, fill=sampleName)) 
g <- g + geom_flat_violin(
  position = position_nudge(x = 0, y = 0),
  adjust = 2, trim = TRUE, scale = "width", color = NA, alpha = 0.7
) 
g <- g + geom_boxplot(
  position = position_nudge(x = -.2, y = 0),
  outlier.shape = NA, width = .2, lwd = .2,
) 
g <- g + theme_light() + ylim(c(0,15000))
g <- g + scale_fill_manual(values=colors_multiomics.l$condition_TP)
g <- g + ggtitle("nCount_RNA") + NoLegend()


## Fig S11B: QC raw ATAC (scMultiome)
toPlot.df=seuratObj_merged@meta.data
toPlot.df$sampleName=factor(toPlot.df$orig.ident, levels=c("WT_0h","HPA_6h","HPA_24h","HPA_72h"))

#Replace y=nCount_ATAC by nFeature_ATAC, TSS.enrichment, and nucleosome_signal accordingly
g <- ggplot(toPlot.df, aes(x=sampleName, y=nCount_ATAC  , fill=sampleName)) #+ geom_violin()
g <- g + geom_flat_violin(
  position = position_nudge(x = 0, y = 0),
  adjust = 2, trim = TRUE, scale = "width", color = NA, alpha = 0.7
) 
g <- g + geom_boxplot(
  position = position_nudge(x = -.2, y = 0),
  outlier.shape = NA, width = .2, lwd = .2,
) 
g <- g + theme_light() + ylim(c(0,75000))
g <- g + scale_fill_manual(values=colors_multiomics.l$condition_TP)
g <- g + ggtitle("nCount_ATAC") + NoLegend()


### Fig S11C, S11D and S11E: UMAPs (RNA, ATAC and combined), colored by cell-type & conditions  (scMultiome)
# colored per CellTypeExt
p_cellTypeExt_rna <- DimPlot(seu_multiomics, reduction="umap.rna" , group.by="CellTypeExt", label=FALSE, repel = TRUE, cols=colors_multiomics.l$CellTypeExt) + ggtitle("RNA", subtitle="CellTypeExt")
p_cellTypeExt_atac <- DimPlot(seu_multiomics, reduction="umap.atac", group.by="CellTypeExt", label=FALSE, repel = TRUE, cols=colors_multiomics.l$CellTypeExt, order=TRUE) + ggtitle("ATAC", subtitle="CellTypeExt")
p_cellTypeExt_rnaAtac <- DimPlot(seu_multiomics, reduction="wnn.umap" , group.by="CellTypeExt", label=FALSE, repel = TRUE, cols=colors_multiomics.l$CellTypeExt) + ggtitle("ATAC + RNA integrated", subtitle="CellTypeExt")

# colored per condition
p_condTP_rna=DimPlot(seu_multiomics, reduction="umap.rna", group.by = "condition_TP", cols=colors_multiomics.l$condition_TP)
p_condTP_atac=DimPlot(seu_multiomics, reduction="umap.atac", group.by = "condition_TP", cols=colors_multiomics.l$condition_TP)
p_condTP_rnaAtac=DimPlot(seu_multiomics, reduction="wnn.umap", group.by = "condition_TP", cols=colors_multiomics.l$condition_TP, order=TRUE)


### Fig S11F:barplot condition x cell-type Extended (scMultiome)
ct=table(seu_multiomics$CellTypeExt,seu_multiomics$condition_TP)
ct.df=data.frame(ct)
colnames(ct.df)=c("cellType","condition","Freq")
order_CellTypeExt = c("CT",paste("CT_",1:18,sep=""),"AER","Surface Ectoderm","Muscle","Blood","Immune","Endothelial","Pericyte","Schwann","Neuron","Erythrocytes")

ct.df$cellType=factor(ct.df$cellType, levels=order_CellTypeExt[-1])
g_cond_cellTypeExt <- ggplot(ct.df, aes(fill=cellType, y=Freq, x=condition)) + 
  geom_bar(position="fill", stat="identity") + theme_light() + scale_fill_manual(values=colors_multiomics.l$CellTypeExt) + ggtitle("Cell-TypesExt per condition",subtitle=seu_multiomics@project.name)
g_cond_cellTypeExt


### Fig S11G: heatmap markers (scMultiome)
I.l=split(rownames(seu_multiomics@meta.data), f=seu_multiomics@meta.data[,"CellTypeExt"])
I.ll=lapply(I.l, function(v){sample(v,size=min(length(v),200))}) # down-sample to 200 cells per cell-type 
seurat_subset = subset(seu_multiomics, cells=unlist(I.ll))
seurat_subset$CellTypeExt = droplevels(seurat_subset$CellTypeExt)

genesToPlot <- c("Csf1r","Epor","Spi1","Cd48","Krt1","Thy1","Sp6","Fgf8","Sp9","Arg1","Gpr34","Cdh4","Msx1","Wnt5a","Epcam","Krt14","Gata3","Krt5","Frem2","Trp63","Mbp","Myl1","Myog","Acta2","Rbm24","Des","Krt10","Pecam1","Dll4","Sox7","Cxcr4","Cd38","Prdm1","Zeb2","Pax3","Myf5","Pax7","Eya2","Runx3","Ptch1","Col8a2","Sox6","Sox5","Sox9","Gdf5","Pdgfra","Prrx1","Pdgfrb","Grem1","Postn","Fgf10","Lox")
DefaultAssay(seurat_subset)="RNA";
seurat_subset=ScaleData(seurat_subset, features=genesToPlot)

colToUse=colors_multiomics.l$CellTypeExt
p <- DoHeatmap(seurat_subset, assay="RNA", features = genesToPlot, group.by = "CellTypeExt", group.colors = colToUse, disp.min = -1, raster = FALSE, size = 3) + scale_fill_gradientn(limits = c(-1,3), colors = CustomPalette(low = "lightblue", high = "orangered3", mid = "snow", k=50), values = scales::rescale(c(-1, 0, 3))  ,na.value = "white") #indianred1 #breaks=seq(-1, 3, length.out = 5), labels = seq(-1, 3, by=1)
p <- p + theme(axis.text.y = element_text(size = 10)) + ggtitle(seurat_subset@project.name, subtitle="CellTypeExt")
print(p + guides(color = "none" ) )


#*******************************#
########## Figure S12 #########
#*******************************#

### Fig S12A: Violin plot DA peaks intact vs. 6 hpa (scMultiome)
pathToData="./multiSpecies_limbRegeneration/data/multiomics_files/" #points to the GitHub folder (https://github.com/BICC-UNIL-EPFL/multiSpecies_limbRegeneration/)
infile=paste0(pathToData,"data_VlnPlot_FigS12A.csv.gz")
da_peaks_filtered <- read.csv(gzfile(infile))

g <- ggplot(da_peaks_filtered , aes(x = cellType, y = avg_log2FC, fill=cellType)) +
  geom_violin() +
  theme_minimal() +
  labs(y = "Log2 Fold Change", x = "") +
  scale_fill_manual(values = setNames(c("cadetblue3","#476EA9","darkolivegreen3","firebrick3"),c("connectiveTissues","AER_and_SurfaceEctoderm","Muscle","Endothelial")))
g <- g + geom_hline(yintercept=0, linetype="dashed", color = "grey33", size=1)
g <- g+ggtitle("Differential Chromatin Accessibility") + NoLegend()

### Fig S12B: Feature Distribution - peaks in all cells (scMultiome)
pathToData="./multiSpecies_limbRegeneration/data/multiomics_files/" #points to the GitHub folder (https://github.com/BICC-UNIL-EPFL/multiSpecies_limbRegeneration/)
rdsfile=paste0(pathToData,"peaks_cleaned_annot_l.rds")
peaks_cleaned_annot.l=readRDS(rdsfile) 

anno.df.l=sapply(names(peaks_cleaned_annot.l), function(curCat){d=peaks_cleaned_annot.l[[curCat]]@annoStat; return(data.frame("Category"=rep(curCat, nrow(d)),d) )}, simplify=FALSE)
anno.df=do.call(rbind, anno.df.l)
allAnnoLevels=c("Promoter","5' UTR","3' UTR","1st Exon","Other Exon","1st Intron","Other Intron","Downstream (<=3kb)","Distal Intergenic")
anno.df$Feature=factor(anno.df$Feature, levels=allAnnoLevels)
anno.df$Category=factor(anno.df$Category, levels=c("0h","6h","24h","72h"))
anno.df$Frequency=as.numeric(anno.df$Frequency)
p <- ggplot(anno.df, aes_string(x = "Category", fill = "Feature", y = "Frequency"))
p <- p + geom_bar(stat="identity") + coord_flip() + theme_bw()
p <- p + ylab("Percentage(%)") + xlab("") + ggtitle("Feature Distribution - peaks all cellTypes")
p <- p + scale_fill_manual(values=rev(col_featDist[1:length(unique(anno.df$Feature))]), guide=guide_legend(reverse=TRUE))


### Fig S12C: average profiles AER & random genes across time points
pathToData="./multiSpecies_limbRegeneration/data/multiomics_files/" #points to the GitHub folder (https://github.com/BICC-UNIL-EPFL/multiSpecies_limbRegeneration/)
infile=paste0(pathToData,"compareSignals_AER_random_data.xlsx")
toPlot.df=readxl::read_xlsx(infile) %>% as.data.frame()
toPlot.df$geneList=factor(toPlot.df$geneList, levels=c("AERgenes","200random"))
toPlot.df$condTP=factor(toPlot.df$condTP, levels=c("WT_0h","HPA_6h","HPA_24h","HPA_72h"))

curLabels=seq(-1000, 1000, length=21)
gg <- ggplot(toPlot.df%>% filter(Tissue=="AER_and_SurfaceEctoderm_signacNorm"), aes(x=pos_n, y=meanNormScore, col=geneList)) + geom_point(size=0.3) + geom_line() 
gg <- gg + geom_vline(xintercept = 21, col="cornflowerblue", linetype="dashed") + theme_light() 
gg <- gg + scale_x_continuous("TSS position", breaks=seq(1,41,by=2), labels=curLabels) 
gg <- gg + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), panel.grid.minor = element_blank(),panel.grid.major.x = element_blank()) 
gg <- gg + facet_grid(Tissue~condTP, scales="free_y") 


### Fig S12D ( (enriched-) Heatmap depicting AER gene set chromatin accessibility at each time point in basal epithelial cells)
pathToData="./multiSpecies_limbRegeneration/data/multiomics_files/" #points to the GitHub folder (https://github.com/BICC-UNIL-EPFL/multiSpecies_limbRegeneration/)
rdsfile=paste(pathToData,"mat_enrichedHeatmap_FigureS12D_l.rds")
mat_EnrichedHeatmap.l=readRDS(rdsfile)
bw_col_fun = circlize::colorRamp2(c(0, 10, 80), c("snow", "grey90", "grey10"))
plots.l=lapply(mat_EnrichedHeatmap.l, function(mat){
  EnrichedHeatmap::EnrichedHeatmap(mat, col = bw_col_fun,  name=NULL, show_row_names = TRUE, use_raster = TRUE)
})
draw(plots.l[[1]]+plots.l[[2]]+plots.l[[3]]+plots.l[[4]])


### Prepare data for Fig S12E & S12F. Enrichment profiles AER gene and random genes at 1 dpa. 
pathToGEO="files_GSE314013/" # points to the data from GEO 
infile=paste0(pathToGEO,"signal_matrix_20bp.gz")
mat <- read.table(gzfile(infile), skip=1) # +/-1kb, 20bp bins (2 x 50 bins)

signal_matrix <- as.matrix(mat[,7:ncol(mat)])
mat_coords <- mat[,1:6]
colnames(mat_coords)[4]="geneName"

AER_genes <- c("Arid3b","Axin2","Bambi","Bmp2","Bmp4","Bmp7","Cd44","Cdh1","Dkk1","Dlx2","Dlx5","Dlx6","En1","Fgf17","Fgf4","Fgf8","Fgf9","Fn1","Gja1","Jag2","Lef1","Lgr5","Msx1","Msx2","Notch1","Notch2","Notch3","Notch4","Rfng","Rspo2","Sp6","Sp8","Sp9","Tcf7","Trp63","Wnt3","Wnt5a")
I_AERgenes <- which(mat[,4] %in% AER_genes)
signal_matrix=signal_matrix[I_AERgenes,]
mat_coords=mat_coords[I_AERgenes,]

bins_per_sample <- 100
n_samples <- 6
sample_matrices <- lapply(1:n_samples, function(i) {
  start_col <- (i - 1) * bins_per_sample + 1
  end_col <- i * bins_per_sample
  signal_matrix[, start_col:end_col]
})
names(sample_matrices) <- c("96w-1", "96w-2", "96w-3","ALI-1","ALI-2","ALI-3")

data_meanProfiles.l=lapply(sample_matrices, function(curMat){
  rownames(curMat)=mat_coords$geneName
  grouped_means <- tapply(1:nrow(curMat), mat_coords$grp, function(rows) {
    colMeans(curMat[rows, , drop = FALSE])})
  grouped_means.df <- do.call(rbind, grouped_means)
  return(grouped_means.df)
})

data_meanProfiles_AER.df <- do.call(rbind,lapply(data_meanProfiles.l, function(d){d["AER",]}))
data_meanProfiles_random200.df <- do.call(rbind,lapply(data_meanProfiles.l, function(d){d["random200",]}))

bins <- seq(-1000, 1000, by = 20)  
bins = bins[1:(length(bins)-1)] # 100 bins
bins[length(bins)]=1000


## Fig S12E:
toPlot=rbind( melt(data_meanProfiles_AER.df), melt(data_meanProfiles_random200.df))
toPlot$geneGrp <- rep(c("AER","random200"), each=nrow(melt(data_meanProfiles_AER.df)))
toPlot$variable=sapply(as.character(toPlot$Var1), function(s){unlist(strsplit(s,"-"))[1]})
toPlot$replicate=sapply(as.character(toPlot$Var1), function(s){unlist(strsplit(s,"-"))[2]})
toPlot$bins=rep(bins, each=6)

g <- ggplot(toPlot, aes(x=bins, y=value, col=geneGrp)) 
g <- g + geom_vline(xintercept=0, lwd=0.7, linetype=2, col="grey")
g <- g + geom_line(lwd=1.2)
g <- g + theme_light() 
g <- g + theme(
  panel.grid.major.x = element_blank(),
  panel.grid.minor.x = element_blank(),
  panel.grid.minor.y = element_blank()
)
g <- g + ggtitle("mean profiles around TSS, per replicates", subtitle="AER genes vs. random genes")
g <- g + facet_grid(variable~replicate) + theme(legend.position = "top")


### Fig S12F: 
tick_positions <- seq(-1000, 1000, by = 500)
tick_indices <- match(tick_positions, bins)
label_row <- rep("", length(bins))
tick_labels <- paste0(tick_positions / 1000, "kb")
tick_labels[which(tick_labels=="0kb")]="TSS"
label_row[tick_indices] <- tick_labels

tmpH.l=list()
for(curSample in names(sample_matrices)){
  curMat=sample_matrices[[curSample]]
  rownames(curMat)=mat_coords$geneName
  print(max_val)
  
  # Create a bottom annotation with custom labels
  bottom_annot <- HeatmapAnnotation(
    Position = anno_text(label_row, rot = 90, gp = gpar(fontsize = 8)),
    annotation_name_side = "left",
    show_annotation_name = FALSE
  )
  
  tmpH=ComplexHeatmap::Heatmap(curMat,
                               name = curSample, 
                               border_gp = gpar(col = "grey", lty = 1),
                               bottom_annotation = bottom_annot,
                               cluster_rows = TRUE,
                               cluster_columns = FALSE,
                               show_row_names = TRUE, 
                               show_column_names = FALSE,
                               col = colorRamp2(c(0, max_val), c("white", "grey35")))
  tmpH.l[[curSample]]=tmpH
} 
draw(tmpH.l$`96w-1`+tmpH.l$`ALI-1`+tmpH.l$`96w-2`+tmpH.l$`ALI-2`+tmpH.l$`96w-3`+tmpH.l$`ALI-3`) 



#*******************************#
########## Figure S15 #########
#*******************************#

## Fig S15B: QC plots scRNA-seq (Mouse 1dpa)
# barplot median gene per cell
pathToData="./multiSpecies_limbRegeneration/data/scRNAseq_files/" #points to the GitHub folder (https://github.com/BICC-UNIL-EPFL/multiSpecies_limbRegeneration/)
infile=paste0(pathToData,"someNumbers.txt")

toPlot.df=t(read.delim(infile)) %>% data.frame()
toPlot.df$sampleName=rownames(toPlot.df)
toPlot.df$sampleName=factor(toPlot.df$sampleName, levels=c("Mouse_96w","Mouse_ALI","Xen_96w","Xen_ALI","Mouse_VehCtrl96w","Mouse_AERfactors96w","Mouse_AERfactorsALI"))
g <- ggplot(toPlot.df %>% filter(sampleName %in% c("Mouse_96w","Mouse_ALI")), aes(x=sampleName, y=median_gene_per_cell, fill=sampleName, label = median_gene_per_cell)) + geom_bar(stat="identity")
g <- g + scale_fill_manual(values=colors_scRNAseq.l$specie_conditions)
g <- g + theme_light() + ggtitle("Median gene per cell")
g + geom_text(aes(label = median_gene_per_cell), vjust = -0.2) + NoLegend()

# violin plots QC metrics
pathToData="./multiSpecies_limbRegeneration/data/scRNAseq_files/" #points to GitHub folder (https://github.com/BICC-UNIL-EPFL/multiSpecies_limbRegeneration/)
filename=paste0(pathToData,"nFeature_nCount.xlsx")

toPlot.df = readxl::read_excel(filename, sheet = "Mouse_96w_ALI") %>% as.data.frame # sheet names are:  "Mouse_AERfactors" "Mouse_96w_ALI"    "Xen_96w_ALI"  
#Replace y=nCount_RNA by nFeature_RNA and percent.mt accordingly
g <- ggplot(toPlot.df, aes(x=sampleName, y=nCount_RNA, fill=sampleName)) 
g <- g + geom_flat_violin(
  position = position_nudge(x = 0, y = 0),
  adjust = 2, trim = TRUE, scale = "width", color = NA, alpha = 0.7
) 
g <- g + geom_boxplot(
  position = position_nudge(x = -.2, y = 0),
  outlier.shape = NA, width = .2, lwd = .2,
) 
g <- g + theme_light() 
g <- g + scale_fill_manual(values=colors_scRNAseq.l$specie_conditions)
g + ggtitle("nCount_RNA") + NoLegend()

## Fig S15C:  umaps (Mouse 1dpa)
seuratObj=seuratObj_Mouse
p_cond <- DimPlot(seuratObj, reduction="umap.int", group.by = "orig.ident", cols=colors_scRNAseq.l$specie_conditions_transparent,label=FALSE, order=TRUE)
p_cond
p_cellTypeExt <- DimPlot(seuratObj, reduction="umap.int", group.by = "CellTypeExt",cols=colors_scRNAseq.l$cellTypeExt_Mouse,order=TRUE) + ggtitle("Mouse ALI_96w", subtitle="CellTypeExt") +  guides(color = guide_legend(override.aes = list(size=4), ncol=1) )
p_cellTypeExt


## Fig S15D: barplot cell types (Mouse 1dpa)
seuratObj=seuratObj_Mouse
ct=table(seuratObj$CellTypeExt,seuratObj$orig.ident)
ct.df=data.frame(ct[which(rowSums(ct)>0),])
colnames(ct.df)=c("cellType","condition","Freq")
ct.df$cellType=factor(ct.df$cellType, levels=order_CellTypeExt_Mouse);
colToUse=colors_scRNAseq.l$cellTypeExt_Mouse

g_cond_cellTypeExt <- ggplot(ct.df, aes(fill=cellType, y=Freq, x=condition)) + 
  geom_bar(position="fill", stat="identity") + theme_light() + scale_fill_manual(values=colToUse) + ylab("Population %") + ggtitle("Cell-Types Ext per condition",subtitle=seuratObj@project.name)
g_cond_cellTypeExt



## ### Fig S15E: heatmap marker genes (Mouse 1dpa)
seuratObj=seuratObj_Mouse
genesToPlot=c("Msx1","Prrx1","Grem1", "Fgf10", "Ptch1" ,"Acta2","Sox5", "Sox6","Sox9","Gdf5","Runx3","Pax7","Pax3","Myog","Des","Rbm24","Myl1","Krt5","Krt14","Epcam","Trp63","Postn","Frem2","Prdm1","Eya2","Sox7","Arg1","Zeb2","Spi1","Cxcr4","Gata3","Dll4","Cdh4","Pecam1")

I.l=split(rownames(seuratObj@meta.data), f=seuratObj@meta.data[,"CellTypeExt"])
I.ll=lapply(I.l, function(v){sample(v,size=min(length(v),200))}) # down-sample to 200 cells per cell-type 
seurat_subset = subset(seuratObj, cells=unlist(I.ll))
seurat_subset$CellTypeExt = droplevels(seurat_subset$CellTypeExt)
DefaultAssay(seurat_subset)="RNA";
seurat_subset=ScaleData(seurat_subset, features=genesToPlot)

colToUse=colors_scRNAseq.l$cellTypeExt_Mouse
p <- DoHeatmap(seurat_subset, assay="RNA", features = genesToPlot, group.by = "CellTypeExt", group.colors = colToUse, 
               disp.min = -1, raster = FALSE, size = 3) 
p <- p + scale_fill_gradientn(limits = c(-1,3), colors = CustomPalette(low = "lightblue", high = "orangered3", mid = "snow", k=50), 
                              values = scales::rescale(c(-1, 0, 3))  ,na.value = "white") 
p <- p + theme(axis.text.y = element_text(size = 10)) + ggtitle(seuratObj@project.name, subtitle="CellTypeExt")
print(p + guides(color = "none" ) )



### Fig S15F: barplot sample+cell-type x Phase (Mouse 1dpa) 
seuratObj=seuratObj_Mouse
seuratObj$condition_CellType=paste(seuratObj$orig.ident, seuratObj$CellType,sep="_")
ct=table(seuratObj$condition_CellType,seuratObj$Phase)
ct.df=data.frame(ct[which(rowSums(ct)>0),])
order_grp=order_grp=c("Mouse_96w","Mouse_ALI")

colnames(ct.df)=c("condition_CellType","Phase","Freq")
ct.df$cellType=sapply(as.character(ct.df$condition_CellType), function(s){unlist(strsplit(s,"_"))[3]})
ct.df$condition=sapply(as.character(ct.df$condition_CellType), function(s){paste(unlist(strsplit(s,"_"))[1:2],collapse="_")})
ct.df$condition=factor(ct.df$condition, levels=order_grp)
ct.df$Phase=factor(ct.df$Phase, levels=rev(c("G1", "G2M", "S")))
N=length(unique(ct.df$cellType))
N_cond=length(order_grp)
ct.df$condition_CellType=factor(ct.df$condition_CellType, levels=paste(rep(order_grp,each=N),rep(unique(ct.df$cellType),N_cond),sep="_"))
g <- ggplot(ct.df %>% filter(cellType %in% names(table(ct.df$cellType))[-1]), aes(fill=Phase, y=Freq, x=condition_CellType)) + 
  geom_bar(position="fill", stat="identity") + scale_fill_manual(values=colors_scRNAseq.l$Phase) 
g <- g + theme_light() + theme(axis.text.x = element_text(color = "black", size = 12, angle = 30, hjust = 1))
g_cellType_Phase <- g +  ggtitle("Phase",subtitle=seuratObj@project.name)
g_cellType_Phase <- g_cellType_Phase + scale_x_discrete(labels = as.character(sapply(order_grp, function(s){unlist(strsplit(s,"_"))[2]})), name = "")
g_cellType_Phase <- g_cellType_Phase + facet_grid(.~cellType, scales = "free_x")



#*******************************#
########## Figure S17 #########
#*******************************#

### Fig S17A: dot plot AER gene expression (Mouse 1dpa) 
seuratObj=seuratObj_Mouse
genesToPlot.l=apply(aer_blastema.l$aer_blastema_Mouse, 2, function(v){intersect(unique(v), rownames(seuratObj@assays$RNA@data))})
dd=DotPlot(seuratObj, assay="RNA", features=genesToPlot.l$aer_genes, group.by="orig.ident", scale=FALSE) 
dd.df=data.frame(dd$data)
g <- ggplot(dd.df, aes(x=features.plot,y=id, size=pct.exp, col=avg.exp)) + geom_point(aes(fill=avg.exp),color='grey75', shape=21)
g <- g + theme_light() + scale_fill_gradient(name = "avg exp",  low = "snow2", high = "#7d0025") 
g <- g + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
g <- g + ylab("")
g <- g + ggtitle("AER genes", subtitle="Mouse")


### Fig S17B: dot plot blastema gene expression (Mouse 1dpa) 
dd=DotPlot(seuratObj, assay="RNA", features=sort(genesToPlot.l$blastema_genes), group.by="orig.ident", scale=FALSE) 
dd.df=data.frame(dd$data)
g <- ggplot(dd.df, aes(x=features.plot,y=id, size=pct.exp, col=avg.exp)) + geom_point(aes(fill=avg.exp),color='grey75', shape=21)
g <- g + theme_light() + scale_fill_gradient(name = "avg exp",  low = "snow2", high = "#7d0025")
g <- g + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
g <- g + ylab("")
g <- g + ggtitle("Blastema genes", subtitle="Mouse")



#*******************************#
########## Figure S20 #########
#*******************************#

## Fig S20B: QC plots scRNA-seq (Xen 1dpa)
# barplot median gene per cell
pathToData="./multiSpecies_limbRegeneration/data/scRNAseq_files/" #points to the GitHub folder (https://github.com/BICC-UNIL-EPFL/multiSpecies_limbRegeneration/)
infile=paste0(pathToData,"someNumbers.txt")

toPlot.df=t(read.delim(infile)) %>% data.frame()
toPlot.df$sampleName=rownames(toPlot.df)
toPlot.df$sampleName=factor(toPlot.df$sampleName, levels=c("Mouse_96w","Mouse_ALI","Xen_96w","Xen_ALI","Mouse_VehCtrl96w","Mouse_AERfactors96w","Mouse_AERfactorsALI"))
g <- ggplot(toPlot.df %>% filter(sampleName %in% c("Xen_96w","Xen_ALI")), aes(x=sampleName, y=median_gene_per_cell, fill=sampleName, label = median_gene_per_cell)) + geom_bar(stat="identity")
g <- g + scale_fill_manual(values=colors_scRNAseq.l$specie_conditions)
g <- g + theme_light() + ggtitle("Median gene per cell")
g + geom_text(aes(label = median_gene_per_cell), vjust = -0.2) + NoLegend()

# violin plots QC metrics
pathToData="./multiSpecies_limbRegeneration/data/scRNAseq_files/" #points to GitHub folder (https://github.com/BICC-UNIL-EPFL/multiSpecies_limbRegeneration/)
filename=paste0(pathToData,"nFeature_nCount.xlsx")

toPlot.df = readxl::read_excel(infile, sheet = "Xen_96w_ALI") %>% as.data.frame # sheet names are:  "Mouse_AERfactors" "Mouse_96w_ALI"    "Xen_96w_ALI"  
#Replace y=nCount_RNA by nFeature_RNA and percent.mt accordingly
g <- ggplot(toPlot.df, aes(x=sampleName, y=nCount_RNA, fill=sampleName)) 
g <- g + geom_flat_violin(
  position = position_nudge(x = 0, y = 0),
  adjust = 2, trim = TRUE, scale = "width", color = NA, alpha = 0.7
) 
g <- g + geom_boxplot(
  position = position_nudge(x = -.2, y = 0),
  outlier.shape = NA, width = .2, lwd = .2,
) 
g <- g + theme_light() 
g <- g + scale_fill_manual(values=colors_scRNAseq.l$specie_conditions)
g + ggtitle("nCount_RNA") + NoLegend()


## Fig S20C: umaps  (Xen 1dpa)
seuratObj=seuratObj_Xen
order_CellTypeExt_Xen=c(paste("CT_",1:15,sep=""),"Basal Ectoderm","Surface Ectoderm","Differentiating Ectoderm","Muscle","Blood","Immune","Endothelial","Pericytes","Schwann","Neuron")
seuratObj$CellTypeExt=factor(seuratObj$CellTypeExt, levels=order_CellTypeExt_Xen)
curRed="umap.int"
p_cond <- DimPlot(seuratObj, reduction=curRed, group.by = "orig.ident", cols=colors_scRNAseq.l$specie_conditions_transparent,label=FALSE, order=TRUE)
p_cond
p_cellTypeExt <- DimPlot(seuratObj, reduction=curRed, group.by = "CellTypeExt",cols=colors_scRNAseq.l$cellTypeExt_Xen,order=TRUE) + ggtitle("Xen ALI_96w", subtitle="CellTypeExt") +  guides(color = guide_legend(override.aes = list(size=4), ncol=1) )
p_cellTypeExt


## Fig S20D: barplot cell types (Xen 1dpa)
seuratObj=seuratObj_Xen
ct=table(seuratObj$CellTypeExt,seuratObj$orig.ident)
ct.df=data.frame(ct[which(rowSums(ct)>0),])
colnames(ct.df)=c("cellType","condition","Freq")
ct.df$cellType=factor(ct.df$cellType, levels=order_CellTypeExt_Xen);
colToUse=colors_scRNAseq.l$cellTypeExt_Xen

g_cond_cellTypeExt <- ggplot(ct.df, aes(fill=cellType, y=Freq, x=condition)) + 
  geom_bar(position="fill", stat="identity") + theme_light() + scale_fill_manual(values=colToUse) + ylab("Population %") + ggtitle("Cell-Types Ext per condition",subtitle=seuratObj@project.name)
g_cond_cellTypeExt


### Fig S20E: heatmap marker genes (Xen 1dpa)
genesToPlot=c("msx1.L","grem1.L","ptch1.S" ,"acta2.L", "thy1.L", "sox5.L", "sox6.S","sox9.S","bmp10.S","col8a2.S","gdf5.L","otos.S", "sox10.S","runx3.S","pax7.L","myl1.S","Myh.L","epcam.S",
              "postn.L","lox.L","tp73.L","frem2.l","sp6.S","sp9.L","prdm1.S","eya2.S","foxa1.L","atp6v1b1.L","znf750.L","sox7.L","muc1.S","arg1.L","Znf50.L","sox7.L","muc1.S","arg1.L","foxi1.L",
              "zeb2.L","gpr34.L","spi1.L","lyz.L","il7r.L","cd38.L","cxcr4.S","ccr7.S","cxcr3.L","cd74.L","cd79a.S","gata3.S","dll4.L","pecam1.L","mbp.S","gata1.L","sox9.L")

seuratObj=seuratObj_Xen
genesForDotplots.l=genesForDotplots_Xen.l
genesForDotplots_in.l=sapply(names(genesForDotplots.l), function(curList){
  intersect(genesForDotplots.l[[curList]], rownames(seuratObj@assays$RNA))
})

I.l=split(rownames(seuratObj@meta.data), f=seuratObj@meta.data[,"CellTypeExt"])
I.ll=lapply(I.l, function(v){sample(v,size=min(length(v),200))}) # down-sample to 200 cells per cell-type 
seurat_subset = subset(seuratObj, cells=unlist(I.ll))
seurat_subset$CellTypeExt = droplevels(seurat_subset$CellTypeExt)
DefaultAssay(seurat_subset)="RNA";
seurat_subset=ScaleData(seurat_subset, features=genesToPlot)

colToUse=colors_scRNAseq.l$cellTypeExt_Xen
p <- DoHeatmap(seurat_subset, assay="RNA", features = genesToPlot, group.by = "CellTypeExt", group.colors = colToUse, 
               disp.min = -1, raster = FALSE, size = 3) 
p <- p + scale_fill_gradientn(limits = c(-1,3), colors = CustomPalette(low = "lightblue", high = "orangered3", mid = "snow", k=50), 
                              values = scales::rescale(c(-1, 0, 3))  ,na.value = "white") 
p <- p + theme(axis.text.y = element_text(size = 10)) + ggtitle(seuratObj@project.name, subtitle="CellTypeExt")
print(p + guides(color = "none" ) )



### Fig S20F barplot sample+cell-type x Phase (Xen 1dpa)
seuratObj=seuratObj_Xen
seuratObj$condition_CellType=paste(seuratObj$orig.ident, seuratObj$CellType,sep="_")
ct=table(seuratObj$condition_CellType,seuratObj$Phase)
ct.df=data.frame(ct[which(rowSums(ct)>0),])
order_grp=order_grp=c("Xen_96w","Xen_ALI")

colnames(ct.df)=c("condition_CellType","Phase","Freq")
ct.df$cellType=sapply(as.character(ct.df$condition_CellType), function(s){unlist(strsplit(s,"_"))[3]})
ct.df$condition=sapply(as.character(ct.df$condition_CellType), function(s){paste(unlist(strsplit(s,"_"))[1:2],collapse="_")})
ct.df$condition=factor(ct.df$condition, levels=order_grp)
ct.df$Phase=factor(ct.df$Phase, levels=rev(c("G1", "G2M", "S")))
N=length(unique(ct.df$cellType))
N_cond=length(order_grp)
ct.df$condition_CellType=factor(ct.df$condition_CellType, levels=paste(rep(order_grp,each=N),rep(unique(ct.df$cellType),N_cond),sep="_"))
g <- ggplot(ct.df %>% filter(cellType %in% c("CT","Ectoderm","Endothelial","Immune","Muscle")), aes(fill=Phase, y=Freq, x=condition_CellType)) + 
  geom_bar(position="fill", stat="identity") + scale_fill_manual(values=colors_scRNAseq.l$Phase) 
g <- g + theme_light() + theme(axis.text.x = element_text(color = "black", size = 12, angle = 30, hjust = 1))
g_cellType_Phase <- g +  ggtitle("Phase",subtitle=seuratObj@project.name)
g_cellType_Phase <- g_cellType_Phase + scale_x_discrete(labels = as.character(sapply(order_grp, function(s){unlist(strsplit(s,"_"))[2]})), name = "")
g_cellType_Phase <- g_cellType_Phase + facet_grid(.~cellType, scales = "free_x")



#*******************************#
########## Figure S21 #########
#*******************************#
seuratObj=seuratObj_Xen
genesToPlot.l=lapply(aer_blastema.l$aer_blastema_Xen, function(v){intersect(unique(v), rownames(seuratObj@assays$RNA@data))})

### Fig S21A: dot plot AER gene expression (Xen 1dpa)
dd=DotPlot(seuratObj, assay="RNA", features=sort(genesToPlot.l$aer_genes), group.by="orig.ident", scale=FALSE) #Warning: Scaling data with a low number of groups may produce misleading results
dd.df=data.frame(dd$data)
g <- ggplot(dd.df, aes(x=features.plot,y=id, size=pct.exp, col=avg.exp)) + geom_point(aes(fill=avg.exp),color='grey75', shape=21)
g <- g + theme_light() + scale_fill_gradient(name = "avg exp",  low = "snow2", high = "#7d0025") #, limits=c(0,10)
g <- g + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
g <- g + ylab("")
g <- g + ggtitle("AER genes", subtitle="Xenopus")

### Fig S21A: dot plot blastema gene expression (Xen 1dpa)
dd=DotPlot(seuratObj, assay="RNA", features=sort(genesToPlot.l$blastema_genes), group.by="orig.ident", scale=FALSE) #Warning: Scaling data with a low number of groups may produce misleading results
dd.df=data.frame(dd$data)
g <- ggplot(dd.df, aes(x=features.plot,y=id, size=pct.exp, col=avg.exp)) + geom_point(aes(fill=avg.exp),color='grey75', shape=21)
g <- g + theme_light() + scale_fill_gradient(name = "avg exp",  low = "snow2", high = "#7d0025") #, limits=c(0,10)
g <- g + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
g <- g + ylab("")
g <- g + ggtitle("Blastema genes", subtitle="Xenopus")



#*******************************#
########## Figure 5 #########
#*******************************#

## NOTE: [Same than Figures S20C & S20D]


#*******************************#
########## Figure 6 #########
#*******************************#

## Fig 6B: dot plot HIF1A & regulators multi-species 
pathToData="./multiSpecies_limbRegeneration/data/scRNAseq_files/" #points to GitHub folder (https://github.com/BICC-UNIL-EPFL/multiSpecies_limbRegeneration/)
filename=paste0(pathToData,"data_combined_dotplot_multispecies_HIF1aGenes_allCells.csv.gz")

toPlot.df=read.csv(gzfile(filename))
g <- ggplot(toPlot.df, aes(x=gene,y=sampleTP, size=pct.exp, col=avg.exp)) + geom_point(aes(fill=avg.exp),color='grey75', shape=21)
g <- g + theme_light() + scale_fill_gradient(name = "avg exp",  low = "snow2", high = "#7d0025") #, limits=c(0,10)
g <- g + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
g <- g + ylab("")
g <- g + ggtitle("HIF1alpha_genes", subtitle="Basal Ectoderm")
g + facet_wrap(.~regenerative, scales = "free_y", ncol=1)


## Fig 6D: RidgePlots Glyco genes (Mouse Veh AER)
geneSetName="glycolytic_process"

seuratObj=seuratObj_Mouse
cells_AUC_Mouse=data_AUCellScores.l$AUCellScores_Mouse
scores_AUC_Mouse=setNames(cells_AUC_Mouse[,geneSetName],cells_AUC_Mouse[,"cellID"])
seuratObj <- AddMetaData(seuratObj, scores_AUC_Mouse, col.name="AUCell_glyco")
d_mouse=seuratObj@meta.data

seuratObj=seuratObj_Xen
cells_AUC_Xen=data_AUCellScores.l$AUCellScores_Xen
scores_AUC_Xen=setNames(cells_AUC_Xen[,geneSetName],cells_AUC_Xen[,"cellID"])
seuratObj <- AddMetaData(seuratObj, scores_AUC_Xen, col.name="AUCell_glyco")
d_xen=seuratObj@meta.data

toPlot=rbind(d_mouse[,c("AUCell_glyco","orig.ident")], d_xen[,c("AUCell_glyco","orig.ident")])
toPlot$Specie=sapply(toPlot$orig.ident, function(s){unlist(strsplit(s,"_"))[1]})
p <- ggplot(toPlot, aes(x = AUCell_glyco, y = orig.ident, fill=orig.ident)) +
  geom_density_ridges(scale = 3, quantile_lines = TRUE, quantiles = c(0.5), vline_color = c("grey33"), vline_linetype="dashed", alpha = 0.7, colour="black"
  ) + theme_light() + ylab("AUCell score")  + facet_grid(Specie~., scales = "free_y") + scale_fill_manual(values=colors_scRNAseq.l$specie_conditions_transparent) + NoLegend() + ggtitle("glycolytic_process") + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + geom_boxplot(width=0.2, outlier.shape = NA) 


## Fig 6E: dot plot FC glyco Mouse Xen
seuratObj=seuratObj_Mouse
dd=DotPlot(seuratObj, assay="RNA", features=genesForAUCell_Mouse.l$glycolytic_process, group.by="orig.ident", scale=FALSE) #Warning: Scaling data with a low number of groups may produce misleading results
dd.df=data.frame(dd$data)
dd.df$geneMouse=dd.df$features.plot

seuratObj=seuratObj_Xen
d=DotPlot(seuratObj, assay="RNA", features=genesForAUCell_Xen.l$glycolytic_process, group.by="orig.ident", scale=FALSE)
d.df=d$data
d.df$geneMouse=sapply(as.character(d.df$features.plot), function(v){vv=unlist(strsplit(v,"\\."))[1];Hmisc::capitalize(tolower(vv))})
d.df=d.df %>% filter(features.plot %in% curMarkers_Xen_glyco.df$gene)

toPlot=rbind(dd.df, d.df)
geneOrder <- c("Aldoa","Gapdh","Eno1","Tpi1","Pgk1","Slc2a1","Pkm","Slc2a3","Pgam1","Eno2","Hk2",
               "Pfkl","Hk1","Hkdc1","Pfkm","Pgam2","Bpgm","Gpi","Pklr","Pgam4","Slc2a2")
toPlot$geneMouse=factor(toPlot$geneMouse, levels=geneOrder) 
toPlot$id=factor(toPlot$id, levels=rev(c("Mouse_ALI","Mouse_96w","Xen_ALI","Xen_96w")))
g <- ggplot(toPlot, aes(x=geneMouse,y=id, size=pct.exp, col=avg.exp)) + geom_point(aes(fill=avg.exp),color='grey75', shape=21)
g <- g + theme_light() + scale_fill_gradient(name = "avg exp",  low = "snow2", high = "#7d0025") #, limits=c(0,10)
g <- g + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
g <- g + ylab("")
g <- g + ggtitle("Glycolysis genes", subtitle="")



#*******************************#
########## Figure S22 #########
#*******************************#

pathToData="./multiSpecies_limbRegeneration/data/scRNAseq_files/" #points to GitHub folder (https://github.com/BICC-UNIL-EPFL/multiSpecies_limbRegeneration/)
filename=paste0(pathToData,"data_combined_dotplot_multispecies_HIF1aGenes_allCells.csv.gz")

## FigS22A: heatmap FC glyco Mouse Xen (scRNA-seq Mouse and Xen)
toPlot.df = read.csv(gzfile(filename))
g <- ggplot(toPlot.df, aes(x=gene,y=sampleTP, size=pct.exp, col=avg.exp)) + geom_point(aes(fill=avg.exp),color='grey75', shape=21)
g <- g + theme_light() + scale_fill_gradient(name = "avg exp",  low = "snow2", high = "#7d0025") #, limits=c(0,10)
g <- g + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
g <- g + ylab("")
g <- g + ggtitle("HIF1alpha_genes", subtitle="all cells")


## FigS22B: violin plots HIF1a Mouse Xen (scRNA-seq Mouse and Xen)
## Mouse
DefaultAssay(seuratObj_Mouse)="RNA"
p = VlnPlot(seuratObj_Mouse, features = "Hif1a", group.by = "CellTypeExt", split.by = "orig.ident", pt.size=0) + ggtitle(curGene) + scale_fill_manual(values=colors_scRNAseq.l$specie_conditions) + NoLegend()
p <- p + theme(axis.title.x = element_text(size = 12),axis.text.x = element_text(size = 12),
               axis.title.y = element_text(size = 12),axis.text.y = element_text(size = 12),
               panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black")) 

## Xenopus
DefaultAssay(seuratObj_Xen)="RNA"
p = VlnPlot(seuratObj_Xen, features = "hif1a.L", group.by = "CellTypeExt", split.by = "orig.ident", pt.size=0) + ggtitle(curGene) + scale_fill_manual(values=colors_scRNAseq.l$specie_conditions) + NoLegend()
p <- p + theme(axis.title.x = element_text(size = 12),axis.text.x = element_text(size = 12),
               axis.title.y = element_text(size = 12),axis.text.y = element_text(size = 12),
               panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black")) 


## FigS22C: heatmap FC glyco Mouse Xen (scRNA-seq Mouse and Xen)
pathToData="./multiSpecies_limbRegeneration/data/scRNAseq_files/" #points to GitHub folder (https://github.com/BICC-UNIL-EPFL/multiSpecies_limbRegeneration/)

curMarkers_Mouse_glyco.df=read.csv(gzip(paste0(pathToData,"markers_Mouse_glyco_96w_vs.ALI.csv.gz")))
curMarkers_Xen_glyco.df=rread.csv(gzip(paste0(pathToData,"markers_Xen_glyco_96w_vs.ALI.csv.gz")))

toPlot=merge(curMarkers_Mouse_glyco.df, curMarkers_Xen_glyco.df, by.x="gene",by.y="geneMouse",all.x=TRUE, all.y=TRUE, suffixes = c("_Mouse","_Xen"))
rownames(toPlot)=toPlot$gene
toPlot=toPlot[order(desc(toPlot$avg_log2FC_Mouse)),]
toPlot_labels=toPlot[rownames(toPlot),c("p_val_adj_Mouse","p_val_adj_Xen")]
toPlot_labels=apply(mat2,2,function(v){v[is.na(v)]=1; return(v)})

ComplexHeatmap::Heatmap(scale(toPlot[,c("avg_log2FC_Mouse","avg_log2FC_Xen")],center=FALSE, scale=FALSE), 
                             name="log2FC 96w vs ALI", column_title = "Glycolitic genes fold change 96-well vs ALI",
                             col = colorRamp2::colorRamp2(c(-1, 0, 2), c("green", "white", "red")),
                             na_col = "grey85",
                             cluster_rows = FALSE, cluster_columns = FALSE,
                             row_names_gp = gpar(fontsize = 12), 
                             column_names_rot = 45, column_labels = setNames(c("Mouse","Xenopus"),c("avg_log2FC_Mouse","avg_log2FC_Xen")),
                             cell_fun = function(j, i, x, y, w, h, fill) {
                               if(toPlot_labels[i, j] < 0.001) {
                                 grid.text("**", x, y)
                               } else if(toPlot_labels[i, j] < 0.05) {
                                 grid.text("*", x, y)
                               }}
                        )


##Fig S22D: violin plots AUCells for glycolytic_process & oxidative_phosphorylation (Mouse)
seuratObj=seuratObj_Mouse
cells_AUC=data_AUCellScores.l$AUCellScores_Mouse
plots.l=list()
for(geneSetName in c("glycolytic_process","oxidative_phosphorylation"))
{
  print(geneSetName)
  scores_AUC=setNames(cells_AUC[,geneSetName],cells_AUC[,"cellID"])
  seuratObj <- AddMetaData(seuratObj, scores_AUC, col.name="cellColor")
  p = VlnPlot(seuratObj, features = "cellColor", group.by = "CellTypeExt", split.by = "orig.ident", pt.size=0) + ggtitle(geneSetName) + scale_fill_manual(values=colors_scRNAseq.l$specie_conditions) + NoLegend()
  plots.l[[geneSetName]] = p + theme(axis.title.x = element_text(size = 12),axis.text.x = element_text(size = 12),
                    axis.title.y = element_text(size = 12),axis.text.y = element_text(size = 12),
                    panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black")) + ylab("AUCell scores")
}
cowplot::plot_grid(plotlist = plots.l, ncol=2)


##Fig S22E: violin plots AUCells for glycolytic_process & oxidative_phosphorylation (Xen)
seuratObj=seuratObj_Xen
cells_AUC=data_AUCellScores.l$AUCellScores_Xen
plots.l=list()
for(geneSetName in c("glycolytic_process","oxidative_phosphorylation"))
{
  print(geneSetName)
  scores_AUC=setNames(cells_AUC[,geneSetName],cells_AUC[,"cellID"])
  seuratObj <- AddMetaData(seuratObj, scores_AUC, col.name="cellColor")
  p = VlnPlot(seuratObj, features = "cellColor", group.by = "CellTypeExt", split.by = "orig.ident", pt.size=0) + ggtitle(geneSetName) + scale_fill_manual(values=colors_scRNAseq.l$specie_conditions) + NoLegend()
  plots.l[[geneSetName]] = p + theme(axis.title.x = element_text(size = 12),axis.text.x = element_text(size = 12),
                                     axis.title.y = element_text(size = 12),axis.text.y = element_text(size = 12),
                                     panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black")) + ylab("AUCell scores")
}
cowplot::plot_grid(plotlist = plots.l, ncol=2)


## Fig S21E: violin plots AUCell scores oxidative_phosphorylation (Mouse & Xen)
geneSetName="oxidative_phosphorylation"
print(geneSetName)

# Mouse
seuratObj=seuratObj_Mouse
cells_AUC=data_AUCellScores.l$AUCellScores_Mouse
scores_AUC=setNames(cells_AUC[,geneSetName],cells_AUC[,"cellID"])
seuratObj <- AddMetaData(seuratObj, scores_AUC, col.name="cellColor")

g_Mouse <- ggplot(seuratObj@meta.data , aes(x=orig.ident,y=cellColor)) +
  geom_violin(aes(fill=orig.ident), scale = "width") + geom_boxplot(width=0.1, color="grey33", alpha=0.2) +
  scale_fill_manual(values=colors_scRNAseq.l$specie_conditions) + 
  theme_light() + xlab("") + ggtitle("", subtitle="") + ylab("AUCell scores") +
  theme(axis.text.x = element_text(angle=45,hjust=1,vjust=1), plot.title = element_text(hjust = 0.5, face="bold")) + NoLegend()
g_Mouse <- g_Mouse + ggtitle(geneSetName)

# Xenopus
seuratObj=seuratObj_Xen
cells_AUC=data_AUCellScores.l$AUCellScores_Xen
scores_AUC=setNames(cells_AUC[,geneSetName],cells_AUC[,"cellID"])
seuratObj <- AddMetaData(seuratObj, scores_AUC, col.name="cellColor")
g_Xen <- ggplot(seuratObj@meta.data , aes(x=orig.ident,y=cellColor)) +
    geom_violin(aes(fill=orig.ident), scale = "width") + geom_boxplot(width=0.1, color="grey33", alpha=0.2) +
    scale_fill_manual(values=colors_scRNAseq.l$specie_conditions) + 
    theme_light() + xlab("") + ggtitle("", subtitle="") + ylab("AUCell scores") +
    theme(axis.text.x = element_text(angle=45,hjust=1,vjust=1), plot.title = element_text(hjust = 0.5, face="bold")) + NoLegend()
g_Xen <- g_Xen + ggtitle(geneSetName)

g <- cowplot::plot_grid(g_Mouse + xlab("Mouse") + ggtitle(""), g_Xen + xlab("Xenopus") + ggtitle(""), ncol=2)
title <- cowplot::ggdraw() + cowplot::draw_label("Oxidative phosphorylation gene set", fontface='bold')
cowplot::plot_grid(title, g, ncol=1, rel_heights=c(0.1, 1)) # rel_heights values control title margins





