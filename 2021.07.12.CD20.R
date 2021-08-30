
options(max.print=200)
#### ------- startup ------- ####
if(F){
	install.packages("devtools")
	install.packages("SeuratObject")
	install.packages( c("gt","cowplot","gtools","patchwork","pheatmap","openxlsx","ggpubr","tidyr","dplyr","ggrepel","ggplot2","Seurat","stringr") )
	if (!requireNamespace("BiocManager",quietly=TRUE)) install.packages("BiocManager")
	BiocManager::install()
	BiocManager::install('limma') # For a more efficient implementation of the Wilcoxon Rank Sum Test,(default method for FindMarkers) please install the limma package 
	BiocManager::install(c("codetools", "KernSmooth", "nlme"))
	devtools::install_github("immunogenomics/harmony")
	devtools::install_github("atakanekiz/CIPR-Package", build_vignettes = TRUE)
}
suppressPackageStartupMessages({
	library(limma) # For a more efficient implementation of the Wilcoxon Rank Sum Test for seurat
	library(Seurat) 
	library(ggplot2) #ggplot
	library(ggrepel) #geom_label_repel() not working with Dimplot?
	library(dplyr) #%>%
	library(tidyr) #separate
	library(ggpubr) #grids
	library(openxlsx) # read.xlsx
	library(patchwork) #wrap_plots
	library(gtools) #mixedsort
	library(cowplot) # help for ggplot2; themes, functions to align plots and arrange them
	library(harmony)
	require(CIPR)
	library(Nebulosa) #plot_density (bioconductor) #https://bioconductor.org/packages/devel/bioc/vignettes/Nebulosa/inst/doc/nebulosa_seurat.html
	library(nichenetr) #devtools::install_github("saeyslab/nichenetr")
})
source("https://raw.githubusercontent.com/GildasLepennetier/R-stuff/master/not_in.R")
source("https://raw.githubusercontent.com/GildasLepennetier/R-stuff/master/NA_padding.R")
#WORKINGDIR="/dss/dsshome1/lxc01/ga94rac2/RSTUDIO/210408.CITEseq.CD20.FINAL/"
WORKINGDIR="/home/ga94rac/WORK/RSTUDIO/210408.CITEseq.CD20.FINAL/"
dir.create(WORKINGDIR,s=F) ; setwd(WORKINGDIR)
INPUTPATH="L:/210228_AGLH_CITEseq/2021.02.28.seq/OUT_CellRangMult/"

#library(future)
#plan("multiprocess", workers = 4)
#plan()
#library(parallel)
#options(Ncpus = 4)

#### create GEX from raw data ####
if(F){
	# RAW DATA = 304.7 GO in /media/ga94rac/toshibackup_3/DATA_SAVE/210228.CITE-seq_AGLH
	SAMPLES=c("H109","H124","H132","I24","I68","I105")
	list_path<-paste0(INPUTPATH,SAMPLES,"/outs/count/filtered_feature_bc_matrix")
	names(list_path)<-SAMPLES
	# Setup Seurat object
	GEX_raw = Read10X(data.dir=list_path,strip.suffix=T)
	rm(list_path)
	#The names are really not amazing, updateing that
	print ( rownames(GEX_raw$`Antibody Capture`) )
	# [1] "TotalSeq-C 0001 amouse CD4 (RM4-5)"       "TotalSeq-C 0002 amouse CD8a (53-6.7)"    
	# [3] "TotalSeq-C 0093 amouse CD19 (6D5)"        "TotalSeq-C 0015 amouse Ly6G (1A8)"       
	# [5] "TotalSeq-C 0118 amouse NK-1.1 (PK136)"    "TotalSeq-C 0301 amouse #1 (M1/42;30-F11)"
	# [7] "TotalSeq-C 0302 amouse #2 (M1/42;30-F11)" "TotalSeq-C0096 anti-mouse CD45 Antibody"
	rownames(GEX_raw$`Antibody Capture`)<-c("CD4","CD8a","CD19","Ly6G","NK","tag1","tag2","CD45")
	GEX <- CreateSeuratObject(counts=GEX_raw$`Gene Expression`,assay="RNA",project="AGLH_CITEseq",min.cells=0,min.features=0) #feature=gene in row
	GEX_Ab <- CreateAssayObject(counts = GEX_raw$`Antibody Capture`)
	GEX[["Ab"]] <- GEX_Ab
	saveRDS(GEX,"01.GEX.Rds")
	# An object of class Seurat 
	# 32293 features across 29992 samples within 2 assays 
	# Active assay: RNA (32285 features, 0 variable features)
	# 1 other assay present: Ab
	rm(GEX_raw,GEX_Ab)
	table(GEX@meta.data$orig.ident)
	# H109  H124  H132  I105   I24   I68 
	# 4072  5163  2290  4482  3923 10062
}
#### Create VDJ from raw data ####
if(F){
	SAMPLES=c("H109","H124","H132","I24","I68","I105") 
	list_path<-paste0(INPUTPATH,SAMPLES,"/outs/vdj_b/filtered_contig_annotations.csv")
	names(list_path)<-SAMPLES
	VDJ<-lapply(list_path, function(x){DF=read.csv(x) ; DF$origin = basename(dirname(dirname(dirname(x)))) ; return(DF)})
	VDJ<-do.call("rbind", VDJ)
	saveRDS(VDJ,"01.VDJ.Rds")
	rm(list_path)
	table(VDJ$origin)
	# H109 H124 H132 I105  I24  I68 
	# 188   12  649 2561  926 2190
	table(VDJ$origin,VDJ$chain)
	# IGH  IGK  IGL
	# H109   64  116    8
	# H124    7    5    0
	# H132   80  524   45
	# I105 1021 1405  135
	# I24    27  867   32
	# I68    70 1981  139
	
	# PREPARE MERGE data VDJ + GEX
	# carefull: one should not add all chains at the same time: otherwise cells will be duplicated
	
	# select useful columns: only few are useless, so I can leave them
	VDJ$barcode <- gsub(pattern="-1$", replacement="", x=VDJ$barcode) #remove trailing -1
	VDJ$sample_ID <- gsub(pattern="[.][[:digit:]]*", replacement="", x=rownames(VDJ) )
	VDJ$sample_barcode <- paste0( VDJ$sample_ID, "_", VDJ$barcode)#"I18_AAACCTGAGAAGAAGC"
	VDJ$barcode<-NULL
	VDJ$is_cell<-NULL#table(VDJ$is_cell)#all TRUE
	VDJ$high_confidence<-NULL#table(VDJ$high_confidence)#all TRUE
	VDJ <- VDJ %>% filter(chain != "Multi")
	VDJ <- VDJ %>% group_by(sample_barcode) %>% mutate( chains_list=paste0(sort(unique(chain)),collapse=","), chains_count=length(sort(unique(chain))))
	VDJ <- VDJ %>% filter(chains_count < 3)
	table(VDJ$chain)
	# IGH  IGK  IGL 
	# 1229 4858  319
	VDJ_split <- split(VDJ, VDJ$chain)
	table(VDJ$raw_consensus_id == "None") # 0 are None, 6406 have clonotype and consensus
	VDJ_split <- lapply(VDJ_split, as.data.frame)
	saveRDS(VDJ_split,"01.VDJ.filter.Rds")
}
#### add annotation ####
if(F){
	#only from one of the chain. Changing the column name to avoid overwriting
	# Problem: duplicated cell barcodes for the same chain (not biologicaly possible to have two heavy or two ligh)
	table( gsub( pattern="_.*", replacement="", x=VDJ_split$IGH$sample_barcode ) )
	# H109 H124 H132 I105  I24  I68 
	# 63    7   80  982   27   70
	VDJ_split <- lapply(VDJ_split, function(x){x %>% filter(v_gene != "None")})
	#VDJ_split$IGH [ VDJ_split$IGH$sample_barcode == "I105_CCTCAGTAGAGTCTGG" , c("reads","chain","length","v_gene","d_gene","j_gene","c_gene","cdr3","raw_clonotype_id","chains_list","sample_barcode")]
	#There are still some line with actually different IGH matches but the same cell id
	# The CDR3 is very different, so those are probably doublets
	# chain length   v_gene d_gene j_gene c_gene          cdr3 raw_clonotype_id chains_list        sample_barcode
	# 578   IGH    640  IGHV8-8         IGHJ2   IGHM CARIPNYYYFDYW     clonotype232     IGH,IGK I105_CCTCAGTAGAGTCTGG
	# 579   IGH    694 IGHV10-3         IGHJ2   IGHM  CVREGWFYFDYW     clonotype232     IGH,IGK I105_CCTCAGTAGAGTCTGG
	#sum(duplicated(VDJ_split$IGH$sample_barcode)) #3
	# Keeping only the first one (semi random choice)
	VDJ_split <- lapply(VDJ_split, function(x){ x=x [ ! duplicated( x$sample_barcode ) , ] ; return(x) })
	VDJ_split <- lapply(VDJ_split, function(x){ rownames(x)<-x$sample_barcode ; return(x) })
	#sum(duplicated($IGH$sample_barcode)) #0

	### IGH ------------
	SELECTION=c("chain","v_gene","d_gene","j_gene","c_gene","productive","cdr3","cdr3_nt","raw_clonotype_id")
	DF=VDJ_split$IGH [ ,SELECTION ]
	names(DF) <- paste0("IGH", "_", names(DF) )
	GEX <- AddMetaData(GEX, metadata=DF )
	### IG_light ------------
	DF=rbind( VDJ_split$IGK [ ,SELECTION ], VDJ_split$IGL [ ,SELECTION ] )
	names(DF) <- paste0("IG_ligh", "_", names(DF) )
	GEX <- AddMetaData(GEX, metadata=DF )
	# We do not have TRA and TRB, but here is the code
	### TRA ------------
	# DF=VDJ_split$TRA [ ,SELECTION ]
	# names(DF) <- paste0("TRA", "_", names(DF) )
	# GEX <- AddMetaData(GEX, metadata=DF )
	# ### TRB ------------
	# DF=VDJ_split$TRB [ ,SELECTION ]
	# names(DF) <- paste0("TRB", "_", names(DF) )
	# GEX <- AddMetaData(GEX, metadata=DF )
	saveRDS(GEX,"02.GEX+VDJ.Rds")
	rm(DF, VDJ, VDJ_split, SELECTION)
}
# GEX <- readRDS("01.GEX.Rds")
# VDJ <- readRDS("01.VDJ.Rds")

#### ANALYSIS part 1 (first filter and clustering) ####
#GEX<-readRDS("02.GEX+VDJ.Rds")
if(F){
	
	table(GEX$orig.ident)
	# H109  H124  H132  I105   I24   I68 
	# 4072  5163  2290  4482  3923 10062
	
	# annotate
	GEX$treatment <- "None"
	GEX$treatment [ GEX$orig.ident %in% c("I24","H109","H124") ] <- "CD20"
	GEX$treatment [ GEX$orig.ident %in% c("I105","I68","H132") ] <- "Isotype"
	
	# Add mito, ribo etc.
	if(T){
		#add mitochondrial gene percent of reads 
		GEX[["percent.mito"]] <- PercentageFeatureSet(GEX, pattern="^[Mm][Tt]-")
		# add hemoglobine gene percent of reads 
		GEX[["percent.hemo"]] <- PercentageFeatureSet(GEX, pattern="^Hb[ba]")
		# add ribosimal gene percent of reads 
		GEX[["percent.ribo"]] <- PercentageFeatureSet(GEX, pattern="^Rp[sl][[:digit:]]" )
	}
	
	nFeature_RNA_min=200
	nFeature_RNA_max=2000
	nCount_RNA_min=100
	nCount_RNA_max=10000
	percent.mito_max=5
	percent.hemo_max=5 #we do not expect cell that express hemoglobin
	
	# PLOTS - quality
	SUBDIR="02_quality_plots";dir.create(SUBDIR)
	if(T){
		PLOT=wrap_plots(
			VlnPlot(GEX, features="nFeature_RNA",pt.size=0) + 
				theme(legend.position="none") + 
				geom_hline(yintercept=nFeature_RNA_min,linetype="dashed",color="red") +
				geom_hline(yintercept=nFeature_RNA_max,linetype="dashed",color="red"),
			VlnPlot(GEX, features="nCount_RNA",pt.size=0) + theme(legend.position="none") + 
				geom_hline(yintercept=nCount_RNA_min,linetype="dashed",color="red")+
				geom_hline(yintercept=nCount_RNA_max,linetype="dashed",color="red"),
			VlnPlot(GEX, features="percent.mito",pt.size=0) + theme(legend.position="none") + 
				geom_hline(yintercept=percent.mito_max,linetype="dashed",color="red"),
			ncol=3)
		OUTNAME=paste0(SUBDIR,"/ViolinPlot.raw1.pdf");pdf(OUTNAME,width=12,height=7,useDingbats=F);print(PLOT);dev.off()
		
		
		PLOT=wrap_plots(
			VlnPlot(GEX, features="percent.hemo",pt.size=0) + theme(legend.position="none") +
				geom_hline(yintercept=percent.hemo_max,linetype="dashed",color="red"),
			VlnPlot(GEX, features="percent.ribo",pt.size=0),
			ncol=2)
		OUTNAME=paste0(SUBDIR,"/ViolinPlot.raw2.pdf");pdf(OUTNAME,width=9,height=7,useDingbats=F);print(PLOT);dev.off()
		
		
		#Pearson correlation displayed
		PLOT <- FeatureScatter(GEX, feature1="nCount_RNA", feature2="nFeature_RNA")
		OUTNAME=paste0(SUBDIR,"/FeatureScatter.readCount_vs_geneCount.pdf");pdf(OUTNAME,width=7,height=7,useDingbats=F);print(PLOT);dev.off()
		
		PLOT <- FeatureScatter(GEX, feature1="nCount_RNA", feature2="percent.mito")
		OUTNAME=paste0(SUBDIR,"/FeatureScatter.readCount_vs_mito.pdf");pdf(OUTNAME,width=7,height=7,useDingbats=F);print(PLOT);dev.off()
		
		PLOT <- FeatureScatter(GEX, feature1="nCount_RNA", feature2="percent.hemo")
		OUTNAME=paste0(SUBDIR,"/FeatureScatter.readCount_vs_hemo.pdf");pdf(OUTNAME,width=7,height=7,useDingbats=F);print(PLOT);dev.off()
		
		PLOT <- FeatureScatter(GEX, feature1="nCount_RNA", feature2="percent.ribo")
		OUTNAME=paste0(SUBDIR,"/FeatureScatter.readCount_vs_ribo.pdf");pdf(OUTNAME,width=7,height=7,useDingbats=F);print(PLOT);dev.off()
		
	}
	
	##### FILTER ##### 
	
	# remove sample I24: Its was of very poor quality of the CD20 B-cell depletion
	GEX<- subset( GEX, orig.ident != "I24" ) # minus 3923 cells (898 B-cells when compared to 128 in H109 and 8 in H124)
	# re-do the factors
	GEX$orig.ident = factor(as.character(GEX$orig.ident))
	
	GEX
	# 32293 features across 26069 samples within 2 assay
	if(T){
		GEX <- subset(GEX, subset = 
					  	nFeature_RNA > nFeature_RNA_min & 
					  	nFeature_RNA < nFeature_RNA_max & 
					  	nCount_RNA > nCount_RNA_min & 
					  	nCount_RNA < nCount_RNA_max & 
					  	percent.mito < percent.mito_max & 
					  	percent.hemo < percent.hemo_max )
	}
	GEX #32293 features across 17226 samples within 2 assays
	# 26069-17226  = 8843 cells-barcodes removed
	# 8843/26069 * 100    #  33.9 % cell removed by filtering
	
	# PLOTS - quality after filter
	if(T){
		PLOT=wrap_plots(
			VlnPlot(GEX, features="nFeature_RNA",pt.size=0) + 
				theme(legend.position="none") + 
				geom_hline(yintercept=nFeature_RNA_min,linetype="dashed",color="red") +
				geom_hline(yintercept=nFeature_RNA_max,linetype="dashed",color="red"),
			VlnPlot(GEX, features="nCount_RNA",pt.size=0) + theme(legend.position="none") + 
				geom_hline(yintercept=nCount_RNA_min,linetype="dashed",color="red")+
				geom_hline(yintercept=nCount_RNA_max,linetype="dashed",color="red"),
			VlnPlot(GEX, features="percent.mito",pt.size=0) + theme(legend.position="none") + 
				geom_hline(yintercept=percent.mito_max,linetype="dashed",color="red"),
			ncol=3)
		OUTNAME=paste0(SUBDIR,"/ViolinPlot.filter1.pdf");pdf(OUTNAME,width=12,height=7,useDingbats=F);print(PLOT);dev.off()
		
		
		PLOT=wrap_plots(
			VlnPlot(GEX, features="percent.hemo",pt.size=0) + theme(legend.position="none") +
				geom_hline(yintercept=percent.hemo_max,linetype="dashed",color="red"),
			VlnPlot(GEX, features="percent.ribo",pt.size=0),
			ncol=2)
		OUTNAME=paste0(SUBDIR,"/ViolinPlot.filter2.pdf");pdf(OUTNAME,width=9,height=7,useDingbats=F);print(PLOT);dev.off()
	}
	
	# PROCESS - GEX
	if(T){
		gc()
		GEX <- NormalizeData(GEX) #, normalization.method="LogNormalize", scale.factor=10000
		GEX <- FindVariableFeatures(GEX) #, selection.method="vst", nfeatures=2000
		GEX <- ScaleData(GEX) # not recommanded any more: vars.to.regress="percent.mt"
		GEX <- RunPCA(GEX,assay="RNA") #, features=VariableFeatures(object=GEX)
		Idents(GEX)<-"orig.ident"
		gc()
	}
	
	#make chain
	if(T){
		GEX@meta.data$chain = paste(GEX@meta.data$IGH_chain,GEX@meta.data$IG_ligh_chain,sep=",")
		GEX@meta.data$chain <- gsub("NA,?","",GEX@meta.data$chain)
		GEX@meta.data$chain <- gsub(",+$","",GEX@meta.data$chain)
		GEX@meta.data$chain <- factor(GEX@meta.data$chain)
		levels(GEX@meta.data$chain)[ levels(GEX@meta.data$chain) == ""] <- "none"
	}
	
	#check if ribosomal genes are on the variable genes: only one on the top 2000 var features #"Rps27l" "Rpl39l"
	if(T){
		SUBDIR="02_quality_plots/"
		topVar <- head(VariableFeatures(GEX), 2000)
		grep("^Rp[sl][[:digit:]]",topVar, value=T) #"Rps27l"
		topVar <- head(VariableFeatures(GEX), 20)
		PLOT=LabelPoints(plot=VariableFeaturePlot(GEX), points=topVar, repel=TRUE, xnudge=0,ynudge=0)
		OUTNAME=paste0(SUBDIR,"/VariableFeaturePlot.Top2000MostVarGenes.pdf");pdf(OUTNAME,width=16,height=8,useDingbats=F);print(PLOT);dev.off()
		
		rm(topVar,PLOT,OUTNAME,SUBDIR)
	}
	#### PLOT - PCA and heatmap
	SUBDIR="03.pca-heatmaps";dir.create(SUBDIR)
	if(T){
		#VizDimLoadings(GEX, dims=1:2, reduction="pca")
		Idents(GEX)<-"orig.ident"
		PLOT=DimPlot(GEX, reduction="pca", dims=c(1,2))
		OUTNAME=paste0(SUBDIR,"/PCA.pdf");pdf(OUTNAME,width=7,height=7,useDingbats=F);print(PLOT);dev.off()
		rm(PLOT,OUTNAME)
		
		OUTNAME=paste0(SUBDIR,"/DimHeatmap.pdf");pdf(OUTNAME,width=7,height=7,useDingbats=F);
		DimHeatmap(GEX, dims=1, cells=500, balanced=TRUE) ; dev.off()
		rm(OUTNAME)
	}
	
	
	#### Cluster
	NUMBER_PC=25
	if(T){
		GEX <- RunUMAP(GEX, dims=1:NUMBER_PC) #default seed is 42
		GEX <- RunTSNE(GEX, dims=1:NUMBER_PC) #default seed is 42
		GEX <- FindNeighbors(GEX, dims=1:NUMBER_PC) # knn graph
		#We find that setting this parameter between 0.4-1.2 typically returns good results for single-cell datasets of around 3K cells.
		for( RESO in c(0.1,0.5,1) ){GEX <- FindClusters(GEX, resolution=RESO, random.seed=2020)} ; rm(RESO)
		for(index in grep("RNA_snn_res",x=names(GEX@meta.data),value=T)){print(paste(index,length(table(GEX@meta.data[,index]))))};rm(index)
		# [1] "RNA_snn_res.0.1 9"
		# [1] "RNA_snn_res.0.5 19"
		# [1] "RNA_snn_res.1 25"
	}
	
	SUBDIR="04.raw_cluster_plot";dir.create(SUBDIR)
	# PLOT chains
	if(T){
		table(GEX@meta.data$chain)
		# none     IGH IGH,IGK IGH,IGL     IGK     IGL 
		# 14030      22     950      64    2011     149
		
		NAMEIDENT<-"chain";Idents(GEX)<-NAMEIDENT
		PLOT=DimPlot( subset(GEX,chain!="none"), reduction="umap", label=TRUE, pt.size=.1, label.size=8, repel=T) + 
			theme_minimal(base_size=16)+
			theme(legend.position='top')+
			coord_fixed()+
			ggtitle(paste(NAMEIDENT))
		OUTNAME=paste0(SUBDIR,"/umap.",NAMEIDENT,".pdf");pdf(OUTNAME,width=8,height=8,useDingbats=F);print(PLOT);dev.off()
	}
	# PLOT clusters
	if(T){
		NAMES_CLU=c( "orig.ident", grep(pattern = "RNA_snn_res",x = names(GEX@meta.data),value = T) )
		for ( NAMEIDENT in NAMES_CLU ){
			Idents(GEX)<-NAMEIDENT
			PLOT=DimPlot(GEX, reduction="umap", label=TRUE, pt.size=.1, label.size=8, repel=T) + 
				theme_minimal(base_size=16)+
				theme(legend.position='top')+
				coord_fixed()+
				ggtitle(paste(NAMEIDENT))
			OUTNAME=paste0(SUBDIR,"/umap.",NAMEIDENT,".pdf");pdf(OUTNAME,width=8,height=8,useDingbats=F);print(PLOT);dev.off()
			rm(PLOT,OUTNAME,NAMEIDENT)
		}
		#	names(GEX@meta.data)
	}
	
	# PREPARE CHAINS INFO
	if(T){
		# Here: removed part concerning the TCR, since no TCR in this experiment
		GEX$immune_chains=paste(GEX$IGH_chain,GEX$IG_ligh_chain,sep=",")
		df<-GEX@meta.data %>% select(IGH_chain,IG_ligh_chain) %>%tidyr::unite(col="immune_chains",IGH_chain,IG_ligh_chain, na.rm=TRUE, sep=",")
		df$immune_chains[ df$immune_chains == ""]<-"no info"
		GEX<-AddMetaData(GEX,df)
		
		KEEP=c("TRA","TRB","TRA,TRB","IGH","IGH,IGK","IGH,IGL","IGK","IGL","no info")
		GEX$immune_chains_clean<-GEX$immune_chains
		GEX$immune_chains_clean[ ! GEX$immune_chains_clean %in% KEEP ] <- "multi"
		
		GEX$immune_chains_clean2=GEX$immune_chains_clean
		GEX$immune_chains_clean2[GEX$immune_chains_clean2 %in% c("TRA","TRB","TRA,TRB") ] <- "T-cell"
		GEX$immune_chains_clean2[GEX$immune_chains_clean2 %in% c("IGH","IGH,IGK","IGH,IGL","IGK","IGL") ] <- "B-cell"
		GEX$immune_chains_clean2[ ! GEX$immune_chains_clean2 %in% c("T-cell","B-cell","multi") ] <- "other"
		rm(KEEP,df)
		
		TABLE=as.data.frame(table(GEX$immune_chains_clean2, GEX$RNA_snn_res.0.1))
		# TABLE = subset(TABLE, Var1 %in% c("B-cell","T-cell"))
		TABLE=TABLE %>% group_by(Var1) %>% mutate(percent = Freq / sum(Freq) * 100 )
		PLOT=ggplot(TABLE)+geom_tile(aes(x=Var1,y=Var2,fill=percent)) + 
			geom_text(aes(x=Var1,y=Var2,label=round(percent,0))) +
			xlab("VDJ information") + ylab("Cluster from GEX") + ggtitle("Percent of VDJ cells in cluster")
		OUTNAME=paste0("04.Percent.VDJ_vs_GEX.pdf");pdf(OUTNAME,width=4,height=4,useDingbats=F);print(PLOT);dev.off()
		
	}
	
	# find doublets ///
	NUMBER_PC=25 #<- REQUIRED
	SUBDIR="05.Doublet";dir.create(SUBDIR)
	if(T){
		# ---- DoubletFinder
		USE_ANNOT_REF="orig.ident" #here expect groups with same doublet proba
		if(T){
			#remotes::install_github('chris-mcginnis-ucsf/DoubletFinder')
			require(DoubletFinder)
			#notes
			if(F){
				#https://github.com/chris-mcginnis-ucsf/DoubletFinder
				#(1) Generate artificial doublets from existing scRNA-seq data
				#(2) Pre-process merged real-artificial data
				#(3) Perform PCA and use the PC distance matrix to find each cell's proportion of artificial k nearest neighbors (pANN)
				#(4) Rank order and threshold pANN values according to the expected number of doublets
				
				# pK identification (principal componant neighborhood size)
				# Optimal pK values can be determined using mean-variance-normalized bimodality coefficient.
				
				## pK Identification (ground-truth=when we know some Doublet already) ------------------------------------------------------------------------------------------
				#sweep.res.list_kidney <- paramSweep_v3(seu_kidney, PCs=1:10, sct=F)
				#gt.calls <- seu_kidney@meta.data[rownames(sweep.res.list_kidney[[1]]), "GT"].   ## GT is a vector containing "Singlet" and "Doublet" calls recorded using sample multiplexing classification and/or in silico geneotyping results 
				#sweep.stats_kidney <- summarizeSweep(sweep.res.list_kidney, GT=T, GT.calls=gt.calls)
				#bcmvn_kidney <- find.pK(sweep.stats_kidney)
				
				## Homotypic Doublet Proportion Estimate -------------------------------------------------------------------------------------
				#homotypic.prop <- modelHomotypic(  GEX@meta.data[,USE_ANNOT_REF]  )           ## ex: annotations <- seu_kidney@meta.data$ClusteringResults
				
			}
			## pK Identification (no ground-truth) ---------------------------------------------------------------------------------------
			sweep.res.list <- paramSweep_v3(GEX, PCs=1:NUMBER_PC, sct=F) #, num.cores=3 !not supported on Windows. mclapply(as.list(1:length(pN)), FUN=parallel_paramSweep_v3
			saveRDS(sweep.res.list,paste0(SUBDIR,"/sweep.res.list.Rds"))
			#sweep.res.list<-readRDS( paste0(SUBDIR,"/sweep.res.list.",USE_ANNOT_REF,".Rds"))
			sweep.stats <- summarizeSweep(sweep.res.list, GT=F)
			bcmvn <- find.pK(sweep.stats)
			#Optimal pK for any scRNA-seq data can be manually discerned as maxima in BCmvn distributions
			pK_choice=as.numeric(paste(bcmvn[which.max(bcmvn$BCmetric),"pK"]))
			
			#This set the threshold: 10%
			nExp_poi <- round(0.1*nrow(GEX@meta.data))  ## Assuming 7.5% doublet formation rate - tailor for your dataset
			#nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
			## Run DoubletFinder with varying classification stringencies ----------------------------------------------------------------
			#using nExp_poi
			GEX <- doubletFinder_v3(GEX,PCs=1:15,pN=0.25,pK=pK_choice,nExp=nExp_poi,reuse.pANN=F,sct=F)
			#pANN_NAME=names(GEX@meta.data)[length(names(GEX@meta.data))-1]
			pANN_NAME=paste0("pANN_0.25_",pK_choice,"_",nExp_poi)
			print(paste("pk=",pK_choice,"pANN_NAME=",pANN_NAME))
			# "pk= 0.11 pANN_NAME= pANN_0.25_0.11_1723"
			#clean
			rm(nExp_poi,sweep.stats,sweep.res.list) #rm(bcmvn,homotypic.prop)
		}
		# DoubletFinder find doublets v2 - scan several % and check count of cell lost in categorie "VDJ-Multi" and "Other (no immune info)"
		#pK_choice=0.2 #0.005
		#pANN_NAME="pANN_0.25_0.005_1292" #names(GEX@meta.data)
		SCAN_VALUES=seq( 0.05, .3, 0.05 ) #0.075
		# make the expected number vary
		if(T){
			# ---- DoubletFinder: 
			homotypic.prop <- modelHomotypic(  GEX@meta.data[,USE_ANNOT_REF]  )   
			for( RATE in SCAN_VALUES ){
				#This set the threshold
				nExp_poi <- round(RATE*nrow(GEX@meta.data)) # changing expected doublet formation rate
				#nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
				#calc
				GEX <- doubletFinder_v3(GEX,PCs=1:10,pN=0.25,pK=pK_choice,nExp=nExp_poi,reuse.pANN=pANN_NAME,sct=F)
				DF.CLASS=names(GEX@meta.data)[length(names(GEX@meta.data))]
				print(paste("pk=",pK_choice,"pANN_NAME=",pANN_NAME, "DF.CLASS0",DF.CLASS, "RATE=",RATE))
			}
			#clean
			rm(nExp_poi,homotypic.prop,DF.CLASS)
		}
	}
	## PLOT Doublet
	if(T){
		NAMES_DF=mixedsort(grep("DF.classifications",x=names(GEX@meta.data), value=T))
		for(NAMEIDENT in NAMES_DF){
			Idents(GEX)<-NAMEIDENT
			TABLE=table(GEX@meta.data[,NAMEIDENT])
			PLOT=DimPlot(GEX, reduction="umap", label=F, pt.size=1, label.size=6, repel=T, order="Doublet") + 
				theme(axis.line=element_line(size=1),
					  text=element_text(size=18),
					  axis.text=element_text(size=18),
					  axis.ticks=element_line(size=1),
					  legend.position='top') + coord_fixed() +
				ggtitle( paste0("Doublet: ",TABLE["Doublet"],", Singlet: ",TABLE["Singlet"], " - ", round(TABLE["Doublet"]/sum(TABLE)*100,2),"%\n",NAMEIDENT))
			OUTNAME=paste0(SUBDIR,"/DoubletFinder.",NAMEIDENT,".pdf");pdf(OUTNAME,width=10,height=10,useDingbats=F) ; print(PLOT) ; dev.off() ; print(OUTNAME)
			rm(PLOT,OUTNAME)
		}
	}
	
	
	# HTO: Ab staining
	if(T){
		# HTO TISSUES
		GEX_AbTissues <- CreateAssayObject(counts = GEX@assays$Ab@counts [c("tag1","tag2"),] )
		GEX[["AbTissues"]] <- GEX_AbTissues ; rm(GEX_AbTissues)
		SUBDIR="05.AbTissues";dir.create(SUBDIR)
		GEX <- NormalizeData(GEX,assay="AbTissues",normalization.method="CLR") #Applies a centered log ratio transformation
		GEX<-HTODemux(GEX,assay = "AbTissues")
		Idents(GEX) <- "AbTissues_classification"  #"Ab_maxID" may be biases because of low expression
		for ( FEAT in rownames(GEX[["AbTissues"]]) ){
			PLOT=RidgePlot(GEX,assay="AbTissues",features=FEAT,ncol=1)
			OUTNAME=paste0(SUBDIR,"/RidgePlot.",FEAT,".pdf");pdf(OUTNAME,width=8,height=8,useDingbats=F);print(PLOT);dev.off()
		}
		table(GEX$AbTissues_classification)
		PLOT=HTOHeatmap(GEX,assay="AbTissues",ncells= dim(GEX)[2] )
		OUTNAME=paste0(SUBDIR,"/HTOHeatmap.pdf");pdf(OUTNAME,width=12,height=8,useDingbats=F);print(PLOT);dev.off()
		GEX$tag_tissue<-GEX$AbTissues_classification
		GEX$tag_tissue [ grep("_",GEX$tag_tissue) ] <- "Doublet"
		GEX$tag_tissue<-factor(GEX$tag_tissue)
		
		levels(GEX$tag_tissue)<-c("Doublet","Negative","CSF","mELT") #"Doublet"  "Negative" "tag1"     "tag2"
		
		TABLE=table(GEX$tag_tissue);TABLE
		# Doublet Negative     CSF     mELT 
		# 1842     9278     5309      797
		print( round(TABLE/ sum(TABLE) * 100))
		# Doublet Negative     CSF     mELT 
		# 11       54       31        5
		
		# gene expression is in log1p() ( log(1+x) ), the reverse operation is exp( X )-1
		
		summary(  as.vector(exp(GEX[["AbTissues"]]["tag1"]) -1 ) )
		# in.  1st Qu.   Median     Mean  3rd Qu.     Max. 
		# 0.000    0.352    0.822    3.419    2.348 3325.934
		summary(  as.vector(exp(GEX[["AbTissues"]]["tag2"]) -1 ) )
		# Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
		# 0.03737  0.64111  0.94003  1.24703  1.52214 35.87251
		THR=1.5
		TAGS1 = ( exp(GEX[["AbTissues"]]["tag1"]) -1 ) > THR * ( exp( GEX[["AbTissues"]]["tag2"]) -1 )
		TAGS2 = ( exp(GEX[["AbTissues"]]["tag2"]) -1 ) > THR * ( exp( GEX[["AbTissues"]]["tag1"]) -1 )
		# relative expression approach
		GEX$tag_tissue_thr1.5 = "none"
		GEX$tag_tissue_thr1.5[ TAGS1 ] <-"CSF"
		GEX$tag_tissue_thr1.5[ TAGS2 ] <-"mELT"
		
		TABLE = table( TAGS1, TAGS2 ) ; TABLE
		print(round(TABLE/ sum(TABLE) * 100))
		rm(TAGS1,TAGS2)
		#		TAGS2
		# TAGS1   FALSE TRUE
		# FALSE  3999 7532
		# TRUE   5695    0
		#			TAGS2				#percent: 77% good
		# TAGS1   FALSE TRUE
		# FALSE    23   44
		# TRUE     33    0
		TABLE=table(GEX$tag_tissue_thr1.5);TABLE
		# CSF mELT none 
		# 5695 7532 3999
		print( round(TABLE/ sum(TABLE) * 100))
		# CSF mELT none 	#percent of cells: 77 properly annotated
		# 33   44   23
		
		
		
		# HTO CELLTYPES
		GEX_AbCelltypes <- CreateAssayObject(counts = GEX@assays$Ab@counts [c("CD4","CD8a","CD19","Ly6G","NK"),] )
		GEX[["AbCelltypes"]] <- GEX_AbCelltypes ; rm(GEX_AbCelltypes)
		SUBDIR="05.AbCelltypes";dir.create(SUBDIR)
		GEX <- NormalizeData(GEX,assay="AbCelltypes",normalization.method="CLR") #Applies a centered log ratio transformation
		GEX<-HTODemux(GEX,assay = "AbCelltypes")
		Idents(GEX) <- "AbCelltypes_classification"  #"Ab_maxID"
		for ( FEAT in rownames(GEX[["AbCelltypes"]]) ){
			PLOT=RidgePlot(GEX,assay="AbCelltypes",features=FEAT,ncol=1)
			OUTNAME=paste0(SUBDIR,"/RidgePlot.",FEAT,".pdf");pdf(OUTNAME,width=8,height=8,useDingbats=F);print(PLOT);dev.off()
		}
		
		PLOT=HTOHeatmap(GEX,assay="AbCelltypes",ncells= dim(GEX)[2] )
		OUTNAME=paste0(SUBDIR,"/HTOHeatmap.pdf");pdf(OUTNAME,width=12,height=8,useDingbats=F);print(PLOT);dev.off()
		GEX$tag_celltype<-factor(GEX$AbCelltypes_classification)
		levels(GEX$tag_celltype) [ grep("_",levels(GEX$tag_celltype)) ] <- "Doublet"
		TABLE = table(GEX$tag_celltype); TABLE
		# CD19  Doublet      CD4     CD8a     Ly6G Negative       NK 
		# 594    10348     2921      317      390     2337      319
		print( round(TABLE/ sum(TABLE) * 100))
		# CD19  Doublet      CD4     CD8a     Ly6G Negative       NK 
		# 3       60       17        2        2       14        2
	}
	
	table(GEX$tag_celltype, GEX$DF.classifications_0.25_0.02_1723)
	#				Doublet Singlet	(DoubletFinder)
	# CD19          24     570
	# Doublet     1141    9207			<- probably wrong annotation of the HTODemux
	# CD4          120    2801
	# CD8a          39     278
	# Ly6G          52     338
	# Negative     315    2022
	# NK            32     287
	#table(GEX$tag_celltype, GEX$DF.classifications_0.25_0.02_5168)
	
	saveRDS(GEX,"03.GEX+VDJ.Rds")
	
	
}

####  ANALYSIS part 2 - load and filter putative doublets #### 
GEX<-readRDS("03.GEX+VDJ.Rds")
# An object of class Seurat 
# 32300 features across 17226 samples within 4 assays 
# Active assay: RNA (32285 features, 2000 variable features)
# 3 other assays present: Ab, AbTissues, AbCelltypes
# 3 dimensional reductions calculated: pca, umap, tsne
GEX <- subset(GEX, DF.classifications_0.25_0.02_1723 == "Singlet")
GEX
# An object of class Seurat 
# 32300 features across 15503 samples within 4 assays		#before 32293 genes, now +7: 2tissues and 5 celltypes in the HTO assay
# Active assay: RNA (32285 features, 2000 variable features)
# 3 other assays present: Ab, AbTissues, AbCelltypes
# 3 dimensional reductions calculated: pca, umap, tsne


#### CELL CYCLE PHASE -- REGRESS OUT EFFECTS + Harmony ####
#library(future)
#plan("multiprocess", workers = 4)
#plan()
#options(future.globals.maxSize = 1000 * 1024^2) # 1 Go instead of 500 Mo RAM limit
if(F){
	s.genes=na.omit(nichenetr::convert_human_to_mouse_symbols(cc.genes$s.genes))
	g2m.genes=na.omit(nichenetr::convert_human_to_mouse_symbols(cc.genes$g2m.genes))
	
	GEX <- CellCycleScoring(GEX,s.genes,g2m.genes)
	rm(s.genes,g2m.genes)
	GEX<-NormalizeData(GEX) # re-normalize after each 
	GEX <- FindVariableFeatures(GEX)
	system.time(GEX <- ScaleData(GEX, vars.to.regress = c("S.Score", "G2M.Score"), features = VariableFeatures(GEX) ) )
	
	#using only 2K variable feature now, really faster and take less place
	GEX <- RunPCA(GEX) #, features=VariableFeatures(object=GEX)
	
	GEX$SAMPLES=GEX$orig.ident
	GEX <- GEX %>% RunHarmony("SAMPLES", plot_convergence = TRUE)
	
	GEX <- RunUMAP(GEX, reduction = "harmony", dims = 1:30)
	GEX <- FindNeighbors(GEX, reduction = "harmony", dims = 1:30)
	NAMES=grep("RNA_snn",names(GEX@meta.data),value=T)
	names(GEX@meta.data)[ names(GEX@meta.data) %in% NAMES ] <- paste0("pre_cc_",NAMES) ; rm(NAMES)
	names(GEX@meta.data)
	
	plan("multiprocess", workers = 1)# this is to avoid a seed problem with parallel 
	for( RESO in c(0.1,0.5,1,1.5,2) ){GEX <- FindClusters(GEX, resolution = RESO, random.seed = 2020)} ; rm(RESO)
	#Louvain algorithm
	
	for(index in grep("RNA_snn_res",x=names(GEX@meta.data),value=T)){print(paste(index,length(table(GEX@meta.data[,index]))))};rm(index)
	# [1] "pre_cc_RNA_snn_res.0.1 9"
	# [1] "pre_cc_RNA_snn_res.0.5 19"
	# [1] "pre_cc_RNA_snn_res.1 25"
	
	# [1] "RNA_snn_res.0.1 11"
	# [1] "RNA_snn_res.0.5 18"
	# [1] "RNA_snn_res.1 23"
	# [1] "RNA_snn_res.1.5 30"
	# [1] "RNA_snn_res.2 33"

}

#saveRDS(GEX,"06.cellcycle.Rds")
#GEX<-readRDS("06.cellcycle.Rds")

#### FindAllMarkers for the clusters: used by CIPR + to check markers ####
SUBDIR="11.Cluster_Harmony";dir.create(SUBDIR)
library(future)
#plan("multiprocess", workers = 6)
#plan()
#options(future.globals.maxSize = 1000 * 1024^2)
if(F){
	#DefaultAssay(GEX)
	#*NAMES_CLU=grep("^RNA_snn_res",names(GEX@meta.data),value=T)
	NAMES_CLU=c( "RNA_snn_res.0.1", "RNA_snn_res.0.5", "RNA_snn_res.1","RNA_snn_res.1.5","RNA_snn_res.2") # dir(SUBDIR) 
	
	for ( NAMEIDENT in NAMES_CLU ){
		RDS_MARKERS = paste0(SUBDIR,"/allMarkers.",NAMEIDENT,".Rds")
		if( ! file.exists( RDS_MARKERS ) | TRUE){
			Idents(GEX)<-NAMEIDENT
			print(paste("DOING:",NAMEIDENT, "/",length(unique(Idents(GEX)))))
			allMarkers <- FindAllMarkers(GEX,verbose=T)
			saveRDS(allMarkers,RDS_MARKERS) # 
		}else{print(paste("already exists:",RDS_MARKERS))}
		if( file.exists( RDS_MARKERS ) | TRUE ){
			print(paste("making top20 + average expr for identity:",NAMEIDENT))
			allMarkers <- readRDS(RDS_MARKERS)
			allMarkers_top20 = allMarkers %>% group_by(cluster) %>% arrange(cluster,p_val_adj)  %>% slice_head(n = 20)
			write.xlsx(allMarkers_top20,paste0(SUBDIR,"/allMarkers+top20.",NAMEIDENT,".xlsx"))
			if(F){
				Idents(GEX)<-NAMEIDENT
				avgexp <- AverageExpression(GEX,verbose=T,assays="RNA")
				avgexp=as.data.frame(avgexp$RNA)
				avgexp$gene = rownames(avgexp)
				allMarkers_top20_Expr = merge.data.frame(x = allMarkers_top20,y = avgexp, by = "gene") %>% arrange(cluster,p_val_adj)
				write.xlsx(allMarkers_top20_Expr,paste0(SUBDIR,"/allMarkers+avgexp+top20.",NAMEIDENT,".xlsx"))
				rm(allMarkers_top20_Expr,avgexp)
			}
			rm(allMarkers,allMarkers_top20)
		}
	}
}
#### annotation celltypes - CIPR - ####
#https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-020-3538-2
SUBDIR="11.Cluster_Harmony";dir.create(SUBDIR,showWarnings=F)
#select: carefull select function also exist in AnnotationDbi::
if(F){
	#NAMES_CLU=grep("^RNA_snn_res",names(GEX@meta.data),value=T)
	NAMES_CLU=c( "RNA_snn_res.0.1", "RNA_snn_res.0.5", "RNA_snn_res.1", "RNA_snn_res.1.5", "RNA_snn_res.2") # dir(SUBDIR)
	for ( NAMEIDENT in NAMES_CLU ) { # c("RNA_snn_res.0.1", "RNA_snn_res.0.2") ) { #
		Idents(GEX) <- NAMEIDENT
		print(paste0("doing: ",NAMEIDENT, " (",length(unique(Idents(GEX)))," clusters)"))
		# # # # CIPR - Cluster Identity PRedictor https://github.com/atakanekiz/CIPR-Package
		if(T){
			#required for CIPR
			
			allMarkers<-readRDS(paste0(SUBDIR,"/allMarkers.",NAMEIDENT,".Rds"))
			#DATABASE="mmrnaseq"
			DATABASE="immgen"
			
			# Plot summarizing top scoring references per cluster (logFC comparison)
			PLOT <- CIPR(input_dat = allMarkers,comp_method="logfc_dot_product",reference=DATABASE,top_num=1,keep_top_var=50,plot_top=F)
			
			ANNOT = CIPR_top_results %>% 
				dplyr::select( cluster,reference_cell_type) %>% 
				dplyr::distinct() %>%
				dplyr::group_by(cluster) %>% dplyr::mutate( annot = paste0(reference_cell_type,collapse="+")) %>%
				dplyr::select(cluster,annot) %>% dplyr::distinct()
			#Deal with unannoted clusters
			NOT_ANNOT = levels(Idents(GEX)) [ ! levels(Idents(GEX)) %in% ANNOT$cluster ]
			if(length(NOT_ANNOT)>0){
				tmp_df=data.frame(cluster=NOT_ANNOT, annot="unknown")
				ANNOT=rbind(ANNOT,tmp_df) ; rm(tmp_df)
			}
			#USE THE ROWNAME TO AVOID ERROR IN CLUSTER NAME ~ LEVELS
			ANNOT=as.data.frame(ANNOT[order(ANNOT$cluster),])
			rownames(ANNOT)=as.character(ANNOT$cluster)
			write.xlsx(ANNOT, paste0(SUBDIR,"/CIPR.",NAMEIDENT,".xlsx"))
			#ADD in the data
			NEWNAME=paste0("annot_",NAMEIDENT,"_CIPR")
			GEX@meta.data[,NEWNAME] <- GEX@meta.data[,NAMEIDENT]
			levels(GEX@meta.data[,NEWNAME]) <- ANNOT[ levels(GEX@meta.data[,NEWNAME]),"annot"]
			rm(CIPR_top_results,CIPR_all_results, top_plots,PLOT,allMarkers,ANNOT)
			rm(DATABASE,NEWNAME,NOT_ANNOT)
			
		} 
		#clusters:
		#12
		#18
		#24
		#26
	}
}
##### MERGE EXCEL FILES: gene markers + CIPR annot ####
SUBDIR="11.Cluster_Harmony/"
if(F){
	NAMES_CLU=grep(pattern = "^RNA_snn_res", x = names(GEX@meta.data), value = T)
	for ( NAMEIDENT in NAMES_CLU ){
		allMarkersExpr <- readRDS(paste0(SUBDIR,"/allMarkers.",NAMEIDENT,".Rds"))
		Annotation <-read.xlsx(paste0(SUBDIR,"/CIPR.",NAMEIDENT,".xlsx"))
		MERGED=merge.data.frame(allMarkersExpr,Annotation,by="cluster")
		write.xlsx(MERGED,paste0(SUBDIR,"/allMarkers.",NAMEIDENT,"_CIPR.xlsx"))
		rm(allMarkersExpr,Annotation,MERGED)
	}
}



#### PLOT clusters + annot ####
SUBDIR="12.Plots";dir.create(SUBDIR)
if(F){
	for(index in 1:4){
		if(index==1){NAMES_CLU=c("RNA_snn_res.2","annot_RNA_snn_res.2_CIPR","tag_tissue_thr1.5","orig.ident") ; ORDER=NULL }
		if(index==2){NAMES_CLU=c("immune_chains_clean2"); ORDER=c("B-cell","other")  }
		if(index==3){NAMES_CLU="tag_celltype"; ORDER=c("NK","Ly6G","CD8a","CD4","CD19","Negative","unknown") }
		if(index==4){NAMES_CLU="Phase"; ORDER=c("G2M","S","G1") }
		
		for ( NAMEIDENT in NAMES_CLU){
			print(paste("doing:",NAMEIDENT))
			Idents(GEX)<-NAMEIDENT
			PLOT=DimPlot(GEX,reduction="umap",label=T,order=ORDER,pt.size=.1,label.size=4,repel=T) + 
				theme_minimal(base_size = 16)+
				theme(legend.position = 'top')+
				ggtitle(paste(NAMEIDENT)) + coord_fixed() #+ geom_hline( yintercept = -9) #+ geom_vline( xintercept = -3.5)
			OUTNAME=paste0(SUBDIR,"/umap.",NAMEIDENT,".pdf");pdf(OUTNAME,width=10,height=8,useDingbats=F);print(PLOT);dev.off()
			rm(PLOT,OUTNAME)
		}
	}
	
	#GEX$split = paste(GEX$treatment, GEX$tag_tissue)
	GEX$split = paste(GEX$treatment, GEX$tag_tissue_thr1.5)
	NAMES_CLU=c("annot_RNA_snn_res.2_CIPR") ; ORDER=NULL # dir(SUBDIR)
	for(index in 1:3){
		if(index==1){SPLIT="tag_tissue_thr1.5" ; WIDTH=20;HEIGHT=8; NCOL=4}
		if(index==2){SPLIT="orig.ident"; WIDTH=20;HEIGHT=15; NCOL=3}
		if(index==3){SPLIT="split" ; WIDTH=20;HEIGHT=12; NCOL=3}
		
		for ( NAMEIDENT in NAMES_CLU){
			print(paste(NAMEIDENT,"split:",SPLIT))
			Idents(GEX)<-NAMEIDENT
			PLOT=DimPlot(GEX, reduction="umap", label=TRUE,split.by=SPLIT,order=ORDER,pt.size=.5,
						 label.size=6,repel=T,ncol=NCOL) + 
				theme_minimal(base_size=20)+
				theme(legend.position='top')+
				ggtitle(paste(NAMEIDENT)) + coord_fixed()
			OUTNAME=paste0(SUBDIR,"/umap.",NAMEIDENT,".Split_",SPLIT,".pdf");pdf(OUTNAME,width=WIDTH,height=HEIGHT,useDingbats=F);print(PLOT);dev.off()
		}
	}
	
	
	rm(NAMES_CLU,NAMEIDENT,WIDTH,HEIGHT,NCOL,ORDER,SPLIT,OUTNAME)
	# axis.line = element_line(size=1),	# text = element_text(size = 12),	# axis.text = element_text(size = 12),	# axis.ticks = element_line(size=1),
}
#### GENE MARKERS: FeaturePlot / plot_density ####
SUBDIR="13.GeneMarkers"; dir.create(SUBDIR)
#problem at index 14 Nos2, index 29 Aicda ---> zero expresion in all cells
if(FALSE){
	Idents(GEX)<-"RNA_snn_res.2"
	VlnPlot(GEX,features = "Aicda")
	VlnPlot(GEX,features = "Nos2")
}
# plot the markers
SUBDIR="13.GeneMarkers"; dir.create(SUBDIR)
if(F){
	#MarkersTable<-read.xlsx("/home/ga94rac/LRZ Sync+Share/PROJECT-GILDAS/MarkersTable.xlsx")
	#	grep("TNFSF13B",x = rownames(GEX),ignore.case = T, value = T)
	#SPLIT="split" #GEX@meta.data
	for( index in c( 1:45 ) ){
		# other cell types
		if(index==1){ FEATURES=c("Flt3");TITLE="Dendritic cells"} #"Itgax","Cd209a","Siglech"
		#if( index==18 ){ FEATURES=c("Siglech","Ly6c1","Ly6c2","Bst2");TITLE="Dendritic cells plasmacytoids"}
		if(index==2){ FEATURES=c("Adgre1","Cd14","Fcgr1");TITLE="Macrophages"}
		if(index==3){ FEATURES=c("Pros1","Siglech");TITLE="Microglia"} #"Mrc2","P2yr12"
		if(index==4){ FEATURES=c("Csf1r","Cd68");TITLE="Macrophages+Monocytes" }
		if(index==5){ FEATURES=c("Slc1a2");TITLE="Glial cell"} # "Gfap"
		if(index==6){ FEATURES=c("Pecam1");TITLE="Platelet-Endothelial"}
		if(index==7){ FEATURES=c("Cd34");TITLE="Mastocytes"}
		if(index==8){ FEATURES=c("Mbp","Cldn11");TITLE="Myelin"}
		#if(index==7){ FEATURES=c("Aldoc");TITLE="Purkinje cells"} #platelets and mast cells (MCs)
		if(index==9){ FEATURES=c("Krt18");TITLE="Epithelial"}
		if(index==10){ FEATURES=c("Xcr1","Clec9a","Sirpa");TITLE="Dendritic cells classical"} # "Itgam" ??
		if(index==11){ FEATURES=c("Rbfox3");TITLE="Neurons"}
		if(index==12){ FEATURES=c("Mog","Plp1");TITLE="Myelin"}
		if(index==13){ FEATURES=c("Hbb-bt","Hbb-bs");TITLE="Erythroid cells"} # Hbb genes # "Ter119" > Ly76   ##CD47? Gypa ??? Humain gene present in mouse, WEIRD
		if(index==14){ FEATURES=c("Tlr2","Tlr4","Il12a","Il18","Ly6c1");TITLE="Macrophages M1"} #
		#"Nos2" creates problems: ignored
		
		#M1 macrophges: Nos2, Tlr2, Tlr4, IL12a, Il18, Ly6c1 (higher than in M2), (Cd86 and CD80, which we tested already for Bcells)
		if(index==15){ FEATURES=c("Cd163","Mrc1","Tgfb1","Tgfb2","Tgfb3","Arg1","Il4ra","Pparg");TITLE="Macrophages M2"}
		#M2 macrophages: cd163, Mrc1 (= cd206), Tgfb1, Tgfb2, Tgfb3, Arg1, Il4r, Pparg
		if(index==16){ FEATURES=c("S100a8","S100a9");TITLE="Granulocyte" }
		if(index==17){ FEATURES=c("Olig1");TITLE="Oligodendrocytes" }
		
		
		# B-cell subset signature
		if(index==18){ FEATURES=c("Cd19","Cd79a","Cd79b","Ly6d","Ms4a1","Bank1","Blk");TITLE="B-cell"} #Ms4a1 = CD20
		if(index==19){ FEATURES=c("Hmmr");TITLE="B-cell mature"}
		if(index==20){ FEATURES=c("Cd38","Cd22");TITLE="B-cell Follicular"}
		if(index==21){ FEATURES=c("Cd1d1","Cd9","Fcer2a","Bcl6");TITLE="B-cell Marginal (Sec. Lymph. Organs)"}
		if(index==22){ FEATURES=c("Cd80","Nt5e","Pdcd1lg2","Cd84","Cd86");TITLE="B-cell memory"} #"Cd80","Nt5e","Pdcd1lg2","Cd84","Cd86"
		if(index==23){ FEATURES=c("Ighm","Cd93","Cd24a","Ebf1");TITLE="B-cell immature"}
		if(index==24){ FEATURES=c("Ighg2b","Ighg3","Cd86","Cd80","Cd9","Cd1d1","Tnfrsf8");TITLE="B-cell activated"} # expression of IgG genes, expression of IgM but no IgD
		if(index==25){ FEATURES=c("Sdc1","Tnfrsf17");TITLE="Plasma B-cell"} #Sdc1=CD138 ... Sdc1","Tnfrsf17","Il6ra","Cxcr4","Cd320","Prdm1","Irf4","Xbp1"
		if(index==26){ FEATURES=c("Ighd","Ighm");TITLE="B-cell mature when ighm+ighd same time"}
		if(index==27){ FEATURES=c("Ighg1","Ighg2c","Ighg2b","Ighg3","Ighe","Igha");TITLE="B-cell antibodies"}
		
		if(index==28){ FEATURES=c("Cd1d1","S1pr3","Tlr3");TITLE="Marginal zone B-cell" }
		if(index==29){ FEATURES=c("Basp1","Fas","Neil1","Plxnb2","Rgs13","S1pr2","Tnfsf9");TITLE="Germinal center B-cell" }
		#"Aicda" creating problems, ignored
		
		if(index==30){ FEATURES=c("Cd5");TITLE="B1 B-cell" }
		#	index=0; FEATURES=c("Tnfsf13b");TITLE="Act. B-cell"
		
		# T-cell subset signature
		if(index==31){ FEATURES=c("Ifngr2","Il21r","Il17ra","Foxp1","Sell","Ccr7");TITLE="T-cell naive"}
		if(index==32){ FEATURES=c("Gja1");TITLE="T-cell mature"} #
		if(index==33){ FEATURES=c("Pdcd1","Il2","Rora");TITLE="T-cell activated"} #"Cd44" "Il2ra" "Tnf" "Cd69"  DOWN REG of CD62L (Sell), CD25 (Il2ra) "Nme1", "Ifit3" 
		if(index==34){ FEATURES=c("Cd3g","Cd3d","Cd3e","Cd28","Cd247"); TITLE="T-cell"}
		if(index==35){ FEATURES=c("Il10");TITLE="B-cell reg"}
		if(index==36){ FEATURES=c("Foxp3");TITLE="T-cell reg"} # Ctla4="CTL4,CD152", "Ikzf2","Nrp1","Fosb","Tnfsf11","Irf8"
		if(index==37){ FEATURES=c("Cd4");TITLE="T-cells Cd4"}
		if(index==38){ FEATURES=c("Cd8a");TITLE="T-cells Cd8"}
		if(index==39){ FEATURES=c("Cxcr6","Il18r1","Tbx21","Ifng");TITLE="Th1"} #"Csf1","Il12rb2","Klrc1","Hopx"
		if(index==40){ FEATURES=c("Il4");TITLE="Th2"} #"Tnfsf13b","Batf","Nfil3","Atf5"
		if(index==41){ FEATURES=c("Il2","Il17f","Il21","Rora");TITLE="Th17"} #"Il17a","Tnfrsf13b","Ptgfrn","Ahr","Irf4","Plagl2"
		if(index==42){ FEATURES=c("Klrk1","Klre1","Klrd1","Klrb1c");TITLE="NK-cell"} #"Ncr1","Nkg7","Klra8","Klra1","Gzma","Klrb1c","Klrb1a"
		if(index==43){ FEATURES=c("Ptprc");TITLE="CD45"}
		if(index==44){ FEATURES=c("Tcrg-C1","Tcrg-C2","Trdc"); TITLE="gd T-cell"}  #gamma delta T cells: "Btnl6"
		if(index==45){ FEATURES=c("Cxcr5","Pdcd1","Ccr7");TITLE="T-cell follicular helper cell - Ccr7 LOW" }
		
		
		#grep("Cd24",x = rownames(GEX),ignore.case = T, value = T)
		#SPLIT=
		for(FEATURE in FEATURES){
			print(paste("index:",index,"current:",TITLE,"->",FEATURE))
			# GENE MARKERS: FeaturePlot
			if(FALSE){
				PLOT=FeaturePlot(GEX,reduction="umap",features=FEATURE,pt.size=1,label=T,repel=T,min.cutoff="q10",ncol=1,coord.fixed=T,order=T)
				#OUTNAME=paste0(SUBDIR,"/umap.",TITLE,".",FEATURE,".pdf");pdf(OUTNAME,width=8,height=8,useDingbats=F);print(PLOT);dev.off()
				OUTNAME=paste0(SUBDIR,"/umap.",TITLE,".",FEATURE,".png");png(OUTNAME,width=500,height=500);print(PLOT);dev.off()
			}
			# GENE MARKERS: FeaturePlot SPLIT
			if(FALSE){
				PLOT=FeaturePlot(GEX,reduction="umap",features=FEATURE,pt.size=1,label=F,split.by=SPLIT,
								 repel=T,min.cutoff="q10",ncol=1,coord.fixed=T,order=T)
				#OUTNAME=paste0(SUBDIR,"/umap.",TITLE,".",FEATURE,".pdf");pdf(OUTNAME,width=8,height=8,useDingbats=F);print(PLOT);dev.off()
				OUTNAME=paste0(SUBDIR,"/umap.","split_",SPLIT,".",TITLE,".",FEATURE,".png");png(OUTNAME,width=1200,height=500);print(PLOT);dev.off()
			}
			
			# GENE MARKERS - nebulosa: plot_density: nicer than FeaturePlot
			if(T){
				PLOT=plot_density(GEX,reduction="umap",features=FEATURE,size=1)
				OUTNAME=paste0(SUBDIR,"/umap.density.",TITLE,".",FEATURE,".png");png(OUTNAME,width=500,height=500);print(PLOT);dev.off()
			}
			# GENE MARKERS: plot_density: SPLIT by tissue - 
			#the plot_density indeed does not plot using a split option
			
		}
	}
	
	#grep("Baff",x = rownames(GEX),ignore.case = T, value = T)
	
	
	system(paste0("mkdir -p ",SUBDIR,"/B-cells","; mv ",SUBDIR,"/umap*B-cell*.png ",SUBDIR,"/B-cells/"))
	system(paste0("mkdir -p ",SUBDIR,"/T-cells",
				  "; mv ",SUBDIR,"/umap*T-cell*.png ",SUBDIR,"/T-cells/",
				  "; mv ",SUBDIR,"/umap*.Th*.png ",SUBDIR,"/T-cells/"))
	
	system(paste0("mkdir -p ",SUBDIR,"/DC ; mv ",SUBDIR,"/umap*Dendritic*.png ",SUBDIR,"/DC/"))
	system(paste0("mkdir -p ",SUBDIR,"/Macrophages ; mv ",SUBDIR,"/umap*Macrophages*.png ",SUBDIR,"/Macrophages/"))
	system(paste0("mkdir -p ",SUBDIR,"/NK ; mv ",SUBDIR,"/umap*NK-cell*.png ",SUBDIR,"/NK/"))
	system(paste0("mkdir -p ",SUBDIR,"/Myelin ; mv ",SUBDIR,"/umap*Myelin*.png ",SUBDIR,"/Myelin/"))
	
	rm(index,FEATURES,PLOT,OUTNAME,FEATURE,TITLE,SUBDIR,SPLIT)
}

#### CELL TYPES ANNOT external tool -- hints ####
#http://webtools.mcdb.ucla.edu/			#checking cell types
if(F){
	Idents(GEX)<-"RNA_snn_res.2"
	#Idents(GEX)<-"CellTypeFinal"
	avgexp <- AverageExpression(GEX,assays="RNA")
	avgexp <- as.data.frame(avgexp$RNA)
	colnames(avgexp)<-paste0("cluster_",colnames(avgexp))
	#avgexp<-avgexp[order(colnames(avgexp)),]
	avgexp$gene <- rownames(avgexp)
	avgexp<-avgexp[, c(ncol(avgexp),1:(ncol(avgexp)-1)) ]
	
	write.table(avgexp,file = "avgexp.res.2.csv",quote=F,sep="\t",row.names=F,col.names=T)
	
	# Mouse Body Atlas		- ADD THE FOXP3 / no B cells apparently
	
	# Tabula muris (not enough cell types)
	# ImmGen (NOT AMAZING)
}
#### CLEAN CLUSTERS - update from gene markers and auto annotation ####
if(T){
	GEX@meta.data$CellTypeFinal <- as.character(GEX@meta.data$annot_RNA_snn_res.2_CIPR)
	#DF <- as.data.frame(GEX@reductions$umap@cell.embeddings)
	GEX$UMAP_1<-GEX@reductions$umap@cell.embeddings[,1]
	GEX$UMAP_2<-GEX@reductions$umap@cell.embeddings[,2]
	
	# NK cell as NK + NKT cell [cluster 17]
	if(T){
		GEX$CellTypeFinal [ GEX$CellTypeFinal == "NK cell" ] <- "NK + NK T cell"
		#table(GEX$CellTypeFinal, GEX$immune_chains_clean)
	}
	
	# 4,6,8,13,14,16		   -> naive T cell CD4+		(based on Ccr7, Foxp1   -   Cd4 /Cd8a expression )
	# cluster 0, 15		-> activated T cell			(based on Il2, Pdcd1   -   Cd4 /Cd8a expression )
	if(T){
		GEX$CellTypeFinal[ GEX$RNA_snn_res.2 %in% c("4","6","8","13","14","16") ] <- "naive T cell CD4+"
		GEX$CellTypeFinal[ GEX$RNA_snn_res.2 %in% c("0") ] <- "act. T cell CD4+"
	}
	# Stromal clu 12>> Oligodendrocytes			(based on Plp1, Cldn11, Mog)
	if(T){
		GEX$CellTypeFinal[ GEX$annot_RNA_snn_res.2_CIPR %in% c("Stromal") ] <- "Oligodendrocytes"
	}
	# Cd8a T cells (~= cluter 23)							(based on Cd8a expression)
	if(T){
		GEX$CellTypeFinal[ GEX$RNA_snn_res.2 %in% c("23") ] <- "T cell CD8+"
	}
	
	# Erythroid cell Cluster 22							(based on Hbb-bt / Hbb-bs markers)      			NOT: pre-B cell, NOT follicular B cell (based on GEDIT: http://webtools.mcdb.ucla.edu/)
	if(T){
		GEX$CellTypeFinal[ GEX$RNA_snn_res.2 %in% c("22") ] <- "Erythroid cell" # Blood cell
	}
	
	#Plasma-cells: clu 32  						(based on Sdc1, Tnfrsf17
	if(T){
		GEX$CellTypeFinal[ GEX$RNA_snn_res.2 %in% c("32") ]  <- "Plasma cell"
	}
	
	#Glial cells: 25, 30							"(based on Slc1a2"
	if(T){
		GEX$CellTypeFinal[ GEX$RNA_snn_res.2 %in% c("25","30") ]  <- "Glial"
	}
	
	
	#rm-1		"pre T-cell"
	if(T){
		GEX$CellTypeFinal[ GEX$RNA_snn_res.2 %in% c("18") ]  <- "rm-1"
		
		
	}
	#rm-2		Stem-Progenitor
	if(T){
		GEX$CellTypeFinal[ GEX$RNA_snn_res.2 %in% c("20","26","30") ]  <- "rm-2"
	}
	
}
#### PLOT clusters + annot ####
SUBDIR="14.Plots.CellTypeFinal";dir.create(SUBDIR)
## Do with and without rm-X clusters -> change also the name manually #### 
#BACK=GEX		GEX=BACK		rm(BACK)# GEX<-subset(GEX, CellTypeFinal %ni% c("rm-1","rm-2","rm-3","rm-4"))
if(F){
	ORDER=c("act. T cell CD4+"="#7cae00",
			"naive T cell CD4+"="#cd9600",
			"T cell CD8+"="#00c19a",
			"Treg"="#00be67",
			"NK + NK T cell"="#0cb702",
			"B cell"="#00bfc4",
			"Plasma cell"="#00b8e7",
			"DC"="#e68613",
			"Erythroid cell"="#ff68a1",
			"Glial"="#00a9ff",
			"Granulocyte"="#ff61cc",
			"Macrophage"="#aba300",
			"Monocyte"="#8494ff",
			"Oligodendrocytes"="#f8766d",
			"rm-2"="gray65","rm-1"="gray70")
	
	
	NAMES_CLU=c("CellTypeFinal") 
	GEX$split = paste(GEX$treatment, GEX$tag_tissue_thr1.5)
	CELLTYPES=unique(GEX$CellTypeFinal)
	CELLTYPES=grep("^rm-",CELLTYPES,value=T,invert=T)
	
	
	
	for ( NAMEIDENT in NAMES_CLU){
		print(paste("doing:",NAMEIDENT))
		Idents(GEX)<-NAMEIDENT
		PLOT=DimPlot(GEX,reduction="umap",label=T, cols=ORDER,pt.size=.1,label.size=4,repel=T) + 
			theme_minimal(base_size = 16)+
			theme(legend.position = 'top')+
			ggtitle(paste(NAMEIDENT)) + coord_fixed()
		OUTNAME=paste0(SUBDIR,"/umap.",NAMEIDENT,".pdf");pdf(OUTNAME,width=10,height=8,useDingbats=F);print(PLOT);dev.off()
		rm(PLOT,OUTNAME)
	}
	
	# remove the rm clusters (ignored)
	#CELLTYPES=unique(GEX$CellTypeFinal)
	#CELLTYPES=grep("^rm-",CELLTYPES,value=T,invert=T) 
	GEX.tmp=subset( GEX , CellTypeFinal %in% CELLTYPES)
	for ( NAMEIDENT in NAMES_CLU){
		print(paste("doing:",NAMEIDENT))
		Idents(GEX.tmp)<-NAMEIDENT
		PLOT=DimPlot(GEX.tmp,reduction="umap",label=T,order=ORDER,pt.size=.1,label.size=3,repel=T, cols=ORDER) + 
			theme_minimal(base_size = 16)+
			theme(legend.position = 'top')+
			ggtitle(paste(NAMEIDENT)) + coord_fixed()
		OUTNAME=paste0(SUBDIR,"/umap.",NAMEIDENT,".noRm.pdf");pdf(OUTNAME,width=10,height=8,useDingbats=F);print(PLOT);dev.off()
		rm(PLOT,OUTNAME)
	}
	
	# split
	for (index in 1:3){
		if(index==1){SPLIT="tag_tissue_thr1.5" ; WIDTH=20;HEIGHT=8; NCOL=4; LABSIZE=3}
		if(index==2){SPLIT="orig.ident"; WIDTH=20;HEIGHT=15; NCOL=3; LABSIZE=3}
		if(index==3){SPLIT="split" ; WIDTH=20;HEIGHT=12; NCOL=3; LABSIZE=3}
		
		for ( NAMEIDENT in NAMES_CLU){
			print(paste(NAMEIDENT,"split:",SPLIT))
			Idents(GEX)<-NAMEIDENT
			PLOT=DimPlot(GEX, reduction="umap", label=TRUE,split.by=SPLIT, cols=ORDER,pt.size=.5,label.size=LABSIZE,repel=T,ncol=NCOL) + 
				theme_minimal(base_size=20)+
				theme(legend.position='top')+
				ggtitle(paste(NAMEIDENT)) + coord_fixed()
			OUTNAME=paste0(SUBDIR,"/umap.",NAMEIDENT,".Split_",SPLIT,".pdf");pdf(OUTNAME,width=WIDTH,height=HEIGHT,useDingbats=F);print(PLOT);dev.off()
		}
	}
	
	# split - remove the rm clusters (ignored)
	#CELLTYPES=unique(GEX$CellTypeFinal)
	#CELLTYPES=grep("^rm-",CELLTYPES,value=T,invert=T) 
	GEX.tmp=subset( GEX , CellTypeFinal %in% CELLTYPES)
	for (index in 1:3){
		if(index==1){SPLIT="tag_tissue_thr1.5" ; WIDTH=20;HEIGHT=8; NCOL=4; LABSIZE=3}
		if(index==2){SPLIT="orig.ident"; WIDTH=20;HEIGHT=15; NCOL=3; LABSIZE=3}
		if(index==3){SPLIT="split" ; WIDTH=20;HEIGHT=12; NCOL=3; LABSIZE=3}
		
		for ( NAMEIDENT in NAMES_CLU){
			print(paste(NAMEIDENT,"split:",SPLIT))
			Idents(GEX.tmp)<-NAMEIDENT
			PLOT=DimPlot(GEX.tmp, reduction="umap", label=TRUE,split.by=SPLIT, cols=ORDER,pt.size=.5,label.size=LABSIZE,repel=T,ncol=NCOL) + 
				theme_minimal(base_size=20)+
				theme(legend.position='top')+
				ggtitle(paste(NAMEIDENT)) + coord_fixed()
			OUTNAME=paste0(SUBDIR,"/umap.",NAMEIDENT,".Split_",SPLIT,".noRm.pdf");pdf(OUTNAME,width=WIDTH,height=HEIGHT,useDingbats=F);print(PLOT);dev.off()
		}
	}
	
	# split even moooooorrrrre 
	if(F){
		CELLTYPES=unique(GEX$CellTypeFinal) # this is with rm clusters
		CELLTYPES=grep("^rm-",CELLTYPES,value=T,invert=T) 
		
		ELEMENTS=unique(GEX$split)
		for( ELEMENT in ELEMENTS ){
			GEX.tmp=subset( GEX , split %in% ELEMENT & CellTypeFinal %in% CELLTYPES)
			
			for ( NAMEIDENT in NAMES_CLU){
				print(paste(NAMEIDENT,"split:",ELEMENT))
				Idents(GEX.tmp)<-NAMEIDENT
				PLOT=DimPlot(GEX.tmp, reduction="umap", label=TRUE,cols=ORDER,pt.size=.5,label.size=3,repel=T,ncol=1) + #split.by=SPLIT, 
					theme_minimal(base_size=20)+
					theme(legend.position='top')+
					ggtitle(paste(NAMEIDENT)) + coord_fixed()
				#OUTNAME=paste0(SUBDIR,"/umap.",NAMEIDENT,".",ELEMENT,".pdf");pdf(OUTNAME,width=12,height=8,useDingbats=F);print(PLOT);dev.off()
				OUTNAME=paste0(SUBDIR,"/umap.",NAMEIDENT,".",ELEMENT,".noRm.pdf");pdf(OUTNAME,width=12,height=8,useDingbats=F);print(PLOT);dev.off()
			}
		}
		
		
		CELLTYPES=unique(GEX$CellTypeFinal) # this is with rm clusters
		CELLTYPES=grep("^rm-",CELLTYPES,value=T,invert=T) 
		ELEMENTS=unique(GEX$tag_tissue_thr1.5)
		for( ELEMENT in ELEMENTS ){
			GEX.tmp=subset( GEX , tag_tissue_thr1.5 %in% ELEMENT & CellTypeFinal %in% CELLTYPES)
			
			for ( NAMEIDENT in NAMES_CLU){
				print(paste(NAMEIDENT,"split:",ELEMENT))
				Idents(GEX.tmp)<-NAMEIDENT
				PLOT=DimPlot(GEX.tmp, reduction="umap", label=TRUE,cols=ORDER,pt.size=.5,label.size=3,repel=T,ncol=1) + #split.by=SPLIT, 
					theme_minimal(base_size=20)+
					theme(legend.position='top')+
					ggtitle(paste(NAMEIDENT)) + coord_fixed()
				#OUTNAME=paste0(SUBDIR,"/umap.",NAMEIDENT,".",ELEMENT,".pdf");pdf(OUTNAME,width=12,height=8,useDingbats=F);print(PLOT);dev.off()
				OUTNAME=paste0(SUBDIR,"/umap.",NAMEIDENT,".",ELEMENT,".noRm.pdf");pdf(OUTNAME,width=12,height=8,useDingbats=F);print(PLOT);dev.off()
			}
		}
		
		CELLTYPES=unique(GEX$CellTypeFinal) # this is with rm clusters
		CELLTYPES=grep("^rm-",CELLTYPES,value=T,invert=T) 
		ELEMENTS=unique(GEX$treatment)
		for( ELEMENT in ELEMENTS ){
			GEX.tmp=subset( GEX , treatment %in% ELEMENT & CellTypeFinal %in% CELLTYPES)
			
			for ( NAMEIDENT in NAMES_CLU){
				print(paste(NAMEIDENT,"split:",ELEMENT))
				Idents(GEX.tmp)<-NAMEIDENT
				PLOT=DimPlot(GEX.tmp, reduction="umap", label=TRUE,cols=ORDER,pt.size=.5,label.size=3,repel=T,ncol=1) + #split.by=SPLIT, 
					theme_minimal(base_size=20)+
					theme(legend.position='top')+
					ggtitle(paste(NAMEIDENT)) + coord_fixed()
				#OUTNAME=paste0(SUBDIR,"/umap.",NAMEIDENT,".",ELEMENT,".pdf");pdf(OUTNAME,width=12,height=8,useDingbats=F);print(PLOT);dev.off()
				OUTNAME=paste0(SUBDIR,"/umap.",NAMEIDENT,".",ELEMENT,".noRm.pdf");pdf(OUTNAME,width=12,height=8,useDingbats=F);print(PLOT);dev.off()
			}
		}
		
		
	
	}
	rm(NAMES_CLU,NAMEIDENT,WIDTH,HEIGHT,NCOL,ORDER,SPLIT,OUTNAME,PLOT)
	# axis.line = element_line(size=1),	# text = element_text(size = 12),	# axis.text = element_text(size = 12),	# axis.ticks = element_line(size=1),
	
	
	
	
}
#### DOTPLOT gene markers ####
SUBDIR="14.Plots.CellTypeFinal/";dir.create(SUBDIR)
if(T){
	#Zhaorong
	# FEATURES=c("Cd3g","Cd3e","Cd3d",#T-cell ok
	# 		   "Cd2", #T-cell / NK cells ???
	# 		   "Cd86","Cd80",#memory B
	# 		   "Cd68", #Monocytes
	# 		   "S100a9","S100a8",#Granulocytes
	# 		   "Itgb3bp","Krt18","Pecam1","Cd34","Cx3cr1","Itgam","Olig1","Cldn11","Mbp","Gja1","Aldoc","Slc1a2","Hmmr")
	
	
	
	FEATURES<-c("Cd4","Ccr7","Il2","Cd8a","Foxp3",
				"Klrk1",
				"Cd19","Ighd","Ighm","Sdc1",
				"Flt3","Cd68","S100a8",
				"Slc1a2","Olig1",
				"Hbb-bt")
	names(FEATURES)=c("CD4","T naive","T act.","CD8","Treg",
					  "NK",
					  "B","B","B","PC",
					  "DC","Mo+Ma","Gran.",
					  "Glial","Oligo",
					  "Eryth.")
	
	NAMEIDENT = "CellTypeFinal"
	Idents(GEX)<-NAMEIDENT
	Idents(GEX) <- factor(GEX$CellTypeFinal,levels=rev( c("naive T cell CD4+","act. T cell CD4+","T cell CD8+","Treg",
														  "NK + NK T cell",
														  "B cell","Plasma cell",
														  "DC","Macrophage","Monocyte","Granulocyte",
														  "Glial","Oligodendrocytes",
														  "Erythroid cell","rm-1","rm-2") ))

	IDENTS = unique(Idents(GEX))
	IDENTS<-IDENTS[IDENTS%ni%c("rm-1","rm-2","rm-3","rm-4")]
	
	PLOT=DotPlot(GEX,features=FEATURES,assay="RNA",idents=IDENTS,dot.scale=10) + 
		theme_minimal(base_size=20) + 
		theme(axis.text.x=element_text(angle=90)) + 
		ggtitle("Gene markers")
	#OUTNAME=paste0(SUBDIR,"/DotPlot_GeneMarkers_with_rm.pdf") ; pdf(OUTNAME,width=20,height=8) ; print(PLOT);dev.off()
	OUTNAME=paste0(SUBDIR,"/DotPlot_GeneMarkers.pdf") ; pdf(OUTNAME,width=20,height=8) ; print(PLOT);dev.off()
	
}

#options(max.print=800)
#table(GEX$RNA_snn_res.2,GEX$CellTypeFinal)

#saveRDS(GEX,"07.final.tmp.Rds")

#### Annotate the HTO + celltype ####
if(F){
	## Previous logic -- to update for simplicity ?
	if(F){
		# GEX$Selecti^
	# table(GEX$CellTypeFinal, GEX$immune_chains_clean3,useNA="ifany")
	}
	
	# New logic: remove when clear tag in the wrong cluster
	if(T){
		GEX$Selection_HTO_Celltype = TRUE
		
		table(GEX$CellTypeFinal)
		table(GEX$immune_chains_clean2)
		
		B_ANNOT_LIST=c("B cell","Plasma cell")
		T_ANNOT_LIST=c("act. T cell CD4+","naive T cell CD4+","T cell CD8+","Treg","NK + NK T cell") #NK T cells have are here and with NK cells
		NK_ANNOT_LIST=c("NK + NK T cell")
		GRA_MONO_ANNOT_LIST=c("Granulocyte","Monocyte")
		
		GEX$Selection_HTO_Celltype [ GEX$immune_chains_clean2 == "B-cell" & ! GEX$CellTypeFinal %in% B_ANNOT_LIST] <- FALSE
		table(GEX$tag_celltype)
		#"CD19"     "Doublet"  "CD4"      "CD8a"     "Ly6G"     "Negative" "NK" 
		#if it is not in the B-cell clusters but had a CD19 HTO, then do not select it.
		GEX$Selection_HTO_Celltype [ ! GEX$CellTypeFinal %in% B_ANNOT_LIST & GEX$tag_celltype %in% c("CD19") ] <- FALSE
		GEX$Selection_HTO_Celltype [ ! GEX$CellTypeFinal %in% T_ANNOT_LIST & GEX$tag_celltype %in% c("CD4","CD8a") ] <- FALSE
		GEX$Selection_HTO_Celltype [ ! GEX$CellTypeFinal %in% NK_ANNOT_LIST & GEX$tag_celltype %in% c("NK") ] <- FALSE
		GEX$Selection_HTO_Celltype [ ! GEX$CellTypeFinal %in% GRA_MONO_ANNOT_LIST & GEX$tag_celltype %in% c("Ly6G") ] <- FALSE
		
		rm(B_ANNOT_LIST,T_ANNOT_LIST,NK_ANNOT_LIST,GRA_MONO_ANNOT_LIST)
	}
	#Idents(GEX)<-"CellTypeFinal"
	#VlnPlot(GEX,features="Ly6G",assay = "AbCelltypes")
	#VlnPlot(GEX,features=c("Csf1r","Cd68"))
	#rownames(GEX@assays$AbCelltypes)
}
#### annotate using VDJ (B-cells only) ####
if(F){
	
	GEX$immune_chains_clean3 <- "other"
	GEX$immune_chains_clean3[GEX$immune_chains_clean %in% c("IGH,IGK","IGH,IGL") ] <- "B-cell"
	
	#TABLE=table(GEX$Selection_HTO_Celltype, GEX$immune_chains_clean3,useNA="ifany"); TABLE
	#TABLE=as.matrix(TABLE) ; TABLE/rowSums(TABLE)*100
	#			B-cell other
	# FALSE     84  1742
	# TRUE     826 12851
	
	#remove if VDJ B-cell present in unexpected cluster
	B_ANNOT_LIST=c("B cell","Plasma cell")
	
	#using the same selection vector, using first NA to check
	GEX$Selection_HTO_Celltype[  GEX$immune_chains_clean3 == "B-cell" & GEX$CellTypeFinal %ni% B_ANNOT_LIST ] <- NA
	
	TABLE=table(GEX$Selection_HTO_Celltype, GEX$CellTypeFinal,useNA="ifany"); TABLE
	
	GEX$FinalSelection <- GEX$Selection_HTO_Celltype
	GEX$FinalSelection[  GEX$immune_chains_clean3 == "B-cell" & GEX$CellTypeFinal %ni% B_ANNOT_LIST ] <- FALSE
	
	TABLE=table(GEX$FinalSelection, GEX$CellTypeFinal,useNA="ifany"); TABLE
	
}

#saveRDS(GEX,"07.final.Rds")

#### ------ FINAL FILE LOAD - please then filter  ------ ####

# system.time( CountPerCell <- apply(X=GEX@assays$RNA[,],MARGIN=2,FUN=function(x){return(sum(x>0))}) )
# dim(GEX@assays$RNA[,])
# GEX@meta.data$GeneExprCount <- CountPerCell
# ggplot( GEX@meta.data ) + geom_boxplot(aes(x=CellTypeFinal,y=GeneExprCount))
# summary(GEX@meta.data$GeneExprCount )

GEX<-readRDS("07.final.Rds")

# An object of class Seurat 
# 32300 features across 15503 samples within 4 assays 
# Active assay: RNA (32285 features, 2000 variable features)
# 3 other assays present: Ab, AbTissues, AbCelltypes
# 4 dimensional reductions calculated: pca, umap, tsne, harmony

GEX <- subset(GEX, FinalSelection)

#32300 features across 13677 samples within 4 assays
#13677/15503*100 = 88.22163 % good cells

# once this selection is done -> do the dotplot again, do the cluster UMAPs 

#### TABLES #### 
SUBDIR="15.tables";dir.create(SUBDIR)
if(F){
	TABLE=table(GEX$CellTypeFinal, GEX$tag_tissue_thr1.5);TABLE
	write.xlsx(TABLE,paste0(SUBDIR,"/table_CellTypeFinal_X_tissue_thr1.5.xlsx"))
	TABLE=table(GEX$CellTypeFinal, paste(GEX$tag_tissue_thr1.5,GEX$treatment));TABLE
	write.xlsx(TABLE,paste0(SUBDIR,"/table_CellTypeFinal_X_tissue_thr1.5+treatment.xlsx"))
	
}

#### do the DEG - mELT versus CSF ############ 
if(F){
	# selection for the DEG
	SUBDIR="16.DEG.mELT_vs_CSF";dir.create(SUBDIR)
	GEX$TEST <- paste(GEX$CellTypeFinal, GEX$split, sep="-") ; table(GEX$TEST) #GEX$split = paste(GEX$treatment, GEX$tag_tissue_thr1.5)
	TREATMENTS = c("CD20","Isotype")#TREATMENTS = unique(GEX$split)
	CELLTYPES=unique(GEX$CellTypeFinal)
	CELLTYPES<-CELLTYPES[CELLTYPES %ni% c( grep("^rm-",CELLTYPES,value=T),"Erythroid cells")]
	##selection for the DEG
	# SUBDIR="16.DEG.SelectionHTO.CSF_vs_mELT";dir.create(SUBDIR)
	# GEX$TEST <- paste(GEX$SelectionHTO, GEX$split, sep="-") ; table(GEX$TEST)
	# TREATMENTS = c("CD20","Isotype")#TREATMENTS = unique(GEX$split)
	# CELLTYPES=unique(GEX$SelectionHTO)
	# CELLTYPES<-CELLTYPES[! is.na(CELLTYPES)]
	## do the DEG - mELT versus CSF
	
	
	
	Idents(GEX)<-"TEST"
	# TABLE=table(GEX$TEST);TABLE
	# write.xlsx(TABLE,"table_annot+tissue+treatment.xlsx")
	#only once:
	avgexp <- AverageExpression(GEX,verbose=T,assays="RNA") #,slot='scale.data' -> we are not using that
	avgexp <- as.data.frame(avgexp$RNA)
	avgexp$gene <- rownames(avgexp)

	for( TREATMENT in TREATMENTS ) {
		#for( TISSUE in unique( GEX$Ab_tissue) ) {
		for( CELLTYPE in CELLTYPES ) {
			# HERE
			TREAT = paste(CELLTYPE, paste(TREATMENT,"mELT"), sep="-")
			CONTROL=paste(CELLTYPE, paste(TREATMENT,"CSF"),sep="-")
			
			if( ! sum( GEX$TEST == TREAT) >= 50){print(paste("not enough cells for:",TREAT));next}
			if( ! sum( GEX$TEST == CONTROL) >= 50){print(paste("not enough cells for:",CONTROL));next}
			print(paste(TREAT,"versus",CONTROL))
			Markers <- FindMarkers(GEX,ident.1=TREAT, ident.2=CONTROL)
			Markers$gene <- rownames(Markers)
			Markers$test = paste0(TREAT,"_vs_",CONTROL)
			Markers2 <- Markers %>% filter(p_val < 0.05 & p_val_adj < 0.1) %>% arrange(desc(avg_log2FC))
			
			Markers2_Expr = merge.data.frame(x = Markers2,y = avgexp[,c("gene",TREAT,CONTROL)], by = "gene") %>% arrange(desc(avg_log2FC))
			if( nrow(Markers2_Expr) > 0 ){
				OUTNAME=paste0(SUBDIR,"/Markers_",TREAT,"_vs_",CONTROL,".xlsx")
				write.xlsx(Markers2_Expr,OUTNAME)
				
				# plot some heatmap
				if(F){
					#Markers2_Expr
					SUB=avgexp[avgexp$gene %in% Markers2_Expr$gene, grep(pattern = "B cell",x = colnames(avgexp), value = T) ]
					GENES=SUB$gene
					SUB$gene=NULL
					MAT=as.matrix.data.frame(SUB)
					
					#rownames(MAT)<-GENES
					PLOT=pheatmap(MAT,cellwidth=10,cellheight=10,scale="row",drop_levels=T,cluster_cols=T,legend = F,main = paste0(TREAT," vs ",CONTROL) )
					OUTNAME=paste0("pheat.",TREAT,"_vs_",CONTROL,".pdf");pdf(OUTNAME,width=7,height=20,useDingbats=F);print(PLOT);dev.off()
					
					#remove ribosomal protein/genes
					IGNORE=c(grep(pattern="^Rp[ls]",x=Markers2_Expr$gene,value=T),
							 grep(pattern="^mt-",x=Markers2_Expr$gene,value=T) )
					length(Markers2_Expr$gene)
					length(IGNORE)
					table(Markers2_Expr[!Markers2_Expr$gene%in%IGNORE,"avg_log2FC"] > 0)
					
					Markers2_Expr[ ! Markers2_Expr$gene %in% IGNORE & Markers2_Expr$avg_log2FC < 0 , "gene" ]
					Markers2_Expr[ ! Markers2_Expr$gene %in% IGNORE & Markers2_Expr$avg_log2FC > 0 , "gene" ]
					
					MAT<- MAT[ grep(pattern = "^Rp[ls]",x = rownames(MAT),value = T,invert = T),]
					MAT<- MAT[ grep(pattern = "^mt-",x = rownames(MAT),value = T,invert = T),]
					PLOT=pheatmap(MAT,cellwidth=10,cellheight=10,scale="row",drop_levels=T,cluster_cols=T,legend = F,main = paste0(TREAT," vs ",CONTROL) )
					OUTNAME=paste0("pheat.",TREAT,"_vs_",CONTROL,"_noRiboMt.pdf");pdf(OUTNAME,width=7,height=20,useDingbats=F);print(PLOT);dev.off()
					
				}
				
			}else{
				OUTNAME=paste0(SUBDIR,"/Markers_",TREAT,"_vs_",CONTROL,"_EMPTY.xlsx")
				write.xlsx("",OUTNAME)
			}
			rm(Markers2_Expr,Markers,Markers2,OUTNAME)
			
			
			#The following features were omitted as they were not found in the scale.data slot for the RNA assay:
			# CELLS=rownames( GEX@meta.data [ GEX$TEST %in% c(TREAT,CONTROL) , ] )
			# PLOT=DoHeatmap(GEX,features=Markers2_Expr$gene,slot="data",size=0,combine=T,cells=CELLS) + NoLegend() # + NoAxes()
			# OUTNAME=paste0(SUBDIR,"/heat.",TREAT,"_vs_",CONTROL,".pdf");pdf(OUTNAME,width=7,height=20,useDingbats=F);print(PLOT);dev.off()
		}
	}
	rm(TREATMENT,CELLTYPE,TREAT,CONTROL,avgexp)
	
}
#### do the DEG - CD20 versus Isotype in CSF-CSF / mELT-mELT ############
# using TEST = celltypefinal ~ split   ->    CD20 CSF    CD20 mELT    CD20 none  Isotype CSF Isotype mELT Isotype none
if(F){
	## selection for the DEG
	SUBDIR="17.DEG.CD20_vs_Isotype";dir.create(SUBDIR)
	GEX$TEST <- paste(GEX$CellTypeFinal, GEX$split, sep="-") ; table(GEX$TEST)
	TISSUES = c("CSF","mELT")
	CELLTYPES=unique(GEX$CellTypeFinal)
	CELLTYPES<-CELLTYPES[CELLTYPES %ni% c(grep("^rm-",CELLTYPES,value=T),"Erythroid cells")]
	
	
	### selection for the DEG
	# SUBDIR="17.DEG.SelectionHTO.CD20_vs_Isotype";dir.create(SUBDIR)
	# GEX$TEST <- paste(GEX$SelectionHTO, GEX$split, sep="-") ; table(GEX$TEST)
	# TISSUES = c("CSF","mELT")
	# CELLTYPES=unique(GEX$SelectionHTO)
	# CELLTYPES<-CELLTYPES[! is.na(CELLTYPES)]
	## do the DEG - CD20 versus Isotype - with tissues
	Idents(GEX)<-"TEST"
	# TABLE=table(GEX$TEST);TABLE
	# write.xlsx(TABLE,"table_annot+tissue+treatment.xlsx")
	#only once:
	avgexp <- AverageExpression(GEX,verbose=T,assays="RNA")
	avgexp <- as.data.frame(avgexp$RNA)
	avgexp$gene <- rownames(avgexp)
	
	#for( TREATMENT in TREATMENTS ) {
	for( TISSUE in TISSUES ) { #c("CSF","mELT")
		for( CELLTYPE in CELLTYPES ) {
			# HERE
			TREAT = paste(CELLTYPE, paste("CD20",TISSUE), sep="-")
			CONTROL=paste(CELLTYPE, paste("Isotype",TISSUE),sep="-")
			
			if( ! sum( GEX$TEST == TREAT) >= 50){print(paste("not enough cells for:",TREAT));next}
			if( ! sum( GEX$TEST == CONTROL) >= 50){print(paste("not enough cells for:",CONTROL));next}
			print(paste(TREAT,"versus",CONTROL))
			
			Markers <- FindMarkers(GEX,ident.1=TREAT, ident.2=CONTROL)
			Markers$gene <- rownames(Markers)
			Markers$test = paste0(TREAT,"_vs_",CONTROL)
			Markers2 <- Markers %>% filter(p_val < 0.05 & p_val_adj < 0.1) %>% arrange(desc(avg_log2FC))
			
			Markers2_Expr = merge.data.frame(x = Markers2,y = avgexp[,c("gene",TREAT,CONTROL)], by = "gene") %>% arrange(desc(avg_log2FC))
			if( nrow(Markers2_Expr) > 0 ){
				OUTNAME=paste0(SUBDIR,"/Markers_",TREAT,"_vs_",CONTROL,".xlsx")
				write.xlsx(Markers2_Expr,OUTNAME)
				
				# plot some heatmap
				if(F){
					
					#Markers2_Expr
					
					SUB=avgexp[avgexp$gene %in% Markers2_Expr$gene, grep(pattern = "B cell",x = colnames(avgexp), value = T) ]
					GENES=SUB$gene
					SUB$gene=NULL
					MAT=as.matrix.data.frame(SUB)
					
					#rownames(MAT)<-GENES
					PLOT=pheatmap(MAT,cellwidth=10,cellheight=10,scale="row",drop_levels=T,cluster_cols=T,legend = F,main = paste0(TREAT," vs ",CONTROL) )
					OUTNAME=paste0("pheat.",TREAT,"_vs_",CONTROL,".pdf");pdf(OUTNAME,width=7,height=20,useDingbats=F);print(PLOT);dev.off()
					
					#remove ribosomal protein/genes
					IGNORE=c(grep(pattern="^Rp[ls]",x=Markers2_Expr$gene,value=T),
							 grep(pattern="^mt-",x=Markers2_Expr$gene,value=T) )
					length(Markers2_Expr$gene)
					length(IGNORE)
					table(Markers2_Expr[!Markers2_Expr$gene%in%IGNORE,"avg_log2FC"] > 0)
					
					Markers2_Expr[ ! Markers2_Expr$gene %in% IGNORE & Markers2_Expr$avg_log2FC < 0 , "gene" ]
					Markers2_Expr[ ! Markers2_Expr$gene %in% IGNORE & Markers2_Expr$avg_log2FC > 0 , "gene" ]
					
					MAT<- MAT[ grep(pattern = "^Rp[ls]",x = rownames(MAT),value = T,invert = T),]
					MAT<- MAT[ grep(pattern = "^mt-",x = rownames(MAT),value = T,invert = T),]
					PLOT=pheatmap(MAT,cellwidth=10,cellheight=10,scale="row",drop_levels=T,cluster_cols=T,legend = F,main = paste0(TREAT," vs ",CONTROL) )
					OUTNAME=paste0("pheat.",TREAT,"_vs_",CONTROL,"_noRiboMt.pdf");pdf(OUTNAME,width=7,height=20,useDingbats=F);print(PLOT);dev.off()
					
				}
				
				
			}else{
				OUTNAME=paste0(SUBDIR,"/Markers_",TREAT,"_vs_",CONTROL,"_EMPTY.xlsx")
				write.xlsx("",OUTNAME)
			}
			rm(Markers2_Expr,Markers,Markers2,OUTNAME)
			
			
			#The following features were omitted as they were not found in the scale.data slot for the RNA assay:
			# CELLS=rownames( GEX@meta.data [ GEX$TEST %in% c(TREAT,CONTROL) , ] )
			# PLOT=DoHeatmap(GEX,features=Markers2_Expr$gene,slot="data",size=0,combine=T,cells=CELLS) + NoLegend() # + NoAxes()
			# OUTNAME=paste0(SUBDIR,"/heat.",TREAT,"_vs_",CONTROL,".pdf");pdf(OUTNAME,width=7,height=20,useDingbats=F);print(PLOT);dev.off()
		}
	}
	rm(TISSUE,CELLTYPE,TREAT,CONTROL,avgexp)
}
############ selection for the DEG - allCD20_vs_allIsotype use treatment with tag_tissue_thr1.5  ############
# using TEST = celltypefinal ~ treatment     ->     CD20     Isotype			(we therefore include the "none")
if(F){
	SUBDIR="17.DEG.allCD20_vs_allIsotype";dir.create(SUBDIR)
	GEX$TEST <- paste(GEX$CellTypeFinal, GEX$treatment, sep="-") ; table(GEX$TEST) 
	TISSUES = c("CSF","mELT")
	CELLTYPES=unique(GEX$CellTypeFinal)
	CELLTYPES<-CELLTYPES[CELLTYPES %ni% c(grep("^rm-",CELLTYPES,value=T),"Erythroid cells")]
	## do the DEG - CD20 versus Isotype
	
	Idents(GEX)<-"TEST"
	# TABLE=table(GEX$TEST);TABLE
	# write.xlsx(TABLE,"table_annot+tissue+treatment.xlsx")
	#only once:
	avgexp <- AverageExpression(GEX,verbose=T,assays="RNA")
	avgexp <- as.data.frame(avgexp$RNA)
	avgexp$gene <- rownames(avgexp)

	for( CELLTYPE in CELLTYPES ) {
		# HERE
		TREAT = paste(CELLTYPE, "CD20", sep="-")
		CONTROL=paste(CELLTYPE, "Isotype",sep="-")
		
		if( ! sum( GEX$TEST == TREAT) >= 50){print(paste("not enough cells for:",TREAT));next}
		if( ! sum( GEX$TEST == CONTROL) >= 50){print(paste("not enough cells for:",CONTROL));next}
		print(paste(TREAT,"versus",CONTROL))
		
		Markers <- FindMarkers(GEX,ident.1=TREAT, ident.2=CONTROL)
		Markers$gene <- rownames(Markers)
		Markers$test = paste0(TREAT,"_vs_",CONTROL)
		Markers2 <- Markers %>% filter(p_val < 0.05 & p_val_adj < 0.1) %>% arrange(desc(avg_log2FC))
		
		Markers2_Expr = merge.data.frame(x = Markers2,y = avgexp[,c("gene",TREAT,CONTROL)], by = "gene") %>% arrange(desc(avg_log2FC))
		if( nrow(Markers2_Expr) > 0 ){
			OUTNAME=paste0(SUBDIR,"/Markers_",TREAT,"_vs_",CONTROL,".xlsx")
			write.xlsx(Markers2_Expr,OUTNAME)
			
			# plot some heatmap
			if(F){
				
				#Markers2_Expr
				
				SUB=avgexp[avgexp$gene %in% Markers2_Expr$gene, grep(pattern = "B cell",x = colnames(avgexp), value = T) ]
				GENES=SUB$gene
				SUB$gene=NULL
				MAT=as.matrix.data.frame(SUB)
				
				#rownames(MAT)<-GENES
				PLOT=pheatmap(MAT,cellwidth=10,cellheight=10,scale="row",drop_levels=T,cluster_cols=T,legend = F,main = paste0(TREAT," vs ",CONTROL) )
				OUTNAME=paste0("pheat.",TREAT,"_vs_",CONTROL,".pdf");pdf(OUTNAME,width=7,height=20,useDingbats=F);print(PLOT);dev.off()
				
				#remove ribosomal protein/genes
				IGNORE=c(grep(pattern="^Rp[ls]",x=Markers2_Expr$gene,value=T),
						 grep(pattern="^mt-",x=Markers2_Expr$gene,value=T) )
				length(Markers2_Expr$gene)
				length(IGNORE)
				table(Markers2_Expr[!Markers2_Expr$gene%in%IGNORE,"avg_log2FC"] > 0)
				
				Markers2_Expr[ ! Markers2_Expr$gene %in% IGNORE & Markers2_Expr$avg_log2FC < 0 , "gene" ]
				Markers2_Expr[ ! Markers2_Expr$gene %in% IGNORE & Markers2_Expr$avg_log2FC > 0 , "gene" ]
				
				MAT<- MAT[ grep(pattern = "^Rp[ls]",x = rownames(MAT),value = T,invert = T),]
				MAT<- MAT[ grep(pattern = "^mt-",x = rownames(MAT),value = T,invert = T),]
				PLOT=pheatmap(MAT,cellwidth=10,cellheight=10,scale="row",drop_levels=T,cluster_cols=T,legend = F,main = paste0(TREAT," vs ",CONTROL) )
				OUTNAME=paste0("pheat.",TREAT,"_vs_",CONTROL,"_noRiboMt.pdf");pdf(OUTNAME,width=7,height=20,useDingbats=F);print(PLOT);dev.off()
				
			}
		}else{
			OUTNAME=paste0(SUBDIR,"/Markers_",TREAT,"_vs_",CONTROL,"_EMPTY.xlsx")
			write.xlsx("",OUTNAME)
		}
		rm(Markers2_Expr,Markers,Markers2,OUTNAME)
		
		
		#The following features were omitted as they were not found in the scale.data slot for the RNA assay:
		# CELLS=rownames( GEX@meta.data [ GEX$TEST %in% c(TREAT,CONTROL) , ] )
		# PLOT=DoHeatmap(GEX,features=Markers2_Expr$gene,slot="data",size=0,combine=T,cells=CELLS) + NoLegend() # + NoAxes()
		# OUTNAME=paste0(SUBDIR,"/heat.",TREAT,"_vs_",CONTROL,".pdf");pdf(OUTNAME,width=7,height=20,useDingbats=F);print(PLOT);dev.off()
	}
	rm(TISSUE,CELLTYPE,TREAT,CONTROL,avgexp)
}



#### Welch’s t-test - cell type enrichment per tissues #### 
#plotting fold change (log10) against p value (−log10) based on beta-binomial regression
SUBDIR="15.Welch-t-test";dir.create(SUBDIR)
#all cells - mELT reference versus tissues ~ FOR EACH CELLTYPEFINAL
OUTNAME_DATA="Welch_t_tests_CellType"
if(F){
	GEX$mouse <- factor(GEX$orig.ident)
	GEX$tissue<-factor(GEX$tag_tissue_thr1.5)
	GEX$CellTypeFinal_Welch <- factor(GEX$CellTypeFinal)
	
	
	#table(GEX$CellTypeFinal_Welch)
	# # # # # # calculate percent of celltype per tissus for X mice - T-tests
	# n = count of cell per mouse / tissue / celltype
	# we want the percent of celltype in tissue / mouse
	COUNTS <- GEX@meta.data %>% 
		filter( CellTypeFinal_Welch %ni% c("Erythroid cells",grep("^rm-",unique(GEX$CellTypeFinal_Welch),value=T)) ) %>%
		filter( tissue %ni% c("none") ) %>%
		group_by( mouse, tissue, CellTypeFinal_Welch ) %>% count() %>% 
		group_by( mouse,tissue ) %>% mutate( mouse.tissue.count=sum(n)) %>%
		group_by( tissue, CellTypeFinal_Welch ) %>% mutate( tissue.CellType.count=sum(n)) %>%
		mutate( mouse.count=length(unique(mouse)), percent.celltype.in.tissue.per.mouse = n / mouse.tissue.count ) %>% arrange(tissue)
	
	COUNTS=as.data.frame(COUNTS)
	OUTNAME=paste0(SUBDIR,"/table.",OUTNAME_DATA,".xlsx")
	write.xlsx(COUNTS,OUTNAME)
	#COUNTS<-read.xlsx( paste0(SUBDIR,"/table_counts.xlsx"))
	
	# # # # # # Calculate T-tests
	TISSUES<-c("CSF" ) # # "CSF", unique( COUNTS$tissue )
	CELLTYPES<-unique(COUNTS$CellTypeFinal_Welch)
	
	BONFERRONI=length(CELLTYPES) #compared to mELT 
	OUTNAME=paste0(SUBDIR,"/table.",OUTNAME_DATA,".txt")
	
	cat(c("target","tissue","celltype","treat_count","ref_count","ratio","Pvalue",paste0("Padj_bonferroni",BONFERRONI),"Tscore\n"),file=OUTNAME,append=F,sep="\t")
	for(CELLTYPE in CELLTYPES){
		for(TISSUE in TISSUES){
			D.mELT <- COUNTS[ COUNTS$tissue == "mELT" & COUNTS$CellTypeFinal_Welch == CELLTYPE, ]
			D.other <- COUNTS[ COUNTS$tissue == TISSUE & COUNTS$CellTypeFinal_Welch == CELLTYPE, ]
			REFERENCE = D.mELT$percent.celltype.in.tissue.per.mouse
			TREATMENT = D.other$percent.celltype.in.tissue.per.mouse
			RATIO = mean(TREATMENT) / mean(REFERENCE)
			# test normality, should be at least 3 obs.
			if( length(REFERENCE) >2 & length(TREATMENT) >2 ){
				SHA_TEST = shapiro.test( c(REFERENCE,TREATMENT) )#plot(density( c(REFERENCE,TREATMENT) )) ; #library("ggpubr") ; #ggqqplot( c(REFERENCE,TREATMENT) )
				if(SHA_TEST$p.value<0.05){
					PValue=1
					Statistic="not norm. distri."
				}else{
					TEST = t.test(x = REFERENCE,y = TREATMENT) # at least 3 observation AND normale distribution
					PValue=TEST$p.value
					Statistic=TEST$statistic
					rm(TEST)
				}
				rm(SHA_TEST)
			}else{
				PValue=1
				Statistic="< 3 obs."
			}
			
			#print(paste("mELT",TISSUE,CELLTYPE,PValue,Statistic,sep=";"))
			cat( c("mELT",TISSUE,CELLTYPE,length(TREATMENT),length(REFERENCE),RATIO,PValue,PValue*BONFERRONI, paste0(Statistic,"\n")),file=OUTNAME,append=T,sep="\t")
		}
	}
	
	TTEST=read.table(OUTNAME,h=T,sep="\t")
	TTEST$log10FC <- log10(TTEST$ratio)
	TTEST$log10Pval_m <- -log10(TTEST$Pvalue) # minus log10 : high values = low p-value
	write.xlsx(TTEST, paste0(OUTNAME,".xlsx"))
	
	LAB_SIZE=.4
	NUDGE_X=0 #-.1
	NUDGE_Y=.1
	
	BONFNAME=grep("Padj_bonferroni",names(TTEST),value=T)
	TTEST$celltype2 <- ifelse( TTEST[,BONFNAME] < 0.05, paste0(TTEST$celltype," *"),TTEST$celltype )
	
	PLOT=ggplot(TTEST) + 
		geom_hline(yintercept=-log10(0.05),linetype="dashed",col="red") + 
		geom_hline(yintercept=-log10(0.05/BONFERRONI),linetype="dashed",col="orange") +
		geom_vline(xintercept=0,linetype="dashed",col="black") + 
		geom_point(aes(log10FC,log10Pval_m,size=abs(log10FC)*log10Pval_m ,col=celltype),show.legend=F) + #to remove color: here
		geom_label_repel(aes(log10FC,log10Pval_m,label=celltype2,size=LAB_SIZE,col=celltype),max.overlaps = 40,min.segment.length=0,nudge_x=NUDGE_X,nudge_y=NUDGE_Y,show.legend=F) + 
		scale_x_continuous(limits = c(-1,1)) +
		theme_bw(base_size=22) + 
		ggtitle( label = paste("Celltypes % (reference: mELT)"),
				 subtitle = paste("Bonferroni correction on",BONFERRONI,"tests")) + 
		facet_wrap(~tissue)
	OUTNAME=paste0(SUBDIR,"/Figure_Welch.",OUTNAME_DATA,".pdf") ; pdf(OUTNAME,width=10,height=10) ; print(PLOT);dev.off()
	print(paste("done:",OUTNAME_DATA))
	
	rm(COUNTS,D.mELT,D.other,TTEST,PLOT,avgexp)
}
#all cells - mELT reference versus tissues ~ FOR EACH CELLTYPEFINAL -- general (B-cell, T-cells...)
OUTNAME_DATA="Welch_t_tests_GeneralCellType"
if(F){
	GEX$mouse <- factor(GEX$orig.ident)
	GEX$tissue<-factor(GEX$tag_tissue_thr1.5)
	GEX$CellTypeFinal_Welch <- factor(GEX$CellTypeFinal)
	
	#grep("B",levels(GEX$CellTypeFinal_Welch),value=T) ; grep("B",levels(GEX$CellTypeFinal_Welch),value=F)
	levels(GEX$CellTypeFinal_Welch) [ grep("B",levels(GEX$CellTypeFinal_Welch),value=F) ] <- "B cell"
	levels(GEX$CellTypeFinal_Welch) [ levels(GEX$CellTypeFinal_Welch) %in% c("Plasma cell") ] <- "B cell"
	levels(GEX$CellTypeFinal_Welch) [ grep("T",levels(GEX$CellTypeFinal_Welch),value=F) ] <- "T cell"
	levels(GEX$CellTypeFinal_Welch) [ levels(GEX$CellTypeFinal_Welch) %in% c("DC","Monocyte","Glial") ] <- "Myeloid cells"
	
	levels(GEX$CellTypeFinal_Welch) [ levels(GEX$CellTypeFinal_Welch) %in% c("NK cell","NK T cell") ] <- "NK cell"
	#?dplyr::recode() # case_when() # recode_factor()
	
	# # # # # # calculate percent of celltype per tissus for X mice - T-tests
	# n = count of cell per mouse / tissue / celltype
	# we want the percent of celltype in tissue / mouse
	COUNTS <- GEX@meta.data %>% 
		filter( CellTypeFinal_Welch %ni% c("Erythroid cell",grep("^rm-",unique(GEX$CellTypeFinal_Welch),value=T) ) ) %>%
		filter( tissue %ni% c("none") ) %>%
		group_by( mouse, tissue, CellTypeFinal_Welch ) %>% count() %>% 
		group_by( mouse,tissue ) %>% mutate( mouse.tissue.count=sum(n)) %>%
		group_by( tissue, CellTypeFinal_Welch ) %>% mutate( tissue.CellType.count=sum(n)) %>%
		mutate( mouse.count=length(unique(mouse)), percent.celltype.in.tissue.per.mouse = n / mouse.tissue.count ) %>% arrange(tissue)
	
	COUNTS=as.data.frame(COUNTS)
	OUTNAME=paste0(SUBDIR,"/table.",OUTNAME_DATA,".xlsx")
	write.xlsx(COUNTS,OUTNAME)
	#COUNTS<-read.xlsx( paste0(SUBDIR,"/table_counts.xlsx"))
	
	# # # # # # Calculate T-tests
	TISSUES<-c("CSF" ) # # "CSF", unique( COUNTS$tissue )
	CELLTYPES<-unique(COUNTS$CellTypeFinal_Welch)
	BONFERRONI=length(CELLTYPES) #compared to mELT 
	OUTNAME=paste0(SUBDIR,"/table.",OUTNAME_DATA,".txt")
	
	cat(c("target","tissue","celltype","treat_count","ref_count","ratio","Pvalue",paste0("Padj_bonferroni",BONFERRONI),"Tscore\n"),file=OUTNAME,append=F,sep="\t")
	for(CELLTYPE in CELLTYPES){
		for(TISSUE in TISSUES){
			D.mELT <- COUNTS[ COUNTS$tissue == "mELT" & COUNTS$CellTypeFinal_Welch == CELLTYPE, ]
			D.other <- COUNTS[ COUNTS$tissue == TISSUE & COUNTS$CellTypeFinal_Welch == CELLTYPE, ]
			REFERENCE = D.mELT$percent.celltype.in.tissue.per.mouse
			TREATMENT = D.other$percent.celltype.in.tissue.per.mouse
			RATIO = mean(TREATMENT) / mean(REFERENCE)
			# test normality, should be at least 3 obs.
			if( length(REFERENCE) >2 & length(TREATMENT) >2 ){
				SHA_TEST = shapiro.test( c(REFERENCE,TREATMENT) )#plot(density( c(REFERENCE,TREATMENT) )) ; #library("ggpubr") ; #ggqqplot( c(REFERENCE,TREATMENT) )
				if(SHA_TEST$p.value<0.05){
					PValue=1
					Statistic="not norm. distri."
				}else{
					TEST = t.test(x = REFERENCE,y = TREATMENT) # at least 3 observation AND normale distribution
					PValue=TEST$p.value
					Statistic=TEST$statistic
					rm(TEST)
				}
				rm(SHA_TEST)
			}else{
				PValue=1
				Statistic="< 3 obs."
			}
			
			#print(paste("mELT",TISSUE,CELLTYPE,PValue,Statistic,sep=";"))
			cat( c("mELT",TISSUE,CELLTYPE,length(TREATMENT),length(REFERENCE),RATIO,PValue,PValue*BONFERRONI, paste0(Statistic,"\n")),file=OUTNAME,append=T,sep="\t")
		}
	}
	
	TTEST=read.table(OUTNAME,h=T,sep="\t")
	TTEST$log10FC <- log10(TTEST$ratio)
	TTEST$log10Pval_m <- -log10(TTEST$Pvalue) # minus log10 : high values = low p-value
	write.xlsx(TTEST, paste0(OUTNAME,".xlsx"))
	
	LAB_SIZE=.4
	NUDGE_X=0 #-.1
	NUDGE_Y=.1
	
	BONFNAME=grep("Padj_bonferroni",names(TTEST),value=T)
	TTEST$celltype2 <- ifelse( TTEST[,BONFNAME] < 0.05, paste0(TTEST$celltype," *"),TTEST$celltype )
	
	PLOT=ggplot(TTEST) + 
		geom_hline(yintercept=-log10(0.05),linetype="dashed",col="red") + 
		geom_hline(yintercept=-log10(0.05/BONFERRONI),linetype="dashed",col="orange") +
		geom_vline(xintercept=0,linetype="dashed",col="black") + 
		geom_point(aes(log10FC,log10Pval_m,size=abs(log10FC)*log10Pval_m ,col=celltype),show.legend=F) + #to remove color: here
		geom_label_repel(aes(log10FC,log10Pval_m,label=celltype2,size=LAB_SIZE,col=celltype),min.segment.length=0,nudge_x=NUDGE_X,nudge_y=NUDGE_Y,show.legend=F) + 
		scale_x_continuous(limits = c(-1,1)) +
		theme_bw(base_size=22) + 
		ggtitle( label = paste("Celltypes % (reference: mELT)"),
				 subtitle = paste("Bonferroni correction on",BONFERRONI,"tests")) + 
		facet_wrap(~tissue)
	OUTNAME=paste0(SUBDIR,"/Figure_Welch.",OUTNAME_DATA,".pdf") ; pdf(OUTNAME,width=10,height=10) ; print(PLOT);dev.off()
	
}
#         mELT reference versus tissues ~ Only B cell subsets
OUTNAME_DATA="Welch_t_tests_Bcell"
if(F){
	GEX$mouse <- factor(GEX$orig.ident)
	GEX$tissue<-factor(GEX$tag_tissue_thr1.5)
	GEX$CellTypeFinal_Welch <- factor(GEX$CellTypeFinal)
	
	# # # # # # calculate percent of celltype per tissus for X mice - T-tests
	# n = count of cell per mouse / tissue / celltype
	# we want the percent of celltype in tissue / mouse
	grep("B",levels(GEX$CellTypeFinal_Welch),value=T)
	COUNTS <- GEX@meta.data %>% 
		filter( CellTypeFinal_Welch %in% c("Plasma cell","B cell") ) %>%
		filter( tissue %ni% c("none") ) %>%
		group_by( mouse, tissue, CellTypeFinal_Welch ) %>% count() %>% 
		group_by( mouse,tissue ) %>% mutate( mouse.tissue.count=sum(n)) %>%
		group_by( tissue, CellTypeFinal_Welch ) %>% mutate( tissue.CellType.count=sum(n)) %>%
		mutate( mouse.count=length(unique(mouse)), percent.celltype.in.tissue.per.mouse = n / mouse.tissue.count ) %>% arrange(tissue)
	
	COUNTS=as.data.frame(COUNTS)
	OUTNAME=paste0(SUBDIR,"/table.",OUTNAME_DATA,".xlsx")
	write.xlsx(COUNTS,OUTNAME)
	#COUNTS<-read.xlsx( paste0(SUBDIR,"/table_counts.xlsx"))
	
	# # # # # # Calculate T-tests
	TISSUES<-c("CSF" ) # # "CSF", unique( COUNTS$tissue )
	CELLTYPES<-unique(COUNTS$CellTypeFinal_Welch)
	
	BONFERRONI=length(CELLTYPES) #compared to mELT 
	OUTNAME=paste0(SUBDIR,"/table.",OUTNAME_DATA,".txt")
	
	cat(c("target","tissue","celltype","treat_count","ref_count","ratio","Pvalue",paste0("Padj_bonferroni",BONFERRONI),"Tscore\n"),file=OUTNAME,append=F,sep="\t")
	for(CELLTYPE in CELLTYPES){
		for(TISSUE in TISSUES){
			D.mELT <- COUNTS[ COUNTS$tissue == "mELT" & COUNTS$CellTypeFinal_Welch == CELLTYPE, ]
			D.other <- COUNTS[ COUNTS$tissue == TISSUE & COUNTS$CellTypeFinal_Welch == CELLTYPE, ]
			REFERENCE = D.mELT$percent.celltype.in.tissue.per.mouse
			TREATMENT = D.other$percent.celltype.in.tissue.per.mouse
			RATIO = mean(TREATMENT) / mean(REFERENCE)
			# test normality, should be at least 3 obs.
			if( length(REFERENCE) >2 & length(TREATMENT) >2 ){
				SHA_TEST = shapiro.test( c(REFERENCE,TREATMENT) )#plot(density( c(REFERENCE,TREATMENT) )) ; #library("ggpubr") ; #ggqqplot( c(REFERENCE,TREATMENT) )
				if(SHA_TEST$p.value<0.05){
					PValue=1
					Statistic="not norm. distri."
				}else{
					TEST = t.test(x = REFERENCE,y = TREATMENT) # at least 3 observation AND normale distribution
					PValue=TEST$p.value
					Statistic=TEST$statistic
					rm(TEST)
				}
				rm(SHA_TEST)
			}else{
				PValue=1
				Statistic="< 3 obs."
			}
			
			#print(paste("mELT",TISSUE,CELLTYPE,PValue,Statistic,sep=";"))
			cat( c("mELT",TISSUE,CELLTYPE,length(TREATMENT),length(REFERENCE),RATIO,PValue,PValue*BONFERRONI, paste0(Statistic,"\n")),file=OUTNAME,append=T,sep="\t")
		}
	}
	
	TTEST=read.table(OUTNAME,h=T,sep="\t")
	TTEST$log10FC <- log10(TTEST$ratio)
	TTEST$log10Pval_m <- -log10(TTEST$Pvalue) # minus log10 : high values = low p-value
	write.xlsx(TTEST, paste0(OUTNAME,".xlsx"))
	
	LAB_SIZE=.1
	NUDGE_X=0
	NUDGE_Y=.1
	BONFNAME=grep("Padj_bonferroni",names(TTEST),value=T)
	TTEST$celltype2 <- ifelse( TTEST[,BONFNAME] < 0.05, paste0(TTEST$celltype," *"),TTEST$celltype )
	
	PLOT=ggplot(TTEST) + 
		geom_hline(yintercept=-log10(0.05),linetype="dashed",col="red") + 
		geom_hline(yintercept=-log10(0.05/BONFERRONI),linetype="dashed",col="orange") +
		geom_vline(xintercept=0,linetype="dashed",col="black") + 
		geom_point(aes(log10FC,log10Pval_m,size=abs(log10FC)*log10Pval_m ,col=celltype),show.legend=F) + #to remove color: here
		geom_label_repel(aes(log10FC,log10Pval_m,label=celltype2,size=LAB_SIZE,col=celltype),force_pull=10,force=10,
						 min.segment.length=0,nudge_x=NUDGE_X,nudge_y=NUDGE_Y,show.legend=F) + 
		scale_x_continuous(limits = c(-1,1)) +
		theme_bw(base_size=22) + 
		ggtitle( label = paste("Celltypes % (reference: mELT)"),
				 subtitle = paste("Bonferroni correction on",BONFERRONI,"tests")) + 
		facet_wrap(~tissue)
	OUTNAME=paste0(SUBDIR,"/Figure_Welch.",OUTNAME_DATA,".pdf") ; pdf(OUTNAME,width=10,height=10) ; print(PLOT);dev.off()
}
#         mELT reference versus tissues ~ Only T cell subsets
OUTNAME_DATA="Welch_t_tests_Tcell"
if(F){
	GEX$mouse <- factor(GEX$orig.ident)
	GEX$tissue<-factor(GEX$tag_tissue_thr1.5)
	GEX$CellTypeFinal_Welch <- factor(GEX$CellTypeFinal)
	
	# # # # # # calculate percent of celltype per tissus for X mice - T-tests
	# n = count of cell per mouse / tissue / celltype
	# we want the percent of celltype in tissue / mouse
	grep("T",levels(GEX$CellTypeFinal_Welch),value=T)
	COUNTS <- GEX@meta.data %>% 
		filter( CellTypeFinal_Welch %in% c("act. T cell CD4+","naive T cell CD4+","T cell CD8+","Treg") ) %>%
		filter( tissue %ni% c("none") ) %>%
		group_by( mouse, tissue, CellTypeFinal_Welch ) %>% count() %>% 
		group_by( mouse,tissue ) %>% mutate( mouse.tissue.count=sum(n)) %>%
		group_by( tissue, CellTypeFinal_Welch ) %>% mutate( tissue.CellType.count=sum(n)) %>%
		mutate( mouse.count=length(unique(mouse)), percent.celltype.in.tissue.per.mouse = n / mouse.tissue.count ) %>% arrange(tissue)
	
	COUNTS=as.data.frame(COUNTS)
	OUTNAME=paste0(SUBDIR,"/table.",OUTNAME_DATA,".xlsx")
	write.xlsx(COUNTS,OUTNAME)
	#COUNTS<-read.xlsx( paste0(SUBDIR,"/table_counts.xlsx"))
	
	# # # # # # Calculate T-tests
	TISSUES<-c("CSF" ) # # "CSF", unique( COUNTS$tissue )
	CELLTYPES<-unique(COUNTS$CellTypeFinal_Welch)
	
	BONFERRONI=length(CELLTYPES) #compared to mELT 
	OUTNAME=paste0(SUBDIR,"/table.",OUTNAME_DATA,".txt")
	
	cat(c("target","tissue","celltype","treat_count","ref_count","ratio","Pvalue",paste0("Padj_bonferroni",BONFERRONI),"Tscore\n"),file=OUTNAME,append=F,sep="\t")
	for(CELLTYPE in CELLTYPES){
		for(TISSUE in TISSUES){
			D.mELT <- COUNTS[ COUNTS$tissue == "mELT" & COUNTS$CellTypeFinal_Welch == CELLTYPE, ]
			D.other <- COUNTS[ COUNTS$tissue == TISSUE & COUNTS$CellTypeFinal_Welch == CELLTYPE, ]
			REFERENCE = D.mELT$percent.celltype.in.tissue.per.mouse
			TREATMENT = D.other$percent.celltype.in.tissue.per.mouse
			RATIO = mean(TREATMENT) / mean(REFERENCE)
			# test normality, should be at least 3 obs.
			if( length(REFERENCE) >2 & length(TREATMENT) >2 ){
				SHA_TEST = shapiro.test( c(REFERENCE,TREATMENT) )#plot(density( c(REFERENCE,TREATMENT) )) ; #library("ggpubr") ; #ggqqplot( c(REFERENCE,TREATMENT) )
				if(SHA_TEST$p.value<0.05){
					PValue=1
					Statistic="not norm. distri."
				}else{
					TEST = t.test(x = REFERENCE,y = TREATMENT) # at least 3 observation AND normale distribution
					PValue=TEST$p.value
					Statistic=TEST$statistic
					rm(TEST)
				}
				rm(SHA_TEST)
			}else{
				PValue=1
				Statistic="< 3 obs."
			}
			
			#print(paste("mELT",TISSUE,CELLTYPE,PValue,Statistic,sep=";"))
			cat( c("mELT",TISSUE,CELLTYPE,length(TREATMENT),length(REFERENCE),RATIO,PValue,PValue*BONFERRONI, paste0(Statistic,"\n")),file=OUTNAME,append=T,sep="\t")
		}
	}
	
	TTEST=read.table(OUTNAME,h=T,sep="\t")
	TTEST$log10FC <- log10(TTEST$ratio)
	TTEST$log10Pval_m <- -log10(TTEST$Pvalue) # minus log10 : high values = low p-value
	write.xlsx(TTEST, paste0(OUTNAME,".xlsx"))
	
	LAB_SIZE=.4
	NUDGE_X=0 #-.1
	NUDGE_Y=.1
	BONFNAME=grep("Padj_bonferroni",names(TTEST),value=T)
	TTEST$celltype2 <- ifelse( TTEST[,BONFNAME] < 0.05, paste0(TTEST$celltype," *"),TTEST$celltype )
	
	PLOT=ggplot(TTEST) + 
		geom_hline(yintercept=-log10(0.05),linetype="dashed",col="red") + 
		geom_hline(yintercept=-log10(0.05/BONFERRONI),linetype="dashed",col="orange") +
		geom_vline(xintercept=0,linetype="dashed",col="black") + 
		geom_point(aes(log10FC,log10Pval_m,size=abs(log10FC)*log10Pval_m ,col=celltype),show.legend=F) + #to remove color: here
		geom_label_repel(aes(log10FC,log10Pval_m,label=celltype2,size=LAB_SIZE,col=celltype),min.segment.length=0,nudge_x=NUDGE_X,nudge_y=NUDGE_Y,show.legend=F) + 
		scale_x_continuous(limits = c(-1,1)) +
		theme_bw(base_size=22) + 
		ggtitle( label = paste("Celltypes % (reference: mELT)"),
				 subtitle = paste("Bonferroni correction on",BONFERRONI,"tests")) + 
		facet_wrap(~tissue)
	OUTNAME=paste0(SUBDIR,"/Figure_Welch.",OUTNAME_DATA,".pdf") ; pdf(OUTNAME,width=10,height=10) ; print(PLOT);dev.off()
}
#         mELT reference versus tissues ~ Only T cell subsets
OUTNAME_DATA="Welch_t_tests_Myeloid"
if(F){
	GEX$mouse <- factor(GEX$orig.ident)
	GEX$tissue<-factor(GEX$tag_tissue_thr1.5)
	GEX$CellTypeFinal_Welch <- factor(GEX$CellTypeFinal)
	
	# # # # # # calculate percent of celltype per tissus for X mice - T-tests
	# n = count of cell per mouse / tissue / celltype
	# we want the percent of celltype in tissue / mouse
	grep("T",levels(GEX$CellTypeFinal_Welch),value=T)
	COUNTS <- GEX@meta.data %>% 
		filter( CellTypeFinal_Welch %in% c("DC","Monocyte","Glial") ) %>%
		filter( tissue %ni% c("none") ) %>%
		group_by( mouse, tissue, CellTypeFinal_Welch ) %>% count() %>% 
		group_by( mouse,tissue ) %>% mutate( mouse.tissue.count=sum(n)) %>%
		group_by( tissue, CellTypeFinal_Welch ) %>% mutate( tissue.CellType.count=sum(n)) %>%
		mutate( mouse.count=length(unique(mouse)), percent.celltype.in.tissue.per.mouse = n / mouse.tissue.count ) %>% arrange(tissue)
	
	COUNTS=as.data.frame(COUNTS)
	OUTNAME=paste0(SUBDIR,"/table.",OUTNAME_DATA,".xlsx")
	write.xlsx(COUNTS,OUTNAME)
	#COUNTS<-read.xlsx( paste0(SUBDIR,"/table_counts.xlsx"))
	
	# # # # # # Calculate T-tests
	TISSUES<-c("CSF" ) # # "CSF", unique( COUNTS$tissue )
	CELLTYPES<-unique(COUNTS$CellTypeFinal_Welch)
	
	BONFERRONI=length(CELLTYPES) #compared to mELT 
	OUTNAME=paste0(SUBDIR,"/table.",OUTNAME_DATA,".txt")
	
	cat(c("target","tissue","celltype","treat_count","ref_count","ratio","Pvalue",paste0("Padj_bonferroni",BONFERRONI),"Tscore\n"),file=OUTNAME,append=F,sep="\t")
	for(CELLTYPE in CELLTYPES){
		for(TISSUE in TISSUES){
			D.mELT <- COUNTS[ COUNTS$tissue == "mELT" & COUNTS$CellTypeFinal_Welch == CELLTYPE, ]
			D.other <- COUNTS[ COUNTS$tissue == TISSUE & COUNTS$CellTypeFinal_Welch == CELLTYPE, ]
			REFERENCE = D.mELT$percent.celltype.in.tissue.per.mouse
			TREATMENT = D.other$percent.celltype.in.tissue.per.mouse
			RATIO = mean(TREATMENT) / mean(REFERENCE)
			# test normality, should be at least 3 obs.
			if( length(REFERENCE) >2 & length(TREATMENT) >2 ){
				SHA_TEST = shapiro.test( c(REFERENCE,TREATMENT) )#plot(density( c(REFERENCE,TREATMENT) )) ; #library("ggpubr") ; #ggqqplot( c(REFERENCE,TREATMENT) )
				if(SHA_TEST$p.value<0.05){
					PValue=1
					Statistic="not norm. distri."
				}else{
					TEST = t.test(x = REFERENCE,y = TREATMENT) # at least 3 observation AND normale distribution
					PValue=TEST$p.value
					Statistic=TEST$statistic
					rm(TEST)
				}
				rm(SHA_TEST)
			}else{
				PValue=1
				Statistic="< 3 obs."
			}
			
			#print(paste("mELT",TISSUE,CELLTYPE,PValue,Statistic,sep=";"))
			cat( c("mELT",TISSUE,CELLTYPE,length(TREATMENT),length(REFERENCE),RATIO,PValue,PValue*BONFERRONI, paste0(Statistic,"\n")),file=OUTNAME,append=T,sep="\t")
		}
	}
	
	TTEST=read.table(OUTNAME,h=T,sep="\t")
	TTEST$log10FC <- log10(TTEST$ratio)
	TTEST$log10Pval_m <- -log10(TTEST$Pvalue) # minus log10 : high values = low p-value
	write.xlsx(TTEST, paste0(OUTNAME,".xlsx"))
	
	LAB_SIZE=.4
	NUDGE_X=0 #-.1
	NUDGE_Y=.1
	BONFNAME=grep("Padj_bonferroni",names(TTEST),value=T)
	TTEST$celltype2 <- ifelse( TTEST[,BONFNAME] < 0.05, paste0(TTEST$celltype," *"),TTEST$celltype )
	
	PLOT=ggplot(TTEST) + 
		geom_hline(yintercept=-log10(0.05),linetype="dashed",col="red") + 
		geom_hline(yintercept=-log10(0.05/BONFERRONI),linetype="dashed",col="orange") +
		geom_vline(xintercept=0,linetype="dashed",col="black") + 
		geom_point(aes(log10FC,log10Pval_m,size=abs(log10FC)*log10Pval_m ,col=celltype),show.legend=F) + #to remove color: here
		geom_label_repel(aes(log10FC,log10Pval_m,label=celltype2,size=LAB_SIZE,col=celltype),min.segment.length=0,nudge_x=NUDGE_X,nudge_y=NUDGE_Y,show.legend=F) + 
		scale_x_continuous(limits = c(-1,1)) +
		theme_bw(base_size=22) + 
		ggtitle( label = paste("Celltypes % (reference: mELT)"),
				 subtitle = paste("Bonferroni correction on",BONFERRONI,"tests")) + 
		facet_wrap(~tissue)
	OUTNAME=paste0(SUBDIR,"/Figure_Welch.",OUTNAME_DATA,".pdf") ; pdf(OUTNAME,width=10,height=10) ; print(PLOT);dev.off()
}

rm(PLOT,TTEST,BONFNAME,NUDGE_X,NUDGE_Y,LAB_SIZE,OUTNAME)


#https://metascape.org



#### special klaus talk ####
SUBDIR="20.Plots.Klaus-Hertie-talk";dir.create(SUBDIR)
if(F){
	options(max.print=200)
	#names(GEX@meta.data)
	
	
	# 1 - no filter (all) with CD20 versus isotype
	GEX<-readRDS("06.cc.Rds")
	# remove the "rm clusters
	GEX<-subset(GEX, CellTypeFinal %ni% c("rm-1","rm-2","rm-3","rm-4"))
	table(GEX$treatment)
	NAMES_CLU=c("CellTypeFinal") ; ORDER=NULL # dir(SUBDIR)
	SPLIT="treatment" ; WIDTH=14;HEIGHT=10; NCOL=2; LABSIZE=4
	for ( NAMEIDENT in NAMES_CLU){
		print(paste(NAMEIDENT,"split:",SPLIT))
		Idents(GEX)<-NAMEIDENT
		PLOT=DimPlot(GEX, reduction="umap", label=TRUE,split.by=SPLIT,order=ORDER,pt.size=.5,
					 label.size=LABSIZE,repel=T,ncol=NCOL) + 
			theme_minimal(base_size=20)+
			theme(legend.position='top')+
			ggtitle("All cells") + coord_fixed()
		OUTNAME=paste0(SUBDIR,"/all cell - CD20 and Isotype.umap.pdf");pdf(OUTNAME,width=WIDTH,height=HEIGHT,useDingbats=F);print(PLOT);dev.off()
	}
	TABLE=table(GEX$CellTypeFinal,GEX$treatment)
	write.xlsx(TABLE,file = paste0(SUBDIR,"/all cell - CD20 and Isotype.table.xlsx"))
	
	
	# For the following I assume Klaus wants the HTO + VDJ filtered cells
	
	# 2 ) CSF isotype vs. mELT isotype
	GEX<-readRDS("07.final.Rds")
	GEX<-subset(GEX, CellTypeFinal %ni% c("rm-1","rm-2","rm-3","rm-4"))
	table(GEX$split)
	GEX<-subset(GEX,split %in% c("Isotype CSF","Isotype mELT"))
	
	NAMES_CLU=c("CellTypeFinal") ; ORDER=NULL # dir(SUBDIR)
	SPLIT="split" ; WIDTH=14;HEIGHT=10; NCOL=2; LABSIZE=4
	for ( NAMEIDENT in NAMES_CLU){
		print(paste(NAMEIDENT,"split:",SPLIT))
		Idents(GEX)<-NAMEIDENT
		PLOT=DimPlot(GEX, reduction="umap", label=TRUE,split.by=SPLIT,order=ORDER,pt.size=.5,
					 label.size=LABSIZE,repel=T,ncol=NCOL) + 
			theme_minimal(base_size=20)+
			theme(legend.position='top')+
			ggtitle("HTO+VDJ-filter") + coord_fixed()
		OUTNAME=paste0(SUBDIR,"/filter HTO+VDJ - CSF isotype vs. mELT isotype.umap.pdf");pdf(OUTNAME,width=WIDTH,height=HEIGHT,useDingbats=F);print(PLOT);dev.off()
	}
	TABLE=table(GEX$CellTypeFinal,GEX$split)
	write.xlsx(TABLE,file = paste0(SUBDIR,"/filter HTO+VDJ - CSF isotype vs. mELT isotype.table.xlsx"))
	
	
	
	# 3 ) CSF CD20 vs. CSF isotype
	GEX<-readRDS("07.final.Rds")
	GEX<-subset(GEX, CellTypeFinal %ni% c("rm-1","rm-2","rm-3","rm-4"))
	table(GEX$split)
	GEX<-subset(GEX,split %in% c("CD20 CSF","Isotype CSF"))
	
	NAMES_CLU=c("CellTypeFinal") ; ORDER=NULL # dir(SUBDIR)
	SPLIT="split" ; WIDTH=14;HEIGHT=10; NCOL=2; LABSIZE=4
	for ( NAMEIDENT in NAMES_CLU){
		print(paste(NAMEIDENT,"split:",SPLIT))
		Idents(GEX)<-NAMEIDENT
		PLOT=DimPlot(GEX, reduction="umap", label=TRUE,split.by=SPLIT,order=ORDER,pt.size=.5,
					 label.size=LABSIZE,repel=T,ncol=NCOL) + 
			theme_minimal(base_size=20)+
			theme(legend.position='top')+
			ggtitle("HTO+VDJ-filter") + coord_fixed()
		OUTNAME=paste0(SUBDIR,"/filter HTO+VDJ - CD20 CSF vs. Isotype CSF.umap.pdf");pdf(OUTNAME,width=WIDTH,height=HEIGHT,useDingbats=F);print(PLOT);dev.off()
	}
	TABLE=table(GEX$CellTypeFinal,GEX$split)
	write.xlsx(TABLE,file = paste0(SUBDIR,"/filter HTO+VDJ - CD20 CSF vs. Isotype CSF.table.xlsx"))
	
	
	
}







#### check some genes ####
if(F){
	#GEX$split <- paste(GEX$treatment,GEX$tag_tissue_thr1.5)

	VlnPlot( subset(GEX, tag_tissue_thr1.5 %in% c("CSF","mELT") & 
						CellTypeFinal %ni% grep("^rm-",unique(GEX$CellTypeFinal),value=T) ),
			 features="Mbp",combine=T,
			 group.by = "CellTypeFinal",split.by="tag_tissue_thr1.5")
	
	
	GENE<-"Trdc"
	GENE_bool<-GetAssayData(object=GEX,assay="RNA",slot="data")[GENE,] > 0
	table(GENE_bool,GEX$RNA_snn_res.2)
	
	#grep("Tcrg-C",rownames(GEX),value=T)
	
	VlnPlot( subset(GEX, tag_tissue_thr1.5 %in% c("CSF","mELT") & 
						RNA_snn_res.2 %in% c("13","8","4","6","16","14","0","10","23")),
			 features=c("Tcrg-C1","Tcrg-C2","Tcrg-C4"),combine=T,
			 group.by = "RNA_snn_res.2")
	
	
}
##### Background gene list ####
SUBDIR="16.background";dir.create(SUBDIR)
if(T){
	BACKGROUND_GENES=rownames(GEX)
	OUTNAME=paste0(SUBDIR,"/background_all_genes.txt")
	write.table(BACKGROUND_GENES,OUTNAME,quote=F,sep="\t",row.names=F,col.names=F)
	
	CELLTYPES<-unique(GEX$CellTypeFinal)
	CELLTYPES<-CELLTYPES[ ! CELLTYPES %in% grep("^rm-",CELLTYPES,value=T)]
	NAMEIDENT<-"CellTypeFinal"
	Idents(GEX)<-NAMEIDENT
	
	avgexp <- AverageExpression(GEX,verbose=T,assays="RNA")
	avgexp <- as.data.frame(avgexp$RNA)
	

	#gene taken if at least expression in 10 % of the cells in the cluster
	for( MinPercent in c(0,1,5,10)){
		print(paste0("-------------------gt",MinPercent,"-------------------"))
		X=sapply( CELLTYPES , function(identName){
			CELLS_index=which( GEX@meta.data[ , NAMEIDENT] == identName )
			THRESHOLD= ceiling(length(CELLS_index)/100*MinPercent)
			RSUM = data.frame( x=rowSums(GEX@assays$RNA@counts[,CELLS_index] > 0) )
			GENELIST=rownames(RSUM[ RSUM$x > THRESHOLD , , drop=F])
			print(paste(identName,":",length(GENELIST),"genes"))
			OUTNAME=paste0(SUBDIR,"/background_gt",MinPercent,"percent_",identName,".txt")
			write.table(x=GENELIST,file=OUTNAME,quote=F,row.names=F,col.names=F)
		} )
	}
	
	system( paste0("mkdir -p ",SUBDIR,"/gt0 ; mv ",SUBDIR,"/*gt0*.txt ",SUBDIR,"/gt0"))
	system( paste0("mkdir -p ",SUBDIR,"/gt10 ; mv ",SUBDIR,"/*gt10*.txt ",SUBDIR,"/gt10"))
	system( paste0("mkdir -p ",SUBDIR,"/gt1 ; mv ",SUBDIR,"/*gt1*.txt ",SUBDIR,"/gt1"))
	system( paste0("mkdir -p ",SUBDIR,"/gt5 ; mv ",SUBDIR,"/*gt5*.txt ",SUBDIR,"/gt5"))
	
}

#### add background to gene list - metascape #### 
if(F){
	# Load background genes
	BACKfiles<-list.files(paste0( "16.background", "/", c("gt0","gt1","gt5","gt10")),full.names=T)
	BACKG<-lapply(BACKfiles, function(x){ TXT<-read.table(x) ; return(TXT$V1) } )
	names(BACKG)<-basename(BACKfiles)
	
	#DEGfiles<-list.files(paste0( "16.DEG.mELT_vs_CSF"),pattern=".xlsx",full.names=T)
	#DEGfiles<-list.files(paste0( "17.DEG.allCD20_vs_allIsotype/"),pattern=".xlsx",full.names=T)
	DEGfiles<-list.files(paste0( "17.DEG.CD20_vs_Isotype/"),pattern=".xlsx",full.names=T)
	
	DEGfiles<-grep("gt[0-9]+",x = DEGfiles, invert=T,value=T)# in case we already have some there
	DEGG<-lapply(DEGfiles, read.xlsx )
	names(DEGG)<-basename(DEGfiles)
	
	lapply( DEGfiles ,function( FILENAME ){
		CELLTYPE<-gsub( ".xlsx$", "", basename(FILENAME) )
		CELLTYPE<-gsub("Markers_","", strsplit(CELLTYPE,"-",fixed=T)[[1]][1])
		
		GenesList <- DEGG[[ basename(FILENAME) ]] %>% filter(p_val_adj < 0.05) %>% pull(gene)
		
		for( SUB in c("gt0","gt1","gt5","gt10")){
			Background <- BACKG[[ paste0("background_",SUB,"percent_",CELLTYPE,".txt") ]]
			DF<-data.frame( NA_padding( list( GenesList=GenesList, Background=Background )  ) )
			names(DF)<-c(CELLTYPE,"_BACKGROUND")
			OUTNAME<-paste0( gsub( ".xlsx$", "", FILENAME ), "_", SUB, ".xlsx" )
			write.xlsx( DF, OUTNAME)
		}
		
	})
	rm(BACKfiles,BACKG,DEGfiles,DEGG)
}
# save all DEG for each cell type: 
if(F){
	
	# pick one
	DEGfiles<-list.files(paste0( "16.DEG.mELT_vs_CSF"),pattern=".xlsx",full.names=T) ; OUTNAME="16.DEG.list.mELT_vs_CSF.xlsx"
	DEGfiles<-list.files(paste0( "17.DEG.allCD20_vs_allIsotype/"),pattern=".xlsx",full.names=T) ; OUTNAME="17.DEG.all.allCD20_vs_allIsotype.xlsx"
	DEGfiles<-list.files(paste0( "17.DEG.CD20_vs_Isotype/"),pattern=".xlsx",full.names=T); OUTNAME="17.DEG.list.CD20_vs_Isotype.xlsx"
	
	#filter is already done
	DEGfiles<-grep("gt[0-9]+",x = DEGfiles, invert=T,value=T)# in case we already have some there
	DEGG<-lapply(DEGfiles, function(x){ DF<-read.xlsx(x) ; DF <- DF %>% filter(p_val_adj < 0.05); return(DF$gene) })
	CELLTYPE<-gsub("Markers_","",  gsub( ".xlsx$", "", basename(DEGfiles) ) )
	CELLTYPE<-unlist(lapply( strsplit(CELLTYPE,"_vs_",fixed=T) , "[[", 1))
	names(DEGG)<-1:length(DEGG)
	DF<-as.data.frame( NA_padding( DEGG ) )
	names(DF)<-CELLTYPE
	write.xlsx(DF,OUTNAME)
}


#### Metascape post-analysis - dotplot - mainly manual process ####
if(F){
	
	
	# load all the GO_FINAL -> SAVE
	if(F){
		# do the following code for each of the directory
		setwd("/home/ga94rac/LRZ Sync+Share/AG_Lehmann-Horn/Gildas_storage/2021.02.16.CITEseq/experiment_CD20/19.GO_Enrichment_Metascape.auto/Mus_to_Mus/CD20_allCD20_vs_allIsotype")
		
		setwd("/home/ga94rac/LRZ Sync+Share/AG_Lehmann-Horn/Gildas_storage/2021.02.16.CITEseq/experiment_CD20/19.GO_Enrichment_Metascape.auto/Mus_to_Mus/CD20_CD20_mELT_versus_Isotype_mELT/")
		
		setwd("/home/ga94rac/LRZ Sync+Share/AG_Lehmann-Horn/Gildas_storage/2021.02.16.CITEseq/experiment_CD20/19.GO_Enrichment_Metascape.auto/Mus_to_Mus/CD20_CD20_vs_Isotype")
		setwd("/home/ga94rac/LRZ Sync+Share/AG_Lehmann-Horn/Gildas_storage/2021.02.16.CITEseq/experiment_CD20/19.GO_Enrichment_Metascape.auto/Mus_to_Mus/CD20_mELT_vs_CSF")
		
		LIST<-list.files(".",pattern="_FINAL_GO.csv",full.names=T,recursive=T)
		#LIST2<-grep(pattern="mELT",x=LIST,value=T)
		LIST2<-LIST
		FILES<-lapply(LIST2, function(x){df<-read.csv(x) ; df$source=x ; return(df) })
		names(FILES) <- gsub("Markers_","",basename(dirname(dirname(LIST2))))
		# remove column with specific names
		FILES2<-lapply(FILES, function(x){  x$tissue=gsub("X_MEMBER_","", names(x)[1]) ; names(x)[1:2]<-c("X_MEMBER","X_LogP");return(x)  })
		FILES2<-lapply(FILES2, function(x){ x$percent = x$X.GeneInGOAndHitList / x$X.GeneInGO * 100 ; return(x) })
		#SIGNIF_SUBSET<-lapply(FILES2, function(x){ subset( x , x$Log.q.value. < LOG_Q_VAL_LIMIT ) })
		
		
		#lapply(FILES2, ncol)
		#NORM=names(FILES2$`act._T_cell_CD4+-Isotype_mELT_vs_act._T_cell_CD4+-Isotype_CSF_gt0`)
		#WEIRD=names(FILES2$`Macrophage-CD20_mELT_vs_Macrophage-CD20_CSF_gt0`)
		#NORM [ ! NORM %in% WEIRD ]
		# I had to add the following in the file:
		# Macrophage-CD20_mELT_vs_Macrophage-CD20_CSF_gt0
		#"BestLogPInGroup"       "BestEnrichmentInGroup" "URL" 
		
		DF <- do.call(rbind.data.frame, FILES2)
		rm(FILES,FILES2,LIST,LIST2)
		
		#Mus_to_Mus/GO_MusMus_CD20_allCD20_vs_allIsotype
		if(T){
			dir.create("../GO_MusMus_CD20_allCD20_vs_allIsotype")
			openxlsx::write.xlsx(DF,file = paste0("../GO_MusMus_CD20_allCD20_vs_allIsotype/GO_MusMus_CD20_allCD20_vs_allIsotype.xlsx"))
			saveRDS(DF, paste0("../GO_MusMus_CD20_allCD20_vs_allIsotype/GO_MusMus_CD20_allCD20_vs_allIsotype.Rds"))
		}
		
		#Mus_to_Mus/GO_MusMus_CD20_CD20_mELT_versus_Isotype_mELT
		if(T){
			dir.create("../GO_MusMus_CD20_CD20_mELT_versus_Isotype_mELT")
			openxlsx::write.xlsx(DF,file = paste0("../GO_MusMus_CD20_CD20_mELT_versus_Isotype_mELT/GO_MusMus_CD20_CD20_mELT_versus_Isotype_mELT.xlsx"))
			saveRDS(DF, paste0("../GO_MusMus_CD20_CD20_mELT_versus_Isotype_mELT/GO_MusMus_CD20_CD20_mELT_versus_Isotype_mELT.Rds"))
		}
		
		#Mus_to_Mus/GO_MusMus_CD20_CD20_vs_Isotype
		if(T){
			dir.create("../GO_MusMus_CD20_CD20_vs_Isotype")
			openxlsx::write.xlsx(DF,file = paste0("../GO_MusMus_CD20_CD20_vs_Isotype/GO_MusMus_CD20_CD20_vs_Isotype.xlsx"))
			saveRDS(DF, paste0("../GO_MusMus_CD20_CD20_vs_Isotype/GO_MusMus_CD20_CD20_vs_Isotype.Rds"))
		}
		
		#Mus_to_Mus/GO_MusMus_CD20_mELT_vs_CSF
		if(T){
			dir.create("../GO_MusMus_CD20_mELT_vs_CSF")
			openxlsx::write.xlsx(DF,file = paste0("../GO_MusMus_CD20_mELT_vs_CSF/GO_MusMus_CD20_mELT_vs_CSF.xlsx"))
			saveRDS(DF, paste0("../GO_MusMus_CD20_mELT_vs_CSF/GO_MusMus_CD20_mELT_vs_CSF.Rds"))
		}
		
		
	}
	
	setwd("/home/ga94rac/LRZ Sync+Share/AG_Lehmann-Horn/Gildas_storage/2021.02.16.CITEseq/experiment_CD20/19.GO_Enrichment_Metascape.auto/Mus_to_Mus/")
	
	
	# process using the top for each category
	LOG_Q_VAL_LIMIT=log(0.05)
	if(F){
		# choose one
		TOP=30  ; WIDTH=12 ; HEIGHT=8
		
		# useless here, since max 16 signif diff
		#TOP=100 ; WIDTH=12 ; HEIGHT=12
		
		
		for( index_GO_CAT in 1:3 ){
			
			if( index_GO_CAT == 1) { GO_CATEGORY="GO Biological Processes" ; GO_CAT_NAME="GO.Bio.Pro" }
			if( index_GO_CAT == 2) { GO_CATEGORY="Reactome Gene Sets" ; GO_CAT_NAME="GO.Reacto" }
			if( index_GO_CAT == 3) { GO_CATEGORY="CORUM" ; GO_CAT_NAME="GO.CORUM" }
			
			for( index_DF in 1:4 ){
				if( index_DF == 1) {
					DF <- readRDS("GO_MusMus_CD20_allCD20_vs_allIsotype/GO_MusMus_CD20_allCD20_vs_allIsotype.Rds") ; OUTNAME=paste0("GO_allTogether_CD20_vs_Isotype.",GO_CATEGORY,".top",TOP,".pdf")
					TITLENAME=paste("Top",TOP,"pathways in","CD20 versus Isotype",paste0("\n(log q-FDR < ",round(LOG_Q_VAL_LIMIT,2),")"))
				}
				
				if( index_DF == 2) {
					DF <- readRDS("GO_MusMus_CD20_CD20_mELT_versus_Isotype_mELT/GO_MusMus_CD20_CD20_mELT_versus_Isotype_mELT.Rds") ; OUTNAME=paste0("GO_CD20_mELT_vs_Isotype_mELT.",GO_CATEGORY,".top",TOP,".pdf")
					TITLENAME=paste("Top",TOP,"pathways in","CD20 mELT versus Isotype mELT",paste0("\n(log q-FDR < ",round(LOG_Q_VAL_LIMIT,2),")"))
				}
				
				if( index_DF == 3) {
					DF <- readRDS("GO_MusMus_CD20_mELT_vs_CSF/GO_MusMus_CD20_mELT_vs_CSF.Rds") ; OUTNAME=paste0("GO_MusMus_mELT_vs_CSF.",GO_CATEGORY,".top",TOP,".pdf")
					TITLENAME=paste("Top",TOP, "pathways in","mELT versus CSF",paste0("\n(log q-FDR < ",round(LOG_Q_VAL_LIMIT,2),")"))
				}
				if( index_DF ==4) {
					DF <- readRDS("GO_MusMus_CD20_CD20_vs_Isotype/GO_MusMus_CD20_CD20_vs_Isotype.Rds") ; OUTNAME=paste0("GO_MusMus_CD20_vs_Isotype.",GO_CATEGORY,".top",TOP,".pdf")
					TITLENAME=paste("Top",TOP, "pathways in","CD20 versus Isotype",paste0("\n(log q-FDR < ",round(LOG_Q_VAL_LIMIT,2),")"))
				}
				
				# select the top category
				DF <- DF %>% 
					filter( Category == GO_CATEGORY) %>% filter( FirstInGroupByEnrichment == 1) %>% 
					filter( Log.q.value. < LOG_Q_VAL_LIMIT) %>% top_n( n = TOP, wt = desc( Log.q.value. ) ) # carefull with those negative logs
				
				# if we want to overlap
				if(F){
					WIDTH=12 ; HEIGHT=6
					D1 = SIGN_1$`act._T_cell_CD4+-mELT_vs_act._T_cell_CD4+-LN_lumbar_gt0`
					D2 = SIGN_2$`act._T_cell_CD4+-mELT_vs_act._T_cell_CD4+-LN_inguinal_gt0`
					D3 = SIGN_3$`act._T_cell_CD4+-mELT_vs_act._T_cell_CD4+-Spleen_gt0`
					LI1 = D1[ D1$Category == GO_CATEGORY, "Description" ]
					LI2 = D2[ D2$Category == GO_CATEGORY, "Description" ]
					LI3 = D3[ D3$Category == GO_CATEGORY, "Description" ]
					#paste(LI1,sep = "", collapse = ","); 	paste(LI2,sep = "", collapse = ",") ;	paste(LI3,sep = "", collapse = ",")
					OVERLAP_GO=Reduce( intersect, list(LI1,LI2,LI3) ) ; rm(LI1,LI2,LI3, D1,D2,D3)
					GO_selection_order = rev(gtools::mixedsort(OVERLAP_GO)) ; rm(OVERLAP_GO)
				}
				
				# CELLTYPES_order=c("naive T cell CD4+", "act. T cell CD4+", "T cell CD8+", "gd T cell", "Treg" ,
				# 				  "NK + NK T cell",
				# 				  "fol. zone B cell","mar. zone B cell", "ger. center B cell","Plasma cell",
				# 				  "DC", "Granulocyte", "Macrophage", "Monocyte",
				# 				  "Erythroid cell")
				#CELLTYPES_order = unique( c(DF1$tissue,DF2$tissue,DF2$tissue))
				#CELLTYPES_order = unique( DF$tissue)
				
				CELLTYPES_order = c("naive T cell CD4+","act. T cell CD4+","Treg", "T cell CD8+","NK + NK T cell",
									"B cell", "Plasma cell",
									"DC","Monocyte","Macrophage","Granulocyte",
									"Oligodendrocytes","Glial")
				#Erythroid cell
				
				OVERLAP_GO=unique(DF$Description)
				GO_selection_order = rev(gtools::mixedsort(OVERLAP_GO)) ; rm(OVERLAP_GO)
				
				# plot
				if(T){
					DF = subset( DF, DF$Category %in% GO_CATEGORY)
					
					
					#rename stuff
					DF$Cell_type <- rownames(DF)
					DF$Cell_type <- gsub("_gt0[.][[:digit:]]+$", "", DF$Cell_type)
					DF$Cell_type <- gsub("[-].*", "", DF$Cell_type)
					DF$Cell_type <- gsub("_", " ", DF$Cell_type)
					DF$Cell_type <- factor(DF$Cell_type)
					
					DF$Pathway <- DF$Description
					DF$Genes <- DF$X.GeneInGOAndHitList
					DF$log_qFDR <- - DF$Log.q.value.
					
					# add empty stuff for the celltypes that have no matches
					DF.plot <- subset(DF, DF$Description %in% GO_selection_order)
					
					MISSING=CELLTYPES_order [ ! CELLTYPES_order %in% DF.plot$Cell_type ]
					if( length(MISSING)>0 ){
						DF.plot.tmp=data.frame(Cell_type=unlist(lapply(MISSING, rep, length(GO_selection_order))),Pathway=rep(GO_selection_order,length(MISSING)))
						DF.plot <- plyr::rbind.fill( DF.plot, DF.plot.tmp)
						rm(DF.plot.tmp)
					}
					
					# play with the levels
					
					DF.plot$Cell_type=factor(DF.plot$Cell_type, levels=CELLTYPES_order)
					DF.plot$Pathway = factor(DF.plot$Pathway, levels=GO_selection_order)
					
					PLOT<-ggplot( DF.plot )+
						geom_point(aes(x=Cell_type, y=Pathway, colour=percent, size=log_qFDR))+ # shape=19 #, shape=Category
						scale_colour_gradient(high='darkred',low='gray80')+ #, limits=c(0,100)
						#scale_size_continuous(breaks=c(10, 50, 80), limits=c(0,100)) +
						ylab('')+xlab('')+ggtitle( TITLENAME )+theme_bw()+theme(axis.text.x = element_text(angle=45,hjust=1)) + #coord_fixed(ratio=0.7) + 
						labs(colour="Genes %",size="- log( q-FDR )")
					
					ggsave(OUTNAME, plot=PLOT,width=WIDTH,height=HEIGHT)
					rm(PLOT,OUTNAME,TITLENAME,DF.plot,DF,MISSING)
				}
			}
		}
		
	}
	
	# process using the selection from Jolien
	if(F){
		SELECTION_J=c("Eukaryotic Translation Initiation","Glycolysis","ATP metabolic process","antigen processing and presentation of exogenous antigen",
					  "adaptive immune response","inflammatory response","leukocyte cell-cell adhesion","regulation of leukocyte activation",
					  "leukocyte proliferation","leukocyte differentiation","regulation of cytokine production","T cell activation","T-helper cell differentiation",
					  "T cell selection","T cell receptor signaling pathway","T cell mediated cytotoxicity","lymph node development","response to interferon-gamma",
					  "tumor necrosis factor production","interleukin-6 production","apoptotic signaling pathway","myeloid leukocyte activation",
					  "myeloid leukocyte migration","Neutrophil degranulation","chemotaxis","phagocytosis","positive regulation of MAPK cascade",
					  "cellular response to oxidative stress","cell killing","Classical antibody-mediated complement activation")
		length(SELECTION_J)
		
		CELLTYPES_order=c("naive T cell CD4+", "act. T cell CD4+", "T cell CD8+", "gd T cell", "Treg" ,
						  "NK + NK T cell",
						  "fol. zone B cell","mar. zone B cell", "ger. center B cell","Plasma cell",
						  "DC", "Granulocyte", "Macrophage", "Monocyte",
						  "Erythroid cell")
		# split the categories
		WIDTH=12 ; HEIGHT=8
		if(F){
			for( index_GO_CAT in 1:3 ){
				
				if( index_GO_CAT == 1) { GO_CATEGORY="GO Biological Processes" ; GO_CAT_NAME="GO.Bio.Pro" }
				if( index_GO_CAT == 2) { GO_CATEGORY="Reactome Gene Sets" ; GO_CAT_NAME="GO.Reacto" }
				if( index_GO_CAT == 3) { GO_CATEGORY="CORUM" ; GO_CAT_NAME="GO.CORUM" }
				
				for( index_DF in 1:3 ){
					if( index_DF == 1) {
						DF <- readRDS("GO_MusMus_LN_lum/GO_MusMus_LN_lum.Rds") ; OUTNAME=paste0("GO_MusMus_LN_lum.",GO_CATEGORY,".pdf")
						TITLENAME=paste("mELT versus LN lumbar")
					}
					if( index_DF == 2) {
						DF <- readRDS("GO_MusMus_LN_ing/GO_MusMus_LN_ing.Rds") ; OUTNAME=paste0("GO_MusMus_LN_ing.",GO_CATEGORY,".pdf")
						TITLENAME=paste("mELT versus LN inguinal")
					}
					if( index_DF == 3) {
						DF <- readRDS("GO_MusMus_Spleen/GO_MusMus_Spleen.Rds") ; OUTNAME=paste0("GO_MusMus_Spleen.",GO_CATEGORY,".pdf")
						TITLENAME=paste("mELT versus Spleen")
					}
					#rename stuff
					DF$Cell_type <- rownames(DF)
					DF$Cell_type <- gsub("_gt0[.][[:digit:]]+$", "", DF$Cell_type)
					DF$Cell_type <- gsub("[-].*", "", DF$Cell_type)
					DF$Cell_type <- gsub("_", " ", DF$Cell_type)
					DF$Cell_type <- factor(DF$Cell_type)
					
					
					GO_selection_order = rev(gtools::mixedsort(SELECTION_J))
					
					# plot - split data origin
					# split GO_Category
					
					DF = subset( DF, DF$Category %in% GO_CATEGORY)
					
					
					DF$Pathway <- DF$Description
					DF$Genes <- DF$X.GeneInGOAndHitList
					DF$log_qFDR <- - DF$Log.q.value.
					
					# add empty stuff for the celltypes that have no matches
					DF.plot <- subset(DF, DF$Description %in% GO_selection_order)
					
					MISSING=CELLTYPES_order [ ! CELLTYPES_order %in% DF.plot$Cell_type ]
					if( length(MISSING)>0 ){
						DF.plot.tmp=data.frame(Cell_type=unlist(lapply(MISSING, rep, length(GO_selection_order))),Pathway=rep(GO_selection_order,length(MISSING)))
						DF.plot <- plyr::rbind.fill( DF.plot, DF.plot.tmp)
						rm(DF.plot.tmp)
					}
					
					# play with the levels
					
					DF.plot$Cell_type=factor(DF.plot$Cell_type, levels=CELLTYPES_order)
					DF.plot$Pathway = factor(DF.plot$Pathway, levels=GO_selection_order)
					
					PLOT<-ggplot( DF.plot )+
						geom_point(aes(x=Cell_type, y=Pathway, colour=percent, size=log_qFDR))+ # colour=log_qFDR, size=percent
						scale_colour_gradient(high='darkred',low='gray80')+ #, limits=c(0,100)
						#scale_size_continuous(breaks=c(10, 50, 80), limits=c(0,100)) +
						ylab('')+xlab('')+ggtitle( TITLENAME )+theme_bw()+theme(axis.text.x = element_text(angle=45,hjust=1)) + #coord_fixed(ratio=0.7) + 
						labs(colour="Genes %",size="- log( q-FDR )")
					
					ggsave(OUTNAME, plot=PLOT,width=WIDTH,height=HEIGHT)
					rm(PLOT,OUTNAME,TITLENAME,DF.plot,DF,MISSING)
				}
			}
		}
		
		#do not split the categories
		WIDTH=12 ; HEIGHT=12
		if(F){
			for( index_DF in 1:3 ){
				if( index_DF == 1) {
					DF <- readRDS("GO_MusMus_LN_lum/GO_MusMus_LN_lum.Rds") ; OUTNAME=paste0("GO_MusMus_LN_lum.pdf")
					TITLENAME=paste("mELT versus LN lumbar")
				}
				if( index_DF == 2) {
					DF <- readRDS("GO_MusMus_LN_ing/GO_MusMus_LN_ing.Rds") ; OUTNAME=paste0("GO_MusMus_LN_ing.pdf")
					TITLENAME=paste("mELT versus LN inguinal")
				}
				if( index_DF == 3) {
					DF <- readRDS("GO_MusMus_Spleen/GO_MusMus_Spleen.Rds") ; OUTNAME=paste0("GO_MusMus_Spleen.pdf")
					TITLENAME=paste("mELT versus Spleen")
				}
				#rename stuff
				DF$Cell_type <- rownames(DF)
				DF$Cell_type <- gsub("_gt0[.][[:digit:]]+$", "", DF$Cell_type)
				DF$Cell_type <- gsub("[-].*", "", DF$Cell_type)
				DF$Cell_type <- gsub("_", " ", DF$Cell_type)
				DF$Cell_type <- factor(DF$Cell_type)
				
				GO_selection_order = rev(gtools::mixedsort(SELECTION_J))
				
				# plot - split data origin
				if(T){
					DF$Pathway <- DF$Description
					DF$Genes <- DF$X.GeneInGOAndHitList
					DF$log_qFDR <- - DF$Log.q.value.
					
					# add empty stuff for the celltypes that have no matches
					DF.plot <- subset(DF, DF$Description %in% GO_selection_order)
					
					MISSING=CELLTYPES_order [ ! CELLTYPES_order %in% DF.plot$Cell_type ]
					if( length(MISSING)>0 ){
						DF.plot.tmp=data.frame(Cell_type=unlist(lapply(MISSING, rep, length(GO_selection_order))),Pathway=rep(GO_selection_order,length(MISSING)))
						DF.plot <- plyr::rbind.fill( DF.plot, DF.plot.tmp)
						rm(DF.plot.tmp)
					}
					
					# play with the levels
					
					DF.plot$Cell_type=factor(DF.plot$Cell_type, levels=CELLTYPES_order)
					DF.plot$Pathway = factor(DF.plot$Pathway, levels=GO_selection_order)
					
					#paste the origine of the pathway:
					DF.plot$Pathway_origin=factor(DF.plot$Category)
					levels(DF.plot$Pathway_origin) [ levels(DF.plot$Pathway_origin) == "GO Biological Processes" ] <- "GO-BP"
					levels(DF.plot$Pathway_origin) [ levels(DF.plot$Pathway_origin) == "Reactome Gene Sets" ] <- "Reactome"
					
					DF.plot <- subset( DF.plot, ! is.na(DF.plot$Pathway_origin) )
					
					DF.plot$Pathway <- paste(DF.plot$Pathway, "-", DF.plot$Pathway_origin)
					
					
					PLOT<-ggplot( DF.plot )+
						geom_point(aes(x=Cell_type, y=Pathway, colour=percent, size=log_qFDR, shape=Category))+ # colour=log_qFDR, size=percent, 
						scale_colour_gradient(high='darkred',low='gray80')+ #, limits=c(0,100)
						#scale_size_continuous(breaks=c(10, 50, 80), limits=c(0,100)) +
						ylab('')+xlab('')+ggtitle( TITLENAME )+theme_bw()+theme(axis.text.x = element_text(angle=45,hjust=1)) + #coord_fixed(ratio=0.7) + 
						labs(colour="Genes %",size="- log( q-FDR )")
					
					ggsave(OUTNAME, plot=PLOT,width=WIDTH,height=HEIGHT)
					rm(PLOT,OUTNAME,TITLENAME,DF.plot,DF,MISSING)
				}
			}
		}
	}
}

#### Cellchat ####
if(F){
	#source("https://raw.githubusercontent.com/GildasLepennetier/SCr/main/CChat.R")
	library(CellChat)
	GEX$Condition=GEX$tag_tissue_thr1.5
	# subset to avoid including unwanted comparisons & noise
	GEX = subset( GEX , CellTypeFinal %ni% c( grep(pattern="^rm-",x = unique(GEX$CellTypeFinal),value = T ) ) )
	GEX = subset( GEX , tag_tissue_thr1.5 %ni% c("none") )
	
	GEX.backup=GEX
	#CChat()
	
	#CONDITIONS <- unique( GEX.backup$Condition )
	
	CONDITIONS <- c("mELT") # unique(GEX$tag_tissue_thr1.5)
	
	TYPE_OF_INTER=c("Secreted Signaling","Cell-Cell Contact","ECM-Receptor") #c("Secreted Signaling","Cell-Cell Contact","ECM-Receptor")

	if(T){
		options(stringsAsFactors = FALSE) # !!!!
		NAMEIDENT="CellTypeFinal" #"annot_RNA_snn_res.0.1_CIPR" 
		for( ORIGIN in CONDITIONS ) {
			GEX <- subset( GEX.backup, Condition == ORIGIN)
			for( SELECT in TYPE_OF_INTER){
				OUTDIR=paste0('15_cellchat_',ORIGIN,"_",gsub(" ","",SELECT)) ; dir.create(OUTDIR,showWarnings = F)
				print(paste(">>>>>>>>> Currently:",OUTDIR))
				### setup
				DATABASE="CellChatDB.mouse"
				PPI=PPI.mouse
				# setup object
				# in the future, when no bugs
				#cellchat <- createCellChat(object = GEX)
				Idents(GEX)<-NAMEIDENT #DefaultAssay(GEX) # RNA
				cellchat <- createCellChat(object = GetAssayData(GEX,assay="RNA",slot="scale.data"),
										   meta=GEX@meta.data,group.by=NAMEIDENT)
				cellchat@DB <- subsetDB(CellChatDB.mouse,search=SELECT) # use Secreted Signaling for cell-cell communication analysis
				
				# process (warning: use the proper input data ~ species: PPI.mouse)
				if(T){
					### We first identify over-expressed ligands or receptors in one cell group, 
					#and then project gene expression data onto protein-protein interaction (PPI) network.
					#The over-expressed ligand-receptor interactions are identified if either the ligand or receptor is over-expressed.
					
					cellchat <- subsetData(cellchat) # subset the expression data of signaling genes for saving computation cost
					#future::plan("multiprocess", workers = 4) # do parallel
					cellchat <- identifyOverExpressedGenes(cellchat)	
					cellchat <- identifyOverExpressedInteractions(cellchat)
					cellchat <- projectData(cellchat,PPI)
					
					#Inference of cell-cell communication network
					cellchat <- computeCommunProb(cellchat,raw.use=F,do.fast=F) #,raw.use=F,trim=0.1,type="truncatedMean",population.size=T)
					
					# Filter out the cell-cell communication if there are only few number of cells in certain cell groups
					cellchat <- filterCommunication(cellchat,min.cells=10)
					
					#Infer the cell-cell communication at a signaling pathway level
					cellchat <- computeCommunProbPathway(cellchat)
					
					#Calculate the aggregated cell-cell communication network
					cellchat <- aggregateNet(cellchat)
					
					saveRDS(cellchat, paste0(OUTDIR,"/cellchat.mouse.Rds"))
				}
				#extract signif interaction - save in excel
				if(T){
					#cellchat@meta
					#returns a data frame consisting of all the inferred cell-cell communications at the level of 
					#ligands/receptors. Set slot.name = "netP" to access the the inferred communications at the 
					#level of signaling pathways
					rm(df.net)
					try({ df.net <- subsetCommunication(cellchat) })
					if( exists("df.net") ) { write.xlsx(df.net,paste0(OUTDIR,"/","df.net.xlsx")) }else{print("no df.net")}
					
				}
				
				rm(df.net,cellchat)
				try({ df.net<-read.xlsx(paste0(OUTDIR,"/","df.net.xlsx")) })
				
				if ( exists("df.net") ){ #& nrow(df.net) > 0
					### If already processed:
					cellchat<-readRDS(paste0(OUTDIR,"/cellchat.mouse.Rds"))
					
					#Visualization and systems analysis of cell-cell communication network
					
					# circular links: not amazing since all at once
					if(T){
						OUTNAME=paste0(OUTDIR,"/",DATABASE,".netVisual_circle.number.png");png(OUTNAME,width=800,height=800)
						netVisual_circle(cellchat@net$count,vertex.weight=as.numeric(table(cellchat@idents)),weight.scale=T,label.edge=F,title.name="Number of interactions")
						dev.off()
						
						OUTNAME=paste0(OUTDIR,"/",DATABASE,".netVisual_circle.strength.png");png(OUTNAME,width=800,height=800)
						netVisual_circle(cellchat@net$weight,vertex.weight = as.numeric(table(cellchat@idents)),weight.scale=T,label.edge=F,title.name="Interaction weights/strength")
						dev.off()
						
						rm(OUTNAME)
					}
					
					# circular links SPLIT per celltype // circle // # This split the previous "weighted" circular plot
					#summarizing the communication probability
					if(T){
						dir.create(paste0(OUTDIR,"/circle/"))
						groupSize <- as.numeric(table(cellchat@idents))
						mat <- cellchat@net$weight
						for (index in 1:nrow(mat)) {
							#create new emtpy matrix
							mat2 <- matrix(0,nrow=nrow(mat),ncol=ncol(mat),dimnames=dimnames(mat))
							#fill the new patrix with only data from the current index/ celltype
							mat2[index, ] <- mat[index, ]
							NAME=rownames(mat)[index]
							NAME=gsub("[/]","-",NAME)
							print(paste("index",index,":",NAME))
							OUTNAME=paste0(OUTDIR,"/circle/",NAME,".pdf");pdf(OUTNAME,width=8,height=8,useDingbats=F)
							netVisual_circle(mat2,vertex.weight=groupSize,
											 weight.scale=T,
											 edge.weight.max=max(mat),
											 title.name=rownames(mat)[index],vertex.label.cex = .8)
							dev.off()
						}
						rm(mat,mat2,NAME,OUTNAME,mat,groupSize,index)
					}
					
					# Compute and visualize the contribution of each (significant)
					# ligand-receptor pair in the overall signaling pathways 
					if(T){
						SIGNALING=cellchat@netP$pathways
						
						WIDTH=10
						HEIGHT=10
						FONT_SIZE=12
						
						OUTNAME=paste0(OUTDIR,"/","L-R pair Contribution.pdf")
						pdf(OUTNAME,width=WIDTH,height=HEIGHT)
						CONTRIB=netAnalysis_contribution(cellchat,signaling=SIGNALING,thresh=0.05,return.data=T,font.size=FONT_SIZE)
						print(CONTRIB$gg.obj)
						dev.off()
						
						#CONTRIB$LR.contribution$name
						CONTRIB$LR.contribution$name = factor(CONTRIB$LR.contribution$name,levels=CONTRIB$LR.contribution$name[order(CONTRIB$LR.contribution$contribution)])
						CONTRIB$LR.contribution <- subset( CONTRIB$LR.contribution[ CONTRIB$LR.contribution$name != "1",])
						PLOT=ggplot(CONTRIB$LR.contribution) +
							geom_bar(aes(x=contribution*100,y=name),stat="identity") +
							theme_minimal(base_size=FONT_SIZE) +
							ggtitle("Contribution of each L-R pair\nin the overall signaling pathways") +
							xlab("Relative contribution (%)")+
							ylab("") 
						
						OUTNAME=paste0(OUTDIR,"/","L-R pair Contribution.2.pdf");pdf(OUTNAME,width=WIDTH,height=HEIGHT);print(PLOT);dev.off()
						#sum(CONTRIB$LR.contribution$contribution)
						rm(CONTRIB,PLOT,OUTNAME)
					}
					
					# Bubble plot
					if(T){
						dir.create( paste0(OUTDIR,"/","BubblePlot/" ) )
						for( SOURCENAME in levels(cellchat@idents) ){
							try({
								#SOURCENAME="T cell" #levels(cellchat@idents)
								print(SOURCENAME)
								SOURCE = which( levels(cellchat@idents) == SOURCENAME)
								TARGET = 1:length(levels(cellchat@idents))
								#TARGET = setdiff( 1:length(unique(cellchat@idents)) ,c(SOURCE) )
								#TARGET = c(2,3,4,6,9,10)
								PLOT=netVisual_bubble(cellchat,sources.use=SOURCE,targets.use=TARGET,remove.isolate=T,font.size=12)
								OUTNAME=paste0(OUTDIR,"/","BubblePlot/",SOURCENAME,".pdf");pdf(OUTNAME,width=7,height=14,useDingbats=F);print(PLOT);dev.off()
							},silent=F)
							rm(PLOT,OUTNAME,SOURCE,TARGET,SOURCENAME)
						}
					}
					
					# circular and hierachical vizualisation 
					if(T){
						# Signif pathway
						SIGNALING=cellchat@netP$pathways
						#same: SIGNALING=unique(df.net [ df.net$pval <=0.05, "pathway_name"])
						dir.create(paste0(OUTDIR,"/pair_LR/"))
						pairLR <- extractEnrichedLR(cellchat, signaling = SIGNALING, geneLR.return = FALSE)
						for( index in 1:nrow(pairLR)){
							LR.show <- pairLR[index,] # show one ligand-receptor pair
							print(LR.show)
							#length(levels(cellchat@idents))
							vertex.receiver = seq(1,4) # a numeric vector: the cell group number TARGET
							PLOT=netVisual_individual(cellchat,signaling=SIGNALING,pairLR.use=LR.show,vertex.receiver=vertex.receiver,layout="hierarchy") #hierarchy, circle, chord
							#PLOT0netVisual_individual(cellchat,signaling=SIGNALING,pairLR.use=LR.show,vertex.receiver=vertex.receiver,layout="circle") #hierarchy, circle, chord
							#df.net[ df.net$interaction_name == LR.show, ]
							
							dev.off()
							
							OUTNAME=paste0(OUTDIR,"/pair_LR/",LR.show,".pdf");pdf(OUTNAME,width=14,height=7,useDingbats=F);print(PLOT);dev.off()
						}
						rm(PLOT,OUTNAME,vertex.receiver,LR.show,pairLR)
					}
				}else{
					print("NO SIGNIF INTERRACTIONS")
				}
				#clean
				rm(cellchat,SIGNALING)
			}
			rm(GEX,ORIGIN,SELECT,OUTDIR,df.net,PPI,DATABASE)
		}
		options(stringsAsFactors = TRUE)
	}
}

#### DATA=sessionInfo() ####
#DATA=devtools::session_info()
writeLines(capture.output(devtools::session_info()), "sessionInfo.txt")



#We are most interested in how mELT is different from the other tissues




#### some figures ####
if(F){
	table(GEX$CellTypeFinal)
	table(GEX$treatment)
	table(GEX$tissue)
	FEATURE="Rpl12" #"Mbp"
	Idents(GEX)<-"split"
	VlnPlot( subset(GEX, tissue %ni% c("none") & 
						CellTypeFinal %ni% c("rm-1","rm-2") & 
						split %ni% c("CD20 mELT")), 
			 features=FEATURE,split.by="CellTypeFinal") + 
		geom_hline(yintercept=c(2,4,6),lty=2)
	
	GEX$TEST <- paste(GEX$treatment, GEX$tissue, GEX$CellTypeFinal)
	Idents(GEX)<-"TEST"
	MARKERS<-FindAllMarkers(GEX, return.thresh = 0.05)
	write.xlsx(MARKERS,file = "MARKERS.xlsx")

}

# velocyto? can not be, needs the reads
#as.loom(GEX,filename="GEX.loom",verbose=T)

#### Mbp isoforms ###
if(F){
	#BiocManager::install("muscle")
	#BiocManager::install("msa")
	library(msa)
	library(Biostrings)
	library(seqinr)
	
	
	#https://www.uniprot.org/uniprot/P04370#sequences
	FASTA<-msa::readAAStringSet("/home/ga94rac/Desktop/Mbp_protein_isoforms.fasta",seqtype="AA",as.string=T)
	
	ALIGN <- muscle(FASTA) # input AAStringSet, output MultipleAlignment
	
	#earch for variant in the reads. Problem is: 
	library(ShortRead)
	library(BiocParallel)
}

