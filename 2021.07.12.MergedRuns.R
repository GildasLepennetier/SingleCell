
#library(future)
#plan("multiprocess", workers = 4)
#plan()
#library(parallel)
#options(Ncpus = 4)

####  ------- startup ------- ####
suppressPackageStartupMessages({
	library(limma) # install.packages("BiocManager",Ncpus=4) ; BiocManager::install("limma") # For a more efficient implementation of the Wilcoxon Rank Sum Test for seurat
	library(Seurat) #install.packages("SeuratObject",Ncpus=4) update to avoid incompatibility with matrix
	library(ggplot2) #ggplot
	library(ggrepel) #geom_label_repel() not working with Dimplot?
	library(dplyr) #%>%
	library(tidyr) #separate
	library(ggpubr) #grids
	library(openxlsx) # read.xlsx
	library(patchwork) #wrap_plots
	library(gtools) #mixedsort
	library(cowplot) # help for ggplot2; themes, functions to align plots and arrange them
	library(harmony) # BiocManager::install("SingleCellExperiment") ; install.packages("devtools",Ncpus=4) ; devtools::install_github("hadley/devtools") ; library(devtools) ; install_github("immunogenomics/harmony")
	require(CIPR) #devtools::install_github("atakanekiz/CIPR-Package", build_vignettes = F)
	library(Nebulosa) #BiocManager::install("Nebulosa") # plot_density (bioconductor) #https://bioconductor.org/packages/devel/bioc/vignettes/Nebulosa/inst/doc/nebulosa_seurat.html
	library(nichenetr) #devtools::install_github("saeyslab/nichenetr")
})
source("https://raw.githubusercontent.com/GildasLepennetier/R-stuff/master/not_in.R")
source("https://raw.githubusercontent.com/GildasLepennetier/R-stuff/master/paste_NA.R")
source("https://raw.githubusercontent.com/GildasLepennetier/R-stuff/master/NA_padding.R")
WORKINGDIR="/home/ga94rac/WORK/RSTUDIO/210511.MergedRuns"
dir.create(WORKINGDIR,s=F) ; setwd(WORKINGDIR)

#library(future)
#plan("multiprocess", workers = 4)
#plan()

#### merge the data from previous experiment - already fildered per run, should be fine to use it as-is ####
if(F){

	#load data 
	#RAW DATA = 191.8 Go in /media/ga94rac/toshibackup_3/DATA_SAVE/190707_AGLH_GEX_Helmholtz_Zentrum
	GEX.201207<-readRDS("/home/ga94rac/WORK/RSTUDIO/201207.AGLH.mELT_LN/09.final.Rds")
	# An object of class Seurat 
	# 32285 features across 18287 samples within 1 assay 
	# Active assay: RNA (32285 features, 2000 variable features)
	# 4 dimensional reductions calculated: pca, umap, tsne, harmony
	
	#load data
	# RAW DATA = 304.7 GO in /media/ga94rac/toshibackup_3/DATA_SAVE/210228.CITE-seq_AGLH
	# Processed cellranger: /media/ga94rac/toshibackup_3/DATA_SAVE/210228.CITE-seq_AGLH_processed/OUT_CellRangMult
	
	GEX.210228<-readRDS("/home/ga94rac/WORK/RSTUDIO/210408.CITEseq.Tissues.FINAL/07.final.Rds")
	# 32300 features across 61243 samples within 4 assays 
	# Active assay: RNA (32285 features, 2000 variable features)
	# 3 other assays present: Ab, AbTissues, AbCelltypes
	# 3 dimensional reductions calculated: pca, harmony, umap
	
	GEX<-merge(x = GEX.201207, y = GEX.210228, add.cell.ids = c("Run201207","Run210228"))
	rm(GEX.201207,GEX.210228) ; gc()
	saveRDS(GEX,"01.Rds")
	# 32300 features across 79530 samples within 4 assays 
	# Active assay: RNA (32285 features, 0 variable features)
	# 3 other assays present: Ab, AbTissues, AbCelltypes
}
GEX<-readRDS("01.Rds")

#library(future)
#plan("multiprocess", workers = 2) # not more than 4 CPU in RStudio session
#plan()
#The total size of the 15 globals that need to be exported for the future expression (‘FUN()’) is 1.20 GiB.
#his exceeds the maximum allowed size of 500.00 MiB (option 'future.globals.maxSize')
#options(future.globals.maxSize= 1300*1024^2)


#### PROCESS - GEX #### 
# After filtering and removing cells, one should run the Norm & scale again | no need to re-normalize in theory | https://github.com/satijalab/seurat/issues/4249
if(F){
	#GEX <- NormalizeData(GEX) #, normalization.method="LogNormalize", scale.factor=10000
	GEX <- FindVariableFeatures(GEX) #, selection.method="vst", nfeatures=2000
	
	#scaling using only the variable features
	GEX <- ScaleData(GEX, vars.to.regress = c("S.Score", "G2M.Score") )
	
	GEX <- RunPCA(GEX) #, features=VariableFeatures(object=GEX)
	
	GEX$SAMPLES=GEX$orig.ident
	table(GEX$SAMPLES)
	GEX <- RunHarmony(GEX, "SAMPLES")
	
	GEX <- RunUMAP(GEX, reduction = "harmony", dims = 1:30)
	
	#	saveRDS(GEX,"02.scaled.Rds")
	#	GEX<-readRDS("02.scaled.Rds")
	
	#remove annot that we will not use
	for ( NAMES in grep("RNA_snn",names(GEX@meta.data),value=T) ) { GEX@meta.data[ , NAMES ] <- NULL }
	for ( NAMES in grep("DF.classifications",names(GEX@meta.data),value=T) ) { GEX@meta.data[ , NAMES ] <- NULL }
	
	GEX <- FindNeighbors(GEX, reduction = "harmony", dims = 1:30)
	
	for( RESO in c( 0.5, 1, 1.5, 2 ) ){GEX <- FindClusters(GEX, resolution = RESO, random.seed = 2020)} ; rm(RESO)
	
	for(index in grep("RNA_snn_res",x=names(GEX@meta.data),value=T)){print(paste(index,length(table(GEX@meta.data[,index]))))};rm(index)
	# [1] "RNA_snn_res.0.5 25"
	# [1] "RNA_snn_res.1 30"
	# [1] "RNA_snn_res.1.5 37"
	# [1] "RNA_snn_res.2 43"
	
}

# library(future)
# plan("multiprocess", workers = 4) ; plan()

#### FindAllMarkers #### 
SUBDIR="11.Cluster_Harmony";dir.create(SUBDIR)
if(F){
	#DefaultAssay(GEX)
	NAMES_CLU=grep("^RNA_snn_res",names(GEX@meta.data),value=T)
	#NAMES_CLU=c("RNA_snn_res.0.1", "RNA_snn_res.0.5", "RNA_snn_res.1","RNA_snn_res.1.5") # "RNA_snn_res.0.1", "RNA_snn_res.0.5", "RNA_snn_res.1","RNA_snn_res.1.5", "RNA_snn_res.2"
	NAMES_CLU="RNA_snn_res.1.5" # "RNA_snn_res.2" 
	NAMES_CLU="RNA_snn_res.2" # "RNA_snn_res.2" 
	
	#options(future.globals.maxSize= 600*1024^2 )
	
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
#### annotation celltypes #### 
SUBDIR="11.Cluster_Harmony";dir.create(SUBDIR,showWarnings=F)
if(F){
	NAMES_CLU=grep("^RNA_snn_res",names(GEX@meta.data),value=T)
	#NAMES_CLU=c( "RNA_snn_res.2") # "RNA_snn_res.0.1", "RNA_snn_res.0.5", "RNA_snn_res.1", "RNA_snn_res.1.5",
	for ( NAMEIDENT in NAMES_CLU ) { # c("RNA_snn_res.0.1", "RNA_snn_res.0.2") ) { #
		Idents(GEX) <- NAMEIDENT
		print(paste0("doing: ",NAMEIDENT, " (",length(unique(Idents(GEX)))," clusters)"))
		# # # # CIPR - Cluster Identity PRedictor https://github.com/atakanekiz/CIPR-Package
		if(T){
			#required for CIPR
			
			allMarkers<-readRDS(paste0(SUBDIR,"/allMarkers.",NAMEIDENT,".Rds"))
			#allMarkers$logfc <- allMarkers$avg_log2FC #log( allMarkers$avg_log2FC^2 , base = 10) #Data frame must have columns named as 'gene', 'logfc', and 'cluster'
			#allMarkers$cluster<-as.character(allMarkers$cluster)
			#DATABASE="mmrnaseq"
			DATABASE="immgen"

			
			# Plot summarizing top scoring references per cluster (logFC comparison)
			PLOT <- CIPR(input_dat = allMarkers,comp_method="logfc_dot_product",reference=DATABASE,top_num=1,keep_top_var=50,plot_ind = F, plot_top = F) #+ grids(axis="x",color="grey92",size=NULL,linetype=NULL) #"dotted");OUTNAME=paste0("7.CIPR.",NAMEIDENT,"_",DATABASE,".pdf");pdf(OUTNAME,width=14,height=7,useDingbats=F);print(PLOT);dev.off()
			ANNOT = CIPR_top_results %>% 
				select( cluster,reference_cell_type) %>% 
				distinct() %>%
				group_by(cluster) %>% mutate( annot = paste0(reference_cell_type,collapse="+")) %>%
				select(cluster,annot) %>% distinct()
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
		
	}
}
#### merge gene markers + CIPR annot #### 
SUBDIR="11.Cluster_Harmony"
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

#### add run info ####
if(F){
	GEX$tissue_merged=paste_NA(GEX$tissue,GEX$tag_tissue_thr1.5)
	GEX$tissue_merged<-gsub("^LN$","LN_lumbar",GEX$tissue_merged) #replace ol LN annotation by new one
	table(GEX$tissue_merged,useNA="ifany")
	
	GEX$run = GEX$orig.ident
	GEX$run[ GEX$orig.ident %in% c("Z3.E21.LN","Z3.E21.mELT","Z3.E72.LN","Z3.E72.mELT","Z3.F21.LN","Z3.F21.mELT")] <-"run1"
	GEX$run[ GEX$orig.ident %in% c("I18A","I18B","J119A","J119B","J119C","J9A","J9B","J9C")] <-"run2"
}
#check run effect on T cell -> fine each cell type has from eah run evenly
#TABLE=as.data.frame(table(GEX$annot_RNA_snn_res.2_CIPR,GEX$run,useNA="ifany"))
#TABLE<-TABLE %>% group_by(Var2) %>% mutate(percent=Freq/sum(Freq)*100) %>% select(Var1,Var2,percent) 
#ggplot(TABLE) + geom_bar(aes(Var1,percent,fill=Var2),stat="identity")

#	saveRDS(GEX,"03.Clu.Rds")
#	GEX<-readRDS("03.Clu.Rds")

#### PLOT clusters + annot #### 
SUBDIR="12.Plots";dir.create(SUBDIR)
if(F){
	NAMES_CLU=c("RNA_snn_res.2","annot_RNA_snn_res.2_CIPR") ; ORDER=NULL # dir(SUBDIR)
	NAMES_CLU=c("immune_chains_clean2"); ORDER=c("B-cell","T-cell","other") # table(GEX$immune_chains_clean2)
	NAMES_CLU="orig.ident"; ORDER=NULL 
	NAMES_CLU="Phase"; ORDER=c("G2M","S","G1")
	
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
	
	
	NAMES_CLU=c("annot_RNA_snn_res.2_CIPR") ; ORDER=NULL # dir(SUBDIR)
	
	SPLIT="tissue_merged" ; WIDTH=16;HEIGHT=12; NCOL=4
	SPLIT="run"; WIDTH=10;HEIGHT=8; NCOL=2

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
	
	
	rm(NAMES_CLU,NAMEIDENT,WIDTH,HEIGHT,NCOL,ORDER,SPLIT,OUTNAME)
	# axis.line = element_line(size=1),	# text = element_text(size = 12),	# axis.text = element_text(size = 12),	# axis.ticks = element_line(size=1),
}

#### GENE MARKERS: FeaturePlot / plot_density : >>> also redo after filtering <<< #### 
SUBDIR="13.GeneMarkers"; dir.create(SUBDIR)
if(F){
	#MarkersTable<-read.xlsx("/home/ga94rac/LRZ Sync+Share/PROJECT-GILDAS/MarkersTable.xlsx")
	#	grep("TNFSF13B",x = rownames(GEX),ignore.case = T, value = T)
	SPLIT="tissue_merged" #GEX@meta.data
	
	#SPLIT_ELEMENTS=unique(GEX@meta.data[,SPLIT])
	#"LN_lumbar"   "mELT"        "none"        "CSF"         "LN_inguinal" "Spleen"      "Blood"
	SPLIT_ELEMENTS=c("LN_lumbar","mELT","CSF","LN_inguinal","Spleen","Blood")
	
	for( index in c(1:46) ){ 
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
		if(index==14){ FEATURES=c("Nos2","Tlr2","Tlr4","Il12a","Il18","Ly6c1");TITLE="Macrophages M1"} #
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
		
		if(index==28){ FEATURES=c("Cd1d1","S1pr3","Tlr3","Mzb1");TITLE="Marginal zone B-cell" }
		if(index==29){ FEATURES=c("Aicda","Basp1","Fas","Neil1","Plxnb2","Rgs13","S1pr2","Tnfsf9");TITLE="Germinal center B-cell" }
		
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
		if(index==46){ FEATURES=c("Lyz2");TITLE="Myeloid cell" }
		
		
		#		grep("Lyz2",x = rownames(GEX),ignore.case = T, value = T)
		
		for(FEATURE in FEATURES){
			
			print(paste("index:",index,"current:",TITLE,"->",FEATURE))
			
			# GENE MARKERS: FeaturePlot
			print(paste("doing:","FeaturePlot"))
			if(T){
				PLOT=FeaturePlot(GEX,reduction="umap",features=FEATURE,pt.size=1,label=F,repel=T,min.cutoff="q10",ncol=1,coord.fixed=T,order=T)
				OUTNAME=paste0(SUBDIR,"/umap.",TITLE,".",FEATURE,".pdf");pdf(OUTNAME,width=8,height=8,useDingbats=F);print(PLOT);dev.off()
				#OUTNAME=paste0(SUBDIR,"/umap.",TITLE,".",FEATURE,".png");png(OUTNAME,width=500,height=500);print(PLOT);dev.off()
			}
			
			
			# GENE MARKERS: FeaturePlot SPLIT : not good if too many to split (>3 = not nice), and only one feature.
			print(paste("doing:","FeaturePlot split"))
			if(FALSE){
				PLOT=FeaturePlot(GEX,reduction="umap",features=FEATURE,pt.size=1,label=F,split.by=SPLIT,
								 repel=T,min.cutoff="q10",ncol=3,coord.fixed=T,order=T)
				OUTNAME=paste0(SUBDIR,"/umap.","split_",SPLIT,".",TITLE,".",FEATURE,".pdf");pdf(OUTNAME,width=16,height=8,useDingbats=F);print(PLOT);dev.off()
				#OUTNAME=paste0(SUBDIR,"/umap.","split_",SPLIT,".",TITLE,".",FEATURE,".png");png(OUTNAME,width=1200,height=500);print(PLOT);dev.off()
			}
			
			
			# GENE MARKERS: FeaturePlot SPLIT using subset
			print(paste("doing:","FeaturePlot split using subset"))
			if(T){
				for( SPLIT_ELEMENT in SPLIT_ELEMENTS ){
					print(paste("    doing:",SPLIT_ELEMENT))
					try({
						PLOT=FeaturePlot(subset(GEX,tissue_merged == SPLIT_ELEMENT),reduction="umap",features=FEATURE,pt.size=1,label=F, #split.by=SPLIT, using the split apparently also remove the legend
										 repel=T,min.cutoff="q10",ncol=1,coord.fixed=T,order=T)
						OUTNAME=paste0(SUBDIR,"/umap.","split_",SPLIT_ELEMENT,".",TITLE,".",FEATURE,".pdf");pdf(OUTNAME,width=8,height=8,useDingbats=F);print(PLOT);dev.off()
						#OUTNAME=paste0(SUBDIR,"/umap.","split_",SPLIT,".",TITLE,".",FEATURE,".png");png(OUTNAME,width=1200,height=500);print(PLOT);dev.off()
					})
				}
			}
			
			
			# GENE MARKERS - nebulosa: plot_density: nicer than FeaturePlot
			print(paste("doing:","nebulosa"))
			if(T){
				PLOT=plot_density(GEX,reduction="umap",features=FEATURE,size=1)
				OUTNAME=paste0(SUBDIR,"/umap.density.",TITLE,".",FEATURE,".pdf");pdf(OUTNAME,width=8,height=8);print(PLOT);dev.off()
				#OUTNAME=paste0(SUBDIR,"/umap.density.",TITLE,".",FEATURE,".png");png(OUTNAME,width=500,height=500);print(PLOT);dev.off()
			}
			
			
			# GENE MARKERS - nebulosa: plot_density: nicer than FeaturePlot, bias if few cells
			# split using subset
			print(paste("doing:","nebulosa split"))
			if(T){
				
				for( SPLIT_ELEMENT in SPLIT_ELEMENTS ){
					print(paste("    doing:",SPLIT_ELEMENT))
					#to avoid crash, some expression shouldexists
					# !!! tissue_merged may have to be set to something else
					try({
						PLOT=plot_density( subset(GEX,tissue_merged == SPLIT_ELEMENT),reduction="umap",features=FEATURE,size=1)
						OUTNAME=paste0(SUBDIR,"/umap.density.","split_",SPLIT_ELEMENT,".",TITLE,".",FEATURE,".pdf");pdf(OUTNAME,width=8,height=8,useDingbats=F);print(PLOT);dev.off()
						#OUTNAME=paste0(SUBDIR,"/umap.density.","split_",SPLIT_ELEMENT,".",TITLE,".",FEATURE,".png");png(OUTNAME,width=500,height=500);print(PLOT);dev.off()
					})
				}
				
			}
			
			# GENE MARKERS: plot_density: SPLIT by tissue - 
			#the plot_density indeed does not plot using a split option
			
		}
	}
	
	#grep("Baff",x = rownames(GEX),ignore.case = T, value = T)
	
	
	system(paste0("mkdir -p ",SUBDIR,"/B-cells",
				  "; mv ",SUBDIR,"/umap*B-cell*.p[nd][gf] ",SUBDIR,"/B-cells/"))
	system(paste0("mkdir -p ",SUBDIR,"/T-cells",
				  "; mv ",SUBDIR,"/umap*T-cell*.p[nd][gf] ",SUBDIR,"/T-cells/",
				  "; mv ",SUBDIR,"/umap*.Th*.p[nd][gf] ",SUBDIR,"/T-cells/"))
	system(paste0("mkdir -p ",SUBDIR,"/DC ; mv ",SUBDIR,"/umap*Dendritic*.p[nd][gf] ",SUBDIR,"/DC/"))
	system(paste0("mkdir -p ",SUBDIR,"/Macrophages ; mv ",SUBDIR,"/umap*Macrophages*.p[nd][gf] ",SUBDIR,"/Macrophages/"))
	system(paste0("mkdir -p ",SUBDIR,"/NK ; mv ",SUBDIR,"/umap*NK-cell*.p[nd][gf] ",SUBDIR,"/NK/"))
	system(paste0("mkdir -p ",SUBDIR,"/Myelin ; mv ",SUBDIR,"/umap*Myelin*.p[nd][gf] ",SUBDIR,"/Myelin/"))
	
	rm(index,FEATURES,PLOT,OUTNAME,FEATURE,TITLE,SUBDIR,SPLIT)
}


#### to check? for humans the HCL-blast / human cell landscape | to annotate cell types  ####
if(F){
	# ONLINE
	#
	# paper: https://jhoonline.biomedcentral.com/articles/10.1186/s13045-020-00941-y#Sec13
	#
	#http://bis.zju.edu.cn/HCL/index.html
	#http://bis.zju.edu.cn/HCL/blast.html
	#require a .txt or .csv file tab/or coma -sep. Fine with no space of special char
	# Gene symbol as rownames, cell as colnames
	
	#Rpackage
	# require ggplot2/reshape2/plotly/shiny/shinythemes/shiny
	#install.packages("shinythemes")
	#require(ggplot2);require(reshape2);require(plotly);require(shiny);require(shinythemes)
	#devtools::install_github("ggjlab/scHCL")
	library(scHCL)
	
	TABLE_before<-table(GEX$RNA_snn_res.2)
	GEX_sampl <- GEX[, sample(colnames(GEX), size =10000, replace=F)]
	TABLE_after<-table(GEX_sampl$RNA_snn_res.2)
	if(any(TABLE_after==0)) print("please do it again")
	#T1=round(TABLE_before/sum(TABLE_before)*100,1)
	#T2=round(TABLE_after/sum(TABLE_after)*100,1)
	#T2-T1
	
	# scHCL has two parameters , single cell expression matrix(scdata) and 
	# the number of most similar cell types
	MAT <- as.matrix(GEX_sampl@assays$RNA[,])
	
	#colnames(MAT) <- GEX@meta.data[ colnames(MAT), "CellTypeFinal" ] 
	CELLnamesBackup=colnames(MAT)
	colnames(MAT) <- as.character( seq(1,length(colnames(MAT))) )
	MAT[1:5,1:5]
	
	hcl_result <- scHCL::scHCL(scdata = MAT, numbers_plot = 20)
	
	table(unlist(hcl_result$scHCL))
	sum(table(unlist(hcl_result$scHCL))) # 327 / 10000 *100 = 3% annotated
	
	
	scHCL::scHCL_vis(hcl_result)

}


#### CLEAN CLUSTERS - update from gene markers and auto annotation #### 
if(F){
	
	GEX$immune_chains_clean4 <- GEX$immune_chains_clean
	GEX$immune_chains_clean4[GEX$immune_chains_clean %in% c("IGH,IGK","IGH,IGL") ] <- "B-cell-paired"
	GEX$immune_chains_clean4[GEX$immune_chains_clean %in% c("TRA,TRB") ] <- "T-cell-paired"
	
	GEX$CellTypeFinal <- as.character(GEX$annot_RNA_snn_res.2_CIPR)
	GEX$UMAP_1<-GEX@reductions$umap@cell.embeddings[,1]
	GEX$UMAP_2<-GEX@reductions$umap@cell.embeddings[,2]
	#DF <- as.data.frame(GEX@reductions$umap@cell.embeddings) # UMAP_1 & UMAP_2
	
	#cells to remove because pre-T cell 
	if(T){
		GEX@meta.data[ GEX$annot_RNA_snn_res.2_CIPR == "Pre-T cell" ,"CellTypeFinal"] <- "rm-1"
	}
	#cells to remove because weird: annot T-cell but contain top gene of Blood + T = putative extra doublets
	if(T){
		GEX@meta.data[ GEX$RNA_snn_res.2 == "39" ,"CellTypeFinal"] <- "rm-2"
	}
	#cells to remove because expression of Mbp in cluster 11 of the B-cells
	if(T){
		GEX@meta.data[ GEX$RNA_snn_res.2 == "11" ,"CellTypeFinal"] <- "rm-3"
	}
	
	# Cd8a T cells 					(based on Cd8a expression)
	if(T){
		GEX@meta.data[ GEX$RNA_snn_res.2 %in% c("21") ,"CellTypeFinal"] <- "T cell CD8+"
	}
	# activated T cell CD4+			(based on Il2, Pdcd1   -   Cd4 /Cd8a expression )
	if(T){
		GEX@meta.data[ GEX$RNA_snn_res.2 %in% c("15") ,"CellTypeFinal"] <- "act. T cell CD4+"
	}
	# naive T cell CD4+				(based on Ccr7, Foxp1   -   Cd4 /Cd8a expression )
	if(T){
		GEX@meta.data[ GEX$RNA_snn_res.2 %in% c("7","22","37","14","3","2","0","36") ,"CellTypeFinal"] <- "naive T cell CD4+"
	}
	# gd T cell CD4+				(based on Tcrg-C1,Tcrg-C2, Trdc expression )
	if(T){
		GEX@meta.data[ GEX$RNA_snn_res.2 %in% c("9") ,"CellTypeFinal"] <- "gd T cell"
	}
	
	#Treg: done auto
	
	# Erythroid cells				(based on Hbb-bt / Hbb-bs markers)      			NOT: pre-B cell, NOT follicular B cell (based on GEDIT: http://webtools.mcdb.ucla.edu/)
	if(T){
		GEX@meta.data[ GEX$RNA_snn_res.2 %in% c("27","28","35") ,"CellTypeFinal"] <- "Erythroid cell" # Blood cell
	}
	# marginal zone B cell:  			(based on Cd1d1, S1pr3, Tlr3	# B-cells activated was: Cd1d1, cd86?, Ighg3
	if(T){
		GEX@meta.data[ GEX$RNA_snn_res.2 == "6" ,"CellTypeFinal"] <- "mar. zone B cell" 
	}
	# fol. zone B cell  			(based on NOT activated + Fcer2a ((ex-naive group))
	if(T){
		GEX@meta.data[ GEX$RNA_snn_res.2 %in% c("20","24","1","4","12","10","5","16","33") ,"CellTypeFinal"] <- "fol. zone B cell"
	}
	# pre-B cell > germinal center Bcells
	if(T){
		GEX$CellTypeFinal [ GEX$annot_RNA_snn_res.2_CIPR == "Pre-B cell" ]<- "ger. center B cell"
	}
	# Plasma-cells: par of the pre-B cells  				(based on Sdc1, Tnfrsf17
	if(T){
		GEX$CellTypeFinal [ GEX$RNA_snn_res.2 == "38" ] <- "Plasma cell"
	}
	
	# NK + NK T cells: because we observed several T-cell markers in the NK cells 
	if(T){
		# NK T cell : in the NK cluster but having CDJ data for T cells
		#mixing instead of split	#GEX@meta.data[ GEX$CellTypeFinal == "NK cell" & GEX$immune_chains_clean4 == "T-cell-paired","CellTypeFinal"] <- "NK T cell"
		
		GEX$CellTypeFinal [ GEX$RNA_snn_res.2 %in% c("13","17","42") ] <- "NK + NK T cell"
		
	}
	
	#cells to remove because expression of Mbp in cluster 11 of the B-cells
	if(T){
		GEX$CellTypeFinal [ GEX$RNA_snn_res.2 == "41" ] <- "rm-4"
		#GEX@meta.data[ GEX$RNA_snn_res.2 == "41" ,"CellTypeFinal"] <- "rm-4"
	}
	
	#removeing cluster 12 because some T-cell markers in top genes
	if(T){
		GEX$CellTypeFinal [ GEX$RNA_snn_res.2 == "12" ] <- "rm-5"
	}
}

#### prepare variable to clean this data using VDJ / HTO ####
if(F){
	#### MEMO concerning the data cleaning using VDJ ###
	# Run1 was not cleaned
	# Run2 was cleaned using HTO (only clear positive removed when in unexpected celltype) + VDJ (paired only)
	if(FALSE){
		# FILTER THE B-CELLS FROM VDJ
		#
		# notes concerning the names
		# "tissue_merged" -> run1 + run2 merged names.
		# "run" 
		# # The first run (mELT/LN + VDJ T cells + B cells)
		# # The second run (Tissues + CITE-seq + VDJ B cells only)
		# # cells already cleaned from the sedonc in theory
		
		options(max.print=250)
		table(GEX$immune_chains_clean, GEX$CellTypeFinal)
		
		
		# PREPARE CHAINS INFO
		# VDJ - SAME FOR run 1 & run 2: MEMO: was done this way
		# NOTE: takes the light chain as "B-cell" but may not be true
		if(FALSE){
			GEX@meta.data$immune_chains=paste(GEX@meta.data$IGH_chain,GEX@meta.data$TRA_chain,GEX@meta.data$TRB_chain,GEX@meta.data$IG_ligh_chain,sep=",")
			df<-GEX@meta.data %>% 
				select(IGH_chain,TRA_chain,TRB_chain,IG_ligh_chain) %>%
				tidyr::unite(col="immune_chains",IGH_chain,TRA_chain,TRB_chain,IG_ligh_chain, na.rm=TRUE, sep=",")
			df$immune_chains[ df$immune_chains == ""]<-"no info"
			GEX<-AddMetaData(GEX,df)
			KEEP=c("TRA","TRB","TRA,TRB","IGH","IGH,IGK","IGH,IGL","IGK","IGL","no info")
			GEX@meta.data$immune_chains_clean<-GEX@meta.data$immune_chains
			GEX@meta.data$immune_chains_clean[ ! GEX@meta.data$immune_chains_clean %in% KEEP ] <- "multi"
			GEX@meta.data$immune_chains_clean2=GEX@meta.data$immune_chains_clean
			GEX@meta.data$immune_chains_clean2[GEX@meta.data$immune_chains_clean2 %in% c("TRA","TRB","TRA,TRB") ] <- "T-cell"
			GEX@meta.data$immune_chains_clean2[GEX@meta.data$immune_chains_clean2 %in% c("IGH","IGH,IGK","IGH,IGL","IGK","IGL") ] <- "B-cell"
			GEX@meta.data$immune_chains_clean2[ ! GEX@meta.data$immune_chains_clean2 %in% c("T-cell","B-cell","multi") ] <- "other"
			rm(KEEP,df)
			table(GEX$immune_chains_clean2, GEX$CellTypeFinal)
			table(GEX$immune_chains_clean2, GEX$run, useNA="ifany")
			# 			run1	run2
			# B-cell	3273	15599		sum=18872
			# other		10059	45644		sum=55703
			# T-cell	4955	0
		}
		# VDJ clean3 : only in RUN 2 ; only B-cells, without taking the light chain into account
		if(FALSE){
			GEX$immune_chains_clean3 <- "other"
			GEX$immune_chains_clean3[GEX$immune_chains_clean %in% c("IGH,IGK","IGH,IGL") ] <- "B-cell"
			table(GEX$immune_chains_clean3, GEX$CellTypeFinal)
			table(GEX$immune_chains_clean3, GEX$run, useNA="ifany") # if "IGH,IGK","IGH,IGL") then "B-cell"
			# 			run1	run2
			# B-cell	0		10676
			# other		0		50567
			# <NA>		18287	0
		}
		# paired - unpaired data in VDJ
		if(FALSE){
			TABLE=table(GEX$immune_chains_clean, GEX$run, useNA="ifany");TABLE
			#			run1  run2
			# IGH,IGK  1554  9684				total B paired = 12414 (  66 % )
			# IGH,IGL   184   992
			# 
			# IGH        49    94				total unpaired = 6458 ( 34 % )
			# IGK      1392  4400
			# IGL        94   429
			# 
			# no info 10059 45644
			# 
			# TRA        47     0
			# TRA,TRB  3958     0				3958 / (3958+47+950) # 0.79
			# TRB       950     0				(47+950) / (3958+47+950) # 0.20
		}
		
	}
	
	GEX$immune_chains_clean4_selection = TRUE
	B_ANNOT_LIST=c("fol. zone B cell","ger. center B cell","mar. zone B cell","Plasma cell")
	T_ANNOT_LIST=c("act. T cell CD4+","gd T cell","naive T cell CD4+","T cell CD8+","Treg","NK T cell") #NK T cells have are here and with NK cells
	NK_ANNOT_LIST=c("NK + NK T cell")
	GRA_MONO_ANNOT_LIST=c("Granulocyte","Monocyte")
	
	GEX$immune_chains_clean4_selection [ GEX$immune_chains_clean4 == "B-cell-paired" & ! GEX$CellTypeFinal %in% B_ANNOT_LIST] <- FALSE
	GEX$immune_chains_clean4_selection [ GEX$immune_chains_clean4 == "T-cell-paired" & ! GEX$CellTypeFinal %in% T_ANNOT_LIST] <- FALSE
	
	GEX$immune_chains_clean4_selection [ ! GEX$CellTypeFinal %in% B_ANNOT_LIST & GEX$tag_celltype %in% c("CD19") ] <- FALSE
	GEX$immune_chains_clean4_selection [ ! GEX$CellTypeFinal %in% T_ANNOT_LIST & GEX$tag_celltype %in% c("CD4","CD8a")] <- FALSE
	GEX$immune_chains_clean4_selection [ ! GEX$CellTypeFinal %in% NK_ANNOT_LIST & GEX$tag_celltype %in% c("NK") ] <- FALSE
	GEX$immune_chains_clean4_selection [ ! GEX$CellTypeFinal %in% GRA_MONO_ANNOT_LIST & GEX$tag_celltype %in% c("Ly6G") ] <- FALSE
	
	rm(B_ANNOT_LIST,T_ANNOT_LIST,NK_ANNOT_LIST,GRA_MONO_ANNOT_LIST)

}

#### UMAPs -- PLOT clusters + annot CellTypeFinal >>> also redo after filtering <<< #### 
SUBDIR="14.Plots.CellTypeFinal";dir.create(SUBDIR)
if(F){
	#table(GEX$CellTypeFinal) ; unique(GEX$CellTypeFinal)
	#unique(GEX$CellTypeFinal)
	NAMES_CLU=c("CellTypeFinal") ; ORDER=c("naive T cell CD4+"="#00bf74",
										   "act. T cell CD4+"="#b79f00",
										   "gd T cell"="#619cff",
										   "T cell CD8+"="#ff699c",
										   "Treg"="#db72fb",
										   
										   "NK + NK T cell"="#00ba38",
										   
										   "fol. zone B cell"="#00adfa",
										   "ger. center B cell"="#d39200",
										   "Plasma cell"="#ae87ff",
										   "mar. zone B cell"="#93aa00",
										   
										   "Erythroid cell"="#f8766d",
										   "Macrophage"="#f564e3",
										   "Granulocyte"="#00b9e3",
										   "DC"="#5eb300",
										   "Monocyte"="#e88526",
										   "rm-5"="gray50","rm-4"="gray55","rm-3"="gray60","rm-2"="gray65","rm-1"="gray70")#first are plotted last, then are in top
	#f8766d, e88526, d39200, b79f00, 93aa00, 5eb300, 00ba38, 00bf74, 00c19f, 00bfc4, 00b9e3, 00adfa, 619cff, ae87ff, db72fb, f564e3, ff61c3, ff699c
	
	
	#rainbow(15)
	#RColorBrewer::display.brewer.all()
	#RColorBrewer::brewer.pal(15, "Dark2")
	#DiscretePalette(15, palette = NULL)
	
	# normal plot with rm
	for ( NAMEIDENT in NAMES_CLU){
		print(paste("doing:",NAMEIDENT))
		Idents(GEX)<-NAMEIDENT
		PLOT=DimPlot(GEX,reduction="umap",label=T,pt.size=.1,label.size=3,repel=T, cols=ORDER) + #order=ORDER,
			theme_minimal(base_size = 16)+
			theme(legend.position = 'top')+
			ggtitle(paste(NAMEIDENT)) + coord_fixed()
		OUTNAME=paste0(SUBDIR,"/umap.",NAMEIDENT,".pdf");pdf(OUTNAME,width=10,height=8,useDingbats=F);print(PLOT);dev.off()
		rm(PLOT,OUTNAME)
	}
	
	# remove the rm clusters (ignored)
	CELLTYPES=unique(GEX$CellTypeFinal)
	CELLTYPES=grep("^rm-",CELLTYPES,value=T,invert=T) 
	for ( NAMEIDENT in NAMES_CLU){
		print(paste("doing:",NAMEIDENT))
		GEX.tmp=subset( GEX , CellTypeFinal %in% CELLTYPES)
		Idents(GEX.tmp)<-NAMEIDENT
		PLOT=DimPlot(GEX.tmp,reduction="umap",label=T,order=ORDER,pt.size=.1,label.size=3,repel=T, cols=ORDER) + 
			theme_minimal(base_size = 16)+
			theme(legend.position = 'top')+
			ggtitle(paste(NAMEIDENT)) + coord_fixed()
		OUTNAME=paste0(SUBDIR,"/umap.",NAMEIDENT,".noRm.pdf");pdf(OUTNAME,width=10,height=8,useDingbats=F);print(PLOT);dev.off()
		rm(PLOT,OUTNAME)
	}
	
	# normal plot for each single tissue
	TISSUES=unique(GEX$tissue_merged) # c("LN_lumbar","mELT","LN_inguinal","Spleen") #unique(GEX$tissue_merged)
	for( TISSUE in TISSUES ){
		GEX.tmp=subset( GEX , tissue_merged == TISSUE)
		for ( NAMEIDENT in NAMES_CLU){
			print(paste("doing:",NAMEIDENT))
			Idents(GEX.tmp)<-NAMEIDENT
			PLOT=DimPlot(GEX.tmp,reduction="umap",label=T, cols=ORDER,pt.size=.1,label.size=3,repel=T) + 
				theme_minimal(base_size = 16)+
				theme(legend.position = 'top')+
				ggtitle(paste(NAMEIDENT)) + coord_fixed()
			OUTNAME=paste0(SUBDIR,"/umap.",NAMEIDENT,".",TISSUE,".pdf");pdf(OUTNAME,width=10,height=8,useDingbats=F);print(PLOT);dev.off()
			rm(PLOT,OUTNAME)
		}
	}
	
	# remove the rm clusters (ignored)
	CELLTYPES=unique(GEX$CellTypeFinal)
	CELLTYPES=grep("^rm-",CELLTYPES,value=T,invert=T) 
	for( TISSUE in TISSUES ){
		GEX.tmp=subset( GEX , tissue_merged == TISSUE & CellTypeFinal %in% CELLTYPES)
		for ( NAMEIDENT in NAMES_CLU){
			print(paste("doing:",NAMEIDENT))
			Idents(GEX.tmp)<-NAMEIDENT
			PLOT=DimPlot(GEX.tmp,reduction="umap",label=T, cols=ORDER,pt.size=.1,label.size=3,repel=T) + 
				theme_minimal(base_size = 16)+
				theme(legend.position = 'top')+
				ggtitle(paste(NAMEIDENT)) + coord_fixed()
			OUTNAME=paste0(SUBDIR,"/umap.",NAMEIDENT,".",TISSUE,".noRm.pdf");pdf(OUTNAME,width=10,height=8,useDingbats=F);print(PLOT);dev.off()
			rm(PLOT,OUTNAME)
		}
	}
	
	# split plot 
	#GEX$split = paste(GEX$treatment, GEX$tag_tissue_thr1.5)
	for( index in c(1,2) ){
		if( index == 1) { SPLIT="tissue_merged" ; WIDTH=20;HEIGHT=16; NCOL=4 ; LABSIZE=3 }
		if( index == 2) { SPLIT="orig.ident"; WIDTH=20;HEIGHT=16; NCOL=3; LABSIZE=3 }
		#SPLIT="split" ; WIDTH=20;HEIGHT=12; NCOL=3
		
		for ( NAMEIDENT in NAMES_CLU){
			print(paste(NAMEIDENT,"split:",SPLIT))
			Idents(GEX)<-NAMEIDENT
			PLOT=DimPlot(GEX, reduction="umap", label=TRUE,split.by=SPLIT, cols=ORDER,pt.size=.5,
						 label.size=LABSIZE,repel=T,ncol=NCOL) + 
				theme_minimal(base_size=20)+
				theme(legend.position='top')+
				ggtitle(paste(NAMEIDENT)) + coord_fixed()
			OUTNAME=paste0(SUBDIR,"/umap.",NAMEIDENT,".Split_",SPLIT,".pdf");pdf(OUTNAME,width=WIDTH,height=HEIGHT,useDingbats=F);print(PLOT);dev.off()
		}
	}
	
	# split plot - without the rm clusters (ignored)
	CELLTYPES=unique(GEX$CellTypeFinal)
	CELLTYPES=grep("^rm-",CELLTYPES,value=T,invert=T) 
	GEX.tmp=subset( GEX , CellTypeFinal %in% CELLTYPES)
	for( index in c(1,2) ){
		if( index == 1) { SPLIT="tissue_merged" ; WIDTH=20;HEIGHT=16; NCOL=4 ; LABSIZE=3 }
		if( index == 2) { SPLIT="orig.ident"; WIDTH=20;HEIGHT=16; NCOL=3; LABSIZE=3 }
		#SPLIT="split" ; WIDTH=20;HEIGHT=12; NCOL=3
		
		for ( NAMEIDENT in NAMES_CLU){
			print(paste(NAMEIDENT,"split:",SPLIT))
			Idents(GEX.tmp)<-NAMEIDENT
			PLOT=DimPlot(GEX.tmp, reduction="umap", label=TRUE,split.by=SPLIT, cols=ORDER,pt.size=.5,
						 label.size=LABSIZE,repel=T,ncol=NCOL) + 
				theme_minimal(base_size=20)+
				theme(legend.position='top')+
				ggtitle(paste(NAMEIDENT)) + coord_fixed()
			OUTNAME=paste0(SUBDIR,"/umap.",NAMEIDENT,".Split_",SPLIT,".noRm.pdf");pdf(OUTNAME,width=WIDTH,height=HEIGHT,useDingbats=F);print(PLOT);dev.off()
		}
	}
	
	rm(NAMES_CLU,NAMEIDENT,WIDTH,HEIGHT,NCOL,ORDER,SPLIT,OUTNAME,PLOT,LABSIZE,index) 
	# axis.line = element_line(size=1),	# text = element_text(size = 12),	# axis.text = element_text(size = 12),	# axis.ticks = element_line(size=1),
}
#### DOTPLOT gene markers >>> also redo after filtering <<< #### 
SUBDIR="14.Plots.CellTypeFinal/";dir.create(SUBDIR)
if(F){
	#Zhaorong
	#FEATURES=c("Cd3g","Cd3e","Cd3d","Cd2","Cd86","Cd80","Cd68", "S100a9","S100a8","Itgb3bp","Krt18","Pecam1","Cd34","Cx3cr1","Itgam","Olig1","Cldn11","Mbp","Gja1","Aldoc","Slc1a2","Hmmr")
	
	#"Il2", 			#"act. T",
	#"S1pr3", "Tlr3",	#"Marg. B","Marg. B",
	#FEATURES=c("Hbb-bt","Sdc1","Cd68","S100a8","Rora","Foxp3","Cd8a","Cd4","Trdc","Klrk1","Flt3","Cd1d1","Fcer2a","Cd19","Ighm","Ighd","Ccr7","Cd3e") #"C1qa",#,"Olig1"
	#names(FEATURES)<-c("Eryth.","PC","Mo+Ma","Gran.","act. T","Treg","CD8","CD4","gdT","NK","DC","Marg. B","Fol. B","B-cell","B-cell","B-cell","T-naive","T-cell") #"Macro",#,"Oligo"
	
	FEATURES=c("Cd3e","Ccr7","Rora","Cd4","Cd8a","Trdc","Foxp3","Klrb1c",
			   "Cd79a","Ighd","Fcer2a","Cd1d1","Aicda","Sdc1",
			   "Lyz2",
			   "Flt3","Cd68","S100a8","Hbb-bt")
	names(FEATURES)<-c("T-cell","T-naive","act. T","CD4","CD8","gdT","Treg","NK+NKT",
					   "B-cell","B-cell","Fol. B","Marg. B","Germ. B","PC",
					   "Myeloid",
					   "DC","Ma+Mo","Gran.","Eryth.")
	
	
	NAMEIDENT = "CellTypeFinal"
	Idents(GEX)<-NAMEIDENT
	
	Idents(GEX) <- factor(GEX$CellTypeFinal,levels=rev(c("naive T cell CD4+","act. T cell CD4+","T cell CD8+","gd T cell","Treg",
													 "NK + NK T cell",
													 "fol. zone B cell","mar. zone B cell","ger. center B cell","Plasma cell",
													 "DC","Macrophage","Monocyte","Granulocyte",
													 "Erythroid cell",
													 "rm-1","rm-2","rm-3","rm-4","rm-5")))
	IDENTS = unique(Idents(GEX))
	
	PLOT=DotPlot(GEX,features=FEATURES,assay="RNA",idents=IDENTS,dot.scale=10) + 
		theme_minimal(base_size=20) + 
		theme(axis.text.x=element_text(angle=90)) + 
		ggtitle("Gene markers")
	OUTNAME=paste0(SUBDIR,"/DotPlot_GeneMarkers_with_rm.pdf") ; pdf(OUTNAME,width=25,height=8) ; print(PLOT);dev.off()
	
	
	IDENTS<-IDENTS[IDENTS %ni% grep("^rm-",IDENTS,value=T)]
	
	PLOT=DotPlot(GEX,features=FEATURES,assay="RNA",idents=IDENTS,dot.scale=10) + 
		theme_minimal(base_size=20) + 
		theme(axis.text.x=element_text(angle=90)) + 
		ggtitle("Gene markers")
	OUTNAME=paste0(SUBDIR,"/DotPlot_GeneMarkers.pdf") ; pdf(OUTNAME,width=25,height=8) ; print(PLOT);dev.off()
	
	rm(FEATURES,IDENTS,PLOT,OUTNAME)
}




#### SAVE BEFORE THE SUBSET ####
if(F){
	GEX$mouse <- factor(GEX$orig.ident)
	levels(GEX$mouse) <-c("I18","I18","J119","J119","J119","J9","J9","J9","E21","E21","E72","E72","F21","F21") #"I18A","I18B","J119A","J119B","J119C","J9A","J9B","J9C","Z3.E21.LN"   "Z3.E21.mELT" "Z3.E72.LN"   "Z3.E72.mELT" "Z3.F21.LN"   "Z3.F21.mELT")
	saveRDS(GEX,"04.final.Rds")
}
#### ------- FINAL FILE LOAD and filter ------- ####

GEX<-readRDS("04.final.Rds")
GEX
# 32300 features across 79530 samples within 4 assays 
# Active assay: RNA (32285 features, 2000 variable features)
# 3 other assays present: Ab, AbTissues, AbCelltypes
# 3 dimensional reductions calculated: pca, harmony, umap
#table(GEX$immune_chains_clean4_selection, GEX$CellTypeFinal) # check that the FALSE is not in majority: we want to filter HTO in wrong celltype
#table(GEX$immune_chains_clean4_selection) #to remove: 6626		72904
# GEX <- subset(GEX,subset = immune_chains_clean4_selection )
# GEX
# 32300 features across 72904 samples within 4 assays
#79530-72904 = 6626 removed because of cleaning -> match previous tables

# once this selection is done -> do the dotplot again, do the cluster UMAPs 

#### TABLES #### 
SUBDIR="15.tables";dir.create(SUBDIR)
if(F){
	TABLE=table(GEX$CellTypeFinal, GEX$tissue_merged);TABLE
	write.xlsx(TABLE,paste0(SUBDIR,"/table_CellTypeFinal_X_tag_tissue_merged.xlsx"))

}

#### Welch’s t-test - cell type enrichment per tissues #### 
#plotting fold change (log10) against p value (−log10) based on beta-binomial regression
SUBDIR="15.Welch-t-test";dir.create(SUBDIR)
# NOW: removing Blood from the figures
# TODO : fix size of the points / make it comparable between figures
if(T){
	
	#for the legend: show.legend = T
	SHOW_LEGEND=F
	
	LAB_SIZE=5
	NUDGE_X=.1 #-.1
	NUDGE_Y=.4
	
	#all cells - mELT reference versus tissues ~ FOR EACH CELLTYPEFINAL
	OUTNAME_DATA="Welch_t_tests_CellType"
	if(T){
		
		#GEX$mouse
		GEX$tissue<-factor(GEX$tissue_merged)
		GEX$CellTypeFinal_Welch <- factor(GEX$CellTypeFinal)
		# # # # # # calculate percent of celltype per tissus for X mice - T-tests
		# n = count of cell per mouse / tissue / celltype
		# we want the percent of celltype in tissue / mouse
		COUNTS <- GEX@meta.data %>% filter( CellTypeFinal_Welch %ni% c("Erythroid cell",grep("^rm-",unique(GEX$CellTypeFinal_Welch),value=T) ) ) %>%
			group_by( mouse, tissue, CellTypeFinal_Welch ) %>% count() %>% 
			group_by( mouse,tissue ) %>% mutate( mouse.tissue.count=sum(n)) %>%
			group_by( tissue, CellTypeFinal_Welch ) %>% mutate( tissue.CellType.count=sum(n)) %>%
			mutate( mouse.count=length(unique(mouse)), percent.celltype.in.tissue.per.mouse = n / mouse.tissue.count ) %>% arrange(tissue)
		
		COUNTS=as.data.frame(COUNTS)
		OUTNAME=paste0(SUBDIR,"/table.",OUTNAME_DATA,".xlsx")
		write.xlsx(COUNTS,OUTNAME,overwrite=T)
		#COUNTS<-read.xlsx( paste0(SUBDIR,"/table_counts.xlsx"))
		
		# # # # # # Calculate T-tests
		TISSUES<-c("LN_inguinal","LN_lumbar","Spleen" ) # # "Blood", "CSF", unique( COUNTS$tissue )
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
		
		rm(D.mELT,D.other)
		
		TTEST=read.table(OUTNAME,h=T,sep="\t")
		TTEST$log10FC <- log10(TTEST$ratio)
		TTEST$log10Pval_m <- -log10(TTEST$Pvalue) # minus log10 : high values = low p-value
		write.xlsx(TTEST, paste0(OUTNAME,".xlsx"),overwrite=T)
		
		
		BONFNAME=grep("Padj_bonferroni",names(TTEST),value=T)
		TTEST$celltype2 <- ifelse( TTEST[,BONFNAME] < 0.05, paste0(TTEST$celltype," *"),TTEST$celltype )
		
		PLOT=ggplot(TTEST) + 
			geom_hline(yintercept=-log10(0.05),linetype="dashed",col="red") + 
			geom_hline(yintercept=-log10(0.05/BONFERRONI),linetype="dashed",col="orange") +
			geom_vline(xintercept=0,linetype="dashed",col="black") + 
			geom_point(aes(log10FC,log10Pval_m,size=abs(log10FC)*log10Pval_m ,col=celltype),show.legend=SHOW_LEGEND) + #to remove color: here
			geom_text_repel(aes(log10FC,log10Pval_m,label=celltype2,size=LAB_SIZE,col=celltype),
							 min.segment.length=0,nudge_x=NUDGE_X,nudge_y=NUDGE_Y,show.legend=F,force=100,max.overlaps=100) + 
			scale_x_continuous(limits = c(-1,1)) +
			theme_bw(base_size=22) + 
			ggtitle( label = paste("Celltypes % (reference: mELT)"),
					 subtitle = paste("Bonferroni correction on",BONFERRONI,"tests")) + 
			facet_wrap(~tissue) # + coord_equal()
		OUTNAME=paste0(SUBDIR,"/Figure_Welch.",OUTNAME_DATA,".pdf") ; pdf(OUTNAME,width=18,height=8) ; print(PLOT);dev.off()
		print(paste("done:",OUTNAME_DATA))
	}
	#all cells - mELT reference versus tissues ~ FOR EACH CELLTYPEFINAL -- general (B-cell, T-cells...)
	OUTNAME_DATA="Welch_t_tests_GeneralCellType"
	if(T){
		GEX$tissue<-factor(GEX$tissue_merged)
		GEX$CellTypeFinal_Welch <- factor(GEX$CellTypeFinal)
		#grep("B",levels(GEX$CellTypeFinal_Welch),value=T) ; grep("B",levels(GEX$CellTypeFinal_Welch),value=F)
		levels(GEX$CellTypeFinal_Welch) [ grep("B",levels(GEX$CellTypeFinal_Welch),value=F) ] <- "B cell"
		#"fol. zone B cell"   "ger. center B cell" "mar. zone B cell"
		levels(GEX$CellTypeFinal_Welch) [ levels(GEX$CellTypeFinal_Welch) %in% c("Plasma cell") ] <- "B cell"
		levels(GEX$CellTypeFinal_Welch) [ grep("T",levels(GEX$CellTypeFinal_Welch),value=F) ] <- "T cell"
		
		levels(GEX$CellTypeFinal_Welch) [ levels(GEX$CellTypeFinal_Welch) %in% c("Granulocyte","Macrophage","DC","Monocyte") ] <- "Myeloid cells" #(macrophages, monocytes, DC and granulocytes)
		
		levels(GEX$CellTypeFinal_Welch) [ levels(GEX$CellTypeFinal_Welch) %in% c("NK cell","NK T cell") ] <- "NK cell"
		#?dplyr::recode() # case_when() # recode_factor()
		
		# # # # # # calculate percent of celltype per tissus for X mice - T-tests
		# n = count of cell per mouse / tissue / celltype
		# we want the percent of celltype in tissue / mouse
		COUNTS <- GEX@meta.data %>% filter( CellTypeFinal_Welch %ni% c("Erythroid cell",grep("^rm-",unique(GEX$CellTypeFinal_Welch),value=T)) ) %>%
			group_by( mouse, tissue, CellTypeFinal_Welch ) %>% count() %>% 
			group_by( mouse,tissue ) %>% mutate( mouse.tissue.count=sum(n)) %>%
			group_by( tissue, CellTypeFinal_Welch ) %>% mutate( tissue.CellType.count=sum(n)) %>%
			mutate( mouse.count=length(unique(mouse)), percent.celltype.in.tissue.per.mouse = n / mouse.tissue.count ) %>% arrange(tissue)
		
		COUNTS=as.data.frame(COUNTS)
		OUTNAME=paste0(SUBDIR,"/table.",OUTNAME_DATA,".xlsx")
		write.xlsx(COUNTS,OUTNAME,overwrite=T)
		#COUNTS<-read.xlsx( paste0(SUBDIR,"/table_counts.xlsx"))
		
		# # # # # # Calculate T-tests
		TISSUES<-c("LN_inguinal","LN_lumbar","Spleen" ) # # "Blood","CSF", unique( COUNTS$tissue )
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
		write.xlsx(TTEST, paste0(OUTNAME,".xlsx"),overwrite=T)
		
		# LAB_SIZE=1
		# NUDGE_X=0 #-.1
		# NUDGE_Y=.1
		
		BONFNAME=grep("Padj_bonferroni",names(TTEST),value=T)
		TTEST$celltype2 <- ifelse( TTEST[,BONFNAME] < 0.05, paste0(TTEST$celltype," *"),TTEST$celltype )
		
		PLOT=ggplot(TTEST) + 
			geom_hline(yintercept=-log10(0.05),linetype="dashed",col="red") + 
			geom_hline(yintercept=-log10(0.05/BONFERRONI),linetype="dashed",col="orange") +
			geom_vline(xintercept=0,linetype="dashed",col="black") + 
			geom_point(aes(log10FC,log10Pval_m,size=abs(log10FC)*log10Pval_m ,col=celltype),show.legend=SHOW_LEGEND) + #to remove color: here
			geom_text_repel(aes(log10FC,log10Pval_m,label=celltype2,size=LAB_SIZE,col=celltype),min.segment.length=0,
							 nudge_x=NUDGE_X,nudge_y=NUDGE_Y,show.legend=F) + 
			scale_x_continuous(limits = c(-1,1)) +
			theme_bw(base_size=22) + 
			ggtitle( label = paste("Celltypes % (reference: mELT)"),
					 subtitle = paste("Bonferroni correction on",BONFERRONI,"tests")) + 
			facet_wrap(~tissue) # + coord_equal()
		OUTNAME=paste0(SUBDIR,"/Figure_Welch.",OUTNAME_DATA,".pdf") ; pdf(OUTNAME,width=18,height=8) ; print(PLOT);dev.off()
		
	}
	#         mELT reference versus tissues ~ Only B cell subsets
	OUTNAME_DATA="Welch_t_tests_Bcell"
	if(T){
		GEX$tissue<-factor(GEX$tissue_merged)
		GEX$CellTypeFinal_Welch <- factor(GEX$CellTypeFinal)
		
		# # # # # # calculate percent of celltype per tissus for X mice - T-tests
		# n = count of cell per mouse / tissue / celltype
		# we want the percent of celltype in tissue / mouse
		grep("B",levels(GEX$CellTypeFinal_Welch),value=T)
		COUNTS <- GEX@meta.data %>% filter( CellTypeFinal_Welch %in% c("fol. zone B cell","ger. center B cell","mar. zone B cell", "Plasma cell") ) %>%
			group_by( mouse, tissue, CellTypeFinal_Welch ) %>% count() %>% 
			group_by( mouse,tissue ) %>% mutate( mouse.tissue.count=sum(n)) %>%
			group_by( tissue, CellTypeFinal_Welch ) %>% mutate( tissue.CellType.count=sum(n)) %>%
			mutate( mouse.count=length(unique(mouse)), percent.celltype.in.tissue.per.mouse = n / mouse.tissue.count ) %>% arrange(tissue)
		
		COUNTS=as.data.frame(COUNTS)
		OUTNAME=paste0(SUBDIR,"/table.",OUTNAME_DATA,".xlsx")
		write.xlsx(COUNTS,OUTNAME,overwrite=T)
		#COUNTS<-read.xlsx( paste0(SUBDIR,"/table_counts.xlsx"))
		
		# # # # # # Calculate T-tests
		TISSUES<-c("LN_inguinal","LN_lumbar","Spleen" ) # "Blood", "CSF", unique( COUNTS$tissue )
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
		write.xlsx(TTEST, paste0(OUTNAME,".xlsx"),overwrite=T)
		
		# LAB_SIZE=1
		# NUDGE_X=0
		# NUDGE_Y=.1
		BONFNAME=grep("Padj_bonferroni",names(TTEST),value=T)
		TTEST$celltype2 <- ifelse( TTEST[,BONFNAME] < 0.05, paste0(TTEST$celltype," *"),TTEST$celltype )
		
		PLOT=ggplot(TTEST) + 
			geom_hline(yintercept=-log10(0.05),linetype="dashed",col="red") + 
			geom_hline(yintercept=-log10(0.05/BONFERRONI),linetype="dashed",col="orange") +
			geom_vline(xintercept=0,linetype="dashed",col="black") + 
			geom_point(aes(log10FC,log10Pval_m,size=abs(log10FC)*log10Pval_m ,col=celltype),show.legend=SHOW_LEGEND) + #to remove color: here
			geom_text_repel(aes(log10FC,log10Pval_m,label=celltype2,size=LAB_SIZE,col=celltype),force_pull=10,force=10,
							 min.segment.length=0,nudge_x=NUDGE_X,nudge_y=NUDGE_Y,show.legend=F) + 
			#scale_x_continuous(limits = c(-1,1)) +
			scale_x_continuous(limits = c(-1.2,1.2)) +
			theme_bw(base_size=22) + 
			ggtitle( label = paste("Celltypes % (reference: mELT)"),
					 subtitle = paste("Bonferroni correction on",BONFERRONI,"tests")) + 
			facet_wrap(~tissue) # + coord_equal()
		OUTNAME=paste0(SUBDIR,"/Figure_Welch.",OUTNAME_DATA,".pdf") ; pdf(OUTNAME,width=18,height=8) ; print(PLOT);dev.off()
	}
	#         mELT reference versus tissues ~ Only T cell subsets
	OUTNAME_DATA="Welch_t_tests_Tcell"
	if(T){
		GEX$tissue<-factor(GEX$tissue_merged)
		GEX$CellTypeFinal_Welch <- factor(GEX$CellTypeFinal)
		
		# # # # # # calculate percent of celltype per tissus for X mice - T-tests
		# n = count of cell per mouse / tissue / celltype
		# we want the percent of celltype in tissue / mouse
		grep("T",levels(GEX$CellTypeFinal_Welch),value=T)
		COUNTS <- GEX@meta.data %>% filter( CellTypeFinal_Welch %in% c("act. T cell CD4+","gd T cell","naive T cell CD4+","T cell CD8+","Treg") ) %>%
			group_by( mouse, tissue, CellTypeFinal_Welch ) %>% count() %>% 
			group_by( mouse,tissue ) %>% mutate( mouse.tissue.count=sum(n)) %>%
			group_by( tissue, CellTypeFinal_Welch ) %>% mutate( tissue.CellType.count=sum(n)) %>%
			mutate( mouse.count=length(unique(mouse)), percent.celltype.in.tissue.per.mouse = n / mouse.tissue.count ) %>% arrange(tissue)
		
		COUNTS=as.data.frame(COUNTS)
		OUTNAME=paste0(SUBDIR,"/table.",OUTNAME_DATA,".xlsx")
		write.xlsx(COUNTS,OUTNAME,overwrite=T)
		#COUNTS<-read.xlsx( paste0(SUBDIR,"/table_counts.xlsx"))
		
		# # # # # # Calculate T-tests
		TISSUES<-c("LN_inguinal","LN_lumbar","Spleen" ) # "Blood", "CSF", unique( COUNTS$tissue )
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
				rm(RATIO,PValue,Statistic)
			}
		}
		rm(TISSUE,TISSUES,CELLTYPE,CELLTYPES,TREATMENT,REFERENCE)
		
		TTEST=read.table(OUTNAME,h=T,sep="\t")
		TTEST$log10FC <- log10(TTEST$ratio)
		TTEST$log10Pval_m <- -log10(TTEST$Pvalue) # minus log10 : high values = low p-value
		write.xlsx(TTEST, paste0(OUTNAME,".xlsx"),overwrite=T)
		
		# LAB_SIZE=1
		# NUDGE_X=0 #-.1
		# NUDGE_Y=.1
		BONFNAME=grep("Padj_bonferroni",names(TTEST),value=T)
		TTEST$celltype2 <- ifelse( TTEST[,BONFNAME] < 0.05, paste0(TTEST$celltype," *"),TTEST$celltype )
		
		PLOT=ggplot(TTEST) + 
			geom_hline(yintercept=-log10(0.05),linetype="dashed",col="red") + 
			geom_hline(yintercept=-log10(0.05/BONFERRONI),linetype="dashed",col="orange") +
			geom_vline(xintercept=0,linetype="dashed",col="black") + 
			geom_point(aes(log10FC,log10Pval_m,size=abs(log10FC)*log10Pval_m ,col=celltype),show.legend=SHOW_LEGEND) + #to remove color: here
			geom_text_repel(aes(log10FC,log10Pval_m,label=celltype2,size=LAB_SIZE,col=celltype),min.segment.length=0,
							 nudge_x=NUDGE_X,nudge_y=NUDGE_Y,show.legend=F) + 
			scale_x_continuous(limits = c(-1,1)) +
			theme_bw(base_size=22) + 
			ggtitle( label = paste("Celltypes % (reference: mELT)"),
					 subtitle = paste("Bonferroni correction on",BONFERRONI,"tests")) + 
			facet_wrap(~tissue) # + coord_equal()
		OUTNAME=paste0(SUBDIR,"/Figure_Welch.",OUTNAME_DATA,".pdf") ; pdf(OUTNAME,width=18,height=8) ; print(PLOT);dev.off()
		
	}
	#         mELT reference versus tissues ~ Only Myeloid cell subsets
	OUTNAME_DATA="Welch_t_tests_Myeloid_cells"
	if(T){
		GEX$tissue<-factor(GEX$tissue_merged)
		GEX$CellTypeFinal_Welch <- factor(GEX$CellTypeFinal)
		
		# # # # # # calculate percent of celltype per tissus for X mice - T-tests
		# n = count of cell per mouse / tissue / celltype
		# we want the percent of celltype in tissue / mouse
		
		COUNTS <- GEX@meta.data %>% filter( CellTypeFinal_Welch %in% c("Granulocyte","Macrophage","DC","Monocyte") ) %>%
			group_by( mouse, tissue, CellTypeFinal_Welch ) %>% count() %>% 
			group_by( mouse,tissue ) %>% mutate( mouse.tissue.count=sum(n)) %>%
			group_by( tissue, CellTypeFinal_Welch ) %>% mutate( tissue.CellType.count=sum(n)) %>%
			mutate( mouse.count=length(unique(mouse)), percent.celltype.in.tissue.per.mouse = n / mouse.tissue.count ) %>% arrange(tissue)
		
		COUNTS=as.data.frame(COUNTS)
		OUTNAME=paste0(SUBDIR,"/table.",OUTNAME_DATA,".xlsx")
		write.xlsx(COUNTS,OUTNAME,overwrite=T)
		#COUNTS<-read.xlsx( paste0(SUBDIR,"/table_counts.xlsx"))
		
		# # # # # # Calculate T-tests
		TISSUES<-c("LN_inguinal","LN_lumbar","Spleen" ) # "Blood", "CSF", unique( COUNTS$tissue )
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
				rm(RATIO,PValue,Statistic)
			}
		}
		rm(TISSUE,TISSUES,CELLTYPE,CELLTYPES,TREATMENT,REFERENCE)
		
		TTEST=read.table(OUTNAME,h=T,sep="\t")
		TTEST$log10FC <- log10(TTEST$ratio)
		TTEST$log10Pval_m <- -log10(TTEST$Pvalue) # minus log10 : high values = low p-value
		write.xlsx(TTEST, paste0(OUTNAME,".xlsx"),overwrite=T)
		
		# LAB_SIZE=1
		# NUDGE_X=0 #-.1
		# NUDGE_Y=.1
		BONFNAME=grep("Padj_bonferroni",names(TTEST),value=T)
		TTEST$celltype2 <- ifelse( TTEST[,BONFNAME] < 0.05, paste0(TTEST$celltype," *"),TTEST$celltype )
		
		PLOT=ggplot(TTEST) + 
			geom_hline(yintercept=-log10(0.05),linetype="dashed",col="red") + 
			geom_hline(yintercept=-log10(0.05/BONFERRONI),linetype="dashed",col="orange") +
			geom_vline(xintercept=0,linetype="dashed",col="black") + 
			geom_point(aes(log10FC,log10Pval_m,size=abs(log10FC)*log10Pval_m ,col=celltype),show.legend=SHOW_LEGEND) + #to remove color: here
			geom_text_repel(aes(log10FC,log10Pval_m,label=celltype2,size=LAB_SIZE,col=celltype),min.segment.length=0,
							 nudge_x=NUDGE_X,nudge_y=NUDGE_Y,show.legend=F) + 
			scale_x_continuous(limits = c(-1,1)) +
			theme_bw(base_size=22) + 
			ggtitle( label = paste("Celltypes % (reference: mELT)"),
					 subtitle = paste("Bonferroni correction on",BONFERRONI,"tests")) + 
			facet_wrap(~tissue) # + coord_equal()
		OUTNAME=paste0(SUBDIR,"/Figure_Welch.",OUTNAME_DATA,".pdf") ; pdf(OUTNAME,width=18,height=8) ; print(PLOT);dev.off()
		
	}
	
	rm(PLOT,TTEST,BONFNAME,NUDGE_X,NUDGE_Y,LAB_SIZE,OUTNAME, BONFERRONI)
}
# code to plot with the Blood
if(F){
	#all cells - mELT reference versus tissues ~ FOR EACH CELLTYPEFINAL
	OUTNAME_DATA="Welch_t_tests_CellType"
	if(F){
		#GEX$mouse
		GEX$tissue<-factor(GEX$tissue_merged)
		GEX$CellTypeFinal_Welch <- factor(GEX$CellTypeFinal)
		# # # # # # calculate percent of celltype per tissus for X mice - T-tests
		# n = count of cell per mouse / tissue / celltype
		# we want the percent of celltype in tissue / mouse
		COUNTS <- GEX@meta.data %>% filter( CellTypeFinal_Welch %ni% c("Erythroid cell",grep("^rm-",unique(GEX$CellTypeFinal_Welch),value=T) ) ) %>%
			group_by( mouse, tissue, CellTypeFinal_Welch ) %>% count() %>% 
			group_by( mouse,tissue ) %>% mutate( mouse.tissue.count=sum(n)) %>%
			group_by( tissue, CellTypeFinal_Welch ) %>% mutate( tissue.CellType.count=sum(n)) %>%
			mutate( mouse.count=length(unique(mouse)), percent.celltype.in.tissue.per.mouse = n / mouse.tissue.count ) %>% arrange(tissue)
		
		COUNTS=as.data.frame(COUNTS)
		OUTNAME=paste0(SUBDIR,"/table.",OUTNAME_DATA,".xlsx")
		write.xlsx(COUNTS,OUTNAME)
		#COUNTS<-read.xlsx( paste0(SUBDIR,"/table_counts.xlsx"))
		
		# # # # # # Calculate T-tests
		TISSUES<-c("Blood","LN_inguinal","LN_lumbar","Spleen" ) # # "CSF", unique( COUNTS$tissue )
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
		
		rm(D.mELT,D.other)
		
		TTEST=read.table(OUTNAME,h=T,sep="\t")
		TTEST$log10FC <- log10(TTEST$ratio)
		TTEST$log10Pval_m <- -log10(TTEST$Pvalue) # minus log10 : high values = low p-value
		write.xlsx(TTEST, paste0(OUTNAME,".xlsx"))
		
		LAB_SIZE=.4
		NUDGE_X=.1 #-.1
		NUDGE_Y=.4
		
		BONFNAME=grep("Padj_bonferroni",names(TTEST),value=T)
		TTEST$celltype2 <- ifelse( TTEST[,BONFNAME] < 0.05, paste0(TTEST$celltype," *"),TTEST$celltype )
		
		PLOT=ggplot(TTEST) + 
			geom_hline(yintercept=-log10(0.05),linetype="dashed",col="red") + 
			geom_hline(yintercept=-log10(0.05/BONFERRONI),linetype="dashed",col="orange") +
			geom_vline(xintercept=0,linetype="dashed",col="black") + 
			geom_point(aes(log10FC,log10Pval_m,size=abs(log10FC)*log10Pval_m ,col=celltype),show.legend=F) + #to remove color: here
			geom_label_repel(aes(log10FC,log10Pval_m,label=celltype2,size=LAB_SIZE,col=celltype),
							 min.segment.length=0,nudge_x=NUDGE_X,nudge_y=NUDGE_Y,show.legend=F,force=100,max.overlaps=100) + 
			scale_x_continuous(limits = c(-1,1)) +
			theme_bw(base_size=22) + 
			ggtitle( label = paste("Celltypes % (reference: mELT)"),
					 subtitle = paste("Bonferroni correction on",BONFERRONI,"tests")) + 
			facet_wrap(~tissue)
		OUTNAME=paste0(SUBDIR,"/Figure_Welch.",OUTNAME_DATA,".pdf") ; pdf(OUTNAME,width=10,height=10) ; print(PLOT);dev.off()
		print(paste("done:",OUTNAME_DATA))
	}
	#all cells - mELT reference versus tissues ~ FOR EACH CELLTYPEFINAL -- general (B-cell, T-cells...)
	OUTNAME_DATA="Welch_t_tests_GeneralCellType"
	if(F){
		GEX$tissue<-factor(GEX$tissue_merged)
		GEX$CellTypeFinal_Welch <- factor(GEX$CellTypeFinal)
		#grep("B",levels(GEX$CellTypeFinal_Welch),value=T) ; grep("B",levels(GEX$CellTypeFinal_Welch),value=F)
		levels(GEX$CellTypeFinal_Welch) [ grep("B",levels(GEX$CellTypeFinal_Welch),value=F) ] <- "B cell"
		#"fol. zone B cell"   "ger. center B cell" "mar. zone B cell"
		levels(GEX$CellTypeFinal_Welch) [ levels(GEX$CellTypeFinal_Welch) %in% c("Plasma cell") ] <- "B cell"
		levels(GEX$CellTypeFinal_Welch) [ grep("T",levels(GEX$CellTypeFinal_Welch),value=F) ] <- "T cell"
		
		levels(GEX$CellTypeFinal_Welch) [ levels(GEX$CellTypeFinal_Welch) %in% c("Granulocyte","Macrophage","DC","Monocyte") ] <- "Myeloid cells" #(macrophages, monocytes, DC and granulocytes)
		
		levels(GEX$CellTypeFinal_Welch) [ levels(GEX$CellTypeFinal_Welch) %in% c("NK cell","NK T cell") ] <- "NK cell"
		#?dplyr::recode() # case_when() # recode_factor()
		
		# # # # # # calculate percent of celltype per tissus for X mice - T-tests
		# n = count of cell per mouse / tissue / celltype
		# we want the percent of celltype in tissue / mouse
		COUNTS <- GEX@meta.data %>% filter( CellTypeFinal_Welch %ni% c("Erythroid cell",grep("^rm-",unique(GEX$CellTypeFinal_Welch),value=T)) ) %>%
			group_by( mouse, tissue, CellTypeFinal_Welch ) %>% count() %>% 
			group_by( mouse,tissue ) %>% mutate( mouse.tissue.count=sum(n)) %>%
			group_by( tissue, CellTypeFinal_Welch ) %>% mutate( tissue.CellType.count=sum(n)) %>%
			mutate( mouse.count=length(unique(mouse)), percent.celltype.in.tissue.per.mouse = n / mouse.tissue.count ) %>% arrange(tissue)
		
		COUNTS=as.data.frame(COUNTS)
		OUTNAME=paste0(SUBDIR,"/table.",OUTNAME_DATA,".xlsx")
		write.xlsx(COUNTS,OUTNAME)
		#COUNTS<-read.xlsx( paste0(SUBDIR,"/table_counts.xlsx"))
		
		# # # # # # Calculate T-tests
		TISSUES<-c("Blood","LN_inguinal","LN_lumbar","Spleen" ) # # "CSF", unique( COUNTS$tissue )
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
		GEX$tissue<-factor(GEX$tissue_merged)
		GEX$CellTypeFinal_Welch <- factor(GEX$CellTypeFinal)
		
		# # # # # # calculate percent of celltype per tissus for X mice - T-tests
		# n = count of cell per mouse / tissue / celltype
		# we want the percent of celltype in tissue / mouse
		grep("B",levels(GEX$CellTypeFinal_Welch),value=T)
		COUNTS <- GEX@meta.data %>% filter( CellTypeFinal_Welch %in% c("fol. zone B cell","ger. center B cell","mar. zone B cell", "Plasma cell") ) %>%
			group_by( mouse, tissue, CellTypeFinal_Welch ) %>% count() %>% 
			group_by( mouse,tissue ) %>% mutate( mouse.tissue.count=sum(n)) %>%
			group_by( tissue, CellTypeFinal_Welch ) %>% mutate( tissue.CellType.count=sum(n)) %>%
			mutate( mouse.count=length(unique(mouse)), percent.celltype.in.tissue.per.mouse = n / mouse.tissue.count ) %>% arrange(tissue)
		
		COUNTS=as.data.frame(COUNTS)
		OUTNAME=paste0(SUBDIR,"/table.",OUTNAME_DATA,".xlsx")
		write.xlsx(COUNTS,OUTNAME)
		#COUNTS<-read.xlsx( paste0(SUBDIR,"/table_counts.xlsx"))
		
		# # # # # # Calculate T-tests
		TISSUES<-c("Blood","LN_inguinal","LN_lumbar","Spleen" ) # # "CSF", unique( COUNTS$tissue )
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
		OUTNAME=paste0(SUBDIR,"/Figure_Welch.",OUTNAME_DATA,".pdf") ; pdf(OUTNAME,width=12,height=12) ; print(PLOT);dev.off()
	}
	#         mELT reference versus tissues ~ Only T cell subsets
	OUTNAME_DATA="Welch_t_tests_Tcell"
	if(T){
		GEX$tissue<-factor(GEX$tissue_merged)
		GEX$CellTypeFinal_Welch <- factor(GEX$CellTypeFinal)
		
		# # # # # # calculate percent of celltype per tissus for X mice - T-tests
		# n = count of cell per mouse / tissue / celltype
		# we want the percent of celltype in tissue / mouse
		grep("T",levels(GEX$CellTypeFinal_Welch),value=T)
		COUNTS <- GEX@meta.data %>% filter( CellTypeFinal_Welch %in% c("act. T cell CD4+","gd T cell","naive T cell CD4+","T cell CD8+","Treg") ) %>%
			group_by( mouse, tissue, CellTypeFinal_Welch ) %>% count() %>% 
			group_by( mouse,tissue ) %>% mutate( mouse.tissue.count=sum(n)) %>%
			group_by( tissue, CellTypeFinal_Welch ) %>% mutate( tissue.CellType.count=sum(n)) %>%
			mutate( mouse.count=length(unique(mouse)), percent.celltype.in.tissue.per.mouse = n / mouse.tissue.count ) %>% arrange(tissue)
		
		COUNTS=as.data.frame(COUNTS)
		OUTNAME=paste0(SUBDIR,"/table.",OUTNAME_DATA,".xlsx")
		write.xlsx(COUNTS,OUTNAME)
		#COUNTS<-read.xlsx( paste0(SUBDIR,"/table_counts.xlsx"))
		
		# # # # # # Calculate T-tests
		TISSUES<-c("Blood","LN_inguinal","LN_lumbar","Spleen" ) # # "CSF", unique( COUNTS$tissue )
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
				rm(RATIO,PValue,BONFERRONI,Statistic)
			}
		}
		rm(TISSUE,TISSUES,CELLTYPE,CELLTYPES,TREATMENT,REFERENCE,OUTNAME)
		
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
	OUTNAME_DATA="Welch_t_tests_Myeloid_cells"
	if(T){
		GEX$tissue<-factor(GEX$tissue_merged)
		GEX$CellTypeFinal_Welch <- factor(GEX$CellTypeFinal)
		
		# # # # # # calculate percent of celltype per tissus for X mice - T-tests
		# n = count of cell per mouse / tissue / celltype
		# we want the percent of celltype in tissue / mouse
	
		COUNTS <- GEX@meta.data %>% filter( CellTypeFinal_Welch %in% c("Granulocyte","Macrophage","DC","Monocyte") ) %>%
			group_by( mouse, tissue, CellTypeFinal_Welch ) %>% count() %>% 
			group_by( mouse,tissue ) %>% mutate( mouse.tissue.count=sum(n)) %>%
			group_by( tissue, CellTypeFinal_Welch ) %>% mutate( tissue.CellType.count=sum(n)) %>%
			mutate( mouse.count=length(unique(mouse)), percent.celltype.in.tissue.per.mouse = n / mouse.tissue.count ) %>% arrange(tissue)
		
		COUNTS=as.data.frame(COUNTS)
		OUTNAME=paste0(SUBDIR,"/table.",OUTNAME_DATA,".xlsx")
		write.xlsx(COUNTS,OUTNAME)
		#COUNTS<-read.xlsx( paste0(SUBDIR,"/table_counts.xlsx"))
		
		# # # # # # Calculate T-tests
		TISSUES<-c("Blood","LN_inguinal","LN_lumbar","Spleen" ) # # "CSF", unique( COUNTS$tissue )
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
				rm(RATIO,PValue,Statistic)
			}
		}
		rm(TISSUE,TISSUES,CELLTYPE,CELLTYPES,TREATMENT,REFERENCE)
		
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
}
# This is old
#all cells - mELT reference versus tissues
TAG="LNing+LNlumb+Spleen"
if(FALSE){
	GEX$mouse <- factor(GEX$orig.ident)
	levels(GEX$mouse) <-c("I18","I18",
						  "J119","J119","J119",
						  "J9","J9","J9",
						  "E21","E21",
						  "E72","E72",
						  "F21","F21") #"I18A","I18B","J119A","J119B","J119C","J9A","J9B","J9C","Z3.E21.LN"   "Z3.E21.mELT" "Z3.E72.LN"   "Z3.E72.mELT" "Z3.F21.LN"   "Z3.F21.mELT")
	GEX$tissue<-factor(GEX$tissue_merged)
	
	#change the specific tissues to the same: effectively groupping them
	table(GEX$tissue)
	levels(GEX$tissue) <- c("Blood","CSF","LNing+LNlumb+Spleen","LNing+LNlumb+Spleen","mELT","none","LNing+LNlumb+Spleen") 
	
	# # # # # # calculate percent of celltype per tissus for X mice - T-tests
	# n = count of cell per mouse / tissue / celltype
	# we want the percent of celltype in tissue / mouse
	COUNTS <- GEX@meta.data %>% 
		group_by( mouse, tissue, CellTypeFinal ) %>% count() %>% 
		group_by( mouse,tissue ) %>% mutate( mouse.tissue.count=sum(n)) %>%
		group_by( tissue, CellTypeFinal ) %>% mutate( tissue.CellType.count=sum(n)) %>%
		mutate( mouse.count=length(unique(mouse)), percent.celltype.in.tissue.per.mouse = n / mouse.tissue.count ) %>% arrange(tissue)
	
	COUNTS=as.data.frame(COUNTS)
	write.xlsx(COUNTS, paste0(SUBDIR,"/table_counts.",TAG,".xlsx"))
	COUNTS<-read.xlsx( paste0(SUBDIR,"/table_counts.",TAG,".xlsx"))
	
	# # # # # # Calculate T-tests
	TISSUES<-c("LNing+LNlumb+Spleen") #unique( COUNTS$tissue ) #"Blood","CSF",
	CELLTYPES<-unique(COUNTS$CellTypeFinal)
	CELLTYPES<-CELLTYPES[ ! CELLTYPES %in% grep("^rm-",CELLTYPES,value=T) ]
	
	BONFERRONI=length(CELLTYPES) #compared to mELT 
	
	cat(c("target","tissue","celltype","treat_count","ref_count","ratio","Pvalue",paste0("Padj_bonferroni",BONFERRONI),"Tscore\n"),file=paste0(SUBDIR,"/table_t-test",TAG,".txt"),append=F,sep="\t")
	for(CELLTYPE in CELLTYPES){
		for(TISSUE in TISSUES){
			D.mELT <- COUNTS[ COUNTS$tissue == "mELT" & COUNTS$CellTypeFinal == CELLTYPE, ]
			D.other <- COUNTS[ COUNTS$tissue == TISSUE & COUNTS$CellTypeFinal == CELLTYPE, ]
			REFERENCE = D.mELT$percent.celltype.in.tissue.per.mouse
			TREATMENT = D.other$percent.celltype.in.tissue.per.mouse
			RATIO = mean(TREATMENT) / mean(REFERENCE)
			
			SHA_TEST = shapiro.test( c(REFERENCE,TREATMENT) )
			#plot(density( c(REFERENCE,TREATMENT) ))
			#library("ggpubr")
			#ggqqplot( c(REFERENCE,TREATMENT) )
			
			#at least X observations ??
			
			DOTEST=TRUE
			if( length(REFERENCE)<2 | length(TREATMENT)<2 ){
				PValue=1
				Statistic="< 2 obs."
				DOTEST=F
			}
			if( SHA_TEST$p.value<0.05 ){
				PValue=1
				Statistic="not norm. distri."
				DOTEST=F
			}
			if(DOTEST){
				TEST = t.test(REFERENCE,TREATMENT) 
				PValue=TEST$p.value
				Statistic=TEST$statistic
				rm(TEST)
			}
			#print(paste("mELT",TISSUE,CELLTYPE,PValue,Statistic,sep=";"))
			cat( c("mELT",TISSUE,CELLTYPE,length(TREATMENT),length(REFERENCE),RATIO,PValue, PValue*BONFERRONI, paste0(Statistic,"\n")),file=paste0(SUBDIR,"/table_t-test",TAG,".txt"),append=T,sep="\t")
		}}
	
	TTEST=read.table(paste0(SUBDIR,"/table_t-test",TAG,".txt"),h=T,sep="\t")
	TTEST$log10FC <- log10(TTEST$ratio)
	TTEST$log10Pval_m <- -log10(TTEST$Pvalue) # minus log10 : high values = low p-value
	
	write.xlsx(TTEST, paste0(SUBDIR,"/table_figure_Welch_t_tests",TAG,".xlsx"))
	
	TTEST<-subset(TTEST, tissue != "CSF")
	LAB_SIZE=1
	NUDGE_X=0 #-.1
	NUDGE_Y=.1
	
	PLOT=ggplot(TTEST) + 
		geom_hline(yintercept=-log10(0.05),linetype="dashed",col="red") + 
		geom_hline(yintercept=-log10(0.05/BONFERRONI),linetype="dashed",col="orange") + 
		geom_vline(xintercept=0,linetype="dashed",col="black") + 
		geom_point(aes(log10FC,log10Pval_m,size=abs(log10FC)*log10Pval_m ),show.legend=F) +
		geom_label_repel(aes(log10FC,log10Pval_m,label=paste0(celltype),size=LAB_SIZE),min.segment.length=0,nudge_x=NUDGE_X,nudge_y=NUDGE_Y,show.legend=F) + 
		theme_bw(base_size=22) + 
		ggtitle( label = paste("Celltypes % (reference: mELT)"),
				 subtitle = paste("Bonferroni correction on",BONFERRONI,"tests")) + 
		facet_wrap(~tissue) #+ coord_fixed() 
	OUTNAME=paste0(SUBDIR,"/Figure_facet.",TAG,".pdf") ; pdf(OUTNAME,width=8,height=8) ; print(PLOT);dev.off()
	
	PLOT=ggplot(TTEST) + 
		geom_hline(yintercept=-log10(0.05),linetype="dashed",col="red") + 
		geom_hline(yintercept=-log10(0.05/BONFERRONI),linetype="dashed",col="orange") + 
		geom_vline(xintercept=0,linetype="dashed",col="black") + 
		geom_point(aes(log10FC,log10Pval_m,size=abs(log10FC)*log10Pval_m ,col=celltype),show.legend=F) +
		geom_label_repel(aes(log10FC,log10Pval_m,label=paste0(celltype),size=LAB_SIZE,col=celltype),min.segment.length=0,nudge_x=NUDGE_X,nudge_y=NUDGE_Y,show.legend=F) + 
		theme_bw(base_size=22) + 
		ggtitle( label = paste("Celltypes % (reference: mELT)"),
				 subtitle = paste("Bonferroni correction on",BONFERRONI,"tests")) + 
		facet_wrap(~tissue)
	OUTNAME=paste0(SUBDIR,"/Figure_facet_color.",TAG,".pdf") ; pdf(OUTNAME,width=8,height=8) ; print(PLOT);dev.off()
	
}



#### do the DEG: depending on the comparison we want to do ####
if(F){
	
	############ selection for the DEG
	SUBDIR="16.DEG";dir.create(SUBDIR)
	GEX$TEST <- paste(GEX$CellTypeFinal, GEX$tissue_merged, sep="-") ; table(GEX$TEST) #NOT tag_tissue_thr1.5 any more
	TISSUES<-c("mELT","CSF","LN_lumbar","LN_inguinal","Spleen","Blood") #unique(GEX$tag_tissue_thr1.5)
	CELLTYPES <- unique(GEX$CellTypeFinal)
	CELLTYPES<-CELLTYPES[ ! CELLTYPES %in% c( grep("rm-",CELLTYPES,value=T),"Erytroid cell")]
	
	#TABLE<-as.data.frame(table(GEX$TEST))
	#TABLE[ TABLE$Var1 %in% grep("Monocyte",TABLE$Var1,value=T), ]
	
	############ selection for the DEG - mELT versus SELECTION
	SELECTION="LNl+LNi+Spleen"
	SUBDIR="17.DEG_mELT_vs_LNl+LNi+Spleen";dir.create(SUBDIR)
	GEX$tissue_group = GEX$tissue_merged
	GEX$tissue_group <- gsub("LN_inguinal",SELECTION,GEX$tissue_group)
	GEX$tissue_group <- gsub("LN_lumbar",SELECTION,GEX$tissue_group)
	GEX$tissue_group <- gsub("Spleen",SELECTION,GEX$tissue_group)
	GEX$TEST <- paste(GEX$CellTypeFinal, GEX$tissue_group, sep="-") ; table(GEX$TEST)
	TISSUES<-c("mELT",SELECTION) #unique(GEX$tag_tissue_thr1.5)
	CELLTYPES <- unique(GEX$CellTypeFinal)
	CELLTYPES<-CELLTYPES[ ! CELLTYPES %in% c( grep("^rm-",CELLTYPES,value=T),"Erytroid cell")]
	
	
	
	Idents(GEX)<-"TEST"
	# TABLE=table(GEX$TEST);TABLE
	# write.xlsx(TABLE,"table_annot+tissue+treatment.xlsx")
	#only once:
	avgexp <- AverageExpression(GEX,verbose=T,assays="RNA")
	avgexp <- as.data.frame(avgexp$RNA)
	avgexp$gene <- rownames(avgexp)
	
	for( TISSUE1 in TISSUES ) {
		for( TISSUE2 in TISSUES[ 2:length(TISSUES) ] ) {
			if(TISSUE1 == TISSUE2){next}
			for( CELLTYPE in CELLTYPES ) {
				TREAT = paste(CELLTYPE, TISSUE1, sep="-")
				CONTROL=paste(CELLTYPE, TISSUE2, sep="-")
				OUTNAME=paste0(SUBDIR,"/Markers_",TREAT,"_vs_",CONTROL,".xlsx")
				if(file.exists(OUTNAME)){print(paste("file already exists:",OUTNAME));next} #the comparizon should be the same, make some 15 th decimal difference in numbers
				if( ! sum( GEX$TEST == TREAT) >= 50){print(paste("not enough cells for treatment:",TREAT));next}
				if( ! sum( GEX$TEST == CONTROL) >= 50){print(paste("not enough cells for control:",CONTROL));next}
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
						IGNORE=c(grep(pattern="^Rp[ls]",x=Markers2_Expr$gene,value=T),grep(pattern="^mt-",x=Markers2_Expr$gene,value=T) )
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
	}
	
	# system( paste0("mkdir -p ",SUBDIR,"/gt0 ; mv ",SUBDIR,"/*gt0*.txt ",SUBDIR,"/gt0"))
	# system(paste0("mkdir -p ",SUBDIR,"/ref-LN_lumbar ",SUBDIR,"/ref-LN_inguinal ",SUBDIR,"/ref-Spleen ",SUBDIR,"/ref-Blood"))
	# #system(paste0("mv ",SUBDIR,"/*CSF.xlsx ",SUBDIR,"/ref-CSF"))
	# system(paste0("mv ",SUBDIR,"/*LN_lumbar.xlsx ",SUBDIR,"/ref-LN_lumbar"))
	# system(paste0("mv ",SUBDIR,"/*LN_inguinal.xlsx ",SUBDIR,"/ref-LN_inguinal"))
	# system(paste0("mv ",SUBDIR,"/*Spleen.xlsx ",SUBDIR,"/ref-Spleen"))
	# system(paste0("mv ",SUBDIR,"/*Blood.xlsx ",SUBDIR,"/ref-Blood"))
	
	rm(TISSUE1,TISSUE2,CELLTYPE,TREAT,CONTROL,avgexp,CELLTYPES,SELECTION,TAG)
}
#### do the DEG: gene marker in cluster versus all rest ####
if(F){
	############ selection for the DEG
	SUBDIR="16.DEG.specific";dir.create(SUBDIR)
	
	GEX_SUB<-subset(GEX,  CellTypeFinal %ni% grep("^rm-",unique(GEX$CellTypeFinal),value=T) & 
						tissue_merged %ni% "none")
	
	GEX_SUB$TEST <- paste(GEX_SUB$CellTypeFinal, GEX_SUB$tissue_merged, sep="-") ; table(GEX_SUB$TEST) #NOT tag_tissue_thr1.5 any more
	Idents(GEX_SUB)<-"TEST"
	AllMarkers <- FindAllMarkers(GEX_SUB)
	OUTNAME=paste0("/home/ga94rac/LRZ Sync+Share/AG_Lehmann-Horn/Gildas_storage/2021.02.16.CITEseq/experiment_Merged_Tissues+mELT_LN","/AllMarkers.xlsx")
	write.xlsx(AllMarkers,OUTNAME)
}

#### check WTF we have so many Ribosomal genes in our DEG ####
if(F){
	
	# Ribosome: mRNA from DNA
	# enzyme of 4 rRNA + ~100 Ribo proteins
	# increased ribo biogenesis associated with increased proliferative activity.
	
	options(max.print=250)
	
	#RiboGenes1 = grep("^Rp",rownames(GEX),value=T)
	RiboGenes = grep("^Rp[sl][[:digit:]]",rownames(GEX),value=T)
	MitoGenes = grep("^[Mm][Tt]-",rownames(GEX),value=T)
	
	#length(RiboGenes1) # 144 #this include some extra genes: not only ribo
	#length(RiboGenes) # 97 #This seem better
	#RiboGenes1 [ RiboGenes1 %ni% RiboGenes ]
	#paste(RiboGenes1,sep=",",collapse=",")
	#paste(RiboGenes,sep=",",collapse=",")
	
	length(RiboGenes) # 97
	length(MitoGenes) # 13
	
	
	Idents(GEX)<-"RNA_snn_res.2"
	
	PLOT<-DoHeatmap(object = GEX, features = "Foxp3")
	DimHeatmap(object = GEX, reduction = "pca", cells = 200)
	
}



#### VLN plot check some genes VIOLIN plot ####
if(F){
	SUBDIR="18.Vln";dir.create(SUBDIR)
	#tissue_merged %in% c("CSF","mELT") & 
	DF<-read.xlsx("/home/ga94rac/WORK/RSTUDIO/210511.MergedRuns/Gene_Selection2.xlsx")
	#LIST=c("Mbp","S100a8","S100a9")
	LIST=na.omit( c( DF[,1],DF[,2],DF[,3],DF[,4],DF[,5]) )
	
	# add some more :
	LIST<-c(LIST, "Mzb1","Cd86") #grep("Cd86",rownames(GEX), value=T)
	LIST<-c("Cd86","Aicda","Cr2","Ccl5","Ccr7","Il7r","Il6ra","Cxcr6","Il2rb",
	"Tnfrsf18","Tnfrsf4","Tnfrsf9","Tnfsf11","Tgfb1",
	"Ifngr1","Igfbp7","Tgfbr2","Il21r","Cxcr4","Ifngr2","Ccr6")
	
	#write.table(LIST,"Gene_Selection2.txt",quote=F,row.names=F,col.names=F)
	
	GEX_SUB<-subset(GEX, CellTypeFinal %ni% grep("^rm-",unique(GEX$CellTypeFinal),value=T) & 
						tissue_merged %ni% c("Blood", "none","CSF") )
	#table(GEX_SUB$tissue_merged)
	for( GENE in LIST){
		PLOT<-VlnPlot(GEX_SUB,sort=T,log=F,features=GENE,combine=T,group.by="CellTypeFinal",split.by="tissue_merged",pt.size=0)
		OUTNAME=paste0(SUBDIR,"/",GENE,".pdf");pdf(OUTNAME,width=30,height=5);print(PLOT);dev.off();print(OUTNAME)
	}
	#table( as.matrix(GEX_SUB@assays$RNA["Gapdh", GEX_SUB$CellTypeFinal == "Macrophage" & GEX_SUB$tissue_merged == "Blood"] ) > 0 )
}

##### Background gene list ####
SUBDIR="16.background";dir.create(SUBDIR)
if(F){
	BACKGROUND_GENES=rownames(GEX)
	OUTNAME=paste0(SUBDIR,"/background_all_genes.txt")
	write.table(BACKGROUND_GENES,OUTNAME,quote=F,sep="\t",row.names=F,col.names=F)
	
	CELLTYPES<-unique(GEX$CellTypeFinal)
	CELLTYPES<-CELLTYPES[ ! CELLTYPES %in% grep("^rm-",CELLTYPES,value=T)]
	
	
	NAMEIDENT<-"CellTypeFinal"
	Idents(GEX)<-NAMEIDENT
	
	avgexp <- AverageExpression(GEX,verbose=T,assays="RNA")
	avgexp <- as.data.frame(avgexp$RNA)
	
	 # gt0
	#just take the gene is expression in the tissue
	# print("-------------------gt0-------------------")
	# X=sapply( CELLTYPES , function(identName){
	# 	GENELIST=rownames( avgexp[ which( avgexp[ , identName] > 0) , ] )
	# 	print(paste(identName,":",length(GENELIST),"genes"))
	# 	OUTNAME=paste0(SUBDIR,"/background_gt0_",identName,".txt")
	# 	write.table(x=GENELIST,file=OUTNAME,quote=F,row.names=F,col.names=F)
	# } )
	
	
	#gene taken if at least expression in 10 % of the cells in the cluster
	for( MinPercent in c(0, 1,5,10)){
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
	
	rm(BACKGROUND_GENES,OUTNAME,CELLTYPES,MinPercent,avgexp)
}

#### add background to gene list - metascape #### 
if(F){
	
	# Load background genes
	BACKfiles<-list.files(paste0( "16.background", "/", c("gt0","gt1","gt5","gt10")),full.names=T)
	BACKG<-lapply(BACKfiles, function(x){ TXT<-read.table(x) ; return(TXT$V1) } )
	names(BACKG)<-basename(BACKfiles)
	
	#DEGfiles<-list.files(paste0( "16.DEG", "/",  c("ref-Blood","ref-LN_inguinal","ref-LN_lumbar","ref-Spleen")),pattern=".xlsx",full.names=T)
	DEGfiles<-list.files(paste0( "16.DEG_mELT_vs_LNl+LNi+Spleen"),pattern=".xlsx",full.names=T)
	
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
	DEGfiles<-list.files(paste0( "16.DEG", "/",  "ref-Blood"),pattern=".xlsx",full.names=T); OUTNAME="16.DEG.list.ref-Blood.xlsx"
	DEGfiles<-list.files(paste0( "16.DEG", "/",  "ref-LN_inguinal"),pattern=".xlsx",full.names=T); OUTNAME="16.DEG.list.ref-LN_inguinal.xlsx"
	DEGfiles<-list.files(paste0( "16.DEG", "/",  "ref-LN_lumbar"),pattern=".xlsx",full.names=T); OUTNAME="16.DEG.list.LN_lumbar.xlsx"
	DEGfiles<-list.files(paste0( "16.DEG", "/",  "ref-Spleen"),pattern=".xlsx",full.names=T); OUTNAME="16.DEG.list.ref-Spleen.xlsx"
	
	DEGfiles<-list.files(paste0( "17.DEG_mELT_vs_LNl+LNi+Spleen"),pattern=".xlsx",full.names=T); OUTNAME="17.DEG.list.mELT_vs_LNl+LNi+Spleen.xlsx"
	
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
		setwd("/home/ga94rac/LRZ Sync+Share/AG_Lehmann-Horn/Gildas_storage/2021.02.16.CITEseq/experiment_Merged_Tissues+mELT_LN/19.GO_Enrichment_Metascape.auto/Mus_to_Mus/Tissues_ref-LN_lumbar")
		setwd("/home/ga94rac/LRZ Sync+Share/AG_Lehmann-Horn/Gildas_storage/2021.02.16.CITEseq/experiment_Merged_Tissues+mELT_LN/19.GO_Enrichment_Metascape.auto/Mus_to_Mus/Tissues_ref-LN_inguinal/")
		setwd("/home/ga94rac/LRZ Sync+Share/AG_Lehmann-Horn/Gildas_storage/2021.02.16.CITEseq/experiment_Merged_Tissues+mELT_LN/19.GO_Enrichment_Metascape.auto/Mus_to_Mus/Tissues_ref-Spleen/")
		
		
		LIST<-list.files(".",pattern="_FINAL_GO.csv",full.names=T,recursive=T)
		LIST2<-grep(pattern="mELT",x=LIST,value=T)
		FILES<-lapply(LIST2, function(x){df<-read.csv(x) ; df$source=x ; return(df) })
		names(FILES) <- gsub("Markers_","",basename(dirname(dirname(LIST2))))
		# remove column with specific names
		FILES2<-lapply(FILES, function(x){  x$tissue=gsub("X_MEMBER_","", names(x)[1]) ; names(x)[1:2]<-c("X_MEMBER","X_LogP");return(x)  })
		FILES2<-lapply(FILES2, function(x){ x$percent = x$X.GeneInGOAndHitList / x$X.GeneInGO * 100 ; return(x) })
		#SIGNIF_SUBSET<-lapply(FILES2, function(x){ subset( x , x$Log.q.value. < LOG_Q_VAL_LIMIT ) })
		DF <- do.call(rbind.data.frame, FILES2)
		rm(FILES,FILES2,LIST,LIST2)
		
		#Mus_to_Mus/Tissues_ref-LN_lumbar
		if(T){
			dir.create("../GO_MusMus_LN_lum")
			openxlsx::write.xlsx(DF,file = paste0("../GO_MusMus_LN_lum/GO_MusMus_LN_lum.xlsx"))
			saveRDS(DF, paste0("../GO_MusMus_LN_lum/GO_MusMus_LN_lum.Rds"))
		}
		
		#Mus_to_Mus/Tissues_ref-LN_inguinal
		if(T){
			dir.create("../GO_MusMus_LN_ing")
			openxlsx::write.xlsx(DF,file = paste0("../GO_MusMus_LN_ing/GO_MusMus_LN_ing.xlsx"))
			saveRDS(DF, paste0("../GO_MusMus_LN_ing/GO_MusMus_LN_ing.Rds"))
		}
		
		#Mus_to_Mus/Tissues_ref-_Splee
		if(T){
			dir.create("../GO_MusMus_Spleen")
			openxlsx::write.xlsx(DF,file = paste0("../GO_MusMus_Spleen/GO_MusMus_Spleen.xlsx"))
			saveRDS(DF, paste0("../GO_MusMus_Spleen/GO_MusMus_Spleen.Rds"))
		}
	}
	
	setwd("/home/ga94rac/LRZ Sync+Share/AG_Lehmann-Horn/Gildas_storage/2021.02.16.CITEseq/experiment_Merged_Tissues+mELT_LN/19.GO_Enrichment_Metascape.auto/Mus_to_Mus/")
	
	
	# process using the top for each category
	LOG_Q_VAL_LIMIT=log(0.05)
	if(F){
		# choose one
		TOP=30  ; WIDTH=12 ; HEIGHT=8
		TOP=100 ; WIDTH=12 ; HEIGHT=12
		
		
		for( index_GO_CAT in 1:3 ){
		
			if( index_GO_CAT == 1) { GO_CATEGORY="GO Biological Processes" ; GO_CAT_NAME="GO.Bio.Pro" }
			if( index_GO_CAT == 2) { GO_CATEGORY="Reactome Gene Sets" ; GO_CAT_NAME="GO.Reacto" }
			if( index_GO_CAT == 3) { GO_CATEGORY="CORUM" ; GO_CAT_NAME="GO.CORUM" }
		
			for( index_DF in 1:3 ){
				if( index_DF == 1) {
					DF <- readRDS("GO_MusMus_LN_lum/GO_MusMus_LN_lum.Rds") ; OUTNAME=paste0("GO_MusMus_LN_lum.",GO_CATEGORY,".top",TOP,".pdf")
					TITLENAME=paste("Top",TOP,"pathways in","mELT versus LN lumbar",paste0("\n(log q-FDR < ",round(LOG_Q_VAL_LIMIT,2),")"))
				}
				if( index_DF == 2) {
					DF <- readRDS("GO_MusMus_LN_ing/GO_MusMus_LN_ing.Rds") ; OUTNAME=paste0("GO_MusMus_LN_ing.",GO_CATEGORY,".top",TOP,".pdf")
					TITLENAME=paste("Top",TOP, "pathways in","mELT versus LN inguinal",paste0("\n(log q-FDR < ",round(LOG_Q_VAL_LIMIT,2),")"))
				}
				if( index_DF == 3) {
					DF <- readRDS("GO_MusMus_Spleen/GO_MusMus_Spleen.Rds") ; OUTNAME=paste0("GO_MusMus_Spleen.",GO_CATEGORY,".top",TOP,".pdf")
					TITLENAME=paste("Top",TOP, "pathways in","mELT versus Spleen",paste0("\n(log q-FDR < ",round(LOG_Q_VAL_LIMIT,2),")"))
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
				
				#ALLCELLTYPES=names(c(SIGN_1,SIGN_2,SIGN_3)) ; 	#ALLCELLTYPES <- gsub("_gt0[.][[:digit:]]+$", "", ALLCELLTYPES)#ALLCELLTYPES <- gsub("[-].*", "", ALLCELLTYPES) ; 	#ALLCELLTYPES <- gsub("_", " ", ALLCELLTYPES)#ALLCELLTYPES<-unique(ALLCELLTYPES) ; 	#rm(ALLCELLTYPES)
				CELLTYPES_order=c("naive T cell CD4+", "act. T cell CD4+", "T cell CD8+", "gd T cell", "Treg" ,
								  "NK + NK T cell",
								  "fol. zone B cell","mar. zone B cell", "ger. center B cell","Plasma cell",
								  "DC", "Granulocyte", "Macrophage", "Monocyte",
								  "Erythroid cell")
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

#### check some specific markers about NK NK T cell ####
##  does not filter based on immune_chains_clean4 since we are not sure about the NK cells and NKTcells
if(F){
	NAMEIDENT<-"RNA_snn_res.2"
	Idents(GEX)<-NAMEIDENT
	CLUSTERS=c("13","17","42") #"41",
	SPLIT="tissue_merged"
	# Zoom on the NK cells
	GEX@meta.data$UMAP_1 <- as.data.frame(GEX@reductions$umap@cell.embeddings)$UMAP_1
	GEX@meta.data$UMAP_2 <- as.data.frame(GEX@reductions$umap@cell.embeddings)$UMAP_2
	GEX_SUBSET <- subset(GEX, RNA_snn_res.2 %in% CLUSTERS& UMAP_1 < -4 & UMAP_2 > 4)
	
	GEX_SUBSET$immune_chains_clean5<-factor(GEX_SUBSET$immune_chains_clean4)
	levels(GEX_SUBSET$immune_chains_clean5) <- c("","","","","","T-pair-","","")
	#"B-cell-paired" "IGH" "IGK" "IGL" "no info" "T-cell-paired" "TRA" "TRB"
	GEX_SUBSET$immune_chains_clean5 <- paste0(GEX_SUBSET$immune_chains_clean5 , "clu", GEX_SUBSET$RNA_snn_res.2)
	
	
	
	for (FEATURE in c("Klrb1c","Cd3e","Cd8a") ){
		
		PLOT=FeaturePlot( subset(GEX_SUBSET,tissue_merged != "none"),
						  reduction="umap",features=FEATURE,pt.size=1,label=F,split.by=SPLIT,
						  repel=T,min.cutoff="q10",ncol=3,coord.fixed=T,order=T)
		OUTNAME=paste0(SUBDIR,"/test.",FEATURE,".Expr.split_",SPLIT,".png");png(OUTNAME,width=1200,height=600);print(PLOT);dev.off()
		
		PLOT=plot_density( subset(GEX_SUBSET,tissue_merged != "none"),
						   reduction="umap",features=FEATURE,size=1)
		OUTNAME=paste0(SUBDIR,"/test.",FEATURE,".Density.png");png(OUTNAME,width=400,height=400);print(PLOT);dev.off()
		
		PLOT<-VlnPlot( subset(GEX_SUBSET,tissue_merged != "none"),
					   features=FEATURE,split.by=SPLIT,pt.size=1)
		OUTNAME=paste0(SUBDIR,"/test.",FEATURE,".Violin.png");png(OUTNAME,width=800,height=400);print(PLOT);dev.off()
	}
	
	#dgCMatrix: compressed sparse column matrix
	GENE="Klrb1c" ; M<-as.matrix(GEX_SUBSET@assays$RNA[GENE,]) ; GEX_SUBSET@meta.data[, paste0("positive_",GENE) ] <- M[,] > 0
	GENE="Cd3e" ; M<-as.matrix(GEX_SUBSET@assays$RNA[GENE,]) ; GEX_SUBSET@meta.data[, paste0("positive_",GENE) ] <- M[,] > 0
	GENE="Cd8a" ; M<-as.matrix(GEX_SUBSET@assays$RNA[GENE,]) ; GEX_SUBSET@meta.data[, paste0("positive_",GENE) ] <- M[,] > 0
	
	
	TMP<-subset(GEX_SUBSET, tissue_merged != "none")
	
	TMP$RNA_snn_res.2 <- as.character(TMP$RNA_snn_res.2)
	
	table( cluster=  paste(TMP$RNA_snn_res.2 , TMP$immune_chains_clean5) , positive_Cd3e=TMP$positive_Cd3e )
	table( cluster=  paste(TMP$RNA_snn_res.2 , TMP$immune_chains_clean5) , positive_Klrb1c=TMP$positive_Klrb1c )
	table( cluster=  paste(TMP$RNA_snn_res.2 , TMP$immune_chains_clean5) , positive_Cd8a=TMP$positive_Cd8a )
	
	TABLE <- table( cluster=  paste(TMP$RNA_snn_res.2 , TMP$immune_chains_clean5), 
		   combined_bin_expr= TMP$positive_Klrb1c + TMP$positive_Cd3e - 10*TMP$positive_Cd8a )
	TABLE ; round(TABLE / rowSums(TABLE) * 100 )
	rm(TMP)
	
	
	TABLE <- table( as.character(GEX_SUBSET$RNA_snn_res.2), as.character(GEX_SUBSET$immune_chains_clean5 ) )
	TABLE ; round(TABLE / rowSums(TABLE) * 100 )
	
	
	
	Idents(GEX_SUBSET)<-"immune_chains_clean5"
	table(GEX_SUBSET$immune_chains_clean5)
	PLOT=DimPlot( subset(GEX_SUBSET, tissue_merged %in% c("mELT","LN_lumbar") & 
						 	immune_chains_clean5 %in% c("T-pair-clu13","T-pair-clu17","T-pair-clu42")),
				  reduction="umap",label=T,
				  order=NULL,pt.size=1,label.size=5,repel=T,split.by=SPLIT) + 
		theme_minimal(base_size = 16)+
		theme(legend.position = 'top')+
		ggtitle(paste(NAMEIDENT)) + coord_fixed()
	PLOT
	
	PLOT=DimPlot( subset(GEX_SUBSET, tissue_merged %ni% c("mELT","LN_lumbar","none") ),
				  reduction="umap",label=T,
				  order=NULL,pt.size=1,label.size=5,repel=T,split.by=SPLIT) + 
		theme_minimal(base_size = 16)+
		theme(legend.position = 'top')+
		ggtitle(paste(NAMEIDENT)) + coord_fixed()
	PLOT

	
	
	
	
	#table(GEX$tag_celltype)
	#Nk1.1
	#Nkp46
	#positive Cd3e : T-cell
	
	TABLE <- table( as.character(GEX_SUBSET$RNA_snn_res.2), as.character(GEX_SUBSET$tag_celltype ) )
	TABLE ; round(TABLE / rowSums(TABLE) * 100 )
	
	TABLE <- table( as.character(GEX_SUBSET$immune_chains_clean5), as.character(GEX_SUBSET$tag_celltype ) )
	TABLE ; round(TABLE / rowSums(TABLE) * 100 )
	
	Idents(GEX_SUBSET)<-"tag_celltype"
	table(GEX_SUBSET$tag_celltype)
	PLOT=DimPlot( subset( GEX_SUBSET, tissue_merged %in% c("LN_inguinal","LN_lumbar","mELT","Spleen") &
						  tag_celltype %in% c("NK","CD4","CD8a") ),
				  reduction="umap",label=T,
				  order=NULL,pt.size=2,label.size=5,repel=T,split.by=SPLIT) + 
		theme_minimal(base_size = 16)+
		theme(legend.position = 'top')+
		ggtitle(paste(NAMEIDENT)) + coord_fixed()
	PLOT
	
	
	GEX$mouse <- factor(GEX$orig.ident)
	#"I18A","I18B",		"J119A","J119B","J119C",	"J9A",J9B","J9C",
	#""Z3.E21.LN", "Z3.E21.mELT", "Z3.E72.LN" , "Z3.E72.mELT", "Z3.F21.LN", "Z3.F21.mELT")
	levels(GEX$mouse) <-c("I18","I18","J119","J119","J119","J9","J9","J9","E21","E21","E72","E72","F21","F21") 
	
	table(GEX$mouse, GEX$immune_chains) #F21 mice has no information about TCR?
	
	TABLE<-table( paste(GEX$mouse,GEX$tissue_merged), GEX$CellTypeFinal)
	TABLE[  c( grep("mELT",rownames(TABLE),value=T),grep("LN_lumbar",rownames(TABLE),value=T) ) , c("NK T cell","NK cell")]
	
}





#### Cellchat ####
if(F){
	#source("https://raw.githubusercontent.com/GildasLepennetier/SCr/main/CChat.R")
	library(CellChat)
	GEX$Condition=GEX$tissue_merged
	# subset to avoid including unwanted comparisons & noise
	GEX = subset( GEX , CellTypeFinal %ni% c( grep(pattern="^rm-",x = unique(GEX$CellTypeFinal),value = T ) ) )
	GEX = subset( GEX , tissue_merged %ni% c("none") )
	
	CONDITIONS <- c("CSF","mELT","LN_lumbar","LN_inguinal","Spleen","Blood") # c("CSF","mELT","LN_lumbar","LN_inguinal","Spleen","Blood")
	TYPE_OF_INTER=c("Secreted Signaling","Cell-Cell Contact","ECM-Receptor") #c("Secreted Signaling","Cell-Cell Contact","ECM-Receptor")

	GEX.backup=GEX
	#CChat()
	
	#CONDITIONS <- unique( GEX.backup$Condition )
	
	#those are apparently not gene symbols
	#H2-Q8 H2-T9 H2-T18 H2-Q9 H2-L H2-BI H2-D H60a H2-Ea-ps
	
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
				
				Idents(GEX)<-NAMEIDENT
				#DefaultAssay(GEX) # RNA
				
				cellchat <- createCellChat(object = GetAssayData(GEX,assay="RNA",slot="data"),meta=GEX@meta.data,group.by=NAMEIDENT)
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
#### Overlap of the VDJ repertoire ####
if(F){
	
}
#### DATA=sessionInfo() ####
#DATA=devtools::session_info()
writeLines(capture.output(devtools::session_info()), "sessionInfo.txt")
