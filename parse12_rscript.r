#This Rscript was written to analyze the first run of samples using the Parse Bio mini kit. 
#There are 12 samples, 3 batches of 4 lines, with 2 sporadic parkinson's lines and 2 controls.
#This is my first time analyzing data like this 

#load libraries

library(Seurat)
library(dplyr)
library(Matrix)
library(ggplot2)
library(harmony)
library(DESeq2)

#Load the data from the Parse bio pipeline output

parse12 <- ReadParseBio("C:/Users/tgolds3/Documents/scRNA_seq/12_samp/DGE_filtered")

#From Parse workflow, checking if there are empty gene names, and if so a name is added

table(rownames(parse12) == "")
rownames(parse12)[rownames(parse12) == ""] <- "unknown"

#read in cell metadata 

cell_meta <- read.csv(paste0("C:/Users/tgolds3/Documents/scRNA_seq/12_samp/DGE_filtered", "/cell_metadata.csv"), row.names = 1)

#Create Seurat Object with raw non-normalized data

parse12.seurat.obj <- CreateSeuratObject(parse12, min_genes = 250, min_cells = 10,
names.feild = 0, meta.data = cell_meta)

#From parse workflow, remove idents assigned initially, setting intial cell class to a single type, this will change after clustering

parse12.seurat.obj@meta.data$orig.ident <- factor(rep("parse12.seurat.obj", nrow(parse12.seurat.obj@meta.data)))
Idents(parse12.seurat.obj) <- parse12.seurat.obj@meta.data$orig.ident

#Create a vector for the healthy samples 
# healthy samples: "x3448B3", "x3448B1", "x3448B2", "TD22B1", "TD22B2", "TD22B3"
# disease samples: "x2965B3", "x2965B1", "x2965B2", "TD07B1", "TD07B2", "TD07B3"

healthy <- c("x3448B3", "x3448B1", "x3448B2", "TD22B1", "TD22B2", "TD22B3")

#Create column in metadata to reflect disease state

cell_meta_3 <- mutate(cell_meta, "disease_state" = ifelse(sample %in% healthy, "healthy", "disease"))

#assign groups to the sample batches

x3448 <- c("x3448B1", "x3448B2", "x3448B3")
x2965 <- c("x2965B1", "x2965B2", "x2965B3")
TD07 <- c("TD07B1", "TD07B2", "TD07B3")
TD22 <- c("TD22B1", "TD22B2", "TD22B3")

#assign in the metadata

cell_meta_3s <- mutate(cell_meta_3, "comb_samp" = ifelse(sample %in% x3448, "3448", ifelse(sample %in% x2965, "2965", ifelse(sample %in% TD07, "TD07", ifelse(sample %in% TD22, "TD22", "")))))

#Added the edited meta data to the metadata of the Seurat object. In the future I would just mutate the actual metadata of the Seurat object as this wasn't necessary. 

parse12.seurat.obj <- AddMetaData(parse12.seurat.obj, cell_meta_3s)
head(parse12.seurat.obj@meta.data)

#QC steps and filtering

parse12.seurat.obj[["percent.mt"]] <- PercentageFeatureSet(parse12.seurat.obj, pattern = "^MT-")
View(parse12.seurat.obj@meta.data)

VlnPlot <- VlnPlot(parse12.seurat.obj, pt.size = 0.10, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
QCScatter <- FeatureScatter(parse12.seurat.obj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") + geom_smooth(method = 'lm')

plot(VlnPlot)

parse12.seurat.obj <- subset(parse12.seurat.obj, subset = nFeature_RNA < 12000 &  nFeature_RNA >300 & nCount_RNA < 30000 & percent.mt < 15)

#Normalize the data

parse12.seurat.obj <- NormalizeData(parse12.seurat.obj)
str(parse12.seurat.obj)

# Identify Highly Variable Features

parse12.seurat.obj <- FindVariableFeatures(parse12.seurat.obj, selection.method = "vst", nfeatures = 2800)

top10 <- head(VariableFeatures(parse12.seurat.obj), 10)

varfeat <- VariableFeaturePlot(parse12.seurat.obj)
LabelPoints(plot = varfeat, points = top10, repel = TRUE)

#Scale the data

all_genes <- rownames(parse12.seurat.obj)
parse12.seurat.obj <- ScaleData(parse12.seurat.obj, features = all_genes)

#Run PCA

parse12.seurat.obj <- RunPCA(parse12.seurat.obj, features = VariableFeatures(object = parse12.seurat.obj))

print(parse12.seurat.obj[["pca"]], dims = 1:5, nfeatures = 5)
DimHeatmap(parse12.seurat.obj, dims = 1, cells = 500, balanced = TRUE)

#Elbow Plot to determine number of PCs to work with 
#20 were chosen

ElbowPlot(parse12.seurat.obj)

#Find neighbours

parse12.seurat.obj <- FindNeighbors(parse12.seurat.obj, dims = 1:20)

#Choosing a resolution for clustering
#initally chose 0.1

parse12.seurat.obj <- FindClusters(parse12.seurat.obj, resolution = c(0.1, 0.3, 0.5, 0.7, 1))
View(parse12.seurat.obj@meta.data)

DimPlot(parse12.seurat.obj, group.by = "RNA_snn_res.0.1", label = TRUE)

#Setting identities

Idents(parse12.seurat.obj)
Idents(parse12.seurat.obj) <- "RNA_snn_res_0.1"

#UMAP

parse12.seurat.obj <- RunUMAP(parse12.seurat.obj, dims = 1:20)
umap1 <- DimPlot(parse12.seurat.obj, reduction = "umap", group.by = "comb_samp")
umap1
umap2 <- DimPlot(parse12.seurat.obj, reduction = "umap", group.by = "disease_state")
umap2

#Using Harmony to alleviate the batch effects by grouping by disease state

harmony1 <- parse12.seurat.obj %>%
  RunHarmony(group.by.vars = "disease_state", plot_convergence = FALSE)

harmony1@reductions

harmony1.embed <- Embeddings(harmony1, "harmony")
harmony1.embed[1:10, 1:10]

#Running the UMAP and clustering using harmony embeddings instead of PCA

harmony1 <- harmony1 %>%
  RunUMAP(reduction = "harmony", dims = 1:20) %>%
  FindNeighbors(reduction = "harmony", dims = 1:20) %>%
  FindClusters(resolution = 0.3)

umap3 <- DimPlot(harmony1, reduction = "umap", group.by = "disease_state")

umap2|umap3

#Add cluster info to harmonized data set by disease state

clusters <- DimPlot(harmony1, reduction = "umap", group.by = "seurat_clusters", label = TRUE)
condition <- DimPlot(harmony1, reduction = "umap", group.by = "disease_state")

#Find all markers on harmonized data set

markers <- FindAllMarkers(harmony, logfc.threshold = 0.25, min.pct = 0.25)

#Save

saveRDS(harmony1, "G:/parse12harmony.rds")

#Renamed one cluster that I knew for sure 

Idents(harmony)
harmony <- RenameIdents(harmony, "11" = "Epithelial")

#Compared Cluster 5 and 2 because they segregate apart from eachother and are disease or control

compcluster5_2 <- FindMarkers(harmony, ident.1 = 5, ident.2 = 2)

compcluster5_2deseq2 <- FindMarkers(harmony, ident.1 = 5, ident.2 = 2, test.use = "DESeq2")

write.csv(compcluster5_2, "C:/Users/tgolds3/Documents/scRNA_seq/12_samp/cluster5_2.csv" )

write.csv(compcluster5_2deseq2, "C:/Users/tgolds3/Documents/scRNA_seq/12_samp/cluster5_2deseq2.csv" )


