#!/usr/bin/env Rscript

is_bioconductor<-require("BiocManager")
if(is_bioconductor){
    install.packages("BiocManager")
}
is_demuxmix<-require("demuxmix")
if(is_bioconductor){
    BiocManager::install("demuxmix")
}
library(Seurat)
library(demuxmix)
library(argparse)

# Create a parser
parser <- ArgumentParser("Parameters for Demuxmix")
parser$add_argument("--seuratObject", help = "Seurat object. Assumes that the hash tag oligo (HTO) data has been added and normalized.")
parser$add_argument("--assay", help='Output directory')
parser$add_argument("--model", help='Output directory')
parser$add_argument("--alpha_demuxmix", help='Output directory')
parser$add_argument("--beta_demuxmix", help='Output directory')
parser$add_argument("--tol_demuxmix", help='Output directory')
parser$add_argument("--maxIter_demuxmix", help='Output directory')
parser$add_argument("--k_hto", help='Output directory')
parser$add_argument("--k_rna", help='Output directory')
parser$add_argument("--outputdir", help='Output directory')
args <- parser$parse_args()

#check if seurat object is correct
if (!endsWith(args$seuratObject, ".rds")){
    seuratObj <- list.files(args$seuratObject, pattern = "\\.rds$", full.names = TRUE)[1]
}else{
    seuratObj <- args$seuratObject
}


Argument <- c("seuratObject", "assay", "model", "alpha_demuxmix", "beta_demuxmix", "tol_demuxmix", "maxIter_demuxmix", "k_hto","k_rna")
Value <- c(seuratObj, args$assay, args$model, args$alpha_demuxmix, args$beta_demuxmix, args$tol_demuxmix, args$maxIter_demuxmix,args$k_hto,args$k_rna)

params <- data.frame(Argument, Value)
# Loading Seurat object
hashtag <-readRDS(seuratObj)


### Demuxmix
hto_counts <- as.matrix(GetAssayData(hashtag[["HTO"]], slot = "counts"))
#Demultiplexing process
demuxmix_demul <- demuxmix(hto_counts, model = model)
#hashtag Classification
demuxmix_classify <- dmmClassify(dmm)
#do we want to transform multiples into doublets?
#do we want to join result for negatives and uncertains
#Transformar una vez que los resultados este en tabla

write.csv(demuxmix_classify, paste0(args$outputdir, "/","_assignment_htodemux.csv"),row.names=FALSE)
#write.csv(hashtag[[paste0(args$assay,"_classification")]], paste0(args$outputdir, "/", args$assignmentOutHTOdemux, "_assignment_htodemux.csv"))
saveRDS(demuxmix_demul, file=paste0(args$outputdir, "/", args$objectOutHTOdemux,".rds"))
