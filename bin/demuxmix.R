#!/usr/bin/env Rscript

if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
is_demuxmix<-require("demuxmix")
if (!require("demuxmix", quietly = TRUE))
    BiocManager::install("demuxmix")

library(Seurat)
library(demuxmix)
library(argparse)
library(data.table)

# Create a parser
parser <- ArgumentParser("Parameters for Demuxmix")
parser$add_argument("--seuratObject", help = "Seurat object. Assumes that the hash tag oligo (HTO) data has been added and normalized.")
parser$add_argument("--pAcpt", help='Acceptance probability that must be reached in order to assign a droplet to a hashtag. ')
parser$add_argument("--assay", help='Assay name')
parser$add_argument("--model", help='A character specifying the type of mixture model to be used. Either "naive", "regpos", "reg" or "auto".')
parser$add_argument("--alpha_demuxmix", help='Threshold defining the left tail of the mixture distribution where droplets should not be classified as "positive". ')
parser$add_argument("--beta_demuxmix", help='Threshold for defining the right tail of the mixture distribution where droplets should not be classified as "negative".')
parser$add_argument("--tol_demuxmix", help='Convergence criterion for the EM algorithm used to fit the mixture models.')
parser$add_argument("--maxIter_demuxmix", help='Maximum number of iterations for the EM algorithm')
parser$add_argument("--k_hto", help='Factor to define outliers in the HTO counts.')
parser$add_argument("--k_rna", help='Factor to define outliers in the numbers of detected genes.')
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

res_dt <- as.data.table(demuxmix_classify)
res_dt[,Classification := Type]
res_dt[Classification == "multiplet", Classification := "doublet"]
res_dt[Classification == "uncertain", Classification := "negative"]


write.csv(res_dt, paste0(args$outputdir, "/","_assignment_htodemux.csv"),row.names=FALSE)
#write.csv(hashtag[[paste0(args$assay,"_classification")]], paste0(args$outputdir, "/", args$assignmentOutHTOdemux, "_assignment_htodemux.csv"))
saveRDS(demuxmix_demul, file=paste0(args$outputdir, "/", args$objectOutHTOdemux,".rds"))
