from scipy.io import mmwrite
import scvelo as scv
import os
import numpy as np
vel_folder = "/gpfs/gss1/work/sdusinglecell/Andreas/human_liver_nuclei/new_velocyto/"
for path, subdirs, files in os.walk(vel_folder):
    for subdir in subdirs:
        if subdir != "expression" :
            if subdir != "velocity":

                file = os.path.join(path, subdir, "velocity", "zUMI_"+subdir+".loom")
                #print(file)
                adata = scv.read(file, cache=False )
                mat = adata.X.T
                outfile = file.replace(".loom",".mtx")
                print(outfile)
                mmwrite(outfile, mat, comment='', field=None, precision=None, symmetry=None)
                np.savetxt(file.replace(".loom",".csv"), adata.obs_names.tolist(), delimiter=",", fmt='%s')









files <- c(#"/work/sduknn/Andreas/xenopus/zUMI_out/scCapSt08/zUMIs_output/expression/zUMI_scCapSt08.dgecounts.rds",
"/work/sduknn/Andreas/xenopus/zUMI_out/scCapSt105/zUMIs_output/expression/zUMI_scCapSt105.dgecounts.rds",
"/work/sduknn/Andreas/xenopus/zUMI_out/scCapSt12/zUMIs_output/expression/zUMI_scCapSt12.dgecounts.rds",
#"/work/sduknn/Andreas/xenopus/zUMI_out/scCapSt13/zUMIs_output/expression/zUMI_scCapSt13.dgecounts.rds",
#"/work/sduknn/Andreas/xenopus/zUMI_out/scCapSt16/zUMIs_output/expression/zUMI_scCapSt16.dgecounts.rds",
#"/work/sduknn/Andreas/xenopus/zUMI_out/scCapSt18/zUMIs_output/expression/zUMI_scCapSt18.dgecounts.rds",
#"/work/sduknn/Andreas/xenopus/zUMI_out/scCapSt20/zUMIs_output/expression/zUMI_scCapSt20.dgecounts.rds",
"/work/sduknn/Andreas/xenopus/zUMI_out/scCapSt22/zUMIs_output/expression/zUMI_scCapSt22.dgecounts.rds",
"/work/sduknn/Andreas/xenopus/zUMI_out/scCapSt22b/zUMIs_output/expression/zUMI_scCapSt22b.dgecounts.rds"#,
#"/work/sduknn/Andreas/xenopus/zUMI_out/scCapSt24/zUMIs_output/expression/zUMI_scCapSt24.dgecounts.rds",
#"/work/sduknn/Andreas/xenopus/zUMI_out/scCapSt27/zUMIs_output/expression/zUMI_scCapSt27.dgecounts.rds"
)

#files <- c("/work/sduknn/Andreas/xenopus/zUMI_out/scCapSt08/zUMIs_output/expression/zUMI_scCapSt08.dgecounts.rds")
library(Matrix)
library(DropletUtils)
library(ggplot2)
library(ggthemes)
library(scater)
library(Seurat)

gene_converter <- function(seu){
    genes <- as.data.frame(rownames(seu))
    names(genes)[1] <- 'gene_id'
    genes <- merge(genes, gtf_df3, by="gene_id", sort = F)
    #rownames(seu) <- genes$gene_name
    #print(head(genes) )
    return(genes$gene_name)
}

saver <- function( SCE, seu, file){
    #Save in a cellranger like format
    #df = data.frame(gene_converter(seu)) #gene_names
    df = data.frame(rownames(seu@assays$RNA))
    df[,2] = df[,1] #gene_names
    dir <- dirname(file)
    write.table(SCE@colData@rownames, file=paste(dir, '/barcodes.tsv', sep = ''), quote=FALSE, sep='\t', row.names = F, col.names=F)
    write.table(df, file=paste(dir, '/genes.tsv', sep = ''), quote=FALSE, sep='\t', row.names = F, col.names=F)
    Matrix::writeMM(Matrix(GetAssayData(object = seu)), paste(dir, '/matrix.mtx', sep = ''))

}

for (file in files){
    set.seed(100)
    #mat <- Matrix::readMM(file)
    mat <- readRDS(file)
    mat <- mat$umicount$inex$all

    filtered_mat <- mat[Matrix::rowSums(mat) > 0, ]
    drops <-emptyDrops(filtered_mat)
    sum(drops$FDR <= 0.01, na.rm = TRUE)
    names <- which(drops$FDR <= 0.01)
    name_file <- gsub(".dgecounts.rds", ".csv", file)
    #name_ <- read.csv(name_file)
    colnames(filtered_mat) <- paste0(strsplit(basename(file), '[.]')[[1]][1], ':', colnames(filtered_mat), 'x')
    name_ <- colnames(filtered_mat)[names]
    write.table(name_, file = name_file ,row.names=FALSE, na="",col.names=FALSE, sep=",", quote = FALSE)

    test <- drops[c('Total', 'FDR')]
    test['Rank']  <- rank(-test$Total)
    test['Cell'] <- test$FDR <= 0.01
    test <- as.data.frame(test)
    test['Cell'] <- replace(test$Cell, is.na(test$Cell), FALSE)

    gg_name = gsub(".dgecounts.rds", "_BCrank.png", file)
    library(ggplot2)
    library(ggthemes)
    cb_palette <- c( "lightgray", "navyblue")
    ggplot(test, aes(x = Rank, y = Total, color=Cell)) +
    geom_point() +
    scale_x_log10() +
    scale_y_log10() +
    theme_bw() +
    scale_color_manual(values = cb_palette)+
    theme(legend.position="top")+
    geom_point(size = 0.01)+
    geom_point(alpha = 1/10)+
    theme(axis.text = element_text(size = 10),
        axis.title = element_text(size = 14),
        title = element_text(size = 16)) +
    labs(x = "Rank", y = "UMI counts", title = "Barcode Rank vs Total UMI counts")

    ggsave(gg_name)

    #Make single cell experiment object
    cells <- filtered_mat[, which(drops$FDR <= 0.01)]
    SCE <- SingleCellExperiment(assays = list(counts = cells))
    SCE <- calculateQCMetrics(SCE)

    test2 <- as.data.frame(SCE$total_features_by_counts)
    test2['total_counts'] <- SCE$total_counts
    colnames(test2) <- c('total_features_by_counts', 'total_counts')

    gg_name = gsub(".dgecounts.rds", "_gene_feature.png", file)
    ggplot(test2, aes(x = total_counts, y = total_features_by_counts)) +
    geom_point() +
    theme_bw() +
    geom_point(size = 0.01)+
    geom_point(alpha = 1/100)+
    theme(axis.text = element_text(size = 10),
        axis.title = element_text(size = 14),
        title = element_text(size = 16)) +
    labs(x = "UMI counts", y = "Features", title = "features vs. counts")
    ggsave(gg_name)

    #make seurat object
    mat_seurat <- CreateSeuratObject(counts = filtered_mat)
    seurat_filtered <- subset(mat_seurat , cells = SCE@colData@rownames )
    #Load GTF
    gtf <- rtracklayer::import('/work/sduknn/Andreas/xenopus/xen_ref/XENLA_10X/XENLA_GCA001663975v1_XBv9p2/genes/genes.gtf')
    gtf_df=as.data.frame(gtf)
    gtf_df2 <- as.data.frame(gtf_df$gene_id)
    names(gtf_df2)[1] <- 'gene_id'
    gtf_df2['gene_name'] <- gtf_df$gene_name
    gtf_df3 <- unique(gtf_df2)

    saver(SCE , seurat_filtered , file)

}
