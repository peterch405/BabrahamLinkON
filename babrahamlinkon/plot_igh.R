#!/usr/bin/Rscript

'''
BabrahamLinkON: Analysis pipeline for VDJ-seq
Copyright (C) 2017  Peter Chovanec

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <https://www.gnu.org/licenses/>.
'''

source("https://bioconductor.org/biocLite.R")
if("GenomicFeatures" %in% rownames(installed.packages()) == FALSE) {
  biocLite("GenomicFeatures")
}else{
    require(GenomicFeatures)
}
if("GenomicFeatures" %in% rownames(installed.packages()) == FALSE) {
  biocLite("GenomicFeatures")
}else{
  require(GenomicFeatures)
}
if("Gviz" %in% rownames(installed.packages()) == FALSE) {
  biocLite("Gviz")
}else{
  require(Gviz)
}
if("argparse" %in% rownames(installed.packages()) == FALSE) {
  install.packages("argparse", repos = "http://cran.us.r-project.org")
}else{
  require(argparse)
}
if("biomaRt" %in% rownames(installed.packages()) == FALSE) {
  biocLite("biomaRt")
}else{
  require(biomaRt)
}
#######################
### Set up argparse ###
#######################

parser <- ArgumentParser(description='Plot BAM alignment coverage of IgH genes')

parser$add_argument("-b", "--bam", dest="input_list", type="character", metavar='N', nargs='+', help="BAM input files")
parser$add_argument("-n", "--names", dest="names_list", type="character", metavar='N', nargs='+', help="Names to be used for plotting")
parser$add_argument("-o", "--out_pdf", dest="output", help="Output pdf file location/name")
parser$add_argument("-r", "--region", dest="region", type="character", default="chr12:113531809-116015193", help="Region to plot (e.g. chr12:113531809-116015193)")
parser$add_argument("--genome", dest="genome", type="character", default="mm10", help="Which genome is being used (e.g. mm10, hg38)")

args <- parser$parse_args()
# parser$print_help()
# args <- parser$parse_args(c("-b", "/media/chovanec/My_Passport/Sync/BabrahamLinkON/tests/UMI_tmp/lane3_WT_FrBC_1_GCCAAT_L003_R1_val_1_40k_V1_GACTCGT_umi_ext.bam",
#                             "/media/chovanec/My_Passport/Sync/BabrahamLinkON/tests/UMI_tmp/lane3_WT_FrBC_1_GCCAAT_L003_R1_val_1_40k_V1_CTGCTCCT_umi_ext_dedup_sorted.bam",
#                             "-n", "Before", "After", "-o", "/media/chovanec/My_Passport/Sync/BabrahamLinkON/tests/igh.pdf"))
#
# print(args$input_list)

if(length(args$input_list) != length(args$names_list)){
    stop("Number of files is not equal to number of names supplied")
}


################
### Plot BAMs###
################

rgn <- unlist(strsplit(args$region, "(:|-)"))
chrm <- rgn[1]
start <- as.integer(rgn[2])
end <- as.integer(rgn[3])
chrm_num <- unlist(strsplit(chrm, "chr"))[2]

tracks <- list()

tracks[["gtrack"]] <- GenomeAxisTrack()

tracks[["itrack"]] <- IdeogramTrack(genome=args$genome,chromosome=chrm)

#Loop through all the input files
l_size <- length(tracks)
for(i in 1:length(args$input_list)){
    f_name <- args$input_list[i]
    tracks[[l_size+i]] <- DataTrack(range=f_name,genome=args$genome, name=args$names_list[i],
                                    chromosome=chrm, type = "histogram",
                                    col.histogram= "#377EB8", fill="#377EB8")

}

if(args$genome == "mm10"){
  dataset_nam <- "mmusculus_gene_ensembl"
} else if(args$genome == "hg38"){
  dataset_nam <- "hsapiens_gene_ensembl"
}

#IgH only annotation dataframe from ensembl
ensembl <- useEnsembl(biomart="ensembl", dataset=dataset_nam)
genes <- getBM(attributes=c('chromosome_name','start_position','end_position', 'strand', 'external_gene_name'),
               filters ='chromosome_name', values=chrm_num, mart = ensembl)
igh_genes <- genes[grep("(Ig|IG)", genes$external_gene_name),]
colnames(igh_genes) <- c("chromosome", "start", "end", "strand", "symbol")



tracks[["grtrack"]] <- GeneRegionTrack(igh_genes, genome=args$genome, chromosome=chrm, name = "IgH")
tracks[["grtrack_id"]] <- GeneRegionTrack(igh_genes, genome=args$genome, chromosome=chrm, name = "IgH with labels", showId=TRUE)


#Write to pdf
pdf(args$output, height=20, paper='A4') #,width=10,height=18 ,
plotTracks(tracks, from = start, to = end) #extend.left=200000
dev.off()
