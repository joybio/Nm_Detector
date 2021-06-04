library(DESeq2)
#library(edgeR)
#library(ggplot2)
#library(GenomicFeatures)
library(tidyr)
library(dplyr)
library(stringr)
argv <- commandArgs(T)
if (length(argv) != 3) {
        cat('\nuseage:\nRscript DEseq2.reads_end.r 3pend_counts_matrix.xls group.csv pre_results.xls\nexample:\nRscript scripts/DEseq2.reads_end.r 3pend_counts_matrix.xls group.csv pre_results.xls\n\n')
        q('no')
}

input.dir <- unlist(str_split(argv[1],"/"))
out.dir <- paste0(str_c(input.dir[-length(input.dir)],collapse ="/"), "/")

print(argv[1:10])
reads_matrix <- read.csv(argv[1],header = T,row.names = 1,sep = "\t",stringsAsFactors = F)
head(reads_matrix)
#group里必须有treat，且样品名称不能带-
group <- read.csv(argv[2],header = T,row.names=1,stringsAsFactors = F)
head(group)
#OED组reads end必须大于20才有意义。
treat_group <- subset(group,group$type == "treat")
treat_id <- row.names(treat_group)
treat_reads <- reads_matrix[,treat_id]
treat_reads_filter <- treat_reads %>%  filter_all(all_vars(. > 50))
#过滤。 
reads_matrix_filt <- reads_matrix[row.names(treat_reads_filter),]
head(reads_matrix_filt)
condition <- factor(group$type, levels = c("control","treat"))
colData <- data.frame(row.names=colnames(reads_matrix), condition)
dds <- DESeqDataSetFromMatrix(reads_matrix, colData, design= ~ condition)
dds <- DESeq(dds)
res= results(dds)
result <- paste0(out.dir, "all_sites_results.csv")
write.csv(res,file=result)
table(res$padj<0.05)
diff_gene_deseq2 <-subset(res, padj < 0.01 & log2FoldChange > 3)
diff.result <- paste0(out.dir, "diff_sites_deseq2.csv")
write.csv(diff_gene_deseq2,file=diff.result)
merge.res <- cbind(reads_matrix,res)
merge.res.filter <- merge.res[row.names(reads_matrix_filt),]
merge.res.filter.sig <- subset(merge.res.filter,merge.res.filter$log2FoldChange >= 3 & merge.res.filter$padj <= 0.01)
merge.result <- paste0(out.dir, "OED_filter_50.sig.xls")
write.csv(merge.res.filter.sig,merge.result)
merge.res.filter.sig$pos <- row.names(merge.res.filter.sig)
head(merge.res.filter.sig)
class(merge.res.filter.sig)
merge.res.filter.sig <- as.data.frame(merge.res.filter.sig)
merge.res.filter.sig.sep <- separate(merge.res.filter.sig,pos,c("chr","start","stop"),sep="__")
#pos.result <- paste0(out.dir, "merge.res.filter20.sig.sep.xls")
write.table(merge.res.filter.sig.sep,argv[3],quote = F,sep="\t",row.names = F)


