setwd("/app/personalfolder/canoes_v2.0/CANOES/output")
gc <- read.table("/app/personalfolder/canoes_v2.0/CANOES/output/CANOES/GCpercent.A.tsv.conv")$V2
canoes.reads <- read.table("/app/personalfolder/canoes_v2.0/CANOES/output/CANOES/all.canoes.A.coverage.tsv.conv.no_header")
sample.names <- c("SGT144117","SGT140392","SGT162459")
names(canoes.reads) <- c("chromosome", "start", "end", sample.names)
target <- seq(1, nrow(canoes.reads))
canoes.reads <- cbind(target, gc, canoes.reads)
source("/app/personalfolder/canoes_v2.0/CANOES/CANOES.v3.R")
xcnv.list <- vector('list', length(sample.names))
for (i in 1:length(sample.names)){
xcnv.list[[i]] <- CallCNVs(sample.names[i], canoes.reads)
}
xcnvs <- do.call('rbind', xcnv.list)
write.csv(xcnvs, "/app/personalfolder/canoes_v2.0/CANOES/output/CANOES/results.A.csv", row.names=FALSE)
