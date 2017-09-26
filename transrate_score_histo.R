# histo of transrate scores

# Rscript transrate_scores_histo.R test

args <- commandArgs(TRUE)
srcFile <- args[1]
# inFile <- paste0("/mnt/nfs/ryan/diel1/completed_assemblies/transrate_QCed/", srcFile, ".transrate_files/", srcFile, ".contigs.csv")
inFile <- "/Users/rgroussman/data/test.contigs.csv"
outFile <- paste0(srcFile, ".scores_histo.png")

data = read.table(file=inFile, sep=',', header=TRUE, row.names=1)
png(outFile)
p1 = hist(data$score,
  main="Histogram for transrate scores",
  xlab="transrate scores",
  border="blue",
  col=rgb(0,0,1,1/4),
  las=1,
  breaks=20)
p2 = hist(data$sCseg,
          main="Histogram for transrate scores",
          xlab="transrate scores",
          border="blue",
          col=rgb(1,0,0,1/4),
          las=1,
          breaks=20,
          add=T)


tr_scores <- as.data.frame(data)
tr_scores <- as.data.frame(data$sCseg, data$scores)
scores = data[data$scores]
sCseg = data[data$sCseg]
sCnuc = data[data$sCnuc]
sCord = data[data$sCord]
sCcov = data[data$sCcov]

tr_stack <- stack(tr_scores)
ggplot(dfs, aes(x=values)) + geom_density()

library(ggplot)
allscores <- rbind(carrots, cukes)
ggplot2(data, aes(data$score)) + geom_density(alpha = 0.2)

geom_density(mapping = NULL, data = data, stat = "density",
             na.rm = FALSE, show.legend = NA,
             inherit.aes = TRUE)
ggplot(data$scores, aes(data$scores)) +
  geom_density(adjust = 1/5)


