# histo of transrate scores

data = read.table(file='/Users/rgroussman/data/test.contigs.csv', sep=',', header=TRUE, row.names=1)
hist(data$score, 
  main="Histogram for transrate scores", 
  xlab="transrate scores", 
  border="blue", 
  col="green",
  las=1, 
  breaks=20)
