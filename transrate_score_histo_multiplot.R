# histo of transrate scores

# Rscript transrate_scores_histo.R test

args <- commandArgs(TRUE)
sample_id <- args[1]
first_file <- paste0(sample_id, ".transrate_files/", sample_id, ".contigs.csv")
Cseg95_file <- paste0(sample_id, ".transrate_files/", sample_id, ".Cseg95.contigs.csv")
Cseg90_file <- paste0(sample_id, ".transrate_files/", sample_id, ".Cseg90.contigs.csv")

outFile_scores <- paste0(sample_id, ".scores.multi_scores_histo.png")
outFile_sCnuc <- paste0(sample_id, ".sCnuc.multi_scores_histo.png")
outFile_sCord <- paste0(sample_id, ".sCord.multi_scores_histo.png")
outFile_sCcov <- paste0(sample_id, ".sCcov.multi_scores_histo.png")


raw_data = read.table(file=first_file, sep=',', header=TRUE, row.names=1)
Cseg95_data = read.table(file=Cseg95_file, sep=',', header=TRUE, row.names=1)
Cseg90_data = read.table(file=Cseg90_file, sep=',', header=TRUE, row.names=1)

# plot scores
# png(outFile_scores)
# p1 = hist(raw_data$score,
#   main="Histogram for transrate scores",
#   xlab="transrate scores",
#   border="blue",
#   col=rgb(0,0,1,1/4),
#   las=1,
#   breaks=20)
# p2 = hist(Cseg90_data$score,
#           col=rgb(1,0,0,1/4),
#           las=1,
#           breaks=20,
#           add=T)
# p3 = hist(Cseg95_data$score,
#           col=rgb(0,1,0,1/4),
#           las=1,
#           breaks=20,
#           add=T)

# plot scores
png(outFile_scores)
title <- paste0(sample_id, " density plot of transrate scores")
plot(density(raw_data$score), col="red", main=title)
lines(density(Cseg90_data$score), col="blue")
lines(density(Cseg95_data$score), col="green")
legend(-0.1, 1.5, legend=c("Raw", "Cseg90", "Cseg95"),
       col=c("red", "blue", "green"), lty=1:2, cex=0.8)

png(outFile_sCnuc)
title <- paste0(sample_id, " density plot of sCnuc scores")
plot(density(raw_data$sCnuc), col="red", main=title)
lines(density(Cseg90_data$sCnuc), col="blue")
lines(density(Cseg95_data$sCnuc), col="green")
legend(-0.1, 1.5, legend=c("Raw", "Cseg90", "Cseg95"),
      col=c("red", "blue", "green"), lty=1:2, cex=0.8)

png(outFile_sCord)
title <- paste0(sample_id, " density plot of sCord scores")
plot(density(raw_data$sCord), col="red", main=title)
lines(density(Cseg90_data$sCord), col="blue")
lines(density(Cseg95_data$sCord), col="green")
legend(-0.1, 1.5, legend=c("Raw", "Cseg90", "Cseg95"),
       col=c("red", "blue", "green"), lty=1:2, cex=0.8)

png(outFile_sCcov)
title <- paste0(sample_id, " density plot of sCcov scores")
plot(density(raw_data$sCcov), col="red", main=title)
lines(density(Cseg90_data$sCcov), col="blue")
lines(density(Cseg95_data$sCcov), col="green")
legend(-0.1, 1.5, legend=c("Raw", "Cseg90", "Cseg95"),
      col=c("red", "blue", "green"), lty=1:2, cex=0.8)
