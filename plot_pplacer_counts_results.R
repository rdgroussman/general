library(ggplot2)
library(scales)

#upload the count file and the latitude file
lat = read.table(file='/Users/rgroussman/data/SCOPE/Gradients1/gradients_latitude.txt', sep='\t', header=TRUE, row.names=1)

# couple little runs here to test:

# ISIP2a:FTN ratio
data = read.table(file='/Users/rgroussman/data/SCOPE/Gradients1/runs/ISIP2a/ISIP2a-ftn.counts_results.csv', sep=',', header=TRUE, row.names=1)

# cobS
data = read.table(file='/Users/rgroussman/data/SCOPE/Gradients1/runs/ISIP2a/ISIP2a-ftn.counts_results.csv', sep=',', header=TRUE, row.names=1)

# 6 RPs
data = read.table(file='/Users/rgroussman/data/ribosomal_proteins/6RP.G1.counts_results.csv', sep=',', header=TRUE, row.names=1)

comb = merge(lat, data, 'row.names')
small = comb[comb$size_fraction==0.2,]
large = comb[comb$size_fraction==3.0,]


multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)

  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)

  numPlots = length(plots)

  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                    ncol = cols, nrow = ceiling(numPlots/cols))
  }

 if (numPlots==1) {
    print(plots[[1]])

  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))

    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))

      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}

#make plots for different taxa
S1 = ggplot() + geom_point(data=small, aes(small$latitude, small$NonDiatomStramenopile), col="purple", shape=10) + geom_smooth(data=small, aes(small$latitude, small$NonDiatomStramenopile), se=FALSE, col="purple", linetype='dotted') + geom_point(data=large, aes(large$latitude, large$NonDiatomStramenopile), col="green3") + geom_smooth(data=large, aes(large$latitude, large$NonDiatomStramenopile), se=FALSE, col="green3") + scale_y_log10() + annotation_logticks(sides='l') + ylab("Non-Diatom Stramenopiles") + xlab("") + ggtitle('ISIP2a:FTN ratios')
S2 = ggplot() + geom_point(data=small, aes(small$latitude, small$Bacillariophyta), col="purple", shape=10) + geom_smooth(data=small, aes(small$latitude, small$Bacillariophyta), se=FALSE, col="purple", linetype='dotted') + geom_point(data=large, aes(large$latitude, large$Bacillariophyta), col="green3") + geom_smooth(data=large, aes(large$latitude, large$Bacillariophyta), se=FALSE, col="green3") + scale_y_log10() + annotation_logticks(sides='l') + ylab("Diatoms") + xlab("")
S3 = ggplot() + geom_point(data=small, aes(small$latitude, small$NonDinoAlveolata), col="purple", shape=10) + geom_smooth(data=small, aes(small$latitude, small$NonDinoAlveolata), se=FALSE, col="purple", linetype='dotted') + geom_point(data=large, aes(large$latitude, large$NonDinoAlveolata), col="green3") + geom_smooth(data=large, aes(large$latitude, large$NonDinoAlveolata), se=FALSE, col="green3") + scale_y_log10() + annotation_logticks(sides='l') + ylab("Non-Dinoflagellate Alveolates") + xlab("")
S4 = ggplot() + geom_point(data=small, aes(small$latitude, small$Cryptophyta), col="purple", shape=10) + geom_smooth(data=small, aes(small$latitude, small$Cryptophyta), se=FALSE, col="purple", linetype='dotted') + geom_point(data=large, aes(large$latitude, large$Cryptophyta), col="green3") + geom_smooth(data=large, aes(large$latitude, large$Cryptophyta), se=FALSE, col="green3") + scale_y_log10() + annotation_logticks(sides='l') + ylab("Cryptophytes") + xlab("") + ggtitle("")
S5 = ggplot() + geom_point(data=small, aes(small$latitude, small$Dinophyceae), col="purple", shape=10) + geom_smooth(data=small, aes(small$latitude, small$Dinophyceae), se=FALSE, col="purple", linetype='dotted') + geom_point(data=large, aes(large$latitude, large$Dinophyceae), col="green3") + geom_smooth(data=large, aes(large$latitude, large$Dinophyceae), se=FALSE, col="green3") + scale_y_log10() + annotation_logticks(sides='l') + ylab("Dinoflagellates") + xlab("")
S6 = ggplot() + geom_point(data=small, aes(small$latitude, small$Haptophyta), col="purple", shape=10) + geom_smooth(data=small, aes(small$latitude, small$Haptophyta), se=FALSE, col="purple", linetype='dotted') + geom_point(data=large, aes(large$latitude, large$Haptophyta), col="green3") + geom_smooth(data=large, aes(large$latitude, large$Haptophyta), se=FALSE, col="green3") + scale_y_log10() + annotation_logticks(sides='l') + ylab("Haptophytes") + xlab("")
S7 = ggplot() + geom_point(data=small, aes(small$latitude, small$Viridiplantae), col="purple", shape=10) + geom_smooth(data=small, aes(small$latitude, small$Viridiplantae), se=FALSE, col="purple", linetype='dotted') + geom_point(data=large, aes(large$latitude, large$Viridiplantae), col="green3") + geom_smooth(data=large, aes(large$latitude, large$Viridiplantae), se=FALSE, col="green3") + scale_y_log10() + annotation_logticks(sides='l') + ylab("Chlorophytes") + xlab("")

multiplot(S1,S2,S3,S4,S5,S7,cols=2)

multiplot(S3,S4,cols=2)

