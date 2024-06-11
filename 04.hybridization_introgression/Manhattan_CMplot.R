#!/public/software/apps/R/4.1.0/bin/R
library(CMplot)
cmplot.data <- read.table('ZJ_GS_WP_localFstats__2500_500.Z.raw.cmplot',header=T,sep='\t')
cmplot.new.data <- cmplot.data[cmplot.data$D >= 0,]
CMplot(cmplot.new.data, plot.type="m", col=c("#FF595E","#FFCA3A","#8AC926","#1982C4","#6A4C93"), LOG10=FALSE, ylim=c(0,0.35), threshold=0.121,
        threshold.lty=1, threshold.lwd=1, threshold.col=c("black"), amplify=F, ylab='D',
        chr.den.col=NULL, signal.col=c("red"), signal.cex=c(1),signal.pch=c(19),
        file="jpg",memo="ZJ_GS_WP_Z",dpi=300,file.output=TRUE,verbose=TRUE)
