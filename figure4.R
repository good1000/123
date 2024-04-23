


library(pheatmap)        #

#
bioRiskPlot=function(inputFile=null, project=null){
  rt=read.table(inputFile, header=T, sep="\t", check.names=F, row.names=1)    #
  rt=rt[order(rt$riskScore),]      #
  
  #
  riskClass=rt[,"risk"]
  lowLength=length(riskClass[riskClass=="low"])
  highLength=length(riskClass[riskClass=="high"])
  lowMax=max(rt$riskScore[riskClass=="low"])
  line=rt[,"riskScore"]
  line[line>10]=10
  tiff(file=paste0(project, ".riskScore.tiff"), width=7, height=4, units="in", res=300)
  plot(line, type="p", pch=20,
       xlab="Patients (increasing risk socre)",
       ylab="Risk score",
       col=c(rep("#33FF00",lowLength),rep("#FF0000",highLength)) )
  abline(h=lowMax,v=lowLength,lty=2)
  legend("topleft", c("High risk","Low Risk"),bty="n",pch=19,col=c("#FF0000","#33FF00"),cex=1.2)
  dev.off()
  
  #
  color=as.vector(rt$fustat)
  color[color==1]="#FF0000"     #
  color[color==0]="#33FF00"     #
  tiff(file=paste0(project, ".survStat.tiff"), width=7, height=4, units="in", res=300)
  plot(rt$futime, pch=19,
       xlab="Patients (increasing risk socre)",
       ylab="PFS time (years)",
       col=color)
  legend("topleft", c("progression","No progression"),bty="n",pch=19,col=c("#FF0000","#33FF00"),cex=1.2)
  abline(v=lowLength,lty=2)
  dev.off()
  
  #
  ann_colors=list()
  bioCol=c("#33FF00", "#FF0000")
  names(bioCol)=c("low", "high")
  ann_colors[["Risk"]]=bioCol
  
  #
  rt1=rt[c(3:(ncol(rt)-2))]
  rt1=t(rt1)
  annotation=data.frame(Risk=rt[,ncol(rt)])
  rownames(annotation)=rownames(rt)
  tiff(file=paste0(project, ".heatmap.tiff"), width=7, height=4, units="in", res=300)
  pheatmap::pheatmap(rt1, 
                     annotation=annotation,
                     annotation_colors = ann_colors, 
                     cluster_cols = FALSE,
                     cluster_rows = FALSE,
                     show_colnames = F,
                     scale="row",
                     color = colorRampPalette(c(rep("#33FF00",3.5), "white", rep("#FF0000",3.5)))(50),
                     fontsize_col=3,
                     fontsize=7,
                     fontsize_row=8)
  dev.off()
}

#
bioRiskPlot(inputFile="risk.train.txt", project="train")
#
bioRiskPlot(inputFile="risk.test.txt", project="test")
#
bioRiskPlot(inputFile="risk.all.txt", project="all")

