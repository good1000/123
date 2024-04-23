#
library(limma)
library(oncoPredict)
library(parallel)
set.seed(12345)

expFile="symbol.txt"    #

rt = read.table(expFile, header=T, sep="\t", check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp), colnames(exp))
data=matrix(as.numeric(as.matrix(exp)), nrow=nrow(exp), dimnames=dimnames)
data=avereps(data)
data=data[rowMeans(data)>0.5,]

#
GDSC2_Expr=readRDS(file='GDSC2_Expr.rds')
GDSC2_Res=readRDS(file = 'GDSC2_Res.rds')
GDSC2_Res=exp(GDSC2_Res) 

#
calcPhenotype(trainingExprData = GDSC2_Expr,    
              trainingPtype = GDSC2_Res,        
              testExprData = data,              
              batchCorrect = 'eb',              
              powerTransformPhenotype = TRUE,
              removeLowVaryingGenes = 0.2,      
              minNumSamples = 10,               
              printOutput = TRUE,               
              removeLowVaringGenesFrom = 'rawData')






#
library(limma)
library(ggplot2)
library(ggpubr)

pFilter=0.001                      
riskFile="risk.all.txt"           
drugFile="DrugPredictions.csv"     

#
risk=read.table(riskFile, header=T, sep="\t", check.names=F, row.names=1)

#
senstivity=read.csv(drugFile, header=T, sep=",", check.names=F, row.names=1)
colnames(senstivity)=gsub("(.*)\\_(\\d+)", "\\1", colnames(senstivity))

Group = ifelse(as.numeric(str_sub(rownames(senstivity),14,15)) < 10,'tumor','normal')
senstivity = senstivity[Group == "tumor",]

truncated_rownames <- stringr::str_sub(rownames(senstivity), 1, 12)
duplicated_rows <- duplicated(truncated_rownames)
senstivity <- senstivity[!duplicated_rows,]
rownames(senstivity) = str_sub(rownames(senstivity),1,12)


#
sameSample=intersect(row.names(risk), row.names(senstivity))
risk=risk[sameSample, "risk",drop=F]
senstivity=senstivity[sameSample,,drop=F]
rt=cbind(risk, senstivity)

#
rt$risk=factor(rt$risk, levels=c("low", "high"))
type=levels(factor(rt[,"risk"]))
comp=combn(type, 2)
my_comparisons=list()
for(i in 1:ncol(comp)){my_comparisons[[i]]<-comp[,i]}

#, 
for(drug in colnames(rt)[2:ncol(rt)]){
  rt1=rt[,c(drug, "risk")]
  colnames(rt1)=c("Drug", "risk")
  rt1=na.omit(rt1)
  rt1$Drug=log2(rt1$Drug+1)
  #
  test=wilcox.test(Drug ~ risk, data=rt1)
  diffPvalue=test$p.value
  #
  if(diffPvalue<pFilter){
    boxplot=ggboxplot(rt1, x="risk", y="Drug", fill="risk",
                      xlab="risk",
                      ylab=paste0(drug, " senstivity"),
                      legend.title="risk",
                      palette=c("green", "red")
    )+ 
      stat_compare_means(comparisons=my_comparisons)
    #
    pdf(file=paste0("drugSenstivity.", drug, ".pdf"), width=5, height=4.5)
    print(boxplot)
    dev.off()
  }
}


