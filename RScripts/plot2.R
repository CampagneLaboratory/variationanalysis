library("AUC")
modelDir <- "/Users/fac2003/IdeaProjects/variationanalysis/tests/1465684020704"
rm(test)
rm(training)
results<-c()
for (kind in c("best","latest")) {
  for (type in c("training","test")) {
    filename<-paste(modelDir,paste(kind,"-",sep="",type),sep="/");
    cat(paste("Processing",filename),sep="\n");

    training <- read.delim(file=filename, header=TRUE, nrows=-1)
    t<-training
    rocValues <-roc(predictions=t$ProbabilityMut,labels=factor(t$mutatedLabel))
    aucValue<-auc(rocValues)
    results <- append(results,c(c(c(aucValue,rocValues, paste(kind, "on", type)))))
    cat(paste(kind, "on", type,"AUC=",aucValue),sep="\n")
    par(new=TRUE)
    color<-sample(colours(), 50)
     plot(rocValues,col=color,  legend(x=0.2,y=0.2,legend=paste(kind, "on", type),col=color))
  }
}
