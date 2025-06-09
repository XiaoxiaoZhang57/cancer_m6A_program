library("qvalue")
library("Matrix")
library("igraph")
library("parallel")
setwd("~/polysel")
source('./R/polysel.R')
project.name<-"24l34s"
data.path<-file.path("./data")
code.path<-"./R"
empfdr.path<-"./empfdr"
results.path<-"./results"
source(file.path(code.path,'polysel.R'))
minsetsize<-10
result<-ReadSetObjTables(in.path=data.path,
                         set.info.file="SetInfo.txt",
                         set.obj.file="SetObj.txt",
                         obj.info.file="ObjInfo.txt",
                         minsetsize=minsetsize,
                         obj.in.set=F,
                         merge.similar.sets=T)
set.info<-result$set.info
obj.info<-result$obj.info
set.obj<-result$set.obj
set.info.lnk<-result$set.info.lnk
cat("Number of sets: ", nrow(set.info), "\n", sep="")
cat("Number of genes: ", nrow(obj.info), "\n", sep="")
print(set.info[1:5,],row.names=F, right=F)
print(obj.info[1:5,],row.names=F, right=F)
print(set.obj[1:5,],row.names=F, right=F)
flds<-c("genelength")
options(digits=2)
for (fld in flds){
  cat("Correlation between objStat and ", fld, ": ", sep="")
  ct<-(cor.test(obj.info$objStat,obj.info[[fld]]))
  cat("p-value: ", ct$p.value, " estimate: ", ct$estimate, 
      " [", ct$conf.int[1],
      ", ", ct$conf.int[2],"]\n",sep="")
  PlotStatField(obj.info, fld=fld, ylab=fld, xlab="objStat", 
                logaxis="y", show.bins=F)
}
obj.stat<-obj.info[,c("objID", "objStat", "objBin")]
save(set.info, obj.info, obj.stat, set.obj, set.info.lnk,
     file=file.path(data.path, "polysel_objects.RData"))


approx.null <- FALSE
use.bins <- FALSE
seq.rnd.sampling <- TRUE
nrand <- 400000
# later change to high value
test <- "highertail"
#(((qvalue.method <- "smoother" )))
qvalue.method <- "bootstrap"

# Test sets
result<-EnrichmentAnalysis(set.info, set.obj, obj.stat,
                           nrand=nrand, approx.null=approx.null, 
                           seq.rnd.sampling=seq.rnd.sampling,
                           use.bins=use.bins, test=test,
                           do.pruning=FALSE, minsetsize=minsetsize,
                           project.txt=project.name, do.emp.fdr=FALSE,
                           qvalue.method=qvalue.method
)

set.scores.prepruning <- result$set.scores.prepruning

print(set.scores.prepruning[1:10,],row.names=F, right=F)
save(set.scores.prepruning, 
     file = file.path(results.path,
                      paste(project.txt=project.name,
                            "_setscores_prepruning_", 
                            formatC(nrand,format="d"),".RData", sep="")))

write.table(set.scores.prepruning, quote=FALSE, sep="\t", row.names=FALSE, 
            file = file.path(results.path, paste(project.txt=project.name, 
                                                 "_setscores_prepruning_", formatC(nrand,format="d"),".txt",
                                                 sep=""))) 

# Test sets
result<-EnrichmentAnalysis(set.info, set.obj, obj.stat,
                           nrand=nrand, approx.null=approx.null, 
                           seq.rnd.sampling=seq.rnd.sampling,
                           use.bins=use.bins, test=test,
                           do.pruning=TRUE, minsetsize=minsetsize,
                           project.txt=project.name, do.emp.fdr=FALSE,
                           qvalue.method=qvalue.method
)

set.scores.postpruning <- result$set.scores.postpruning

print(set.scores.postpruning[1:10,],row.names=F, right=F)
#zhege 200-300
emp.fdr.nruns <- 200
run <- function(r){
  if (file.exists(paste0(empfdr.path,project.name,"_shuf",r,"_",formatC(nrand,format="d"),".RData"))) {
    return()
  }
  TestShuffledSets(set.info, set.obj, obj.stat, 
                   nrand=nrand, minsetsize=minsetsize, 
                   approx.null=approx.null, test=test,
                   seq.rnd.sampling=seq.rnd.sampling,
                   use.bins=use.bins, out.path=empfdr.path, 
                   runnr=r, project.txt=project.name)
}

junk <- mclapply(1:emp.fdr.nruns,run,mc.cores = 50)
emp.fdr.est.m0<-TRUE
set.scores.postpruning <- GetEmpericalFDR(set.scores=set.scores.postpruning,
                                          in.path=empfdr.path, 
                                          nruns=emp.fdr.nruns, 
                                          est.m0=emp.fdr.est.m0,
                                          project.txt=project.name,
                                          nrand=nrand)
print(set.scores.postpruning[1:10,],row.names=F, right=F)
save(set.scores.postpruning, 
     file = file.path(results.path,
                      paste(project.txt=project.name,
                            "_setscores_postpruning_", 
                            formatC(nrand,format="d"),"_Padj.RData", sep="")))

write.table(set.scores.postpruning, quote=FALSE, sep="\t", row.names=FALSE, 
            file = file.path(results.path, paste(project.txt=project.name, 
                                                 "_setscores_postpruning_", formatC(nrand,format="d"),"_Padj.txt",
                                                 sep=""))) 
