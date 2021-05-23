
logger <- function(a, b=NULL, c=NULL, d=NULL, e=NULL, f=NULL){
  if (missing(b)){
    pars <- c(a)
  } else if (missing(c)){
    pars <- c(a, b)
  } else if (missing(d)){
    pars <- c(a, b, c)
  } else if (missing(e)){
    pars <- c(a, b, c , d)
  } else if (missing(f)){
    pars <- c(a, b, c, d, e)
  } else {
    pars <- c(a, b, c, d, e, f)
  }

  str = paste(pars, sep = " ")
  cat(format(Sys.time(), "%a %b %d %X %Y"), str, "\n", sep = ": ")
}

mlinha <- function(m, df){
  if (length(m) == 0){
    return(c())
  }
  if (length(m) > 1){
    ret <- which(rowSums(df[,m]) == length(m))
  } else {
    ret <- which(df[,m] == length(m))
  }
  return(ret)
}


subsets <- function(x){
  combs <- list()
  if (length(x) == 0){
    return(combs)
  }

  for (i in 1:length(x)){
    comb <- combn(x, i, simplify = FALSE)
    combs <- c(combs, comb)
  }
  return(combs)
}



coveringEdges <- function(X, Y, df){
  nc <- length(X)
  M <- seq(1,ncol(df),1)

  keys <- unlist(lapply(X, paste, collapse="-"))
  keys[which(keys == "")] <- "0"
  keys <- unlist(lapply(keys, digest, algo = 'md5'))
  h <- hash(keys=keys, values=1:length(combs))

  #for (i in 1:nc){
  out <- foreach (i = 1:nc, .export=c("mlinha", "has.key", "digest")) %dopar%{
    #cat("i: ", i,"\n")
    counts <- rep(0, nc)

    edgesi <- list()
    m <- setdiff(M, Y[[i]])
    if (length(m) > 0){
      for (j in 1:length(m)){
        #cat("j: ", j,"\n")
        #cat("mj: ", m[[j]],"\n")
        inters <- as.numeric(intersect(X[[i]], mlinha(m[[j]],df)))
        #cat("inters: ", inters,"\n")
        # tt <- tic(c())
        # X_inters <- which(unlist(lapply(X, function(x){identical(x, inters)})))
        # tt <- tic(tt)
        # tac(tt)
        if (length(inters) > 0){
          inters <- digest(paste(inters, collapse = "-"), algo = 'md5')
          if (has.key(inters, h)){
            X_inters <- h[inters]
            X_inters <- as.list(X_inters)[[1]]
            counts[X_inters] <- counts[X_inters] + 1
            #cat("counts: ", counts, "\n")
            #cat("Yi: ", Y[[i]], "\n")
            #cat("Yxin: ", Y[[X_inters]], "\n")
            #cat("counts: ", counts, "\n")
            #print((length(Y[[X_inters]]) - length(Y[[i]])) == counts[X_inters])
            if ((length(Y[[X_inters]]) - length(Y[[i]])) == counts[X_inters]){
              #cat("Edge: ", X[[X_inters]],":", Y[[X_inters]], " -> ", X[[i]],":", Y[[i]], "\n")
              edgesi[[length(edgesi)+1]] <- c(i,X_inters)
            }
          }
        }
      }
    }
    edgesi
  }
  edges <- list()
  for(i in 1:length(out)){
    edges <- c(edges,  out[[i]])
  }
  return(edges)
}


arffH_par <- function(dir, fileArffH, arff, arfft, arffv, labelFirst, labelLast, classstr, inss, combs, root, att_types){
  linesToWrite = c()

  arffAtts <- arff[,1:(labelFirst - 1)]
  arffLabels <- arff[,labelFirst:labelLast]
  arffAttst <- arfft[,1:(labelFirst - 1)]
  arffLabelst <- arfft[,labelFirst:labelLast]
  arffAttsv <- arffv[,1:(labelFirst - 1)]
  arffLabelsv <- arffv[,labelFirst:labelLast]


  # cabecalho @relation
  linesToWrite <- c(linesToWrite, "@relation HierarchicalFrom_Flat")

  # atributos
  for (i in 1:ncol(arffAtts)){
    linesToWrite <- c(linesToWrite, paste("@attribute ", colnames(arffAtts)[i], " ", att_types[i], sep=""))
  }

  # class hierarchical
  linesToWrite <- c(linesToWrite, classstr)

  #dados
  linesToWrite <- c(linesToWrite, "@Data")

  linesToWriteT <- linesToWrite
  linesToWriteV <- linesToWrite


  str <- predm <- foreach (i = 1:nrow(arff)) %dopar%{
    #for (i in 1:nrow(dftr)){
    attsi <- paste(paste(arffAtts[i,], collapse = ","))
    classesi <- which(unlist(lapply(inss, function(x){length(which(i %in% x)) > 0})))

    # teste: pegar apenas o nó mais profundo da DAG - não muda o resultado
    classesi <- classesi[which.max(lengths(combs[classesi]))]

    # instancias sem classe positiva precisam ser associas à classe root (birds e yelp)
    if (length(classesi) == 0){
      classesi <- root
    }
    classesi <- paste(classesi, collapse = "@", sep = "")
    line <- paste(attsi, classesi, sep = ",")
    line
  }
  linesToWrite <- c(linesToWrite, unlist(str))


  str <- predm <- foreach (i = 1:nrow(arfft)) %dopar%{
    #for (i in 1:nrow(dfte)){
    attsi <- paste(paste(arffAttst[i,], collapse = ","))
    classesi <- which(unlist(lapply(inss, function(x){length(which(i %in% x)) > 0})))

    classesi <- classesi[which.max(lengths(combs[classesi]))]

    # instancias sem classe positiva precisam ser associas à classe root (birds e yelp)
    if (length(classesi) == 0){
      classesi <- root
    }


    classesi <- paste(classesi, collapse = "@", sep = "")
    line <- paste(attsi, classesi, sep = ",")
    line
  }
  linesToWriteT <- c(linesToWriteT, unlist(str))


  str <- predm <- foreach (i = 1:nrow(arffv)) %dopar%{
    #for (i in 1:nrow(dfva)){
    attsi <- paste(paste(arffAttsv[i,], collapse = ","))
    classesi <- which(unlist(lapply(inss, function(x){length(which(i %in% x)) > 0})))

    classesi <- classesi[which.max(lengths(combs[classesi]))]

    # instancias sem classe positiva precisam ser associas à classe root (birds e yelp)
    if (length(classesi) == 0){
      classesi <- root
    }

    classesi <- paste(classesi, collapse = "@", sep = "")
    line <- paste(attsi, classesi, sep = ",")
  }
  linesToWriteV <- c(linesToWriteV, unlist(str))


  filename <- paste(dir,gsub(".arff", "_train.arff", fileArffH), sep="")
  filenameT <- paste(dir,gsub(".arff", "_test.arff", fileArffH), sep="")
  filenameV <- paste(dir,gsub(".arff", "_valid.arff", fileArffH), sep="")

  fileConn<-file(filename)
  writeLines(paste(linesToWrite, collapse = "\n"), fileConn)
  close(fileConn)
  fileConn<-file(filenameT)
  writeLines(paste(linesToWriteT, collapse = "\n"), fileConn)
  close(fileConn)
  fileConn<-file(filenameV)
  writeLines(paste(linesToWriteV, collapse = "\n"), fileConn)
  close(fileConn)

  ret <- list()
  return(ret)

}



geraConfClusHMC<- function(sfilename, trainfile, testfile, validfile, WType, WParam, OptimizeErrorMeasure){

  s <- "[Data]
File = XXX1
TestSet = XXX2

%using the valid dataset achieved low performance
%PruneSet = XXX6

[Hierarchical]
Type = DAG
WType = XXX3
WParam = XXX4
HSeparator = /
ClassificationThreshold = [0,2,4,6,8,10,12,14,16,18,20,22,24,26,28,30,32,34,36,38,40,42,44,46,48,50,52,54,56,58,60,62,64,66,68,70,72,74,76,78,80,82,84,86,88,90,92,94,96,98,100]
OptimizeErrorMeasure = XXX5

[Model]
MinimalWeight = 5.0

[Attributes]
ReduceMemoryNominalAttrs = yes

[Tree]
ConvertToRules = No
FTest = [0.001,0.005,0.01,0.05,0.1,0.125]


[Ensemble]
EnsembleMethod = Bagging
Iterations = 10

[Output]
WritePredictions = {Train, Test}"

  print(s)

  s <- gsub("XXX1", trainfile, s)
  s <- gsub("XXX2", testfile, s)
  s <- gsub("XXX3", WType, s)
  s <- gsub("XXX4", WParam, s)
  s <- gsub("XXX5", OptimizeErrorMeasure, s)
  s <- gsub("XXX6", validfile, s)

  fileConn<-file(sfilename)
  writeLines(s, fileConn)
  close(fileConn)

}


concToClass <- function(inss, conclist){
  ret <- unique(sort(unlist(inss[as.numeric(conclist)])))
  return(ret)
}

tic <- function(times, label = "auto"){
  times[length(times) + 1] <- Sys.time()
  names(times)[length(times)] <- label
  logger("Finished", label, times[length(times)] - times[max(length(times) - 1,1)])
  return(times)
}

tac <- function(times){
  timesf <- as.data.frame(matrix(0, ncol = 2, nrow = length(times) + 1))
  timesf[1,] <- c(names(times)[1], 0)

  for(i in 2:length(times)){
    timesf[i,] <- c(names(times)[i], times[i] - times[i -1])
  }
  timesf[nrow(timesf),] <- c("Total time", times[length(times)] - times[1])
  return(timesf)
}


measureHtoFlat <- function(predarff, combs, truem_colnames){

  dfpred <- predarff[,which(startsWith(colnames(predarff), "Original"))]
  dfpred <- dfpred[,-ncol(dfpred)]

  colnames(dfpred) <- as.numeric(gsub(x = colnames(dfpred), pattern = "Original-p-", replacement = ""))

  predm <- matrix(0, nrow = nrow(dfpred), ncol = length(truem_colnames))
  for (i in 1:ncol(dfpred)){
    metai <- unlist(combs[as.numeric(colnames(dfpred)[i])])
    for (m in metai){
      menores <- which((predm[,m] - dfpred[,i]) < 0)
      predm[menores, m] <- dfpred[menores, i]
    }
  }
  colnames(predm) <- truem_colnames

  return(predm)
}




range01 <- function(x){
  (x-min(x))/(max(x)-min(x))
}



compute_predictions_prob <- function(truemte, predarff, combs, id){

  original <- data.frame(truemte)
  original <- data.frame(sapply(original, function(x) as.numeric(as.character(x))))
  mymldr <- NULL
  mymldr <- mldr_from_dataframe(original, labelIndices = seq(1,ncol(truemte)), name = "original")

  # obtem as probabilidades das predicoes
  times <- tic(times, paste("Starting matrix thresholds H->Flat on data ", id, sep = ""))
  predmte <- measureHtoFlat(predarff, combs, colnames(truemte))
  times <- tic(times, paste("Finished matrix thresholds H->Flat on data ", id, sep = ""))

  write.csv(truemte, paste("true-", id, ".csv", sep = ""), row.names = FALSE)
  write.csv(predmte, paste("pred-", id, ".csv", sep = ""), row.names = FALSE)


}


compute_results_t1 <- function(id){

  true <- read.csv(paste("true-", id, ".csv", sep=""))
  true <- data.frame(sapply(true, function(x) as.numeric(as.character(x))))
  truemldr <- mldr_from_dataframe(true, labelIndices = seq(1,ncol(true)), name = "true")

  pred <- read.csv(paste("pred-", id, ".csv", sep=""))
  pred <- data.frame(sapply(pred, function(x) as.numeric(as.character(x))))

  # Thresholding t1 - simple range
  # faz o range
  for (i in 1:nrow(pred)){
    if (sum(pred[i,]) > 0 && var(as.numeric(pred[i,])) != 0){
      pred[i,] <- range01(pred[i,])
    }
  }

  write.csv(pred, paste("pred-", id, "-t1.csv",sep = ""), row.names = FALSE)

  result <- NULL
  result <- multilabel_evaluate(truemldr, pred)
  for (i in 1:length(result)){
    line <- paste(names(result)[i], ": ", result[i])
    logger(paste("RESULT-", fileid, "t1-", id,":", sep=""), line)
  }
  write.csv(result, paste("results", id, dsname, fold, fileid, "t1.csv", sep = "-"))


}


compute_results_t2 <- function(id){

  true <- read.csv(paste("true-", id, ".csv", sep=""))
  true <- data.frame(sapply(true, function(x) as.numeric(as.character(x))))
  truemldr <- mldr_from_dataframe(true, labelIndices = seq(1,ncol(true)), name = "true")

  # precisa do label space de treinamento e validacao para calculo do threshold
  truetr <- read.csv(paste("true-", "tr", ".csv", sep=""))
  trueva <- read.csv(paste("true-", "va", ".csv", sep=""))

  pred <- read.csv(paste("pred-", id, ".csv", sep=""))
  pred <- data.frame(sapply(pred, function(x) as.numeric(as.character(x))))

  # precisa das predicoes da validacao para calibrar o threshold
  predva <- read.csv(paste("pred-", "va", ".csv", sep=""))
  predva <- data.frame(sapply(predva, function(x) as.numeric(as.character(x))))

  # Thresholding t2 - One Threshold
  # calibrating the threshold
  lcard <- (sum(truetr) + sum(trueva))/(nrow(truetr) + nrow(trueva))
  t <- seq(0.01,0.99, 0.01)
  bt <- 0
  bts <- 1
  for (i in 1:length(t)){
    s <- abs(lcard - (length(which(predva >= t[i]))/nrow(predva)))
    if (s < bts){
      bt <- t[i]
      bts <- s
    }
  }

  # end calibration

  pred[pred >= bt] <- 1
  pred[pred < bt] <- 0
  write.csv(pred, paste("pred-", id, "-t2.csv",sep = ""), row.names = FALSE)
  write.csv(bt, paste("pred-", id, "-t2-threshold.csv",sep = ""), row.names = FALSE)

  result <- NULL
  result <- multilabel_evaluate(truemldr, pred)
  for (i in 1:length(result)){
    line <- paste(names(result)[i], ": ", result[i])
    logger(paste("RESULT-", fileid, "t2-", id,":", sep=""), line)
  }
  write.csv(result, paste("results", id, dsname, fold, fileid, "t2.csv", sep = "-"))



}



compute_results_t3 <- function(id){

  true <- read.csv(paste("true-", id, ".csv", sep=""))
  true <- data.frame(sapply(true, function(x) as.numeric(as.character(x))))
  truemldr <- mldr_from_dataframe(true, labelIndices = seq(1,ncol(true)), name = "true")

  # precisa do label space de treinamento e validacao para calculo do threshold
  truetr <- read.csv(paste("true-", "tr", ".csv", sep=""))
  trueva <- read.csv(paste("true-", "va", ".csv", sep=""))

  pred <- read.csv(paste("pred-", id, ".csv", sep=""))
  pred <- data.frame(sapply(pred, function(x) as.numeric(as.character(x))))

  # precisa das predicoes da validacao para calibrar o threshold
  predva <- read.csv(paste("pred-", "va", ".csv", sep=""))
  predva <- data.frame(sapply(predva, function(x) as.numeric(as.character(x))))

  # Thresholding t3 - One Threshold per class
  # calibrating the threshold
  lcardpc <- (colSums(trueva) + colSums(truetr))/(nrow(trueva) + nrow(truetr))

  # thresholding cc
  t <- seq(0.01,0.99, 0.01)
  tf <- c()
  for (c in 1:ncol(true)){
    bt <- 0
    bts <- 1
    for (i in 1:length(t)){
      s <- abs(lcardpc[c] - (length(which(predva[c] >= t[i]))/nrow(predva)))
      if (s < bts){
        bt <- t[i]
        bts <- s
        #print(bt)
        #print(bts)
      }
    }
    tf <- c(tf, bt)
  }
  # end calibration

  for (c in 1:ncol(pred)){
    pos1 <- which(pred[,c] >= tf[c])
    pos0 <- which(pred[,c] < tf[c])
    pred[pos1, c] <- 1
    pred[pos0, c] <- 0
  }

  write.csv(pred, paste("pred-", id, "-t3.csv",sep = ""), row.names = FALSE)
  write.csv(tf, paste("pred-", id, "-t3-threshold.csv",sep = ""), row.names = FALSE)

  result <- NULL
  result <- multilabel_evaluate(truemldr, pred)
  for (i in 1:length(result)){
    line <- paste(names(result)[i], ": ", result[i])
    logger(paste("RESULT-", fileid, "t3-", id,":", sep=""), line)
  }
  write.csv(result, paste("results", id, dsname, fold, fileid, "t3.csv", sep = "-"))
}

findClusJar <- function(){
  for (path in .libPaths()){
    pathclus <- file.path(paste(path, "/F2H/java/MyClus.jar", sep=""))
    if (file.exists(pathclus)){
      return(pathclus)
    }
  }

  # TODO: message
  return(NULL)
}


F2H <- function(
  dsname = "yeast",
  dsdir = tempdir(),
  javaExe = "java",
  javaMem = "-Xmx3g",
  clusJar = findClusJar(),
  minSupportConcetps = 0,
  clusWType = "ExpMaxParentWeight",
  clusWParam = 0.8,
  clusOptimizeErrorMeasure = "WeightedAverageAUPRC",
  threads = 1,
  ensembleClus = 0
){
  # define some input vars
  clusExe <- paste(javaExe, javaMem, clusJar, sep = " ")

  times <- c()
  times <- tic(times, "Start")





}


