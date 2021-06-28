.onAttach <- function(libname, pkgname) {
  options(java.parameters = "-Xmx2g")
}

logger <- function(a, b=NULL, c=NULL, d=NULL, e=NULL, f=NULL){

  pars <- c(a, b, c, d, e, f)

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


coveringEdges <- function(X, Y, df, threads = 1){

  clusters <- parallel::makeCluster(threads)
  doParallel::registerDoParallel(clusters)

  nc <- length(X)
  M <- seq(1,ncol(df),1)

  keys <- unlist(lapply(X, paste, collapse="-"))
  keys[which(keys == "")] <- "0"
  keys <- unlist(lapply(keys, digest, algo = 'md5'))
  h <- hash(keys=keys, values=1:length(Y))

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

  parallel::stopCluster(clusters)

  return(edges)
}


arffH_par <- function(dir, fileArffH, arff, arfft, arffv, labelFirst, labelLast, classstr, inss, combs, root, att_types, threads = 1){

  clusters <- parallel::makeCluster(threads)
  doParallel::registerDoParallel(clusters)

  linesToWrite = c()

  arffAtts <- arff[,1:(labelFirst - 1)]
  arffLabels <- arff[,labelFirst:labelLast]
  arffAttst <- arfft[,1:(labelFirst - 1)]
  arffLabelst <- arfft[,labelFirst:labelLast]
  #arffAttsv <- arffv[,1:(labelFirst - 1)]
  #arffLabelsv <- arffv[,labelFirst:labelLast]


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
  #linesToWriteV <- linesToWrite


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
  filename <- paste(dir,gsub(".arff", "_train.arff", fileArffH), sep="")
  fileConn<-file(filename)
  writeLines(linesToWrite, fileConn)
  close(fileConn)
  rm(linesToWrite)


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
  filenameT <- paste(dir,gsub(".arff", "_test.arff", fileArffH), sep="")
  fileConn<-file(filenameT)
  writeLines(linesToWriteT, fileConn)
  close(fileConn)
  rm(linesToWriteT)

  # str <- predm <- foreach (i = 1:nrow(arffv)) %dopar%{
  #   #for (i in 1:nrow(dfva)){
  #   attsi <- paste(paste(arffAttsv[i,], collapse = ","))
  #   classesi <- which(unlist(lapply(inss, function(x){length(which(i %in% x)) > 0})))
  #
  #   classesi <- classesi[which.max(lengths(combs[classesi]))]
  #
  #   # instancias sem classe positiva precisam ser associas à classe root (birds e yelp)
  #   if (length(classesi) == 0){
  #     classesi <- root
  #   }
  #
  #   classesi <- paste(classesi, collapse = "@", sep = "")
  #   line <- paste(attsi, classesi, sep = ",")
  # }
  # linesToWriteV <- c(linesToWriteV, unlist(str))
  #
  # filenameV <- paste(dir,gsub(".arff", "_valid.arff", fileArffH), sep="")
  # fileConn<-file(filenameV)
  # writeLines(linesToWriteV, fileConn)
  # close(fileConn)
  # rm(linesToWriteV)


  parallel::stopCluster(clusters)

  # TODO: remove this dummy return (kept for compatibility issues)
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


#concToClass <- function(inss, conclist){
#  ret <- unique(sort(unlist(inss[as.numeric(conclist)])))
#  return(ret)
#}

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
  #TODO define a globa var to store times
  #times <- tic(times, paste("Starting matrix thresholds H->Flat on data ", id, sep = ""))
  predmte <- measureHtoFlat(predarff, combs, colnames(truemte))
  #times <- tic(times, paste("Finished matrix thresholds H->Flat on data ", id, sep = ""))

  write.csv(truemte, paste("true-", id, ".csv", sep = ""), row.names = FALSE)
  write.csv(predmte, paste("pred-", id, ".csv", sep = ""), row.names = FALSE)

}

showResults <- function(result, id, fileid,  dsname, thresholdid){
  for (i in 1:length(result)){
    line <- paste(names(result)[i], ": ", result[i])
    logger(paste("RESULT-", fileid, thresholdid, "-", id,":", sep=""), line)
  }
  write.csv(result, paste("results", "-", id, "-", dsname, "-", fileid, "-", thresholdid, ".csv", sep = ""))
}

compute_results_t1 <- function(id, fileid, dsname){

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

  write.csv(as.bipartition(as.mlresult(pred)), paste("pred-", id, "-t1.csv",sep = ""), row.names = FALSE)

  result <- NULL
  result <- multilabel_evaluate(truemldr, pred)

  showResults(result, id, fileid, dsname, "t1")

}


compute_results_t2 <- function(id, fileid, dsname){
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

  showResults(result, id, fileid,  dsname, "t2")

  # for (i in 1:length(result)){
  #   line <- paste(names(result)[i], ": ", result[i])
  #   logger(paste("RESULT-", fileid, "t2-", id,":", sep=""), line)
  # }
  # write.csv(result, paste("results", id, dsname, fileid, "t2.csv", sep = "-"))
  #
  #

}



compute_results_t3 <- function(id, fileid, dsname){

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

  showResults(result, id, fileid,  dsname, "t3")

  # for (i in 1:length(result)){
  #   line <- paste(names(result)[i], ": ", result[i])
  #   logger(paste("RESULT-", fileid, "t3-", id,":", sep=""), line)
  # }
  # write.csv(result, paste("results", id, dsname, fileid, "t3.csv", sep = "-"))
}


findF2HLibPath <- function(){
  for (path in .libPaths()){
    pathF2H <- file.path(paste(path, "/F2H", sep=""))
    if (file.exists(pathF2H)){
      return(pathF2H)
    }
  }

  # TODO: message
  return(NULL)
}


findClusJar <- function(){
  pathclus <- file.path(paste(findF2HLibPath(), "/java/MyClus.jar", sep=""))
  if (file.exists(pathclus)){
    return(pathclus)
  }

  # TODO: message
  return(NULL)
}


readArffR<- function(arfffile){

  conn <- file(arfffile,open="r")
  lines = readLines(conn)
  close(conn)

  attindexexs <- which(startsWith(lines, "@ATT"))
  datastart <- which(startsWith(lines, "@DATA")) + 1

  data <- read.table(text=paste(lines[datastart:length(lines)], collapse='\n'), header = FALSE, stringsAsFactors = FALSE, sep=',')


  attnames = c()
  atttypes = c()

  for (i in attindexexs){
    x = strsplit(lines[i], " ")[[1]]
    x = x[x != ""]
    attnames <- c(attnames,x[2])
    atttypes <- c(atttypes,x[3])
  }

  colnames(data) = attnames
  rownames(data) = NULL

  # setting datatypes
  for (i in 1:ncol(data)){
    if (atttypes[i] == "{1,0}" || atttypes[i] == "{0,1}"|| atttypes[i] == "numeric"){
      data[,i] <- as.numeric(data[,i])
    } else {
      if (atttypes[i] == "string"){
        data[,i] <- as.character(data[,i])
      }
    }
  }
  data
}


majority_voting_bp <- function(f2hout, mdatat, mdatav, m, dsname){

  sumbpt1 <- f2hout[[1]]$predte1
  sumbpt2 <- f2hout[[1]]$predte2
  sumbpt3 <- f2hout[[1]]$predte3
  for(iteration in seq(2:m)){
    sumbpt1 = sumbpt1 + f2hout[[iteration]]$predte1
    sumbpt2 = sumbpt2 + f2hout[[iteration]]$predte2
    sumbpt3 = sumbpt3 + f2hout[[iteration]]$predte3
  }

  sumbpt1[sumbpt1 < ceiling(m/2)] <- 0
  sumbpt1[sumbpt1 >= ceiling(m/2)] <- 1
  sumbpt2[sumbpt2 < ceiling(m/2)] <- 0
  sumbpt2[sumbpt2 >= ceiling(m/2)] <- 1
  sumbpt3[sumbpt3 < ceiling(m/2)] <- 0
  sumbpt3[sumbpt3 >= ceiling(m/2)] <- 1


  write.csv(sumbpt1, paste("pred-", "teb", "-t1.csv",sep = ""), row.names = FALSE)
  result <- multilabel_evaluate(mdatat, sumbpt1)
  showResults(result, "teb", "EF2H",  dsname, "t1")


  write.csv(sumbpt2, paste("pred-", "teb", "-t2.csv",sep = ""), row.names = FALSE)
  result <- multilabel_evaluate(mdatat, sumbpt2)
  showResults(result, "teb", "EF2H",  dsname, "t2")

  write.csv(sumbpt3, paste("pred-", "teb", "-t3.csv",sep = ""), row.names = FALSE)
  result <- multilabel_evaluate(mdatat, sumbpt3)
  showResults(result, "teb", "EF2H",  dsname, "t3")

}




F2H <- function(
  dsname = "birds",
  train_file = file.path(paste(findF2HLibPath(), "/data/birds_train_1", sep="")),
  test_file = file.path(paste(findF2HLibPath(), "/data/birds_test_1", sep="")),
  valid_file = file.path(paste(findF2HLibPath(), "/data/birds_valid_1", sep="")),
  #dsname = "yeast",
  #train_file = file.path(paste(findF2HLibPath(), "/data/yeast_train_1", sep="")),
  #test_file = file.path(paste(findF2HLibPath(), "/data/yeast_test_1", sep="")),
  #valid_file = file.path(paste(findF2HLibPath(), "/data/yeast_valid_1", sep="")),
  dsdir = tempdir(),
  javaExe = "java",
  javaMem = "-Xmx3g",
  clusJar = findClusJar(),
  minSupportConcetps = 0,
  clusWType = "ExpMaxParentWeight",
  clusWParam = 0.8,
  clusOptimizeErrorMeasure = "WeightedAverageAUPRC",
  threads = 1,
  ensembleClus = 0,
  retPredsConfs = TRUE){

  # define some input vars
  clusExe <- paste(javaExe, " ", javaMem, " -jar \"", clusJar, "\"", sep = "")

  print(paste("Java Option: ", getOption("java.parameters"), sep = ""))

  # reading input files
  times <- c()
  times <- tic(times, "Start")

  setwd(dsdir)

  # reading training file
  print(train_file)
  dfmldr <- mldr(train_file, force_read_from_file = T)
  dfx <- dfmldr$dataset


  # salva o nome original das classes
  original_classnames <- names(dfmldr$attributes[dfmldr$labels$index])

  # obtem os tipos dos atributos
  att_types <- dfmldr$attributes

  # ordena os atributos colocando as classes no fim
  attorder <- c(dfmldr$attributesIndexes, dfmldr$labels$index)
  dfx <- dfx[,attorder]
  firstLabel <- length(dfmldr$attributesIndexes) + 1
  rm(dfmldr)
  lastLabel <- ncol(dfx)
  lastAtt <- firstLabel -1

  # ordena a lista por os atributos foram reordenados
  att_types <- att_types[attorder]

  # renomeia atributos e classes para evitar problemas de caracteres especiais
  colnames(dfx) <- c(paste("att", seq(1, lastAtt), sep = ""), paste("class", seq(firstLabel:lastLabel), sep=""))


  # obtendo o arquivo de teste
  dfmldr <- mldr(test_file, force_read_from_file = T)
  dfxt <- dfmldr$dataset

  # ordena os atributos colocando as classes no fim
  attorder <- c(dfmldr$attributesIndexes, dfmldr$labels$index)
  dfxt <- dfxt[,attorder]


  # renomeia atributos e classes para evitar problemas de caracteres especiais
  colnames(dfxt) <- c(paste("att", seq(1, lastAtt), sep = ""), paste("class", seq(firstLabel:lastLabel), sep=""))


  # obtendo o arquivo de validacao
  dfmldr <- mldr(valid_file, force_read_from_file = T)
  dfxv <- dfmldr$dataset

  # ordena os atributos colocando as classes no fim
  attorder <- c(dfmldr$attributesIndexes, dfmldr$labels$index)
  dfxv <- dfxv[,attorder]


  # renomeia atributos e classes para evitar problemas de caracteres especiais
  colnames(dfxv) <- c(paste("att", seq(1, lastAtt), sep = ""), paste("class", seq(firstLabel:lastLabel), sep=""))

  times <- tic(times, "Finished file loading")


  # inicia com o espaco de classes ####
  df <- dfx[,firstLabel:lastLabel]
  rownames(df) <- NULL
  colnames(df) <- seq(1,ncol(df),1)

  if (ensembleClus == 0){
    fileid = "F2H"
  } else {
    fileid = "EF2H"
  }
  dirf <- file.path(paste(dsdir,"/", dsname, fileid, sep = ""))

  unlink(dirf, recursive = TRUE)
  dir.create(dirf)
  setwd(dirf)

  logger("Output file defined", dirf)
  times <- tic(times, "Finished temp directory creation")

  combs = Rpcbo::computeExtents(df, threads = threads, minsupport = minSupportConcetps)
  times <- tic(times, "Finished PCBO")

  inss <- Rpcbo::computeIntents(df, combs, threads = threads)

  logger(paste("Found", length(combs), "formal concepts ...", sep=" "))
  times <- tic(times, "Finished extent computing")


  # calcula os edges entre os conceitos ####
  edges <- coveringEdges(inss, combs, df, threads)
  logger(paste("Found", length(edges), "edges ...", sep=" "))
  times <- tic(times, "Finished covering edges")


  # adiciona o no raiz e seus edges quando necessario
  if (length(which(unlist(lapply(inss, function(x){identical(x, as.numeric(rownames(df)))})))) == 0){
    # quais conceitos nunca aparecem no consequente do edge? estes são os filhos do novo raiz
    filhosr <- setdiff(seq(1,length(combs),1), unique(sort(unlist(lapply(edges, function(x){x[2]})))))
    novo <- length(combs) + 1
    combs[[novo]] <- as.numeric(NULL)
    inss[[novo]] <- as.numeric(rownames(df))
    for (i in 1:length(filhosr)){
      edges[[length(edges) + 1]] <- c(novo, filhosr[i])
    }
  }
  times <- tic(times, "Finished adding root")

  # encontra o R que contem todos os atributos ou o maior numero de atributos (atributo aqui = instancia)
  r <- which.max(lengths(inss))

  # descobrir os filhos de r
  fr <- which(lengths(lapply(edges, function(x){which(x[1] == r)})) == 1)
  edge_order <- c(fr,seq(1,length(edges))[-fr])


  clusters <- parallel::makeCluster(threads)
  doParallel::registerDoParallel(clusters)

  # gera a string com a hierarquia de classes para posteriormente gerar os arquivos arffh
  strpf <- foreach (i = 1:length(edge_order)) %dopar%{
    ed <- edge_order[i]
    strt <- paste(edges[[ed]][1], "/", edges[[ed]][2], sep="")
    strt
  }
  str <- paste(unlist(strpf), collapse = ",")
  str <- paste("@ATTRIBUTE class                                   hierarchical ", str, sep="")

  parallel::stopCluster(clusters)

  times <- tic(times, "Finished montagem da string da hierarquia")



  # append validation set to test set - clus doesn't make predictions to validation set
  lastTest <- nrow(dfxt)
  dfxt <- rbind(dfxt, dfxv)


  dsdata <- arffH_par("", paste(dsname, ".arff", sep = ""), dfx, dfxt, dfxv, firstLabel, lastLabel, str, inss, combs, r, att_types, threads)

  times <- tic(times, "Finished arrfH writting")


  # check
  # real classes on train only
  # lista_conc <- dsdata$list_concTr
  # rc <- dfx[,firstLabel:lastLabel]
  # #rc <- rbind(rc, dfxt[,firstLabel:lastLabel])
  #
  # transc <- matrix(0, nrow = nrow(rc), ncol = (lastLabel - firstLabel + 1))
  # for (i in 1:nrow(transc)){
  #   transci <- unique(sort(unlist(combs[as.numeric(lista_conc[[i]])])))
  #   transc[i, transci] <- 1
  #   #  print(nrow(truem))
  # }
  # colnames(transc) <- colnames(rc)
  #
  # if (length(which((rc - transc) != 0))){
  #   cat ("ERRO: diferença entre as classes dos conceitos e as classes verdadeiras")
  # }
  # times <- tic(times, "Finished label equivalence checking")
  #

  # gerar conf Clus
  clusSfile <- paste(dsname, ".s", sep="")
  clusTrfile <- paste(dsname, "_train.arff", sep="")
  clusTefile <- paste(dsname, "_test.arff", sep="")
  clusVafile <- paste(dsname, "_valid.arff", sep="")
  geraConfClusHMC(clusSfile, clusTrfile, clusTefile, clusVafile, clusWType, clusWParam, clusOptimizeErrorMeasure)

  print(clusSfile)
  print(clusTrfile)
  print(clusTefile)

  # executar Clus ####
  if (ensembleClus == 0){
    cmd <- paste(clusExe, " ", dsname,  sep= "")
  } else {
    cmd <- paste(clusExe, " -forest ", dsname,  sep= "")
  }
  print(cmd)
  clusout <- system(cmd, intern = TRUE)

  cat(clusout)

  times <- tic(times, "Finished ClusHMC")

  # obtem os resultados ####
  #predarfftr <- read.arff(paste(dsname, ".train.1.pred.arff", sep = ""))
  #predarffte <- read.arff(paste(dsname, ".test.pred.arff", sep = ""))

  logger("Starting reading arff results from Clus")
  predarfftr <- readArffR(paste(dsname, ".train.1.pred.arff", sep = ""))
  predarffte <- readArffR(paste(dsname, ".test.pred.arff", sep = ""))

  logger("Finished reading arff results from Clus")


  # split test and validations predictions
  predarffva <- predarffte[(lastTest+1):nrow(predarffte),]
  predarffte <- predarffte[1:lastTest,]
  dfxt <- dfxt[1:lastTest,]

  logger("Starting converting predictions H->Flat ...")
  compute_predictions_prob(dfx[, firstLabel:lastLabel], predarfftr, combs, "tr")
  compute_predictions_prob(dfxv[, firstLabel:lastLabel], predarffva, combs, "va")
  compute_predictions_prob(dfxt[, firstLabel:lastLabel], predarffte, combs, "te")
  logger("Finished converting predictions H->Flat.")

  logger("Starting computing results ...")
  compute_results_t1("tr", fileid, dsname)
  compute_results_t1("va", fileid, dsname)
  compute_results_t1("te", fileid, dsname)

  compute_results_t2("tr", fileid, dsname)
  compute_results_t2("va", fileid, dsname)
  compute_results_t2("te", fileid, dsname)

  compute_results_t3("tr", fileid, dsname)
  compute_results_t3("va", fileid, dsname)
  compute_results_t3("te", fileid, dsname)
  logger("Finished computing results.")

  list.save(inss, "list_inss.rds")
  list.save(combs, "list_combs.rds")
  list.save(original_classnames, "list_original_classnames.rds")

  infos <- c()
  infos[1] <- dsname
  infos[2] <- dsdir
  infos[3] <- minSupportConcetps
  infos[4] <- clusExe
  infos[5] <- threads
  infos[6] <- length(combs)
  infos[7] <- length(edges)
  infos[8] <- clusWType
  infos[9] <- clusWParam
  infos[10] <- clusOptimizeErrorMeasure
  infos[11] <- ensembleClus


  names(infos) <- c("Dataset", "Outdir", "Min. Sup. PCBO", "Clus", "Threads", "Number of concepts", "Number of edges","clusWType","clusWParam","clusOptimizeErrorMeasure", "ensemble")
  write.csv(x = infos, file = paste(dsname, "-metadata.csv", sep = ""))

  print(infos)

  times <- tic(times, "Finished writting results")

  timest <- tac(times)
  apply(timest, 1, logger, "TIMES")
  write.csv(timest, "times.csv")

  ret = NULL
  if (retPredsConfs){
    ret = list();
    ret$ClusConfte <- read.csv("pred-te.csv");
    ret$ClusConfva <- read.csv("pred-va.csv");
    ret$ClusConftr <- read.csv("pred-tr.csv");
    ret$truete <- read.csv("true-te.csv");
    ret$trueva <- read.csv("true-va.csv");
    ret$truetr <- read.csv("true-tr.csv");
    ret$predte1 <- read.csv("pred-te-t1.csv")
    ret$predte2 <- read.csv("pred-te-t2.csv")
    ret$predte3 <- read.csv("pred-te-t3.csv")
    ret$predtr1 <- read.csv("pred-tr-t1.csv")
    ret$predtr2 <- read.csv("pred-tr-t2.csv")
    ret$predtr3 <- read.csv("pred-tr-t3.csv")
    ret$predva1 <- read.csv("pred-va-t1.csv")
    ret$predva2 <- read.csv("pred-va-t2.csv")
    ret$predva3 <- read.csv("pred-va-t3.csv")
    ret$resultstet1 <- read.csv(paste("results", "te", dsname, "F2H", "t1.csv", sep = "-"));
    ret$resultstet2 <- read.csv(paste("results", "te", dsname, "F2H", "t2.csv", sep = "-"));
    ret$resultstet3 <- read.csv(paste("results", "te", dsname, "F2H", "t3.csv", sep = "-"));

  }
  return(ret)
}

EF2H <- function(
  dsname = "birds",
  train_file = file.path(paste(findF2HLibPath(), "/data/birds_train_1", sep="")),
  test_file = file.path(paste(findF2HLibPath(), "/data/birds_test_1", sep="")),
  valid_file = file.path(paste(findF2HLibPath(), "/data/birds_valid_1", sep="")),
  #dsname = "yeast",
  #train_file = file.path(paste(findF2HLibPath(), "/data/yeast_train_1", sep="")),
  #test_file = file.path(paste(findF2HLibPath(), "/data/yeast_test_1", sep="")),
  #valid_file = file.path(paste(findF2HLibPath(), "/data/yeast_valid_1", sep="")),
  dsdire = tempdir(),
  javaExe = "java",
  javaMem = "-Xmx3g",
  clusJar = findClusJar(),
  minSupportConcetps = 0,
  clusWType = "ExpMaxParentWeight",
  clusWParam = 0.8,
  clusOptimizeErrorMeasure = "WeightedAverageAUPRC",
  threads = 1,
  threadsf2h = 1,
  ensembleClus = 0,
  m = 10, subsample = 1, attr.space = 1, replacement = TRUE, seed = NA, retPredsConfs = TRUE){

  # reading input files
  times <- c()
  times <- tic(times, "Start F2H Ensemble")
  logger("Start F2H Ensemble ...")

  logger("Output directory", dsdire)
  setwd(dsdire)

  # reading training file
  logger("Using as train file: ", train_file)
  mdata <- mldr(train_file, force_read_from_file = T)


  nrow <- ceiling(mdata$measures$num.instances * subsample)
  ncol <- ceiling(length(mdata$attributesIndexes) * attr.space)

  if (!anyNA(seed)) {
    set.seed(seed)
  }

  idx <- lapply(seq(m), function(iteration) {
    list(
      rows = sample(mdata$measures$num.instances, nrow, replacement),
      cols = sample(mdata$attributesIndexes, ncol)
    )
  })

  clusters <- parallel::makeCluster(threads)
  doParallel::registerDoParallel(clusters)

  f2hout <- foreach (iteration = 1:m) %dopar%{
  #for(iteration in seq(1:m)){
    print(iteration)
    setwd(dsdire)
    ndata <- create_subset(mdata, idx[[iteration]]$rows, idx[[iteration]]$cols)

    train_file_e <- paste(dsname, "_ens_", iteration,"_train", sep="")
    write_arff(ndata, file=train_file_e, write.xml = T)

    sink(paste(train_file_e, ".out", sep=""), )

    retf2h <- F2H(
      dsname = train_file_e,
      train_file = file.path(paste(dsdire, "/", train_file_e, sep="")),
      test_file = test_file,
      valid_file = valid_file,
      dsdir = dsdire,
      javaExe = javaExe,
      javaMem = javaMem,
      minSupportConcetps = minSupportConcetps,
      threads = threadsf2h,
      ensembleClus = 0
    )
    sink()
    retf2h
  }
  parallel::stopCluster(clusters)

  # reading the test and validation files
  mdatat <- mldr(test_file, force_read_from_file = T)
  mdatav <- mldr(valid_file, force_read_from_file = T)

  logger("Starting majority voting schema by using bipartition results ...")
  majority_voting_bp(f2hout, mdatat, mdatav, m, dsname)
  logger("Finished majority voting schema by using bipartition results.")

  logger("Starting majority voting schema by using probabilities results ...")

  sumprob <- f2hout[[1]]$ClusConfte
  sumprobva <- f2hout[[1]]$ClusConfva
  for(iteration in seq(2:m)){
    sumprob = sumprob + f2hout[[iteration]]$ClusConfte
    sumprobva = sumprobva + f2hout[[iteration]]$ClusConfva
  }

  sumprob = sumprob/m
  sumprobva = sumprobva/m

  # TODO: adapt the compute_results_xx functions from F2H to receive data by params avoiding the use of disk files
  id = "tep"
  write.csv(mdatat$dataset[,mdatat$labels$index], paste("true-", id, ".csv", sep = ""), row.names = FALSE)
  write.csv(sumprob, paste("pred-", id, ".csv", sep = ""), row.names = FALSE)
  write.csv(sumprobva, paste("pred-va.csv", sep = ""), row.names = FALSE)
  write.csv(mdata$dataset[,mdata$labels$index], paste("true-", "tr", ".csv", sep=""), row.names = FALSE)
  write.csv(mdatav$dataset[,mdatav$labels$index], paste("true-", "va", ".csv", sep=""), row.names = FALSE)

  compute_results_t1(id, "EF2H", dsname)
  compute_results_t2(id, "EF2H", dsname)
  compute_results_t3(id, "EF2H", dsname)

  logger("Finished majority voting schema by using probabilities results ...")

  times <- tic(times, "Finished EF2H")
  timest <- tac(times)
  apply(timest, 1, logger, "TIMES")

  infos <- c(dsname,  train_file,  test_file, valid_file, dsdire, javaExe,  javaMem,  clusJar, minSupportConcetps, clusWType, clusWParam, clusOptimizeErrorMeasure, threads, threadsf2h, ensembleClus, m, subsample, attr.space, replacement, seed, retPredsConfs)
  names(infos) <- c("dsname", "train_file", "test_file", "valid_file", "dsdire", "javaExe", "javaMem", "clusJar", "minSupportConcetps", "clusWType", "clusWParam", "clusOptimizeErrorMeasure", "threads", "threadsf2h", "ensembleClus", "m", "subsample", "attr.space", "replacement", "seed", "retPredsConfs")
  write.csv(x = infos, file = paste(dsname, "-metadata.csv", sep = ""))

  ret = NULL
  if (retPredsConfs){
    ret = list();
    ret$truetr <- read.csv("true-tr.csv");
    ret$trueva <- read.csv("true-va.csv");
    ret$truete <- read.csv("true-tep.csv");
    ret$predtebt1 <- read.csv("pred-teb-t1.csv");
    ret$resultstebt1 <- read.csv(paste("results", "teb", dsname, "EF2H", "t1.csv", sep = "-"));
    ret$predtebt2 <- read.csv("pred-teb-t2.csv");
    ret$resultstebt2 <- read.csv(paste("results", "teb", dsname, "EF2H", "t2.csv", sep = "-"));
    ret$predtebt3 <- read.csv("pred-teb-t3.csv");
    ret$resultstebt3 <- read.csv(paste("results", "teb", dsname, "EF2H", "t3.csv", sep = "-"));

    ret$predConfstep <- sumprob
    ret$predConfsvap <- sumprobva
    ret$predtept1 <- read.csv("pred-tep-t1.csv");
    ret$resultstept1 <- read.csv(paste("results", "tep", dsname, "EF2H", "t1.csv", sep = "-"));

    ret$predtept2 <- read.csv("pred-tep-t2.csv");
    ret$resultstept2 <- read.csv(paste("results", "tep", dsname, "EF2H", "t2.csv", sep = "-"));
    ret$thresholdstept2 <- read.csv("pred-tep-t2-threshold.csv")

    ret$predtept3 <- read.csv("pred-tep-t3.csv");
    ret$resultstept3 <- read.csv(paste("results", "tep", dsname, "EF2H", "t3.csv", sep = "-"));
    ret$thresholdstept3 <- read.csv("pred-tep-t3-threshold.csv")

    ret$F2HData <- f2hout
    ret$infos <- infos

  }
  return(ret)

}






# arffile <- "/home/mauri/temp/tmc2007_500_ens_1_train.train.1.pred.arff"
#
# start_time <- Sys.time()
# x <- read.arff(arffile)
# end_time <- Sys.time()
# end_time - start_time
#
# start_time <- Sys.time()
# y <- F2H::readArffR(arffile)
# end_time <- Sys.time()
# end_time - start_time



