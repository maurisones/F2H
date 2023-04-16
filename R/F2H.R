.onAttach <- function(libname, pkgname) {
  #options(java.parameters = "-Xmx2g")
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

findClusHSCPerl <- function(){
  pathclus <- file.path(paste(findF2HLibPath(), "/hsc/run_hsc.pl", sep=""))
  if (file.exists(pathclus)){
    return(pathclus)
  }

  # TODO: message
  return(NULL)
}


logger <- function(a, b=NULL, c=NULL, d=NULL, e=NULL, f=NULL){

  pars <- c(a, b, c, d, e, f)

  str = paste(pars, sep = " ")
  cat(format(Sys.time(), "%a %b %d %X %Y"), str, "\n", sep = ": ")
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


# TODO: reference
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

# TODO: reference
coveringEdges <- function(X, Y, df, threads = 1){

  clusters <- parallel::makeCluster(threads)
  doParallel::registerDoParallel(clusters)

  nc <- length(X)
  M <- seq(1,ncol(df),1)

  keys <- unlist(lapply(X, paste, collapse="-"))
  keys[which(keys == "")] <- "0"
  keys <- unlist(lapply(keys, digest, algo = 'md5'))
  h <- hash(keys=keys, values=1:length(Y))

  out <- foreach (i = 1:nc, .export=c("mlinha", "has.key", "digest")) %dopar%{
    counts <- rep(0, nc)

    edgesi <- list()
    m <- setdiff(M, Y[[i]])
    if (length(m) > 0){
      for (j in 1:length(m)){
        inters <- as.numeric(intersect(X[[i]], mlinha(m[[j]],df)))

        if (length(inters) > 0){
          inters <- digest(paste(inters, collapse = "-"), algo = 'md5')
          if (has.key(inters, h)){
            X_inters <- h[inters]
            X_inters <- as.list(X_inters)[[1]]
            counts[X_inters] <- counts[X_inters] + 1
            if ((length(Y[[X_inters]]) - length(Y[[i]])) == counts[X_inters]){
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


arffH_par <- function(dir, fileArffH, arff, arfft, labelFirst, labelLast, classstr, inss, combs, root, att_types, threads = 1){

  clusters <- parallel::makeCluster(threads)
  doParallel::registerDoParallel(clusters)

  linesToWrite = c()

  arffAtts <- arff[,1:(labelFirst - 1)]
  arffLabels <- arff[,labelFirst:labelLast]
  arffAttst <- arfft[,1:(labelFirst - 1)]
  arffLabelst <- arfft[,labelFirst:labelLast]

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

  str <- foreach (i = 1:nrow(arff)) %dopar%{

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


  str <- foreach (i = 1:nrow(arfft)) %dopar%{
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

  parallel::stopCluster(clusters)

}

geraConfClusHMC<- function(sfilename, trainfile, testfile, WType, WParam, OptimizeErrorMeasure){

  s <- "[Data]
File = XXX1
TestSet = XXX2

%using the valid dataset achieved low performance
%we choose do not use validation file for clus in F2H
%PruneSet =

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

  fileConn<-file(sfilename)
  writeLines(s, fileConn)
  close(fileConn)

}


measureHtoFlat <- function(predarff, combs, truem_colnames){

  dfpred <- predarff[,which(startsWith(colnames(predarff), "Original"))]
  dfpred <- dfpred[,-ncol(dfpred)]

  colnames(dfpred) <- as.numeric(gsub(x = colnames(dfpred), pattern = "Original-p-", replacement = ""))

  predm <- matrix(0, nrow = nrow(dfpred), ncol = length(truem_colnames))  # return predictions
  predc <- matrix(-1, nrow = nrow(dfpred), ncol = length(truem_colnames)) # return the level of node used in the prediction
  for (i in 1:ncol(dfpred)){
    metai <- unlist(combs[as.numeric(colnames(dfpred)[i])])
    for (m in metai){
      menores <- which((predm[,m] - dfpred[,i]) < 0)
      predm[menores, m] <- dfpred[menores, i]

      predc[menores, m] <- as.numeric(colnames(dfpred)[i])

    }
  }
  colnames(predm) <- truem_colnames
  colnames(predc) <- truem_colnames

  ret <- list()
  ret$predm <- predm
  ret$predc <- predc

  return(ret)
}


compute_predictions_prob <- function(truemte, predarff, combs, id){

  original <- data.frame(truemte)
  original <- data.frame(sapply(original, function(x) as.numeric(as.character(x))))
  mymldr <- NULL
  mymldr <- mldr_from_dataframe(original, labelIndices = seq(1,ncol(truemte)), name = "original")

  # obtem as probabilidades das predicoes
  mhtoflat <- measureHtoFlat(predarff, combs, colnames(truemte))
  predmte <- mhtoflat$predm
  predmtec <- mhtoflat$predc

  write.csv(truemte, paste("true-", id, ".csv", sep = ""), row.names = FALSE)
  write.csv(predmte, paste("pred-", id, ".csv", sep = ""), row.names = FALSE)
  write.csv(predmtec, paste("predlevel-", id, ".csv", sep = ""), row.names = FALSE)

}

showResults <- function(result, id, fileid,  dsname, thresholdid){
  for (i in 1:length(result)){
    line <- paste(names(result)[i], ": ", result[i])
    logger(paste("RESULT-", fileid, thresholdid, "-", id,":", sep=""), line)
  }
  write.csv(result, paste("results", "-", id, "-", dsname, "-", fileid, "-", thresholdid, ".csv", sep = ""))
}


compute_results_t0 <- function(id, fileid, dsname, threads = 1){
  # t0 = no threshold: useful to obtain ranking measures
  true <- read.csv(paste("true-", id, ".csv", sep=""))
  true <- data.frame(sapply(true, function(x) as.numeric(as.character(x))))
  truemldr <- mldr_from_dataframe(true, labelIndices = seq(1,ncol(true)), name = "true")

  pred <- read.csv(paste("pred-", id, ".csv", sep=""))
  pred <- data.frame(sapply(pred, function(x) as.numeric(as.character(x))))

  write.csv(as.bipartition(as.mlresult(pred)), paste("pred-", id, "-t0.csv",sep = ""), row.names = FALSE)

  result <- NULL
  result <- multilabel_evaluate(truemldr, pred, label = T)

  # a macro-AUC corrections when there is a class with no instances - utiml return null and in these cases
  # those classes are discarded from the mean to solve this issu
  if(is.na(result$multilabel["macro-AUC"])){
    result$multilabel["macro-AUC"] <- mean(result$labels[,"AUC"], na.rm=TRUE)
  }

  showResults(result$multilabel, id, fileid, dsname, "t0")

}

# implement lcard para selecionar threshold global
# https://search.r-project.org/CRAN/refmans/utiml/html/lcard_threshold.html
compute_results_t2a <- function(id, fileid, dsname, threads = 1){
  # t0 = no threshold: useful to obtain ranking measures
  true <- read.csv(paste("true-", id, ".csv", sep=""))
  true <- data.frame(sapply(true, function(x) as.numeric(as.character(x))))
  truemldr <- mldr_from_dataframe(true, labelIndices = seq(1,ncol(true)), name = "true")

  # necessário para obter a cardinalidade
  truetr <- read.csv(paste("true-", "tr", ".csv", sep=""))
  truetr <- mldr_from_dataframe(truetr, labelIndices = seq(1,ncol(truetr)), name = "true")

  pred <- read.csv(paste("pred-", id, ".csv", sep=""))
  pred <- data.frame(sapply(pred, function(x) as.numeric(as.character(x))))

  # aplica a função para selecionar o melhor threshold de acordo com a cardinalidade
  pred <- lcard_threshold(as.matrix(pred), truetr$measures$cardinality, probability = T)
  write.csv(as.bipartition(as.mlresult(pred)), paste("pred-", id, "-t2a.csv",sep = ""), row.names = FALSE)

  result <- NULL
  result <- multilabel_evaluate(truemldr, pred, label = T)

  # a macro-AUC corrections when there is a class with no instances - utiml return null and in these cases
  # those classes are discarded from the mean to solve this issu
  if(is.na(result$multilabel["macro-AUC"])){
    result$multilabel["macro-AUC"] <- mean(result$labels[,"AUC"], na.rm=TRUE)
  }

  showResults(result$multilabel, id, fileid, dsname, "t2a")

}

compute_results_t2 <- function(id, fileid, dsname){
  true <- read.csv(paste("true-", id, ".csv", sep=""))
  true <- data.frame(sapply(true, function(x) as.numeric(as.character(x))))
  truemldr <- mldr_from_dataframe(true, labelIndices = seq(1,ncol(true)), name = "true")

  # precisa do label space de treinamento e validacao para calculo do threshold
  truetr <- read.csv(paste("true-", "tr", ".csv", sep=""))
  #rv trueva <- read.csv(paste("true-", "va", ".csv", sep=""))

  pred <- read.csv(paste("pred-", id, ".csv", sep=""))
  pred <- data.frame(sapply(pred, function(x) as.numeric(as.character(x))))

  # precisa das predicoes da validacao para calibrar o threshold
  #rv predva <- read.csv(paste("pred-", "va", ".csv", sep=""))
  #rv predva <- data.frame(sapply(predva, function(x) as.numeric(as.character(x))))

  predtr <- read.csv(paste("pred-", "tr", ".csv", sep=""))
  predtr <- data.frame(sapply(predtr, function(x) as.numeric(as.character(x))))

  # Thresholding t2 - One Threshold
  # calibrating the threshold
  #rv lcard <- (sum(truetr) + sum(trueva))/(nrow(truetr) + nrow(trueva))
  lcard <- sum(truetr)/nrow(truetr)
  #t <- seq(0.01,0.99, 0.01)
  t <- sort(unique(unlist(c(read.csv(paste("pred-", id, ".csv", sep=""))))))
  bt <- 0
  bts <- 1
  for (i in 1:length(t)){
    s <- abs(lcard - (length(which(predtr >= t[i]))/nrow(predtr)))
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


}


readArffR<- function(arfffile){

  conn <- file(arfffile,open="r")
  lines = readLines(conn)
  close(conn)

  attindexexs <- which(startsWith(lines, "@ATT"))
  datastart <- which(startsWith(lines, "@DATA"))

  data <- read.table(arfffile, header = FALSE, stringsAsFactors = FALSE, sep=',', skip = datastart)

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

showPredictionLevels <- function(){
  predlevels <- read.csv("predlevel-te.csv")
  combs <- list.load("list_combs.rds")

  predconc <- read.csv("pred-conc-te.csv")
  predconc <- predconc[,which(startsWith(colnames(predconc), "Original"))]
  predconc <- predconc[,-ncol(predconc)]

  colnames(predconc) <- as.numeric(gsub(x = colnames(predconc), pattern = "Original.p.", replacement = ""))

  for (c in 1:ncol(predlevels)){
    conc <- unique(predlevels[,c])
    conc <- conc[which(conc > 0)]
    if (length(conc) < 1){
      logger("LEVELPREDINFO", colnames(predlevels)[c], "Has no predictions")
      next
    }
    totconc <- unlist(lapply(conc, function(x){length(which(predlevels[,c] == x))}))
    for (i in 1:length(conc)){
      cc <- conc[i]
      logger("LEVELPREDINFO", colnames(predlevels)[c], cc, totconc[i], paste(combs[[cc]], collapse = "-"))


      if (length(combs[[cc]]) > 1){
        ccc <- which(unlist(lapply(combs, function(x){sum(x %in% combs[[cc]]) == length(x) && length(x) > 0})))
      } else{
        ccc <- cc
      }

      ccci <- unlist(lapply(ccc, function(x){which(colnames(predconc) == as.character(x))}))

      #pegar as linhas em que o conceito foi usado para a classe
      insts <- which(predlevels[,c] == cc)

      for (ins in insts){
        for (co in 1:length(ccci)){
          ci <- ccci[co]
          ca <- ccc[co]
          logger("\tLEVELPREDINFO - inst: ", ins,  paste("concept: ", ca, paste(combs[[ca]], collapse = "-"), sprintf("%.3f",predconc[ins, ci]), sep = "\t"))
        }
      }

    }
  }

}

F2H <- function(
  dsname = "birds",
  train_file = file.path(paste(findF2HLibPath(), "/data/birds_train_1", sep="")),
  test_file = file.path(paste(findF2HLibPath(), "/data/birds_test_1", sep="")),
  # dsname = "yeast",
  # train_file = file.path(paste(findF2HLibPath(), "/data/yeast_train_1", sep="")),
  # test_file = file.path(paste(findF2HLibPath(), "/data/yeast_test_1", sep="")),
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
  retPredsConfs = TRUE,
  HierApproach = "global"){

  # define some input vars
  clusExe <- paste(javaExe, " ", javaMem, " -jar \"", clusJar, "\"", sep = "")

  print(paste("Java Option (for RJava and not for Clus): ", getOption("java.parameters"), sep = ""))

  times <- c()
  times <- tic(times, "Start")

  # reading training file
  setwd(dsdir)
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

  # generate HMC arffs
  arffH_par("", paste(dsname, ".arff", sep = ""), dfx, dfxt, firstLabel, lastLabel, str, inss, combs, r, att_types, threads)
  times <- tic(times, "Finished arrfH writting")

  # gerar conf Clus
  clusSfile <- paste(dsname, ".s", sep="")
  clusTrfile <- paste(dsname, "_train.arff", sep="")
  clusTefile <- paste(dsname, "_test.arff", sep="")
  geraConfClusHMC(clusSfile, clusTrfile, clusTefile, clusWType, clusWParam, clusOptimizeErrorMeasure)

  print(clusSfile)
  print(clusTrfile)
  print(clusTefile)

  # executar Clus ####

  if (HierApproach == "global"){
    if (ensembleClus == 0){
      cmd <- paste(clusExe, " ", dsname,  sep= "")
    } else {
      cmd <- paste(clusExe, " -forest ", dsname,  sep= "")
    }
    print(cmd)
    clusout <- system(cmd, intern = TRUE)

    cat(clusout)

    times <- tic(times, "Finished ClusHMC")

    logger("Starting reading arff results from Clus")
    predarfftr <- readArffR(paste(dsname, ".train.1.pred.arff", sep = ""))
    predarffte <- readArffR(paste(dsname, ".test.pred.arff", sep = ""))
    logger("Finished reading arff results from Clus")
  } else {
    system(paste("cp ", findClusHSCPerl(), ".", sep=" "))

    cmd <- paste("perl run_hsc.pl ", dsname,  " ", strsplit(clusJar, "MyClus")[[1]][1], sep= "")
    print(cmd)
    clusout <- system(cmd, intern = TRUE)

    cat(clusout)

    times <- tic(times, "Finished ClusHMC")

    logger("Starting reading arff results from Clus")
    predarfftr <- readArffR(paste(dsname, ".hsc.combined.train.1.pred.arff", sep = ""))
    predarffte <- readArffR(paste(dsname, ".hsc.combined.test.pred.arff", sep = ""))

    logger("Finished reading arff results from Clus")

  }

  logger("Starting converting predictions H->Flat ...")
  compute_predictions_prob(dfx[, firstLabel:lastLabel], predarfftr, combs, "tr")
  compute_predictions_prob(dfxt[, firstLabel:lastLabel], predarffte, combs, "te")
  logger("Finished converting predictions H->Flat.")

  logger("Starting computing results ...")
  compute_results_t0("tr", fileid, dsname)
  compute_results_t0("te", fileid, dsname)

  compute_results_t2("tr", fileid, dsname)
  compute_results_t2("te", fileid, dsname)

  compute_results_t2a("tr", fileid, dsname)
  compute_results_t2a("te", fileid, dsname)
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
  # TODO: log new parameters

  names(infos) <- c("Dataset", "Outdir", "Min. Sup. PCBO", "Clus", "Threads", "Number of concepts", "Number of edges","clusWType","clusWParam","clusOptimizeErrorMeasure", "ensemble")
  write.csv(x = infos, file = paste(dsname, "-metadata.csv", sep = ""))

  print(infos)
  times <- tic(times, "Finished writting results")

  timest <- tac(times)
  apply(timest, 1, logger, "TIMES")
  write.csv(timest, "times.csv")

  ret = NULL
  if (retPredsConfs){
    ret = list()
    ret$ClusConfte <- read.csv("pred-te.csv")
    ret$ClusConftr <- read.csv("pred-tr.csv")
    ret$truete <- read.csv("true-te.csv")
    ret$truetr <- read.csv("true-tr.csv")
    ret$predte0 <- read.csv("pred-te-t0.csv")
    ret$predte2 <- read.csv("pred-te-t2.csv")
    ret$predte2a <- read.csv("pred-te-t2a.csv")
    ret$predtr0 <- read.csv("pred-tr-t0.csv")
    ret$predtr2 <- read.csv("pred-tr-t2.csv")
    ret$predtr2a <- read.csv("pred-tr-t2a.csv")
    ret$resultstet0 <- read.csv(paste("results", "te", dsname, fileid, "t0.csv", sep = "-"))
    ret$resultstet2 <- read.csv(paste("results", "te", dsname, fileid, "t2.csv", sep = "-"))
    ret$resultstet2a <- read.csv(paste("results", "te", dsname, fileid, "t2a.csv", sep = "-"))
    ret$predlevelte <- read.csv("predlevel-te.csv")

    #showPredictionLevels()

  }
  return(ret)
}


testeF2Hhsc <- function(){
  x <- F2H(dsname = "yeast", threads = 4, HierApproach = "local",
              train_file = file.path(paste(findF2HLibPath(), "/data/yeast_train_1", sep="")),
              test_file = file.path(paste(findF2HLibPath(), "/data/yeast_test_1", sep=""))
  )

  y <- F2H(dsname = "birds", threads = 4, HierApproach = "local",
              train_file = file.path(paste(findF2HLibPath(), "/data/birds_train_1", sep="")),
              test_file = file.path(paste(findF2HLibPath(), "/data/birds_test_1", sep=""))
  )

  y <- F2H(dsname = "birds", threads = 4, HierApproach = "local",
           train_file = file.path(paste(findF2HLibPath(), "/data/birds_train_1", sep="")),
           test_file = file.path(paste(findF2HLibPath(), "/data/birds_test_1", sep=""))
  )

}

testeF2Hhmc <- function(){
    x <- F2H(dsname = "yeast", threads = 4, HierApproach = "global",
           train_file = file.path(paste(findF2HLibPath(), "/data/yeast_train_1", sep="")),
           test_file = file.path(paste(findF2HLibPath(), "/data/yeast_test_1", sep="")),
           ensembleClus = 1
           )

   y <- F2H(dsname = "birds", threads = 4, HierApproach = "global",
           train_file = file.path(paste(findF2HLibPath(), "/data/birds_train_1", sep="")),
           test_file = file.path(paste(findF2HLibPath(), "/data/birds_test_1", sep=""))
        )

   y <- F2H(dsname = "birds", threads = 4, HierApproach = "global",
            train_file = file.path(paste(findF2HLibPath(), "/data/birds_train_1", sep="")),
            test_file = file.path(paste(findF2HLibPath(), "/data/birds_test_1", sep="")),
            ensembleClus = 1
   )


}

