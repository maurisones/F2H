
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

  compute_results_t0(id, "EF2H", dsname)
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




testeEF2Hhmc <- function(){
  x <- EF2H(dsname = "yeast", threads = 3,
            train_file = file.path(paste(findF2HLibPath(), "/data/yeast_train_1", sep="")),
            test_file = file.path(paste(findF2HLibPath(), "/data/yeast_test_1", sep="")),
            valid_file = file.path(paste(findF2HLibPath(), "/data/yeast_valid_1", sep="")),
            threadsf2h = 2, m = 10, subsample = 1, attr.space = 1
  )
  x <- EF2H(dsname = "birds", threads = 2,
            train_file = file.path(paste(findF2HLibPath(), "/data/birds_train_1", sep="")),
            test_file = file.path(paste(findF2HLibPath(), "/data/birds_test_1", sep="")),
            valid_file = file.path(paste(findF2HLibPath(), "/data/birds_valid_1", sep="")),
            threadsf2h = 2, m = 10, subsample = 1, attr.space = 1
  )
}
