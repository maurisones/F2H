
majority_voting_bp <- function(f2hout, mdatat, m, dsname){

  sumbpt0 <- f2hout[[1]]$predte0
  sumbpt2a <- f2hout[[1]]$predte2a
  for(iteration in seq(2:m)){
    sumbpt0 = sumbpt0 + f2hout[[iteration]]$predte0
    sumbpt2a = sumbpt2a + f2hout[[iteration]]$predte2a
  }

  sumbpt0[sumbpt0 < ceiling(m/2)] <- 0
  sumbpt0[sumbpt0 >= ceiling(m/2)] <- 1
  sumbpt2a[sumbpt2a < ceiling(m/2)] <- 0
  sumbpt2a[sumbpt2a >= ceiling(m/2)] <- 1


  write.csv(sumbpt0, paste("pred-", "teb", "-t0.csv",sep = ""), row.names = FALSE)
  result <- multilabel_evaluate(mdatat, sumbpt0)
  showResults(result, "teb", "EF2H",  dsname, "t0")

  write.csv(sumbpt2a, paste("pred-", "teb", "-t2a.csv",sep = ""), row.names = FALSE)
  result <- multilabel_evaluate(mdatat, sumbpt2a)
  showResults(result, "teb", "EF2H",  dsname, "t2a")

}


EF2H <- function(
    dsname = "birds",
    train_file = file.path(paste(findF2HLibPath(), "/data/birds_train_1", sep="")),
    test_file = file.path(paste(findF2HLibPath(), "/data/birds_test_1", sep="")),
    #dsname = "yeast",
    #train_file = file.path(paste(findF2HLibPath(), "/data/yeast_train_1", sep="")),
    #test_file = file.path(paste(findF2HLibPath(), "/data/yeast_test_1", sep="")),
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
    m = 10, subsample = 1, attr.space = 1, replacement = TRUE, seed = NA, retPredsConfs = TRUE,
    HierApproach = "global"
    ){

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

  f2hout <- foreach (iteration = 1:m ) %dopar%{

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
      dsdir = dsdire,
      javaExe = javaExe,
      javaMem = javaMem,
      minSupportConcetps = minSupportConcetps,
      threads = threadsf2h,
      ensembleClus = ensembleClus,
      HierApproach = HierApproach
    )
    sink()
    retf2h
  }

  parallel::stopCluster(clusters)

  # reading the test and validation files
  mdatat <- mldr(test_file, force_read_from_file = T)

  logger("Starting majority voting schema by using bipartition results ...")
  majority_voting_bp(f2hout, mdatat, m, dsname)
  logger("Finished majority voting schema by using bipartition results.")

  logger("Starting majority voting schema by using probabilities results ...")
  sumprob <- f2hout[[1]]$ClusConfte

  for(iteration in seq(2:m)){
    sumprob = sumprob + f2hout[[iteration]]$ClusConfte
  }

  sumprob = sumprob/m


  # TODO: adapt the compute_results_xx functions from F2H to receive data by params avoiding the use of disk files
  id = "tep"
  write.csv(mdatat$dataset[,mdatat$labels$index], paste("true-", id, ".csv", sep = ""), row.names = FALSE)
  write.csv(sumprob, paste("pred-", id, ".csv", sep = ""), row.names = FALSE)
  write.csv(mdata$dataset[,mdata$labels$index], paste("true-", "tr", ".csv", sep=""), row.names = FALSE)

  compute_results_t0(id, "EF2H", dsname)
  compute_results_t2a(id, "EF2H", dsname)

  logger("Finished majority voting schema by using probabilities results ...")

  times <- tic(times, "Finished EF2H")
  timest <- tac(times)
  apply(timest, 1, logger, "TIMES")

  infos <- c(dsname,  train_file,  test_file, dsdire, javaExe,  javaMem,  clusJar, minSupportConcetps, clusWType, clusWParam, clusOptimizeErrorMeasure, threads, threadsf2h, ensembleClus, m, subsample, attr.space, replacement, seed, retPredsConfs)
  names(infos) <- c("dsname", "train_file", "test_file", "dsdire", "javaExe", "javaMem", "clusJar", "minSupportConcetps", "clusWType", "clusWParam", "clusOptimizeErrorMeasure", "threads", "threadsf2h", "ensembleClus", "m", "subsample", "attr.space", "replacement", "seed", "retPredsConfs")
  write.csv(x = infos, file = paste(dsname, "-metadata.csv", sep = ""))

  ret = NULL
  if (retPredsConfs){
    ret = list();
    ret$truetr <- read.csv("true-tr.csv");
    ret$truete <- read.csv("true-tep.csv");
    ret$predtebt2a <- read.csv("pred-teb-t2a.csv");
    ret$resultstebt2a <- read.csv(paste("results", "teb", dsname, "EF2H", "t2a.csv", sep = "-"));

    ret$predConfstep <- sumprob
    ret$predtept0 <- read.csv("pred-tep-t0.csv");
    ret$resultstept0 <- read.csv(paste("results", "tep", dsname, "EF2H", "t0.csv", sep = "-"));

    ret$predtept2a <- read.csv("pred-tep-t2a.csv");
    ret$resultstept2a <- read.csv(paste("results", "tep", dsname, "EF2H", "t2a.csv", sep = "-"));

    ret$F2HData <- f2hout
    ret$infos <- infos

  }
  return(ret)

}




testeEF2Hhmc <- function(){

  x <- EF2H(dsname = "yeast", threads = 10,
            train_file = file.path(paste(findF2HLibPath(), "/data/yeast_train_1", sep="")),
            test_file = file.path(paste(findF2HLibPath(), "/data/yeast_test_1", sep="")),
            threadsf2h = 2, m = 10, subsample = 1, attr.space = 1
  )

  x <- EF2H(dsname = "birds", threads = 10,
            train_file = file.path(paste(findF2HLibPath(), "/data/birds_train_1", sep="")),
            test_file = file.path(paste(findF2HLibPath(), "/data/birds_test_1", sep="")),
            threadsf2h = 2, m = 10, subsample = 1, attr.space = 1, ensembleClus = 1
  )


  x <- EF2H(dsname = "birds", threads = 10,
            train_file = file.path(paste(findF2HLibPath(), "/data/birds_train_1", sep="")),
            test_file = file.path(paste(findF2HLibPath(), "/data/birds_test_1", sep="")),
            threadsf2h = 2, m = 10, subsample = 1, attr.space = 1, ensembleClus = 0,
            HierApproach = "local"
            )

  x <- EF2H(dsname = "yeast", threads = 10,
            train_file = file.path(paste(findF2HLibPath(), "/data/yeast_train_1", sep="")),
            test_file = file.path(paste(findF2HLibPath(), "/data/yeast_test_1", sep="")),
            threadsf2h = 2, m = 10, subsample = 1, attr.space = 1, ensembleClus = 0,
            HierApproach = "local"
  )

  # paper sample models

  # EF2H1
  resultEF2H1 <- EF2H(dsname = "yeast", threads = 5,
            train_file = file.path(paste(findF2HLibPath(), "/data/birds_train_1", sep="")),
            test_file = file.path(paste(findF2HLibPath(), "/data/birds_test_1", sep="")),
            threadsf2h = 2, m = 10, subsample = 1, attr.space = 1,
            ensembleClus = 1
  )

  # EF2H2G
  resultEF2H2G <- EF2H(dsname = "yeast", threads = 5,
            train_file = file.path(paste(findF2HLibPath(), "/data/birds_train_1", sep="")),
            test_file = file.path(paste(findF2HLibPath(), "/data/birds_test_1", sep="")),
            threadsf2h = 2, m = 10, subsample = 1, attr.space = 1,
            ensembleClus = 0,
            HierApproach = "global"
  )

  # EF2H2L
  resultEF2H2L <- EF2H(dsname = "yeast", threads = 5,
            train_file = file.path(paste(findF2HLibPath(), "/data/birds_train_1", sep="")),
            test_file = file.path(paste(findF2HLibPath(), "/data/birds_test_1", sep="")),
            threadsf2h = 2, m = 10, subsample = 1, attr.space = 1,
            ensembleClus = 0,
            HierApproach = "local"
  )

  # EF2H3
  resultEF2H3 <- EF2H(dsname = "yeast", threads = 5,
            train_file = file.path(paste(findF2HLibPath(), "/data/birds_train_1", sep="")),
            test_file = file.path(paste(findF2HLibPath(), "/data/birds_test_1", sep="")),
            threadsf2h = 2, m = 10, subsample = 1, attr.space = 1,
            ensembleClus = 1,
            HierApproach = "global"
  )

}
