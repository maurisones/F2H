---
title: "F2H - under development"
authors: "Mauri Ferrandin and Ricardo Cerri"
date: "09/04/2023"
bibliography: references.bib
---

# Install from github

```{r setup, include=FALSE}
devtools::install_github("maurisones/F2H")
```

## Additional packages

F2H requires the following additional packages:

-   R (\>= 3.6.3);

-   foreach (\>= 1.5.0);

-   doParallel (\>= 1.0.15);

-   parallel;

-   parallelDist;

-   mldr;

-   rlist;

-   hash;

-   digest;

-   BBmisc;

-   utiml;

    ```{r setup, include=FALSE}
    devtools::install_github("rivolli/utiml")
    ```

-   [Rpcbo](https://github.com/maurisones/Rpcbo)

    ```{r setup, include=FALSE}
    devtools::install_github("maurisones/Rpcbo")
    ```

# Using F2H

## Simple test

Run the pre-defined test using birds dataset (installed with the F2H package)

```{r setup, include=FALSE}
library("F2H")
result <- F2H()
```

The output process will provide a lot of information about the classification process. The *result* variable will give access to some data from the classification process, the 4 most important are:

```{r setup, include=FALSE}

# get classification results from the model on test instances using p-cut threshold, as used in the paper
result$resultstet2a

# get classification results of model on test instances without thresholding
result$resultstet0

# get the predictions of model on test instances using p-cut threshold, as used in the paper
result$predte2a

# get the predictions of model on test instances without thresholding
result$predte2a

```

## Parameters

The following parameters are allowed:

-   **dsname**: A string to identify de model. This string is used in temporary files and data structures:

    -   *default* value dsname = "birds"

-   **train_file**: The full path to the arff file containing the training instances (the multi-label arff file must have a additional XML file with the metadata):

    -   *default value* train_file = file.path(paste(findF2HLibPath(), "/data/birds_train_1", sep=""))

-   **test_file**: The full path to the arff file containing the test instances (the multi-label arff file must have a additional XML file with the metadata):

    -   *default value* test_file = file.path(paste(findF2HLibPath(), "/data/birds_test_1", sep=""))

-   **threads**: The number of threads to use

    -   *default value* threads = 1

-   **HierApproach**: The hierarchical classification approach used in the classification process. Allowed values {"global", "local"}

    -   *default value* HierApproach = "global"

-   **dsdir**: a folder to store files generated in the classification process

    -   *default value* dsdir = tempdir()

-   **javaExe**: the full path to the JVM. F2H uses the [Clus Framework](https://dtai.cs.kuleuven.be/clus/) to create the hierarchical classification model and it is provided as a Java implementation

    -   *default value* javaExe = "java"

-   **javaMem**: how much memory will be available to the JVM in the classification process

    -   *default value* javaMem = "-Xmx3g"

-   **clusJar:** full path to the Clus jar file. The F2H package already provides the clus.jar file obtained from the site of [Clus Framework](https://dtai.cs.kuleuven.be/clus/), use this parameter if you want to provide an alternate or customized clus.jar file

    -   *default value* clusJar = findClusJar()

-   **minSupportConcetps**: the minimum number of instances each FCL must contains represented as an integer value. With minSupportConcetps = 0 all FCLs are used. For more information, check the parameter *Smin-support* at [PCBO website](https://fcalgs.sourceforge.net/pcbo-amai.html)

    -   *default value* minSupportConcetps = 0

-   **clusWType**: an internal Clus weighting parameter. For more about check the [Clus documentation](https://dtai.cs.kuleuven.be/clus/)

    -   *default value* clusWType = "ExpMaxParentWeight"

-   **clusWParam**: an internal Clus weighting parameter. For more about check the [Clus documentation](https://dtai.cs.kuleuven.be/clus/)

    -   *default value* clusWParam = 0.8

-   **clusOptimizeErrorMeasure**: an internal Clus error measure optimization parameter. For more about check the [Clus documentation](https://dtai.cs.kuleuven.be/clus/)

    -   *default value* clusOptimizeErrorMeasure = "WeightedAverageAUPRC",

-   **ensembleClus**: use ensemble version of F2H (**Still under development)**, ensembleClus = 0 is the only option for now

    -   *default value* ensembleClus = 0

-   **retPredsConfs**: when TRUE the F2H function will return the all the information as result of its main function (as showed with the result variable in the "Simple test" section

    -   *default value* retPredsConfs = TRUE

## Multi-label datasets repository

Datasets must be provide in multi-label ARFF format along with their XML description files (See the Cometa Repository: <https://cometa.ujaen.es/>).

## How to use F2H with Global Approach

#### Running F2H with most common parameters and global approach (*default*)

```{r setup, include=FALSE}
library("F2H")
result <- F2H(dsname="yeast", train_file="/tmp/yeast_train", test_file="/tmp/yeast_test")
```

#### Running F2H with most common parameters, global approach (*default*) and minSupportConcetps

```{r setup, include=FALSE}
library("F2H")
result <- F2H(
  dsname="yeast", train_file="/tmp/yeast_train", test_file="/tmp/yeast_test", minSupportConcetps=10
)
```

## How to use F2H with Local Approach

#### Running F2H with most common parameters and local approach (*default*)

To run Clus HSC (local hierarchical version) developers created a *perl* script named [run_hsc.pl](https://dtai.cs.kuleuven.be/clus/hmcdatasets/run_hsc.pl) as documented in the [HMC Datasets section of Clus website](https://dtai.cs.kuleuven.be/clus/hmcdatasets/). The F2H package already provides the *run_hsc.pl* with the package files.

```{r setup, include=FALSE}
library("F2H")
result <- F2H(
  dsname="yeast", train_file="/tmp/yeast_train", test_file="/tmp/yeast_test", HierApproach = "local" 
)
```


## How to use EF2H 

## Ensemble parameters

The following parameters are allowed in EF2H ensemble versions:

-   **m**: The number of models:

    -   *default* value m = 10

-   **subsample**: The proportion of instances in subsamples, where 1 means the number of instances equal to that of the original dataset:

    -   *default* value subsample = 1

-   **attr.space**: The proportion of attributos in subsamples, where 1 means the number of attributes to that of the original dataset:

    -   *default* value attr.space = 1
    
-   **replacement**: Select with (or without) replacement the instances in each submodel:

    -   *default* value replacement = TRUE    

-   **ensembleClus**: Indicates to Clus use a single PCT or bag of PCTs:

    -   *default* value ensembleClus = 0    


#### Running EF2H1 (see EF2H paper)

```{r setup, include=FALSE}
library("F2H")
resultEF2H1 <- EF2H(dsname = "birds_folder_1", threads = 5,
  train_file = file.path(paste(findF2HLibPath(), "/data/birds_train_1", sep="")),
  test_file = file.path(paste(findF2HLibPath(), "/data/birds_test_1", sep="")),
  threadsf2h = 2, m = 10, subsample = 1, attr.space = 1,
  ensembleClus = 1
)
```

#### Running EF2H2G (see EF2H paper)

```{r setup, include=FALSE}
library("F2H")
resultEF2H2G <- EF2H(dsname = "birds_folder_1", threads = 5,
  train_file = file.path(paste(findF2HLibPath(), "/data/birds_train_1", sep="")),
  test_file = file.path(paste(findF2HLibPath(), "/data/birds_test_1", sep="")),
  threadsf2h = 2, m = 10, subsample = 1, attr.space = 1,
  ensembleClus = 0,
  HierApproach = "global"
)
```

#### Running EF2H2L (see EF2H paper)

```{r setup, include=FALSE}
library("F2H")
resultEF2H2L <- EF2H(dsname = "birds_folder_1", threads = 5,
  train_file = file.path(paste(findF2HLibPath(), "/data/birds_train_1", sep="")),
  test_file = file.path(paste(findF2HLibPath(), "/data/birds_test_1", sep="")),
  threadsf2h = 2, m = 10, subsample = 1, attr.space = 1,
  ensembleClus = 0,
  HierApproach = "local"
)
```

#### Running EF2H3 (see EF2H paper)

```{r setup, include=FALSE}
library("F2H")
resultEF2H3 <- EF2H(dsname = "birds_folder_1", threads = 5,
  train_file = file.path(paste(findF2HLibPath(), "/data/birds_train_1", sep="")),
  test_file = file.path(paste(findF2HLibPath(), "/data/birds_test_1", sep="")),
  threadsf2h = 2, m = 10, subsample = 1, attr.space = 1,
  ensembleClus = 1,
  HierApproach = "global"
)
```

#### 

## Limitations of F2H and EF2H implementation

-   All implementations were tested only on Linux systems;

-   The Clus Framework requires Java 1.8 and was not tested with other versions;

-   The local version of F2H requires the standard perl installed in the system;

-   The implementations doesn't allow to create a model and use the same model twice, that is a consequence of a problem with the Clus Framework that generates some errors in the model saving process. This limitation makes imperative always to call the F2H() function with the training and test sets, doesn't allowing to make predictions using the same model to a second set of testing instances. We are working on this issue;

-   F2H with global approach is faster;

-   F2H with local approach is slower but achieves best results;

-   To run F2H with datasets with \|Labels\| \> 20 we recommend the use of the minSupportConcetps = 10 or greater.

## How to cite?
```
@article{F2H2023,
  author = {Ferrandin, Mauri and Cerri, Ricardo},
  title = {{Multi-label classification via closed frequent labelsets and label taxonomies}},
  booktitle = {Soft Computing},
  doi = {10.1007/s00500-023-08048-5},
  publisher = {Springer Berlin Heidelberg},
  url = {https://doi.org/10.1007/s00500-023-08048-5},
  year = {2023}
}

@article{EF2H2025,
  title = {Ensemble multi-label classification using closed frequent labelsets and label taxonomies},
  journal = {Applied Soft Computing},
  volume = {171},
  pages = {112853},
  year = {2025},
  issn = {1568-4946},
  doi = {https://doi.org/10.1016/j.asoc.2025.112853},
  url = {https://www.sciencedirect.com/science/article/pii/S1568494625001644},
  author = {Mauri Ferrandin and Ricardo Cerri},
  keywords = {Multi-label classification, Ensemble multi-label classification, Problem transformation},
}

