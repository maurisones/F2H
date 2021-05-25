---
title: "F2H - under development"
author: "Mauri Ferrandin"
date: "25/05/2021"
---


# Install from github
```{r setup, include=FALSE}
devtools::install_github("maurisones/F2H")
```

## Additional packages

F2H requires many other packages to run, the most part of them can be directed installed from CRAN mirrors. The bellow mentioned are also/only available from github:

Rpcbo:  
```{r setup, include=FALSE}
devtools::install_github("maurisones/Rpcbo")
```
Utiml:
```{r setup, include=FALSE}
devtools::install_github("rivolli/utiml")
```

# Using F2H

## Simple teste 
Run the pre-defined test on a fold of birds dataset (installed with the package)
```{r setup, include=FALSE}
library("F2H")
F2H()
```
## Setting datasets
Datasets must be provide in multi-label ARFF format along with their XML description files (See the Cometa Repositóry: https://cometa.ujaen.es/). 

```{r setup, include=FALSE}
library("F2H")
F2H(
  train_file=<fullpath to your training file without arff extension>,
  valid_file=<fullpath to your validation file without arff extension>,
  test_file=<fullpath to your test file without arff extension>
)
```

## Other (not less important) parameters to F2H() function

```{r setup, include=FALSE}
dsname = <your string that identifies the dataset> 
```

This identifier is used in the name of output folders and files to identify the experiment. Default value is "birds"

```{r setup, include=FALSE}
dsdir = <path to your output folder>
```
The F2H is dependent of file system temporary data input and output. This parameter is used to especify the output location. Default value is the OS default temporary folder obtained by the function tempdir().

```{r setup, include=FALSE}
javaExe = <location of java exec file including file path and file name>
```
The F2H uses the java command via system call and needs the JRE or JDK installed. If not given by user, de default value is just "java" aiming to exec de default java in the system.


```{r setup, include=FALSE}
  javaMem = <amount of memory to java processes>
```
You can set the maximum used memory by a java process. Default value is -Xmx3g, equivalent to 3GB.


```{r setup, include=FALSE}
  clusJar = <path to your Clus.jar file>
```
F2H needs the Clus PCT implementation. If not specified, a Clus.jar file is provided along with the package. This version of Clus was configured to avoid some output logs aiming to reduce the storage space used during its execution. More about Clus and its original Clus.jar file can be found at:  https://dtai.cs.kuleuven.be/clus/


```{r setup, include=FALSE}
  minSupportConcetps = <your value for minimun CFL support>
```
This parameter specifies the minimum support a CFL (Closed Frequent LabelSet) must have to be included in the generated label hierarchy. Default value is 0.

```{r setup, include=FALSE}
  threads = <number of threads>
```
Specifies the number of threads to use. The F2H uses parallel processing to execute but, the Clus has no support to parallelism.



The above 3 parameters are Clus specific configuratiosn. For more about them see the Clus User Manual available at: https://dtai.cs.kuleuven.be/clus/. Defaul values are presented bellow:

```{r setup, include=FALSE}
  clusWType = "ExpMaxParentWeight"
```

```{r setup, include=FALSE}
  clusWParam = 0.8
```

```{r setup, include=FALSE}
  clusOptimizeErrorMeasure = "WeightedAverageAUPRC"
```


