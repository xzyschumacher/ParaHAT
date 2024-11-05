# ParaHAT: A Fast Noisy Long Read Alignment Tool with Multi-Level Parallelism

## Introduction

The birth of SMRT technology has resolved many limitations of second-generation sequencing, but it has also brought about an exponential increase in alignment data, as well as a higher error rate in sequencing. ParaHAT is a noise long-read alignment tool that employs multi-level parallel technology. It fully utilizes vector-level, thread-level, process-level, and heterogeneous parallelism to accelerate noise long-read alignment.

## Summary

ParaHAT integrates vector-level, thread-level, process-level, and heterogeneous parallelism to accelerate long-read sequence alignment. ParaHAT mainly focuses on parallelizing the rHAT algorithm without altering the original computational process. In the results, our parallel code achieves consistent results with the original rHAT code. The initial intention of ParaHAT is to explore the use of various parallel techniques to maximize program performance.

## ParaHAT installation and usage

Current version of ParaHAT needs to be run on Linux operating system.

### Requrements

ParaHAT is tested to work under:

* Ubuntu 18.04
* gcc 7.5.0
* g++ 7.5.0
* openmpi 2.1.6

### Usage

ParaHAT consists of two parts: indexing and alignment.

**Makefile.**

```
$ make clean
$ make
```

**Indexing.**

```c
$ ./ParaHAT-indexer [-k k-merSize] <HashIndexDir> <Reference>
```

**Alignment.**

```c
$ mpirun [-n nodeNumber] ./ParaHAT-aligner [-w windowsHits] [-m candidates] [-k kmerSize] [-a match] [-b mismatch]
[-q gapOpen] [-r gapExtension] [-t threadNumber] <HashIndexDir> <ReadFile> <Reference>
```

### Parameters

The basic parameters remain consistent with the original rHAT algorithm. The numbers within square brackets represent default values.

**./ParaHAT-indexer:**

```
-k: the size of the k-mers used for indexing in reference genome. [13]
```

**./PararHAT-aligner:**

```
-n: the number of nodes used for multi-node running.
-w: the max allowed number of windows hit by a k-mer. [1000]
-m: the max number of candidates for extension. [5]
-k: the size of the k-mer used for generating short token matches in reads. [13]
-a: score of match for the alignments in extension phase. [1]
-b: mismatch penalty for the alignments in extension phase. [5]
-q: gap open penalty for the alignments in extension phase. [2]
-r: gap extension penalty for the alignments in extension phase. [1]
-l: the minimum length of the local matches used for SDP. [11]
-t: the number of threads used for multi-thread running. [1]
```

