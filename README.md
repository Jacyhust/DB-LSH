# DB-LSH: Locality-Sensitive Hashing with Query-based Dynamic Bucketing
-----------------------------------------------------------------------------------------------------------------
## Introduction
This is a source code for the algorithm described in the paper **[DB-LSH: Locality-Sensitive Hashing with Query-based Dynamic Bucketing (submitted to ICDE 2022 Second round)]**. We call it as **DB-LSH** project.

## Compilation

**DB-LSH** project is written by **C++** and is simple and easy to use. It can be complied by **g++** in **Linux** and **MSVC** in **Windows**.

### Installation
#### Windows
We can use **Visual Studio 2019** (Other version of Visual Studio should also work but remains untested) to build the project with importing all the files in the directory `./dbLSH/src/`.

#### Linux
```bash
cd ./dbLSH
make
```
The excutable file is then in dbLSH directory, called as `dblsh`

## Usage

### Command Usage

-------------------------------------------------------------------
> dblsh datasetName c k L K beta R_min(optinal)
-------------------------------------------------------------------
(the first parameter specifies the procedure be executed and change)

### Parameter explanation

- datasetName  : dataset name
- c            : 0-1, a float number, the approximate ratio
- k            : 1-100, an integer, the number of returned points
- L            : a positive integer, the number of indexes
- K            : a positive integer, the dimensionality of an index
- beta         : a float number, the maximum ratio of the number of returned points to the total number of dataset   
```
when `beta>0`     , algorithm finishes query with `beta` and other input parameters
when `-10<beta<=0`, algorithm finishes query with a series of different `beta` to generate results for recall/ratio-time curves
when `beta<-10`   , algorithm finishes query with `beta=0.1` and varying 'k'
```
- R_min        : a float number, the inital radius
-------------------------------------------------------------------


FOR EXAMPLE, YOU CAN RUN THE FOLLOWING CODE IN COMMAND LINE AFTER BUILD ALL THE TOOLS:

```bash
cd ./dbLSH
dblsh audio 1.5 50 5 10 0.1 0.3  # a single query
dblsh audio 1.5 50 5 10 -5 0.3   # vary beta
dblsh audio 1.5 50 5 10 -20 0.3  # vary k
```

## Dataset

In our project, the format of the input file (such as `audio.data_new`, which is in `float` data type) is the same as that in [LSHBOX](https://github.com/RSIA-LIESMARS-WHU/LSHBOX). It is a binary file but not a text file, because binary file has many advantages. The binary file is organized as the following format:

>{Bytes of the data type} {The size of the vectors} {The dimension of the vectors} {All of the binary vector, arranged in turn}

For your application, you should also transform your dataset into this binary format, then rename it as `[datasetName].data_new` and put it in the directory `./dataset`.

A sample dataset `audio.data_new` has been put in the directory `./dataset`.
Also, you can get it, `audio.data`, from [here](http://www.cs.princeton.edu/cass/audio.tar.gz)(if so, rename it as `audio.data_new`). If the link is invalid, you can also get it from [data](https://github.com/RSIA-LIESMARS-WHU/LSHBOX-sample-data).

For other dataset we use, you can get the raw data from following links: [MNIST](http://yann.lecun.com/exdb/mnist/index.html), [Cifar](http://www.cs.toronto.edu/~kriz/cifar.html), [Trevi](http://phototour.cs.washington.edu/patches/default.htm), [NUS](https://pan.baidu.com/share/init?surl=kVKfXFx)(Extraction code: hpxg), [Deep1M](https://www.cse.cuhk.edu.hk/systems/hash/gqr/dataset/deep1M.tar.gz), [GIST](http://corpus-texmex.irisa.fr/), [TinyImages80M](https://hyper.ai/tracker/download?torrent=6552), [SIFT](http://corpus-texmex.irisa.fr/). Next, you should transform your raw dataset into the mentioned binary format, then rename it is `[datasetName].data_new` and put it in the directory `./dataset`.


## Result
The experimental result is saved in the directory `./dataset/ANN` as the file
`DB-LSH_result.csv`


## Acknowledgement
**DB-LSH** project is developed by referring to LSHBOX (https://github.com/RSIA-LIESMARS-WHU/LSHBOX). Great appreciation to the contributors of LSHBOX.

## Reference
**[DB-LSH: Locality-Sensitive Hashing with Query-based Dynamic Bucketing (submitted to ICDE 2022)]**

If you meet any issue on the code or take interest in our work, please feel free to contact me (xizhao@ust.hk) or YaoTian (ytianbc@cse.ust.hk). Thank you.
