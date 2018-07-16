# NGSweep
Pipeline for preprocessing Next-Gen Sequencing data

## Installation
### Linux

NGSweep can be easily installed by creating a [Conda](https://conda.io/docs/) environment.  Run:
```
conda create -n ngsweep -c matnguyen ngsweep python=3.6
```

and then activate it with:
```
source activate ngsweep
```

NGSweep can then be run in the environment by typing:
```
ngsweep -h
```

To test if all the dependencies were correctly installed, you can run:
```
ngsweep --test
```


