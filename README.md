# TM retrieval with RAF 2-1 algorithm

This repository stores the demo codes for our manuscript *Nonconvex optimization approach for optimum retrieval of the transmission matrix of a multimode fiber*, which is currently under review and hopefully will be accepted soon.   

The aim of TM retrieval is to retrieve a complex-valued transmission matrix (TM) from intensity-only speckle measurements of a scattering medium or multimode fiber (MMF). We propose a modified nonconvex optimization approach (*i.e.*, RAF 2-1) for TM retrieval, which is featured with optimum efficiency and fast execution in a reference-less and robust setting. This may gain special attention for many deep-tissue imaging and focusing applications with the usage of MMF.

## RAF 2-1 for TM retrieval
The pseudocode of retrieving one row of TM is described by:
<img src="https://github.com/Ford666/RAF21_TM_retrieval/blob/main/images/RAF 2-1.png " width="600px">



## Simulation
RAF 2-1 is compared with the existing representative methods (*e.g.*, prVBEM [1], GGS 2-1 [2]) for TM retrieval. Note in simulations, algorithms of CPU verison were used for fair comparison, in which all the target rows of TM are retrieved in parallel, including `prVBEM_tm_cpu.m`, `GGS21_tm_cpu.m`, `RAF_tm_cpu.m` and `RAF21_tm_cpu.m` in the directory of "algs". 

- Run the demo code `TM_retrieval_algCompare_default.m` that gives a visual comparison of the focusing efficiency achieved by various algorithms at default settings (*i.e.*, sampling rate=5, without noise, running time=20 s).
<img src="https://github.com/Ford666/RAF21_TM_retrieval/blob/main/images/Algorithm_comparison_image_histogram_default.png " width="1000px">

- Performance comparison in terms of the running time, the sampling rate, different numbers of input mode and different levels of signal-to-noise ratio are provided in the manuscript.


## Experiments
In experiment, we perform single-spot & multi-spot focusing through MMF, focus scanning across the fiber region, and image transmission through MMF. Note in experiments, algorithms of GPU verison were used (the acceleration effect may vary with different algorithms), including `prVBEM_tm.m`, `GGS21_tm.m`, `RAF_tm.m` and `RAF21_tm.m`. 

### Experimental setup for reference-less MMF calibration
<img src="https://github.com/Ford666/RAF21_TM_retrieval/blob/main/images/Fig.3.png" width="600px">


### single-spot & multi-spot focusing
<img src="https://github.com/Ford666/RAF21_TM_retrieval/blob/main/images/Fig.4.png" width="800px">


### focus scanning across the fiber region
<img src="https://github.com/Ford666/RAF21_TM_retrieval/blob/main/images/Fig.5.Offaxis-RAF21-MMF-FocScan.png" width="800px">


### image transmission through MMF.
<img src="https://github.com/Ford666/RAF21_TM_retrieval/blob/main/images/Fig.6.png" width="800px">


If you encounter bugs that are hard to deal with by youself, please feel free to contact me: Shengfu Cheng (ford.scu.edu@gmail.com)

