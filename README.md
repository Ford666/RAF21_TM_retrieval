# TM retrieval with RAF 2-1 algorithm

This repository stores the Matlab codes for a simple demo of our manuscript *Nonconvex optimization for optimum retrieval of the transmission matrix of a multimode fiber*, which has been published in [Advanced Photonics Nexus] (https://doi.org/10.1117/1.APN.2.6.066005).   

The aim of TM retrieval is to retrieve a complex-valued transmission matrix (TM) from intensity-only speckle measurements of a scattering medium or multimode fiber (MMF). We propose a modified nonconvex optimization approach (*i.e.*, RAF 2-1) for TM retrieval, which is featured with optimum efficiency and fast execution in a reference-less and robust setting. This may gain special attention for many deep-tissue imaging and focusing applications with the usage of MMF.

## Abstract
Transmission matrix (TM) allows light control through complex media such as multimode fibers (MMFs), gaining great attention in areas like biophotonics over the past decade. Efforts have been taken to retrieve a complex-valued TM directly from intensity measurements with several representative phase retrieval algorithms, which still see limitations of slow or suboptimum recovery, especially under noisy environments. Here, a modified non-convex optimization approach is proposed. Through numerical evaluations, it shows the optimum focusing efficiency is approached with less running time or sampling ratio. The comparative tests under different signal-to-noise levels further indicate its improved robustness. Experimentally, the superior focusing performance of our algorithm is collectively validated by single- and multi-spot focusing, especially with a sampling ratio of 8 it achieves a 93.6\% efficiency of the gold standard holography method. Based on the recovered TM, image transmission through a MMF is realized with high fidelity. Thanks to parallel operation and GPU acceleration, our nonconvex approach retrieves a 8685×1024 TM (sampling ratio=8) with 42.3 s averagely on a regular computer. The proposed method provides optimum efficiency and fast execution for TM retrieval that avoids the need for an external reference beam, which will facilitate applications of deep-tissue optical imaging, manipulation and treatment.

## Principle of RAF 2-1 for TM retrieval
The pseudocode of retrieving one row of TM is described by:  

<img src="https://github.com/Ford666/RAF21_TM_retrieval/blob/main/images/RAF 2-1.PNG " width="800px">


## Simulation
RAF 2-1 is compared with the existing representative methods (*e.g.*, prVBEM [1], GGS 2-1 [2]) for TM retrieval. Note in simulations, algorithms of CPU verison were used for fair comparison, in which all the target rows of TM are retrieved in parallel, including `prVBEM_tm_cpu.m`, `GGS21_tm_cpu.m`, `RAF_tm_cpu.m` and `RAF21_tm_cpu.m` in the directory of "algs". 

- Run the main file `TM_retrieval_algCompare_default.m` that gives a visual comparison of the focusing efficiency achieved by various algorithms at default settings (*i.e.*, sampling rate=5, without noise, running time=20 s).
<img src="https://github.com/Ford666/RAF21_TM_retrieval/blob/main/images/Algorithm_comparison_image_histogram_default.png " width="1000px">

- Performance comparison in terms of the running time, the sampling rate, different numbers of input mode and different levels of signal-to-noise ratio are provided in the manuscript.


## Experiments
In experiment, we perform single-spot & multi-spot focusing through MMF, focus scanning across the fiber region, and image transmission through MMF. Note in experiments, algorithms of GPU verison were used, including `prVBEM_tm.m`, `GGS21_tm.m`, `RAF_tm.m` and `RAF21_tm.m` in the directory of "algs", and the acceleration effects may vary for different algorithms.


If you encounter bugs that are hard to deal with by youself, please feel free to contact me: Shengfu Cheng (ford.scu.edu@gmail.com)

If you use the codes in your project, please consider to cite our work:
> @article{10.1117/1.APN.2.6.066005,
author = {Shengfu Cheng and Xuyu Zhang and Tianting Zhong and Huanhao Li and Haoran Li and Lei Gong and Honglin Liu and Puxiang Lai},
title = {{Nonconvex optimization for optimum retrieval of the transmission matrix of a multimode fiber}},
volume = {2},
journal = {Advanced Photonics Nexus},
number = {6},
publisher = {SPIE},
pages = {066005},
keywords = {transmission matrix, phase retrieval, multimode fiber imaging, wavefront shaping, Matrices, Multimode fibers, Holography, Image restoration, Photonics, Signal to noise ratio, Speckle, Phase modulation, Phase retrieval, Thulium},
year = {2023},
doi = {10.1117/1.APN.2.6.066005},
URL = {https://doi.org/10.1117/1.APN.2.6.066005}
}
