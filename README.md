# TM retrieval with RAF 2-1 algorithm

This repository stores codes for our manuscript *Nonconvex optimization approach for optimum retrieval of the transmission matrix of a multimode fiber*, whcih is currently under review and hopefully will be accepted soon.   

The aim of TM retrieval is to retrieve a complex-valued transmission matrix (TM) from intensity-only speckle measurements of a scattering medium or multimode fiber (MMF). We propose a modified nonconvex optimization approach (*i.e.*, RAF 2-1) for TM retrieval, which is featured with optimum efficiency and fast execution in a reference-less and robust setting. This may gain special attention for many deep-tissue imaging and focusing applications with the usage of MMF.



## Simulation
RAF 2-1 is compared with the existing representative methods for TM retrieval, in terms of the running time, the sampling rate, various number of input modes and different levels of signal-to-noise levels.

Run `sim_prTM_main_compar.m` to conduct simulations for performance comparison.  


## Experiments
In experiment, we perform single-spot & multi-spot focusing through MMF, focus scanning across the fiber region, and image transmission through MMF.

### Prequisite
- Your computer should install [ALp-4.3](https://www.vialux.de/en/download.html) for controlling the digital micromirror device (DMD, DLP9500, Texas Instruments Inc, USA), remember to change the directories of ALP-4.3, including alp4395.dll and alp.h, to your local ones.   

- We use the mex programming (`./utils/CudaLee_1920_1080.mexw64`) to generate the binary computer-generated-hologram (CGH), which is uploaded to a 1920x1080 DMD for phase modulation. To use the mexw64 function, you may need to compile it first in MATLAB with C/C++ compiler being installed. The process may be complicated, and a very basic MATLAB command is: mex -setup C++.

- The MATLAB toolbox of Image Acquisition Explorer is also required, as many APIs are called during the process of image acquisition in experiment.

### single-spot & multi-spot focusing
Run `exp_prTM_main_focus_MMF_V1.m `.  


### focus scanning across the fiber region
Run `exp_prTM_main_PBRMap_MMF_V1.m `, remember to change the sampling rate (gamma) for comparison.


### image transmission through MMF.
Run `exp_prTM_main_PBRMap_MMF_V1.m `. This is performed on the condition that the emperial TM is acquired (under a certain sampling rate).

### ATTENTIONS:
- Although the codes are written for calibrating the TM of a MMF, they are also applicable to regular scattering media. For that, set `fiberAlign=0`, which means the alignment of the fiber imaging area is unnecessary.  

- If you encounter bugs that are hard to deal with by youself, please feel free to contact me: Shengfu Cheng (ford.scu.edu@gmail.com)

