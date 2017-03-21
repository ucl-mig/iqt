# Image Quality Transfer

### Table of Contents
0. [Introduction](#introduction)
0. [Citation](#citation)
0. [Disclaimer and known issues](#disclaimer-and-known-issues)
0. [How to use IQT framework](#models)
0. [Results](#results)
0. [Third-party re-implementations](#third-party-re-implementations)


### Introduction
This repository contains Matlab codes for the IQT framework proposed in our recent NeuroImage paper "Image quality transfer and applications in diffusion MRI" (http://www.sciencedirect.com/science/article/pii/S1053811917302008). 

The paper introduces a new computational imaging technique called image quality transfer (IQT). IQT uses machine learning to transfer the rich information available from one-off experimental medical imaging devices to the abundant but lower-quality data from routine acquisitions. The procedure uses matched pairs to learn mappings from low-quality to corresponding high-quality images. Once learned, these mappings then augment unseen low quality images, for example by enhancing image resolution or information content. In the papaer, we demonstrate IQT using a simple patch-regression implementation and the uniquely rich diffusion MRI data set from the human connectome project (HCP). Results highlight potential benefits of IQT in both brain connectivity mapping and microstructure imaging. In brain connectivity mapping, IQT reveals, from standard data sets, thin connection pathways that tractography normally requires specialised data to reconstruct. In microstructure imaging, IQT shows potential in estimating, from standard “single-shell” data (one non-zero b-value), maps of microstructural parameters that normally require specialised multi-shell data. Further experiments show strong generalisability, highlighting IQT's benefits even when the training set does not directly represent the application domain. The concept extends naturally to many other imaging modalities and reconstruction problems. For details, please refer to the original paper!

### Disclaimer
0. The current version is limited to super-resolution of diffusion tensor images, although readily extendable to the MAP-MRI super-resolution or parameter mapping applications.

### Citation

If you use this pipeline in your research, please cite:

      @article{alexander2017image,
        title={Image quality transfer and applications in diffusion MRI},
        author={Daniel C. Alexander and Darko Zikic and Aurobrata Ghosh and Ryutaro Tanno and Viktor Wottschel and Jiaying Zhang and Enrico Kaden and Tim B Dyrby and Stamatios N Sotiropoulos and Hui Zhang and Antonio Criminisi},
        journal={NeuroImage},
        volume={152},
        pages={283--298},
        year={2017},
        publisher={Elsevier}
      }
      
### How to use IQT framework


### Results
May be show some key results?

### Third-party re-implementations
There is a third-party Python implementation of IQT. Is it worth pointing people to this?
