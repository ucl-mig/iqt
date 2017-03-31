# Image Quality Transfer

### Table of Contents
0. [Introduction](#introduction)
0. [Data](#data)
0. [How to use IQT framework](#models)
0. [Citation](#citation)
0. [Disclaimer and known issues](#disclaimer-and-known-issues)
0. [Third-party re-implementations](#third-party-re-implementations)


### Introduction
This repository contains Matlab codes for the IQT framework proposed in our recent NeuroImage paper "Image quality transfer and applications in diffusion MRI" (http://www.sciencedirect.com/science/article/pii/S1053811917302008). The current version supports super-resolution of diffusion tensor images (DTIs). 

The paper introduces a new computational imaging technique called image quality transfer (IQT). IQT uses machine learning to transfer the rich information available from one-off experimental medical imaging devices to the abundant but lower-quality data from routine acquisitions. The procedure uses matched pairs to learn mappings from low-quality to corresponding high-quality images. Once learned, these mappings then augment unseen low quality images, for example by enhancing image resolution or information content. In the papaer, we demonstrate IQT using a simple patch-regression implementation and the uniquely rich diffusion MRI data set from the human connectome project (HCP). Results highlight potential benefits of IQT in both brain connectivity mapping and microstructure imaging. In brain connectivity mapping, IQT reveals, from standard data sets, thin connection pathways that tractography normally requires specialised data to reconstruct. In microstructure imaging, IQT shows potential in estimating, from standard “single-shell” data (one non-zero b-value), maps of microstructural parameters that normally require specialised multi-shell data. Further experiments show strong generalisability, highlighting IQT's benefits even when the training set does not directly represent the application domain. The concept extends naturally to many other imaging modalities and reconstruction problems. For details, please refer to the original paper!

### Data
This demostration requires  HCP in-vivo human diffusion MRI dataset, which is freely available at http://www.humanconnectome.org/documentation/MGH-diffusion/index.html. To test the IQT framework on DTI-superresolution, you need to dowload the diffusion data with b-value = 1000. 

### How to use IQT framework
The pipeline works in three stages: 1. preprocessing (`preprocess.m`); 2. training (`train.m`); 3. testing (`test.m`).

<table>
<tr><td>preprocess.m  </td><td> preprocess the raw data downloaded from HCP, and generate training sets.
</td></tr> <tr><td>train.m </td><td> train regression trees on specified data sets.
</td></tr><tr><td>test.m </td><td>  perform super-resolution on a given diffusion tensor image
</td></tr></table>

You need to point to paths on your system appropriately. The default experiment settings perform x2 super-resolution with input patch 5x5x5.

`test.m` automatically visualises the mean diffusivity, fractional anisotropy and colour encoded directional map for the low-res input, the predicted high-res and the ground truth, and saves as a FIG file. For example, 
![DTI_SR_3x3x3_3x3x3](https://cloud.githubusercontent.com/assets/14926992/24544089/e2e18f72-15f9-11e7-8f7c-0488a8b197aa.png)

### Disclaimer
0. The current version is limited to super-resolution of diffusion tensor images, although readily extendable to the MAP-MRI super-resolution or parameter mapping applications.
0. Edge reconstruction can take a long time. It's best to try the pipeline first without this feature (makes sure the parameter settings.edge = 0)

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

### Third-party re-implementations
There is a third-party Python implementation of IQT. Is it worth pointing people to this?
