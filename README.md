# Image Quality Transfer

### Demo
# <a href="http://www.youtube.com/watch?feature=player_embedded&v=_738lFAZSUk
# " target="_blank"><img src="http://img.youtube.com/vi/_738lFAZSUk/0.jpg" 
# alt="IQT Demo" width="240" height="180" border="10" /></a>
[![IQT Demo](https://img.youtube.com/vi/_738lFAZSUk/0.jpg)](https://www.youtube.com/watch?v=_738lFAZSUk&feature=youtu.be)

### Table of Contents
0. [Introduction](#introduction)
0. [Data](#data)
0. [How to use IQT framework](#models)
0. [Citation](#citation)
0. [Disclaimer and known issues](#disclaimer-and-known-issues)


### Introduction
This repository contains Matlab codes for the __Image Quality Transfer (IQT)__ framework proposed in our recent NeuroImage paper "**Image quality transfer and applications in diffusion MRI**" [(link)](http://www.sciencedirect.com/science/article/pii/S1053811917302008). The current version supports super-resolution of diffusion tensor images (DTIs). 

IQT uses machine learning to transfer the rich information available from one-off experimental medical imaging devices to the abundant but lower-quality data from routine acquisitions. The procedure uses matched pairs to learn mappings from low-quality to corresponding high-quality images. Once learned, these mappings then augment unseen low quality images, for example by enhancing image resolution or information content. In the papaer, we demonstrate IQT using a simple patch-regression implementation and the uniquely rich diffusion MRI data set from the human connectome project (HCP). Results highlight potential benefits of IQT in both brain connectivity mapping and microstructure imaging. In brain connectivity mapping, IQT reveals, from standard data sets, thin connection pathways that tractography normally requires specialised data to reconstruct. In microstructure imaging, IQT shows potential in estimating, from standard “single-shell” data (one non-zero b-value), maps of microstructural parameters that normally require specialised multi-shell data. Further experiments show strong generalisability, highlighting IQT's benefits even when the training set does not directly represent the application domain. The concept extends naturally to many other imaging modalities and reconstruction problems. For details, please refer to the original paper!

### Data
This demostration requires the in-vivo human diffusion MRI dataset from **Human Connectome Project (HCP)**, which is freely downloadable from [here](http://www.humanconnectome.org/documentation/MGH-diffusion/index.html). To try the IQT framework for DTI super-resolution, you need to dowload the diffusion data with b-value = 1000. 

### How to use IQT framework
The pipeline works in three stages: 1. preprocessing (`preprocess.m`); 2. training (`train.m`); 3. testing (`test.m`).

<table>
<tr><td>preprocess.m  </td><td> preprocess the raw data downloaded from HCP, and generate training sets.
</td></tr> <tr><td>train.m </td><td> train regression trees on specified data sets.
</td></tr><tr><td>test.m </td><td>  perform super-resolution on a given diffusion tensor image
</td></tr></table>

You need to point to paths on your system appropriately. Some trained trees are provided in this repository and can be found in the directory `trees`. For fast trial of this pipeline, we recommend you to skip the most time-consuming step 2, and perform the reconstruction with these already trained trees. The default experiment settings perform x3 super-resolution with input patch 5x5x5. By default, the reconstruction on the outer boundary of the brain requires a separate reconstruction method and is ignored. If you wish to perform boundary complettion, just set `settings.edge = 1` in `test.m`. Reconstruction takes about 10 minutes/2 hours for each subject with 8 trees without/with boundary completion. 

`test.m` also automatically visualises the mean diffusivity (MD), fractional anisotropy (FA) and colour encoded directional map (CFA) of the predicted high-res DTI, and saves as a FIG file. For example, the figure below illustrates results of x3 super-resolution with 3x3x3 input patch on subject 117324. Note that no boundary completion was performed here. 

![DTI_SR_3x3x3_3x3x3](https://cloud.githubusercontent.com/assets/14926992/24544089/e2e18f72-15f9-11e7-8f7c-0488a8b197aa.png)


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

### Disclaimer
0. The current version is limited to super-resolution of diffusion tensor images. 
0. Boundary reconstruction can take some time (~2 hours per subject). It's best to try the pipeline first without this feature (makes sure the parameter `settings.edge = 0`).
