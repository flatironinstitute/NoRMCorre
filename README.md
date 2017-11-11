[![Join the chat at https://gitter.im/epnev/ca_source_extraction](https://badges.gitter.im/Join%20Chat.svg)](https://gitter.im/epnev/ca_source_extraction?utm_source=badge&utm_medium=badge&utm_campaign=pr-badge&utm_content=badge)

# NoRMCorre: Non-Rigid Motion Correction 
This package provides a Matlab implementation of the NoRMCorre algorithm [[1]](#ref), and can be used for online piecewise rigid motion correction of 2d (planar) or 3d (volumetric) calcium imaging data. 

## Synopsis

The algorithm operates by splitting the field of view into a set of overlapping patches. For each patch and each frame a rigid translation is estimated by aligning the patch against a template using an efficient, FFT based, algorithm for subpixel registration [[2]](#reg). The estimated set of translations is further upsampled to a finer resolution to create a smooth motion field that is applied to a set of smaller overlapping patches. Extra care is taken to avoid smearing caused by interpolating overlapping patches with drastically different motion vectors. The registered frame is used to update the template in an online fashion by calculating a running/mean of past registered frames. The pipeline is summarized in the figure below.

![Alt text](pipeline.png?raw=true "piecewise rigid motion correction pipeline")

## Code details

See the function [```demo.m```](https://github.com/simonsfoundation/NoRMCorre/blob/master/demo.m) for an example of the code. The algorithm is implemented in the function ```normcorre.m```. If you have access to the parallel computing toolbox, then the function ```normcorre_batch.m``` can offer speed gains by enabling within mini-batch parallel processing. The user gives a dataset (either as 3D or 4D tensor loaded in RAM or memory mapped, or a pointer to a .tiff stack or .hdf5 file), and a parameters struct ```options```. Optionally, an initial template can also be given. The algorithm can also be used for motion correction of 1p micro-endoscopic data, by estimating the shifts on high pass spatially filtered version of the data. See the script [```demo_1p.m```](https://github.com/simonsfoundation/NoRMCorre/blob/master/demo_1p.m) for an example.

The algorithm can also be ran using the MotionCorrection object. See [```demo_mc_class.m```](https://github.com/simonsfoundation/NoRMCorre/blob/master/demo_mc_class.m) for an example on how to use the object for 2p and 1p data.

The ```options``` struct can be set either manually or by using the function ```NoRMCorreSetParms.m```. Some parameters of the ```options``` struct are the following:

| Parameter name | Description |
|----------------|-------------|
| ```d1,d2,d3``` | dimensions of field of view |
| ```grid_size``` | size of non-overlapping portion of the grid in each direction (x-y-z)|
| ```overlap_pre```| size of overlapping region in each direction before upsampling  |
| ```mot_uf```    | upsampling factor for smoothing and refinement of motion field |
| ```overlap_post ``` | size of overlapping region in each direction after upsampling |
| ```max_shift``` | maximum allowed shift for rigid translation | 
| ```max_dev``` | maximum deviation of each patch from estimated rigid translation |
| ```upd_template``` | update the template online after registering some frames |
| ```bin_width``` | length of bin over which the registered frames are averaged to update the template |
| ```init_batch``` | number of frames to be taken for computing initial template |
| ```iter``` | number of times to go over the dataset |
| ```output_type``` | type of output registered file |
| ```phase_flag``` | flag for using phase correlation |
| ```correct_bidir``` | check for offset due to bidirectional scanning (default: true) |

The performance of registration can be evaluated using the function ```motion_metrics.m```. The function simply computes the correlation coefficient of each (registered) frame, with the mean (registered) frame across time, the mean registered frame, and its crispness.

## Developers

[Eftychios A. Pnevmatikakis](https://github.com/epnev), Flatiron Institure, Simons Foundation

## External packages

This package includes functions from the following packages
- [Save and load a multipage tiff file](https://www.mathworks.com/matlabcentral/fileexchange/35684-save-and-load-a-multiframe-tiff-image/content/loadtiff.m)
- [Savefast](https://www.mathworks.com/matlabcentral/fileexchange/39721-save-mat-files-more-quickly) for saving (and then loading) MAT files more quickly without compressing their contents. 
- [Eficient subpixel registration by cross-correlation](https://www.mathworks.com/matlabcentral/fileexchange/18401-efficient-subpixel-image-registration-by-cross-correlation) for fast alignment of an image against a template.

## Integrations 

This package will be integrated with the [Matlab code](https://www.github.com/epnev/ca_source_extraction) for source extraction and deconvolution using CNMF.

A python version of this algorithm developed from [Andrea A. Giovannuci](https://github.com/agiovann) is included as part of the [CaImAn](https://github.com/simonsfoundation/CaImAn) package that provides a complete pipeline for calcium imaging data pre-processing.

Although the two implementations give almost identical results for the same input file, there are some slight differences in the way they are called and their capabilities. These differences are highlighted [here.](https://github.com/simonsfoundation/NoRMCorre/wiki/Differences-between-Matlab-and-Python-implementations)

## More details, contact information, and citing NoRMCorre

Check the [wiki](https://github.com/simonsfoundation/NoRMCorre/wiki) for more details and some frequently asked questions. 

Please use the [gitter chat room](https://gitter.im/epnev/ca_source_extraction?utm_source=badge&utm_medium=badge&utm_campaign=pr-badge&utm_content=badge) for questions and comments, and create an issue for any bugs you might encounter.

If you find this package useful please cite the following paper:

Eftychios A. Pnevmatikakis and Andrea Giovannucci, *NoRMCorre: An online algorithm for piecewise rigid motion correction of calcium imaging data*, Journal of Neuroscience Methods, vol. 291, pp 83-94, 2017; doi: [https://doi.org/10.1016/j.jneumeth.2017.07.031](https://doi.org/10.1016/j.jneumeth.2017.07.031)

## Acknowledgements

Example dataset is kindly provided from Andrea Giovannucci, taken at Wang lab (Princeton University).

## References 

<a name="ref"></a>[1] Eftychios A. Pnevmatikakis and Andrea Giovannucci, *NoRMCorre: An online algorithm for piecewise rigid motion correction of calcium imaging data*, Journal of Neuroscience Methods, vol. 291, pp 83-94, 2017; doi: [https://doi.org/10.1016/j.jneumeth.2017.07.031](https://doi.org/10.1016/j.jneumeth.2017.07.031)

<a name="reg"></a>[2] Guizar-Sicairos, M., Thurman, S. T., & Fienup, J. R. (2008). Efficient subpixel image registration algorithms. Optics letters, 33(2), 156-158. Matlab implementation available [here](https://www.mathworks.com/matlabcentral/fileexchange/18401-efficient-subpixel-image-registration-by-cross-correlation).
