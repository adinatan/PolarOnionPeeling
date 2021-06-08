[![View PolarOnionPeeling on File Exchange](https://www.mathworks.com/matlabcentral/images/matlab-file-exchange.svg)](https://www.mathworks.com/matlabcentral/fileexchange/41064-polaronionpeeling)
# PolarOnionPeeling
The code implements the method shown at "Toward real-time charged-particle image reconstruction using polar onion-peeling" by G.M. Roberts et-al, Rev. Sci. Instr. 80, 053104 (2009)  https://doi.org/10.1063/1.3126527

There are several extensions included beyond the original paper:

* Support of higher Legendre orders
* Specific quadrants of the image can be selected in the analysis in case of image asymmetries.
* Limit Legendre order toward the origin to restrict overfitting.
* Support of NaN values for masking unwanted pixels or areas of raw images.
* Supports images of up to 4K x 4K pixels. bld_lut.m was added for radial basis set generation of arbitrary size. 
* Added beta2cart.m to get the 2D and 3D reconstruction based on beta parameteres. 
* Added LDSD.m to apply Legendre decomposition for a single data vector spanning [0,2*pi] 
* Added cart2tripolar.m to allow only cart2polar transform for any quandrant an to the range  [0,2*pi].

Run the script pop_example to see the various outputs of POP.m
  
  ![Fig1](https://github.com/adinatan/PolarOnionPeeling/blob/master/fig1.png)
  
  ![Fig2](https://github.com/adinatan/PolarOnionPeeling/blob/master/fig2.png)
  
  ![Fig3](https://github.com/adinatan/PolarOnionPeeling/blob/master/fig3.png)
  
  ![Fig4](https://github.com/adinatan/PolarOnionPeeling/blob/master/s3.png)
