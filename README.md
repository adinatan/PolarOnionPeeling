[![View PolarOnionPeeling on File Exchange](https://www.mathworks.com/matlabcentral/images/matlab-file-exchange.svg)](https://www.mathworks.com/matlabcentral/fileexchange/41064-polaronionpeeling)
# PolarOnionPeeling
The code implements the method shown at "Toward real-time charged-particle image reconstruction using polar onion-peeling" by G.M. Roberts et-al, Rev. Sci. Instr. 80, 053104 (2009). It supports images of up to 4K and accepts arbitrary (even) beta parameters. Specific quadrants of the image can be selected in the analysis in case of image asymmetries. NaN values can be used for masking unwanted pixels or areas.
For more info regarding the method see: https://doi.org/10.1063/1.3126527

There are several extensions included beyond the original paper:

* Support of higher Legendre orders
* Limit Legendre order toward the origin to reduce overfitting errors
* Support of NaN values for masking the raw images
* The code supports images of up to 4K. bld_lut.m was added for radial basis set generation of arbitrary size. 

Run the script pop_example to see the various outputs of POP.m
  
  ![Fig1](https://github.com/adinatan/PolarOnionPeeling/blob/master/fig1.png)
  
  ![Fig2](https://github.com/adinatan/PolarOnionPeeling/blob/master/fig2.png)
  
  ![Fig3](https://github.com/adinatan/PolarOnionPeeling/blob/master/fig3.png)
  
  ![Fig4](https://github.com/adinatan/PolarOnionPeeling/blob/master/s3.png)
