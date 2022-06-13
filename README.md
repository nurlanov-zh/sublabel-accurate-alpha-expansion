# Efficient and Flexible Sublabel-accurate Energy Minimization

## Usage  

### Tested on
  
``` 
Matlab R2015-b 
cmake 3.10.2 
g++ 7.5.0
Ubuntu 18.04
  
Nvidia GeForce GTX 1050, Driver Version: 440.100
CUDA 10.2
``` 
  
### Quick Try Out
* Install all third-parties into `third-party` folder.

* Run MATLAB examples in `image_denoising.m`

#### (Optinally, requires external solver) Sublabel-accurate refinement:
- If you want to reproduce the same timings, consider adding external 
solvers (Gurobi). Add their paths into Matlab by editing `startup.m`.

- Run the test file `third-party/yalmip/yalmiptest.m` from Matlab to check 
if the external solvers are found.

#### (Optionally, requires GPU) Previous sublabel-accurate methods: 
```sh
cd third-party/prost/
mkdir build
cd build
cmake .. 
make -j 4 
  
cd ../../sublabel_relax/cvpr2016/ 
mex compute_convex_conjugate.cpp
``` 

To run previous sublabel-accurate methods consider starting MatlabR2015-b 
via the following command:
```sh
LD_PRELOAD="/usr/lib/x86_64-linux-gnu/libstdc++.so.6" matlab
```
or by running bash script `run_matlab.sh`.
  
  
## References 

**Important**: Please be aware of the Licences and Copyrights of the used third-party software. Check Readme and/or Licence files in corresponding `third-party` folders. 


Code is based on the following third-parties:

- [GCO-v3.0](https://github.com/nsubtil/gco-v3.0)
- [prost](https://github.com/tum-vision/prost)
- [sublabel-relax](https://github.com/tum-vision/sublabel_relax)
- [yalmip](https://yalmip.github.io/)

And on the following publications:

[1] D. Schlesinger and B. Flach, ["Transforming an arbitrary MinSum problem into a binary one"](http://www1.inf.tu-dresden.de/~ds24/publications/tr_kto2.pdf), Technical report TUD-FI06-01, Dresden University of Technology, April 2006.
       
[2] Yuri Boykov and Vladimir Kolmogorov, "An Experimental Comparison of Min-Cut/Max-Flow Algorithms for Energy Minimization in Computer Vision", TPAMI, September 2004. 
 
[3] Vladimir Kolmogorov and Ramin Zabih, "What Energy Functions can be Minimized via Graph Cuts?", TPAMI, February 2004. 
         
[3] Shai Bagon "Matlab Implementation of Schlezinger and Flach Submodular Optimization", June 2012.

[4] Yuri Boykov, Olga Veksler, and Ramin Zabih, ["Fast approximate energy minimization via graph cuts"](http://www.cs.cornell.edu/rdz/Papers/BVZ-iccv99.pdf), TPAMI, November 2001.

[5] Thomas Mollenhoff, Emanuel Laude, Michael Moeller, Jan Lellmann, and Daniel Cremers, ["Sublabel-accurate relaxation of nonconvex energies"](https://arxiv.org/pdf/1512.01383.pdf), CVPR, 2016.

[6] Thomas Pock, Daniel Cremers, Horst Bischof, and Antonin Chambolle, "Global solutions of variational models with convex regularization", SIAM Journal on Imaging Sciences, 2010.

* See full list of references in our paper
