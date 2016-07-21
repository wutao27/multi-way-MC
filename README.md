# Multi-way Monte Carlo Method for Linear Systems

#### Tao Wu
#### David F. Gleich
------

### Files
* `hhead.m` is the function for computing the tensor \head{H}.
* `varian.m` is the function for computing the random walk variances of linear system x = Hx+b.
* `script.m` is the script code for various experimental tasks described in the paper. Task 1: compute the success ratio. Task 2: compute the speed up ratio.

### Usage
Please refer the the demo code blocks in `script.m`:
* compute success ratio for synthetic experiments.
* plot the success ratio for synthetic experiments.
* compute the speed-up ratio for synthetic experiments.
* compute the average speed-up ratio for the Harwell-Boeing sparse matrix collection.
