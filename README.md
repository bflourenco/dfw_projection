This repository contains the code for the numerical experiments and dual Frank-Wolfe solver of 
the paper *Projection onto hyperbolicity cones and beyond: a dual Frank-Wolfe approach*
by Takayuki Nagano, Bruno F. Lourenço and Akiko Takeda.

## Main files

The dual Frank-Wolfe solvers can be found in files 
*solver/FW_HP.m, solver/FW_HP_exp.m, solver/FW_GCP_exp.m*.
The implementation of Renegar's accelerated gradient method 
for hyperbolicity cones can be found in the file *solver/AGM_HP.m*.
For information on how to invoke these functions, see the files' documentation 
or see the scripts in the *experiments* folder for further examples.

The file *solver/poly_proj.m* contains code specialized in 
projecting onto hyperbolicity cones, see file 
*solver/poly_proj_examples.m* for examples.



## Reproducing the experiments in the paper

### Prerequisites
 1. Install the DDS solver:[https://mehdi-karimi-math.github.io/DDS.html](https://mehdi-karimi-math.github.io/DDS.html)
 2. Install the Mosek solver and the Matlab interface: [https://docs.mosek.com/latest/toolbox/install-interface.html](https://docs.mosek.com/latest/toolbox/install-interface.html)
 3. Open Matlab and add the *solver* folder to the Matlab path

### Reproducing Tables 1,2,3,4,5,6
 1. Run the file *experiments/experiments_batch_hyper.m*. It may take a whole day to complete.
 2. Run the file *experiments/experiments_batch_pcone.m*. It will take a few minutes to complete.
 3. Open a terminal shell, go to the *experiments* folder and run the file *proj_data_sorter.sh*.
This will organize the csv files generated in the previous two steps.

If all the scripts are completed successfully, then the following files will be generated.

  * Files for table 1:
    * proj_hyper_n10_d1_30_tol_high.csv
    * proj_hyper_n10_d2_30_tol_high.csv
    * proj_hyper_n20_d1_30_tol_high.csv
    * proj_hyper_n20_d2_30_tol_high.csv
  * Files for table 2:
    * proj_hyper_n30_d27_30_tol_high.csv
    * proj_hyper_n40_d37_30_tol_high.csv
    * proj_hyper_n50_d47_30_tol_high.csv
  * Files for table 3:
    * proj_hyper_n10_d1_30_tol_low.csv
    * proj_hyper_n10_d2_30_tol_low.csv
    * proj_hyper_n20_d1_30_tol_low.csv
    * proj_hyper_n20_d2_30_tol_low.csv
  * Files for table 4:
    * proj_hyper_n30_d27_30_tol_low.csv
    * proj_hyper_n40_d37_30_tol_low.csv
    * proj_hyper_n50_d47_30_tol_low.csv
  * Files for table 5:
    * proj_pcone_n300_tol_high_all.csv
    * proj_pcone_n500_tol_high_all.csv
    * proj_pcone_n1000_tol_high_all.csv
  * Files for table 6:
    * proj_pcone_n300_tol_low_all.csv
    * proj_pcone_n500_tol_low_all.csv
    * proj_pcone_n1000_tol_low_all.csv

Because in the paper we used relative time taking the DDS and Mosek running times as 
baselines, we expect that even if the experiments 
are run in the different computers the times should be quite close to the ones we 
reported.

### Reproducing Figures 1 and 2
See file *experiments/Ex_proj_hyper_eleSym_only_graph.m*

## License
The files under the *solver* folder are licensed under the GPLv3. The files 
 under the *experiments* folder are licensed under the MIT License.

Copyright (C) 2024 Takayuki Nagano, Bruno F. Lourenço, Akiko Takeda