# pac-imbalanced

Learning Reliable Rules under Class Imbalance (SDM 2021)
========================================================

You can compile the source code on command line with the command:
g++ experiments.cpp

The main function is between the lines 1166 - 1175.

* Run experiments for the algorithm Find-S for PAC learning:
  - Uncomment line 1170
  - Comment line 1172
  - Recompile the program and run the executable
  - Several text files will be generated with results and statistics on these results.

* Run experiments for the Swapping algorithm from the framework of evolvability:
  - Comment line 1170
  - Uncomment line 1172
  - Recompile the program and run the executable
  - Several text files will be generated with results and statistics on these results.

Further refinement on the experiments can be obtained with the instructions provided below.



Find-S
======
In lines 570 and 571 we determine the sample size needed for learning.
* Line 570 corresponds to the traditional PAC learning framework that only cares about low risk.
* Line 571 corresponds to our extended PAC learning framework and the sample size is determined by combining Theorem 3.1 and Theorem 2.1 (as described in the paper).

Therefore, depending on the results that you want to observe, you need to comment/uncomment the above lines. Only one should be active at each time. 

* Expected running times for experiments (2.9 GHz Dual-Core Intel Core i5, 16GB RAM):
  - Traditional PAC learning: Less than 5 minutes
  - Extended PAC learning (our framework): About 6 hours and 20 minutes


* Experimenting with distributions beyond Uniform:
As mentioned in the paper, we can experiment with distributions beyond uniform. 
- Line 20 has the constant PROB_THRESHOLD_PER_THOUSAND, which is currently set to 500, corresponding to the uniform distribution as this is the numerator when the denominator is 1000. Thus, the current value corresponds to the uniform distribution where each variable is satisfied with probability 500/1000 = 0.5.
- Note that in the current version of the code, the above change is enough for Find-S. Regarding the Swapping Algorithm, additional modifications are needed when we want to test in distributions beyond uniform.


Swapping Algorithm
==================
* For the Swapping algorithm, as explained in the paper (end of Section 4.2.2), applying Theorem 3.1, results in the error parameter epsilon being very small and as a consequence the algorithm always finds the ground truth function completely.
* So, the relevant part for execution has to do with the Swapping algorithm being run in order to satisfy the traditional PAC criterion, which subsequently gives the results for Table 2.

* Expected running time for experiments (2.9 GHz Dual-Core Intel Core i5, 16GB RAM):
  - Traditional PAC learning: About 4 minutes
  - Extended PAC learning (our framework): - (No experiments needed.)
