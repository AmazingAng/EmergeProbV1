# EmergeProbV1
Code for paper: Emergence of Probabilistic Representation in Cortical Network of Mice Primary Visual Cortex (under review). The code are all written in [R](https://www.r-project.org/). The code uses parallel computing to speed-up simulation, which takes ~15 min to run 15 repetitions of the simulation on a server using 15 cores. If your computer has less cores, try to set repeat_n (number of repetition) smaller to save time. 

## Code for each figure (Figure 3 & 5 are figures of experimental results)
Figure 1: Fig1_stripe_rearing_plot.R, Fig1_stripe_rearing_simulation.R

Figure 2 & 4: Fig24_cardianl_overrepresentation.R

Figure 6: Fig6_uniform_inh.R

Figure 7: Fig7_phase_plane.R

Figure 8: Fig8_bayesian.R

Useful functions: utils.R

## Dependent library in R:
Multi-processing library: parallel, foreach, doParallel

Other: circular, distr, radial, plotrix, phaseR
