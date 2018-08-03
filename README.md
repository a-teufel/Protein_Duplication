 Simulation, figure rendering, and analysis files for "The many nuanced evolutionary consequences of duplicated genes" by Teufel et. al.
 ====================================================
 
This repository contains the simulation scripts written in Python, analysis scripts used to parse the simulated results, and the figure rendering scripts written in R for the above-mentioned manuscript. A preprint of the manuscript is available at  https://www.biorxiv.org/content/early/2018/07/10/366971. This manuscript examines protein evolution when one memeber of a binding pair is duplicated with the use of an evolutionary simulation. Our simulations are based on those published in "Accelerated simulation of evolutionary trajectories in origin--fixation models" by Ashley I. Teufel and Claus O. Wilke.  Alignments of the protein sequences produced by our simulations can be found in the Fasta folder. The naming works as follows: SIM1 is bind both, SIM2 is bind B and not B', SIM3 is bind max, SIM4 is no bind, and SIM5 is bind B.


# Dependencies

To run the simulation scripts you must have python 2.7, Rosetta, and Pyrosetta installed. Pyrosetta must be included in the path that you are running the scripts from.  Most of the figure scripts depend on the R libraries cowplot, ggplot2, and grid.
