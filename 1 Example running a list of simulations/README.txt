# Running multiple input files

inputlist.csv contains the input files required to run 5 simulations with the same parameters and varying seeds.

sim_list.sh is a bash script which will split inputlist.csv into n=5 folders, generate input files, run the simulation and then analyse the results. 

First compile the simulation to "sim" with the command "make", then run "sim_list.sh" as follows:

make
./sim_list.sh

sim_list.sh can be used in conjunction with a compiled sim to run all the simulation that generate the figures in the paper. 