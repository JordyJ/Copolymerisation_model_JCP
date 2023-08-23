# Setup

If compiling on windows, uncomment //#include <windows.h> in Source.cpp

Source.cpp is configured to run until a fixed number of polymers are produced. 

A simulation which ends at a specific time can be used instead by uncommenting line 53 and commenting out line 54 in Source.cpp
// Simulation termination conditions

Compile Source.cpp and PolymerSystemIO to sim using Makefile, by running "make".
This compiles the simulation as "sim".

# Input file
Input.csv contains the input parameters for the simulation.

int seed: Seed sets the seed used for RNG.
int length_limit: Sets the template length. Max value is 101. (Can be modified in PolymerSystem.h).
double time_limit: Sets the time at which the simulation terminates (if sim is compiled with line 53 uncommented and 54 commented in Source.cpp).
int max_pols: The number of polymers which will be produced before the sim terminates.
double k0: Standard rate. Set to 1.
double k: Polymerisation rate.
double (kact: Activation rate. Set to 0. )
double ConcEff: Effective monomer concentration of units in polymer tails.
double (Gact: Activation energy. Set to 0.)
double Gbb: Backbone bond free energy. <0 is stable
double Conc0: [0] Concentration of type 0 monomers. 
double Conc1: [1] Concentration of type 1 monomers. 
double G0: Free energy of copy-template interaction for type 0 monomers.
double G1: Free energy of copy-template interation for type 1 monomers. 
double Ggen: Free energy of generic bond. 
double Gend: Free energy penalty at last template site.
string outputfilename: File name for writing polymers which are produced. Set to "output.csv"
string outputtrajfilename: File name for writing trajectory information. Set to "traj.csv"
{0,1} print_flag: Set to 1 to print state of the system at each step of the simulation.
{0,1} print_valid_trans_flag: Set to 1 to print the valid transitions at each site at each step of the simulation
{0,1} write_output_flag: Set to 1 to write polymers to output file when released. 
{0,1} write_traj_flag: Set to 1 to write trajectory information to trajectory file at set times. 
double traj_write_interval: Duration between writing trajectories to file.
{0,1} breakable_backbone_flag: Set to 1 to enable depolymerisation of polymers. Set to 0 to simulate irreversible polymerisation.
{0,1} off_rate_discrim_flag: Set to 1 for off rate discrimination. Set to 0 for on rate discrimination.
{0,1} t_activation_flag: Set to 0.
{0,1} local_updates_flag: Set to 1 to update the system locally. Speeds up simulation for long templates. Set to 0 if you have time to kill.
{0,1} use_test_rates_flag: Set to 1 to use default hardcoded rates. Useful for debugging. 

# Running a simulation
Run as "./sim input.csv"

This will generate trajectories in traj.csv and a list of produced polymers in output.csv.

# Analysis
Generate a length distribution from the ouput file by running "python plot_output.py output.csv".
This will generate OuputAnalysis
"plot_output.py" requires:
matplotlib      3.1.3
numpy           1.18.1
pandas          1.0.1



Email j.juritz18@ic.ac.uk or t.ouldridge@ic.ac.uk for contact.




