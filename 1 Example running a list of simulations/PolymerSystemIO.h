#pragma once
class PolymerSystemIO
{

	/* October 2021 */
public:
	long double time_limit = 100;

	static const unsigned int length_limit_max = 101;           // Maximum size of the template length, used to predefine array lengths Think about making this dynamic

	unsigned int length_limit = 10;           // Maximum size of the template length, used to predefine array lengths Think about making this dynamic
	unsigned int pols_limit = length_limit; // Limit of the number of simulatenously growing polymers on the template

	std::string inputfilename = "input.csv"; // Filenames for outputting terminated copies.

	std::ofstream outputfile; // Stream for outputting the terminated copies
	std::string outputfilename = "defaultoutput.csv"; // Filenames for outputting terminated copies.
	std::ofstream outputtrajfile; // Stream for outputting the trajectory (observables such as # copy template bonds, # polymers attached # active template sites
	std::string outputtrajfilename = "defaulttraj.csv";

	// State variables
	int seed = 1;  // Random number seed
	std::mt19937_64 psRandNumGen; // Mersenne twister PRNG
	long double time = 0;

	int site_to_update;                 // Used to pass output from select_site_to_update to  select_transition and apply_transition;
	int transition_to_update;           // Used to pass output from select_transition to apply_transition;

	static const int  t_alph_size = 1;          // Number of distinct monomer types in the template alphabet.
	int  t_seq[length_limit_max];                   // Template sequence. Array contains (monomer type) integer i  at site j describing the monomer at site j. -1 if empty

	static const int c_alph_size = 2;           // Number of distinct monomer types in the copy alphabet.
	int  c_seq[length_limit_max][length_limit_max];     // Polymer sequence. Array contains integer i at site j describing the monomer at site j. -1 if empty

	bool t_c_bond[length_limit_max];                // true if there is a monomer bound to the template at site j.
	int  t_occupied_by[length_limit_max];           // site j takes the value integer i if polymer label i monomer is bound to the template at site j. -1 if empty

	bool t_active[length_limit_max];                // True at j if the template is active at site j.
	bool bb_bond[length_limit_max][length_limit_max];   // True at [k][j] if the the polymer label k monomer at site j shares a bond with neighbour at j+1.

	bool polymer_present[length_limit_max];       // True if there is some info about a polymer with label i in row i of other state variables.
												// Used to allocate space for new polymers to grow in other arrays.
	// Parameters for the rates

	// Negative energies are more stable, stronger bonds
	// Default parameters
	double k0 = 1; // Rate of monomer binding
	double k = 1;  // Rate of polymerisation
	double kact = 1; // Rate of template deactivation
	double ConcEff = 100; // Effective monomer concentration
	double Gact = 0;  // Free energy change when deactivating the template
	double Gbb = -5;   // Free energy of backbone formation
	double Conc[2] = { 1, 1 }; // Monomer concentrations up to binary
	double Gx[2] = { -3,-1 };  // Monomer binding energies to template up to binary
	double Ggen = 0;  // Free energy of generic template bond
	double Gend = 2;   // Free energy that destabilises the final site on the template

	static const unsigned int num_rules = 17 * c_alph_size + 4;

	long double total_rate_per_site[length_limit_max];
	long double rates[num_rules];
	bool        valid_transitions[num_rules][length_limit_max];

	bool print_flag; // True to print state on each iteration. Flase suppresses output
	bool write_output_flag; // True will write terminated polymers to file
	bool write_traj_flag; // True will write trajectories to file
	int max_pols = 10; // Maximum number of polymers to write to file
	int pols_written = 0;
	double traj_write_counter = 0; // Write trajectory every traj_write_counter secs
	double traj_write_interval; // Write trajectory every traj_write_counter secs

	bool breakable_backbone_flag; // True enables backbone to make/break everywhere, false polymerisation can only occur at the tip
	bool t_activation_flag; // True enables template activation/inactivation
	bool local_updates_flag; // True enables local updating of the valid transitions matrix
	int  local_update_window = 2;
	bool print_valid_trans_flag; // True if print state will also display vaild transitions matrix
	bool use_test_rates_flag; // True then will use test rates instead of physical ones. 
	bool off_rate_discrim_flag; // If True, then discrimination is on the monomer unbinding rates. If False, then discrimination is on the ON binding rates.

	PolymerSystemIO(); // Initialises the polymer system class
	void setup();      // Triggers Reads in data from input file and then triggers initialisation
	void update_valid_transitions(); // Updates the valid transition matrix based upon the current state across a limited range .
	void sum_rates();                // Sums valid transitions at each site to calculate total_rate_per_site.

	int select_site_to_update();  // Selects site on the template to be updated and updates time.
	int select_transition();      // Given a site on the template, selects a valid transition
	void apply_transition();      // implements the transition using rules of the model of choice
	//void calculate_rates();    // Will need to be included in models that have explicit template sequences
	void grow();                // Macro for implementing the algorithm

	bool is_monomer(int idx);

	void generate_initial_condition(int start_idx, int end_idx, int attached_by, int temp_row);
	void set_seed(int s);
	void set_rates(); // Initialises the rate matrix
	void init_length(); // Initialises the length of the template

	void print_valid_transitions();
	void print_state();
	void print_state_vertically();
	void print_rate_per_site();
	void print_valid_transitions_per_site(int site_idx);

	void write_copy_output(); // If a polymer terminates, write the output to a file
	void write_traj_output();
	void init_output();        // Create the input and output filestreams
	void read_input();        // read the input file


	char rule_labels[num_rules][30] = { "+0 monomer", "+1 monomer",
								"-0 monomer deact ahead", "-1 monomer deact ahead",
								"-0 monomer no deact", "-1 monomer no deact" ,
								"-0 monomer at end", "-1 monomer at end" ,
								"-0 upstream","-1 upstream",
								"+0 upstream", "+1 upstream",
								"Polymerise", "Depolymerise",
								"Activate template","Deactivate template",
								"Terminate lead 0 deact ahead","Terminate lead 1 deact ahead",
								"Terminate lead 0 no deact","Terminate lead 1 no deact",
								"Terminate lead 0 at end", "Terminate lead 1 at end",
								"Terminate mid 0 deact ahead","Terminate mid 1 deact ahead",
								"Terminate mid 0 no deact","Terminate mid 1 no deact",
								"-0 poly deact ahead", "-1 poly deact ahead",
								"-0 poly no deact", "-1 poly no deact" ,
								"-0 leading poly deact ahead", "-1 leading poly deact ahead",
								"-0 leading poly no deact", "-1 leading poly no deact" ,
								"-0 leading poly at end", "-1 leading poly at end",
								"+0 downstream" ,"+1 downstream" };

};