#include <iostream>
#include <iomanip>
#include <algorithm>
#include <math.h>
#include <random>
#include <fstream>
#include <string>
#include "PolymerSystemIO.h"
using namespace std;

PolymerSystemIO::PolymerSystemIO()
{

	/* February 2021 */
	// Initialise the polymer state

	fill(&t_seq[0], &t_seq[0] + length_limit_max, 0); // Set template sequence to all wrong
	fill(&c_seq[0][0], &c_seq[length_limit_max - 1][length_limit_max], -1); // Set all initial copy sequences to  -1.
	fill(&t_c_bond[0], &t_c_bond[length_limit_max], false); // Initialise with no template copy bonds.
	fill(&bb_bond[0][0], &bb_bond[length_limit_max - 1][length_limit_max], false); // Initialise with no template copy bonds.
	fill(&t_occupied_by[0], &t_occupied_by[length_limit_max], -1); // Initialise with no label to copies bound to template	.
	fill(&valid_transitions[0][0], &valid_transitions[num_rules - 1][length_limit_max], false);  //set all transitions to initially invalid
	fill(&total_rate_per_site[0], &total_rate_per_site[length_limit_max], 0); // initialise the total rate array to 0
	fill(&polymer_present[0], &polymer_present[length_limit_max], false); // Initialise with no polymers present.
	fill(&t_active[0], &t_active[length_limit_max], false); // Initialise with fully active template.
	t_active[0] = true; // Start site should be active
	fill(&rates[0], &rates[num_rules], 0);  // Initialise the rates with 0 values.

};

void PolymerSystemIO::grow()
{
	

	update_valid_transitions(); // Update valid transitions matrix
	select_site_to_update(); // First layer of stochastic sampling selects site to update
	select_transition();  // Second layer of stochastic sampling selects the transition at the pre-selected site

	if (print_flag)
	{
				
		if (print_valid_trans_flag)
		{
			update_valid_transitions();
			print_valid_transitions();
		}

		//print_rate_per_site();
		//print_valid_transitions_per_site(site_to_update);
		print_state_vertically();
	}

	apply_transition(); // Update the state. Applies the selected transition at the selected site.

}

void PolymerSystemIO::update_valid_transitions()
{
	// Updates the valid transition matrix site by site. 
	// If local updates are enabled, then the update only occurs
	// a few indices either side of the last site to be updated. 
	// Otherwise the entire matrix is updated. 

	int start_idx;
	int end_idx;
	if (local_updates_flag == true && time != 0)
	{
		// set the start and end idx for updating to a few indices either side of the last site that was updated
		start_idx = site_to_update - local_update_window;
		end_idx = site_to_update + local_update_window;
	}
	else {
		// If local updates are enabled, at time ==0 we need to initialise the whole matrix.
		// Otherwise, if local updates aren't enabled then we update the whole matrix anyway.
		start_idx = 0;
		end_idx = length_limit;
	}

	// Prevent overflow errors
	if (start_idx < 0) { start_idx = 0; }

	if (start_idx >= length_limit) { start_idx = length_limit - 1; }

	if (end_idx < 0) { end_idx = 0; }

	if (end_idx > length_limit) { end_idx = length_limit; }

	if (end_idx < start_idx)
	{
		cerr << "\nINVALID INDEX RANGE IN update_valid_transitions\nExiting\n\n";
		exit(EXIT_FAILURE);
		return;
	}

	if (local_updates_flag == true && time != 0)
	{
		for (size_t rule_idx = 0; rule_idx < num_rules; rule_idx++)
		{
			fill(&valid_transitions[rule_idx][start_idx], &valid_transitions[rule_idx][end_idx], false);  //set all transitions to invalid
		}
	}
	else {
		// First set all transitions to invalid between start and end index
		fill(&valid_transitions[0][0], &valid_transitions[num_rules - 1][length_limit], false);  //set all transitions to invalid

	}


	/* Transitions where A is the alphabet size
	Index        Transition label
	0:A-1		"+0 monomer", "+1 monomer",
	A:2A-1		"-0 monomer deact ahead", "-1 monomer deact ahead",
	2A:3A-1		"-0 monomer no deact", "-1 monomer no deact" ,
	3A:4A-1		"-0 monomer at end", "-1 monomer at end" ,
	4A:5A-1		"-0 upstream","-1 upstream",
	5A:6A-1		"+0 upstream","+1 upstream",
	6A & 6A+1	"Polymerise", "Depolymerise",
	6A+2 & 6A+3	"Activate template","Deactivate template",
	6A+4:7A+3	"Terminate lead 0 deact ahead","Terminate lead 1 deact ahead",
	7A+4:7A+3	"Terminate lead 0 no deact","Terminate lead 1  no deact",
	8A+4:8A+3	"Terminate lead 0 at end", "Terminate lead 0 at end",
	9A+4:9A+3	"Terminate mid 0 deact ahead","Terminate mid 1 deact ahead",
	10A+4:10A+3	"Terminate mid 0 no deact","Terminate mid 1 no deact",
	11A+4:11A+3	"-0 poly deact ahead", "-1 poly deact ahead",
	12A+4:12A+3	"-0 poly no deact", "-1 poly no deact" ,
	13A+4:13A+3	"-0 leading poly deact ahead", "-1 leading poly deact ahead",
	14A+4:14A+3	"-0 leading poly no deact", "-1 leading poly no deact" ,
	15A+4:15A+3	"-0 leading poly at end", "-1 leading poly at end",
	16A+4:17A+3	"+0 downstream", "+1 downstream"
		*/

	// For each site on the template
	for (size_t site_idx = start_idx; site_idx < end_idx; site_idx++)
	{

		// Monomer BINDING
		// If the site is unoccupied (no template bond)
		// and if the site is active
		// then it is valid to add monomers
		if (t_c_bond[site_idx] == false &&
			t_active[site_idx] == true)
		{
			for (size_t transition_idx = 0; transition_idx < c_alph_size; transition_idx++)
			{
				valid_transitions[transition_idx][site_idx] = true;
			}
		}

		// Monomer UN-BINDING and deactivating template
		// If the site is occupied by a monomer
		// and if there is a site ahead to deactivate
		// then it is valid to remove monomer
		// 
		if (is_monomer(site_idx) == true &&
			site_idx < length_limit - 1 &&
			t_c_bond[site_idx + 1] == false &&
			t_active[site_idx + 1] == true &&
			t_activation_flag == true)
		{
			// Set the transition for removing unpolymerised monomer type c_seq[t_occupied_by[site_idx]][site_idx] to true
			valid_transitions[c_seq[t_occupied_by[site_idx]][site_idx] + c_alph_size][site_idx] = true;
		}

		// Monomer UN-BINDING without deactivating
		// If the site is occupied by a monomer
		// But the site ahead cannot be deactivated due to the presence of another monomer or the end of the template
		// then it is valid to remove monomer 
		if (is_monomer(site_idx) == true &&
			site_idx < length_limit - 1 &&
			((t_c_bond[site_idx + 1] == true || t_active[site_idx + 1] == false) || t_activation_flag == false)
			)
		{
			// Set the transition for removing unpolymerised monomer type c_seq[t_occupied_by[site_idx]][site_idx] to true
			valid_transitions[c_seq[t_occupied_by[site_idx]][site_idx] + 2 * c_alph_size][site_idx] = true;
		}

		// Monomer UN-BINDING at end of template
				// If the site is occupied by a monomer
				// But the site ahead cannot be deactivated due to the end of the template
				// then it is valid to remove monomer 
		if (site_idx == length_limit - 1 &&
			is_monomer(site_idx) == true)
		{
			// Set the transition for removing unpolymerised monomer type c_seq[t_occupied_by[site_idx]][site_idx] to true
			valid_transitions[c_seq[t_occupied_by[site_idx]][site_idx] + 3 * c_alph_size][site_idx] = true;
		}

		// UPSTREAM UN-BINDING
		// If the site is occupied
		// If the site+1 ahead is occupied by same polymer
		// if the site-1 is not bound by the same polymer or is the start of the template
		// and site and site+1 are polymerised
		if (t_c_bond[site_idx] == true &&
			site_idx < length_limit - 1 &&
			t_occupied_by[site_idx + 1] == t_occupied_by[site_idx] &&
			((t_occupied_by[site_idx - 1] != t_occupied_by[site_idx] && site_idx > 0) || site_idx == 0) &&
			bb_bond[t_occupied_by[site_idx]][site_idx] == true)
		{
			// Set the transition for breaking template bond of upstream monomer type c_seq[t_occupied_by[site_idx]][site_idx] to true
			valid_transitions[c_seq[t_occupied_by[site_idx]][site_idx] + 4 * c_alph_size][site_idx] = true;
		}



		// UPSTREAM BINDING
		// If the site is unoccupied (no template bond)
		// If the site+1 ahead is occupied by a polymer
		// and there is a monomer in the tail to bind upstream
		// and if the site is active
		// then it is valid to bind polymers upstream
		if (site_idx < length_limit - 1 &&
			t_c_bond[site_idx] == false &&
			t_c_bond[site_idx + 1] == true &&
			t_active[site_idx] == true &&
			bb_bond[t_occupied_by[site_idx + 1]][site_idx] == true)
		{
			valid_transitions[c_seq[t_occupied_by[site_idx + 1]][site_idx] + 5 * c_alph_size][site_idx] = true;
		}


		if (breakable_backbone_flag == true)
		{
			// POLYMERISATION
		// If the site is occupied
		// If the site+1 ahead is occupied
		// by a different polymer
		// with no backbone
		// 
			if (site_idx < length_limit - 1 &&
				t_c_bond[site_idx] == true &&
				t_c_bond[site_idx + 1] == true &&
				t_occupied_by[site_idx] != t_occupied_by[site_idx + 1] &&
				bb_bond[t_occupied_by[site_idx]][site_idx] == false &&
				bb_bond[t_occupied_by[site_idx + 1]][site_idx] == false
				)
			{
				// Polymerisation is valid here
				valid_transitions[6 * c_alph_size][site_idx] = true;
			}

			// DEPOLYMERISATION
			// If the site is occupied
			// If the site+1 ahead is occupied by the same polymer chain
			// if the site and site+1 do share a backbone bond
			if (site_idx < length_limit - 1 &&
				t_c_bond[site_idx] == true &&
				t_occupied_by[site_idx] == t_occupied_by[site_idx + 1] &&
				bb_bond[t_occupied_by[site_idx]][site_idx] == true)
			{
				// Depolymerisation is valid here
				valid_transitions[6 * c_alph_size + 1][site_idx] = true;
			}
		}
		else
		{	// If breakable backbones are not allowed everywhere then we have a tighter ruleset

			// POLYMERISATION
			// If the site is occupied
			// If the site+1 ahead is occupied by a monomer
			// (by a different polymer)
			// with no backbone
			// 
			if (site_idx < length_limit - 1 &&
				t_c_bond[site_idx] == true &&
				t_c_bond[site_idx + 1] == true &&
				is_monomer(site_idx + 1) == true &&
				bb_bond[t_occupied_by[site_idx]][site_idx] == false &&
				bb_bond[t_occupied_by[site_idx + 1]][site_idx] == false
				)
			{
				// Polymerisation is valid here
				valid_transitions[6 * c_alph_size][site_idx] = true;
			}

			// DEPOLYMERISATION
			// If the site is occupied
			// If the site+1 ahead is occupied by the same polymer chain
			// If the site+2 ahead is unoccupied within the same polymer chain
			// if the site and site+1 do share a backbone bond
			if (t_c_bond[site_idx] == true &&
				t_occupied_by[site_idx + 1] == t_occupied_by[site_idx] &&
				bb_bond[t_occupied_by[site_idx + 1]][site_idx + 1] == false &&
				bb_bond[t_occupied_by[site_idx]][site_idx] == true &&
				site_idx < length_limit - 1)
			{
				// Depolymerisation is valid here
				valid_transitions[6 * c_alph_size + 1][site_idx] = true;
			}
		}

		// TEMPLATE ACTIVATION
		// If site on the template is inactive,
		// then it is valid to activate
		if (t_active[site_idx] == false)
		{
			valid_transitions[6 * c_alph_size + 2][site_idx] = true;
		}
		// TEMPLATE DEACTIVATION
		// If site on the template is active,
		// And the site is unoccupied
		// and the site is not the first on the template and there's nothing in the site -1
		// then it is valid to deactivate
		if (site_idx != 0 &&
			t_active[site_idx] == true &&
			t_c_bond[site_idx] == false &&
			(t_activation_flag == true && t_c_bond[site_idx - 1] == false && site_idx > 0)
			)
		{
			valid_transitions[6 * c_alph_size + 3][site_idx] = true;
		}

		// Termination leading edge un-binding with deactivating template
		// If a polymer is bound to the template
		// By its last bond
		// At the leading edge of the polymer
		// and the site ahead is able to be deactivated

		if (site_idx < length_limit - 1 &&
			site_idx >0 &&
			t_c_bond[site_idx] == true &&
			is_monomer(site_idx) == false &&
			t_active[site_idx + 1] == true &&
			t_c_bond[site_idx + 1] == false &&
			t_occupied_by[site_idx + 1] != t_occupied_by[site_idx] &&
			t_occupied_by[site_idx - 1] != t_occupied_by[site_idx] &&
			t_activation_flag == true &&
			c_seq[t_occupied_by[site_idx]][site_idx + 1] == -1
			)
		{
			valid_transitions[c_seq[t_occupied_by[site_idx]][site_idx] + 6 * c_alph_size + 4][site_idx] = true;
		}

		// Termination leading edge un-binding without deactivating ahead
		// If a polymer is bound to the template
		// By its last bond
		// and the site ahead is unable to be deactivated

		if (site_idx > 0 &&
			site_idx < length_limit - 1 &&
			t_c_bond[site_idx] == true &&
			is_monomer(site_idx) == false &&
			t_occupied_by[site_idx + 1] != t_occupied_by[site_idx] &&
			t_occupied_by[site_idx - 1] != t_occupied_by[site_idx] &&
			c_seq[t_occupied_by[site_idx]][site_idx + 1] == -1 &&
			((t_c_bond[site_idx + 1] == true || t_active[site_idx + 1] == false) || t_activation_flag == false)
			)
		{
			valid_transitions[c_seq[t_occupied_by[site_idx]][site_idx] + 7 * c_alph_size + 4][site_idx] = true;
		}

		// Termination un-binding at the end of template
		// If a polymer is bound to the template
		// By its last bond
		// and is at the end

		if (site_idx == length_limit - 1 &&
			t_c_bond[site_idx] == true &&
			is_monomer(site_idx) == false &&
			t_occupied_by[site_idx - 1] != t_occupied_by[site_idx]
			)
		{
			valid_transitions[c_seq[t_occupied_by[site_idx]][site_idx] + 8 * c_alph_size + 4][site_idx] = true;
		}

		// Termination in mid of polymer un-binding with deactivation
		// If a polymer is bound to the template by its last bond
		// and isn't the leading edge of the polymer
		// and can deactivate the site ahead

		if (site_idx < length_limit - 1 &&
			t_c_bond[site_idx] == true &&
			is_monomer(site_idx) == false &&
			t_active[site_idx + 1] == true &&
			t_c_bond[site_idx + 1] == false &&
			t_occupied_by[site_idx + 1] != t_occupied_by[site_idx] &&
			(t_occupied_by[site_idx - 1] != t_occupied_by[site_idx] || site_idx == 0) &&
			t_activation_flag == true &&
			c_seq[t_occupied_by[site_idx]][site_idx + 1] != -1
			)
		{
			valid_transitions[c_seq[t_occupied_by[site_idx]][site_idx] + 9 * c_alph_size + 4][site_idx] = true;
		}

		// Termination in mid of polymer un-binding without deactivating ahead
		// If a polymer is bound to the template
		// By its last bond
		// and the site ahead is unable to be deactivated

		if (site_idx < length_limit - 1 &&
			t_c_bond[site_idx] == true &&
			is_monomer(site_idx) == false &&
			t_occupied_by[site_idx + 1] != t_occupied_by[site_idx] &&
			(t_occupied_by[site_idx - 1] != t_occupied_by[site_idx] || site_idx == 0) &&
			c_seq[t_occupied_by[site_idx]][site_idx + 1] != -1 &&
			((t_c_bond[site_idx + 1] == true || t_active[site_idx + 1] == false) || t_activation_flag == false)
			)
		{
			valid_transitions[c_seq[t_occupied_by[site_idx]][site_idx] + 10 * c_alph_size + 4][site_idx] = true;
		}

		// Polymer in middle unbind downstream and deactivate ahead
		// If the site ahead can deactivate
		// and if the polymer is not at its last bond
		// and the site is not at the leading edge of the polymer
		if (site_idx < length_limit - 1 &&
			site_idx > 0 &&
			t_c_bond[site_idx] == true &&
			is_monomer(site_idx) == false &&
			t_active[site_idx + 1] == true &&
			t_c_bond[site_idx + 1] == false &&
			t_occupied_by[site_idx - 1] == t_occupied_by[site_idx] &&
			t_occupied_by[site_idx + 1] != t_occupied_by[site_idx] &&
			c_seq[t_occupied_by[site_idx]][site_idx + 1] != -1 &&
			t_activation_flag == true
			)
		{
			valid_transitions[c_seq[t_occupied_by[site_idx]][site_idx] + 11 * c_alph_size + 4][site_idx] = true;
		}


		// Polymer in middle unbind downstream and don't deactivate ahead
		// If the site ahead cannot deactivate
		// and if the polymer is not at its last bond
		// and the site is not at the leading edge of the polymer
		if (site_idx < length_limit - 1 &&
			site_idx > 0 &&
			t_c_bond[site_idx] == true &&
			is_monomer(site_idx) == false &&
			t_occupied_by[site_idx - 1] == t_occupied_by[site_idx] &&
			t_occupied_by[site_idx + 1] != t_occupied_by[site_idx] &&
			c_seq[t_occupied_by[site_idx]][site_idx + 1] != -1 &&
			((t_c_bond[site_idx + 1] == true || t_active[site_idx + 1] == false) || t_activation_flag == false)
			)
		{
			valid_transitions[c_seq[t_occupied_by[site_idx]][site_idx] + 12 * c_alph_size + 4][site_idx] = true;
		}

		// Leading edge of polymer unbind and deactivate ahead
		// If the site ahead can deactivate
		// and if the polymer is not at its last bond
		// and the site is at the leading edge of the polymer
		if (site_idx < length_limit - 1 &&
			site_idx > 0 &&
			t_c_bond[site_idx] == true &&
			is_monomer(site_idx) == false &&
			t_active[site_idx + 1] == true &&
			t_c_bond[site_idx + 1] == false &&
			t_occupied_by[site_idx - 1] == t_occupied_by[site_idx] &&
			t_occupied_by[site_idx + 1] != t_occupied_by[site_idx] &&
			c_seq[t_occupied_by[site_idx]][site_idx + 1] == -1 &&
			t_activation_flag == true
			)
		{
			valid_transitions[c_seq[t_occupied_by[site_idx]][site_idx] + 13 * c_alph_size + 4][site_idx] = true;
		}


		// Leading edge unbind downstream and don't deactivate ahead
		// If the site ahead cannot deactivate
		// and if the polymer is not at its last bond
		// and the site is at the leading edge of the polymer
		if (site_idx < length_limit - 1 &&
			site_idx > 0 &&
			t_c_bond[site_idx] == true &&
			is_monomer(site_idx) == false &&
			t_occupied_by[site_idx - 1] == t_occupied_by[site_idx] &&
			t_occupied_by[site_idx + 1] != t_occupied_by[site_idx] &&
			c_seq[t_occupied_by[site_idx]][site_idx + 1] == -1 &&
			((t_c_bond[site_idx + 1] == true || t_active[site_idx + 1] == false) || t_activation_flag == false)
			)
		{
			valid_transitions[c_seq[t_occupied_by[site_idx]][site_idx] + 14 * c_alph_size + 4][site_idx] = true;
		}

		// Leading edge of polymer unbind at end of template
		if (site_idx == length_limit - 1 &&
			t_c_bond[site_idx] == true &&
			is_monomer(site_idx) == false &&
			t_occupied_by[site_idx - 1] == t_occupied_by[site_idx]
			)
		{
			valid_transitions[c_seq[t_occupied_by[site_idx]][site_idx] + 15 * c_alph_size + 4][site_idx] = true;
		}



		// Polymer bind template downstream
		// If the site is unoccupied (no template bond)
		// If the site-1 behind is occupied by a polymer
		// and there is a monomer in the tail to bind downstream
		// and if the site is active
		// then it is valid to bind monomers upstream
		if (site_idx > 0 &&
			t_c_bond[site_idx] == false &&
			t_active[site_idx] == true &&
			t_c_bond[site_idx - 1] == true &&
			bb_bond[t_occupied_by[site_idx - 1]][site_idx - 1] == true)
		{
			valid_transitions[c_seq[t_occupied_by[site_idx - 1]][site_idx] + 16 * c_alph_size + 4][site_idx] = true;
		}
	}
	return;
}

void PolymerSystemIO::sum_rates()
{	// Sums the transition rates for valid transitions at each site
	
	int start_idx;
	int end_idx;
	// At time == 0 we have to sum all the rates for the first time
	if (local_updates_flag && time != 0)
	{
		// If local updates are enabled we only need to update sum_rates over a small window
		// We sum over a few sites either side of the site just updated.
		start_idx = site_to_update - local_update_window;
		end_idx = site_to_update + local_update_window;
	}
	else {
		start_idx = 0;
		end_idx = length_limit;
	}

	// Prevent overflow errors
	if (start_idx < 0) { start_idx = 0; }

	if (start_idx > length_limit) { start_idx = length_limit; }

	if (end_idx < 0) { end_idx = 0; }

	if (end_idx > length_limit) { end_idx = length_limit; }

	if (end_idx < start_idx)
	{
		cerr << "\nINVALID INDEX RANGE IN sum_rates\nExiting\n\n";
		exit(EXIT_FAILURE);
		return;
	}

	long double tempsum = 0;
	// For each site between start and end index on the template
	for (size_t site_idx = start_idx; site_idx < end_idx; site_idx++)
	{
		// Set the total rate to 0 before starting the conditional sum
		total_rate_per_site[site_idx] = 0;

		for (size_t transition_idx = 0; transition_idx < num_rules; transition_idx++)
		{
			// running over possible transitions
			if (valid_transitions[transition_idx][site_idx])
			{
				total_rate_per_site[site_idx] += rates[transition_idx];
			}
		}
	}
	return;
}

int PolymerSystemIO::select_site_to_update()
{
	sum_rates();
	// First sum all the rates to construct total_rate_per_site.
   // Here we implement the MC Gillespie transition select process

	long double total_rate = 0; // set the total rate across all sites to 0

	for (size_t site_idx = 0; site_idx < length_limit; site_idx++)
	{
		total_rate += total_rate_per_site[site_idx];
	}

	// Generate random number to update time
	long double rand_time = psRandNumGen() / (long double)psRandNumGen.max();


	while (rand_time == 0)
	{
		rand_time = psRandNumGen() / (long double)psRandNumGen.max();
	}
	
	// If we want to write out trajectories, do so every traj_interval
	if (write_traj_flag == true && time + -1. * log(rand_time) / total_rate >= traj_write_counter)
	{
		write_traj_output();
		traj_write_counter += traj_write_interval;
	}

	// Increment the time by sampling from a single decaying exponential function with scale 1/total_rate
	time += -1. * log(rand_time) / total_rate;

	long double p_sum = 0;  /// cumulative probability sum set to 0

	// Generate random number to select a site to update.
	long double r = psRandNumGen() / (long double)psRandNumGen.max();

	while (r == 0)
	{
		r = psRandNumGen() / (long double)psRandNumGen.max();
	}
	//cout << "\n" << r << "\n";
	for (size_t site_idx = 0; site_idx < length_limit; site_idx++)
	{
		// Create cumulative proabbility bins and see if the random number falls within
		if (r > p_sum && r <= p_sum + (total_rate_per_site[site_idx] / total_rate))
		{
			// Given that the random number falls in this bin
			// Then choose this site to be updated.
			site_to_update = site_idx;
			return site_idx;
		}
		p_sum += total_rate_per_site[site_idx] / total_rate;
	}

	cerr << "\nAnother attempt to select update site\n";
	cerr << "Reselect update site\n";

	cerr << "Big ol error in site selection.\nExiting\n\n";
	exit(EXIT_FAILURE);
	return -1;
};

int PolymerSystemIO::select_transition()
{
	// select_site_to_update() must be called before this function
	// to ensure that site_to_update has the correct value.
	// Here we implement the MC Gillespie algo to select the transition at site_to_update

	long double total_rate = total_rate_per_site[site_to_update];

	int attempt = 0;
	long double p_sum;
	long double r;

	while (attempt < 10) {
		p_sum = 0;  /// cumulative probability sum set to 0

		// Generate random number to select a transition.
		long double r = psRandNumGen() / (long double)psRandNumGen.max();
		while (r == 0)
		{
			long double r = psRandNumGen() / (long double)psRandNumGen.max();
		}

		for (size_t transition_idx = 0; transition_idx < num_rules; transition_idx++)
		{
			if (valid_transitions[transition_idx][site_to_update] == true)
			{
				if (r > p_sum && r <= p_sum + (rates[transition_idx] / total_rate))
				{
					transition_to_update = transition_idx;
					return transition_idx;
				}
				p_sum += rates[transition_idx] / total_rate;
			}
		}
		cerr << "Reselect transition\n";

		attempt += 1;

	}
	cerr << "Big ol error in transition selection.\nExiting\n\n";
	exit(EXIT_FAILURE);
	return -1;
}

void PolymerSystemIO::apply_transition()
{
	// There are num_rules in the system, and each has consequences.
	// Given that a site and a transition by site_to_update and transition_to_update
	// then I can use apply_transition to implement the change

	//cout << "Transition " << transition_to_update << " at site " << site_to_update << ".\n";

/* Transitions where A is the alphabet size
	0:A-1		"+0 monomer", "+1 monomer",
	A:2A-1		"-0 monomer deact ahead", "-1 monomer deact ahead",
	2A:3A-1		"-0 monomer no deact", "-1 monomer no deact" ,
	3A:4A-1		"-0 monomer at end", "-1 monomer at end" ,
	4A:5A-1		"-0 upstream","-1 upstream",
	5A:6A-1		"+0 upstream","+1 upstream",
	6A & 6A+1	"Polymerise", "Depolymerise",
	6A+2 & 6A+3	"Activate template","Deactivate template",
	6A+4:7A+3	"Terminate lead 0 deact ahead","Terminate lead 1 deact ahead",
	7A+4:7A+3	"Terminate lead 0 no deact","Terminate lead 1  no deact",
	8A+4:9A+3	"Terminate lead 0 at end", "Terminate lead 0 at end",
	9A+4:10A+3	"Terminate mid 0 deact ahead","Terminate mid 1 deact ahead",
	10A+4:11A+3	"Terminate mid 0 no deact","Terminate mid 1 no deact",
	11A+4:12A+3	"-0 poly deact ahead", "-1 poly deact ahead",
	12A+4:13A+3	"-0 poly no deact", "-1 poly no deact" ,
	13A+4:14A+3	"-0 leading poly deact ahead", "-1 leading poly deact ahead",
	14A+4:15A+3	"-0 leading poly no deact", "-1 leading poly no deact" ,
	15A+4:16A+3	"-0 leading poly at end", "-1 leading poly at end",
	16A+4:17A+3	"+0 downstream", "+1 downstream"


*/


/*
BIND MONOMER DOWNSTREAM
We add a monomer at site_to_update in a free polymer row
*/
	if (transition_to_update < c_alph_size)
	{
		int pol_temp = 0;
		while (pol_temp < pols_limit)
		{
			// In the first available polymer sequence row add this monomer
			if (polymer_present[pol_temp] == false)
			{
				// The row is flagged as occupied by a sequence
				polymer_present[pol_temp] = true;
				// The monomer type at this site is given by the transition index
				c_seq[pol_temp][site_to_update] = (int)transition_to_update;
				// A template copy bond has formed
				t_c_bond[site_to_update] = true;
				// The template is occupied by a polymer in the row pol_temp
				t_occupied_by[site_to_update] = pol_temp;

				// Monomers activative the sites ahead.
				if (site_to_update < length_limit - 1 && t_activation_flag)
				{
					t_active[site_to_update + 1] = true;
				}
				return;

			}
			pol_temp += 1;

		}
		cerr << "Error at time = " << time;
		print_state();

		cerr << "\nNo space in polymer sequence array to add another monomer\nExiting\n";
		exit(EXIT_FAILURE);
	}

	/*
	UNBIND MONOMER DOWNSTREAM and deactivate template ahead
	We remove a monomer at site_to_update
	*/
	if (transition_to_update >= c_alph_size && transition_to_update < 2 * c_alph_size)
	{
		// Remove the sequence info for this monomer type
		c_seq[t_occupied_by[site_to_update]][site_to_update] = -1;
		// No polymer is present at the site any more
		polymer_present[t_occupied_by[site_to_update]] = false;
		// Template polymer bond is broken
		t_c_bond[site_to_update] = false;
		// The site is unoccupied
		t_occupied_by[site_to_update] = -1;
		// The site ahead is deactivated
		t_active[site_to_update + 1] = false;
		return;
	}

	/*
	UNBIND MONOMER DOWNSTREAM and don't deactivate template ahead
	We remove a monomer at site_to_update
	*/
	if (transition_to_update >= 2 * c_alph_size && transition_to_update < 3 * c_alph_size)
	{
		// Remove the sequence info for this monomer type
		c_seq[t_occupied_by[site_to_update]][site_to_update] = -1;
		// No polymer is present at the site any more
		polymer_present[t_occupied_by[site_to_update]] = false;
		// Template polymer bond is broken
		t_c_bond[site_to_update] = false;
		// The site is unoccupied
		t_occupied_by[site_to_update] = -1;
		return;
	}

	/*
	UNBIND MONOMER at the end of template
	We remove a monomer at site_to_update
	*/
	if (transition_to_update >= 3 * c_alph_size && transition_to_update < 4 * c_alph_size)
	{
		// Remove the sequence info for this monomer type
		c_seq[t_occupied_by[site_to_update]][site_to_update] = -1;
		// No polymer is present at the site any more
		polymer_present[t_occupied_by[site_to_update]] = false;
		// Template polymer bond is broken
		t_c_bond[site_to_update] = false;
		// The site is unoccupied
		t_occupied_by[site_to_update] = -1;
		return;
	}

	/*
	UNBIND MONOMER UPSTREAM
	We break a tcbond at site_to update
	*/
	if (transition_to_update >= 4 * c_alph_size && transition_to_update < 5 * c_alph_size)
	{
		t_c_bond[site_to_update] = false;
		t_occupied_by[site_to_update] = -1;
		return;
	}

	/*
	BIND MONOMER UPSTREAM
	We make a tcbond at site_to update
	*/
	if (transition_to_update >= 5 * c_alph_size && transition_to_update < 6 * c_alph_size)
	{
		t_c_bond[site_to_update] = true;
		t_occupied_by[site_to_update] = t_occupied_by[site_to_update + 1];
		return;
	}

	/*
	Polymerise
	We make a bb_bond at site_to update
	connecting the monomers at site and site+1
	The monomers from site+1 are copied in to the sequence of site

	*/
	if (transition_to_update == 6 * c_alph_size)
	{
		int old = t_occupied_by[site_to_update + 1];

		int current = t_occupied_by[site_to_update];

		// Make the backbone bond
		bb_bond[current][site_to_update] = true;

		// Copy info from site+1 polymer to site polymer
		size_t temp_index = site_to_update + 1;

		// Run until we reach the end of the old polymer
		while (c_seq[old][temp_index] > -1 && temp_index < length_limit)
		{
			// Copy sequence information
			// Set new
			c_seq[current][temp_index] = c_seq[old][temp_index];
			// Erase old
			c_seq[old][temp_index] = -1;
			// Set new
			bb_bond[current][temp_index] = bb_bond[old][temp_index];
			// Erase old
			bb_bond[old][temp_index] = false;
			// Set new
			if (t_occupied_by[temp_index] == old) {

				t_occupied_by[temp_index] = current;
			}
			temp_index += 1;
		}
		// Flag that the row of the polymer seq matrix is empty
		polymer_present[old] = false;

		return;
	}


	/*
	Depolymerise
	We break a bb_bond at site_to update
	then move all information from the downstream half of the broken polymer
	to a free row.
	*/
	if (transition_to_update == 6 * c_alph_size + 1)
	{

		int old = t_occupied_by[site_to_update];
		int current = -1;
		size_t pol_temp = 0;
		while (pol_temp < pols_limit)
		{
			// Current is the index of the first available polymer sequence row
			if (polymer_present[pol_temp] == false)
			{
				current = pol_temp;
				break;
			}
			pol_temp += 1;
		}
		// Flag that the row of the new polymer seq matrix is occupied
		polymer_present[current] = true;

		// Break the backbone bond
		bb_bond[old][site_to_update] = false;

		// Copy info from old(site) polymer to current(new) polymer (site+1)
		size_t temp_index = site_to_update + 1;

		// Run until we reach the end of the old polymer
		while (c_seq[old][temp_index] > -1 && temp_index < length_limit)
		{
			// Copy sequence information
			c_seq[current][temp_index] = c_seq[old][temp_index];
			c_seq[old][temp_index] = -1;

			bb_bond[current][temp_index] = bb_bond[old][temp_index];
			bb_bond[old][temp_index] = false;
			if (t_occupied_by[temp_index] == old) {

				t_occupied_by[temp_index] = current;
			}

			temp_index += 1;
		}


		return;
	}

	// TEMPLATE ACTIVATION
	if (transition_to_update == 6 * c_alph_size + 2)
	{
		t_active[site_to_update] = true;
		return;
	}

	// TEMPLATE DEACTIVATION

	if (transition_to_update == 6 * c_alph_size + 3)
	{
		t_active[site_to_update] = false;
		return;
	}

	/*
	Terminate polymer with deactivating ahead

	*/
	if ((transition_to_update >= 6 * c_alph_size + 4 &&
		transition_to_update < 7 * c_alph_size + 4) ||
		(transition_to_update >= 9 * c_alph_size + 4 &&
			transition_to_update < 10 * c_alph_size + 4))
	{

		// Write terminated polymer to file
		if (write_output_flag == true)
		{
			write_copy_output();  // Save the sequence information of the polymer to file
		}
		
		
		int temp_occupied_ptr = t_occupied_by[site_to_update];
		// Remove the sequence info for the terminated polymer
		for (size_t i = 0; i < length_limit; i++)
		{
			c_seq[temp_occupied_ptr][i] = -1;
			bb_bond[temp_occupied_ptr][i] = false;
		}
		// Deactivate the template ahead
		t_active[site_to_update + 1] = false;
		// No polymer is present at the site any more
		polymer_present[temp_occupied_ptr] = false;
		// Template polymer bond is broken
		t_c_bond[site_to_update] = false;
		// The site is unoccupied
		t_occupied_by[site_to_update] = -1;

		return;
	}

	/*
		Terminate leading polymer without deactivating ahead
		*/

	if ((transition_to_update >= 7 * c_alph_size + 4 &&
		transition_to_update < 8 * c_alph_size + 4) ||
		(transition_to_update >= 10 * c_alph_size + 4 &&
			transition_to_update < 11 * c_alph_size + 4))
	{

		// Write terminated polymer to file
		if (write_output_flag == true)
		{
			write_copy_output();  // Save the sequence information of the polymer to file
		}


		int temp_occupied_ptr = t_occupied_by[site_to_update];
		// Remove the sequence info for the terminated polymer
		for (size_t i = 0; i < length_limit; i++)
		{
			c_seq[temp_occupied_ptr][i] = -1;
			bb_bond[temp_occupied_ptr][i] = false;
		}

		// No polymer is present at the site any more
		polymer_present[temp_occupied_ptr] = false;
		// Template polymer bond is broken
		t_c_bond[site_to_update] = false;
		// The site is unoccupied
		t_occupied_by[site_to_update] = -1;

		return;
	}
	/*
		Terminate polymer at end
		*/

	if (transition_to_update >= 8 * c_alph_size + 4 && transition_to_update < 9 * c_alph_size + 4)
	{
		// Write terminated polymer to file
		if (write_output_flag == true)
		{
			write_copy_output();  // Save the sequence information of the polymer to file
		}

		int temp_occupied_ptr = t_occupied_by[site_to_update];
		// Remove the sequence info for the terminated polymer
		for (size_t i = 0; i < length_limit; i++)
		{
			c_seq[temp_occupied_ptr][i] = -1;
			bb_bond[temp_occupied_ptr][i] = false;
		}

		// No polymer is present at the site any more
		polymer_present[temp_occupied_ptr] = false;
		// Template polymer bond is broken
		t_c_bond[site_to_update] = false;
		// The site is unoccupied
		t_occupied_by[site_to_update] = -1;

		return;
	}

	// Polymer unbind and deactivate ahead
	if ((transition_to_update >= 11 * c_alph_size + 4
		&& transition_to_update < 12 * c_alph_size + 4) ||
		(transition_to_update >= 13 * c_alph_size + 4
			&& transition_to_update < 14 * c_alph_size + 4))
	{
		t_c_bond[site_to_update] = false;
		t_occupied_by[site_to_update] = -1;
		t_active[site_to_update + 1] = false;
		return;
	}

	// Polymer unbind and don't deactivate ahead or unbind at the end
	if ((transition_to_update >= 12 * c_alph_size + 4
		&& transition_to_update < 13 * c_alph_size + 4) ||
		(transition_to_update >= 14 * c_alph_size + 4
			&& transition_to_update < 15 * c_alph_size + 4) ||
		(transition_to_update >= 15 * c_alph_size + 4
			&& transition_to_update < 16 * c_alph_size + 4))
	{
		t_c_bond[site_to_update] = false;
		t_occupied_by[site_to_update] = -1;
		return;
	}

	// Polymer bind downstream
	if (transition_to_update >= 16 * c_alph_size + 4
		&& transition_to_update < 17 * c_alph_size + 4)
	{
		t_c_bond[site_to_update] = true;
		t_occupied_by[site_to_update] = t_occupied_by[site_to_update - 1];
		if (site_to_update < length_limit - 1 && t_activation_flag) {

			t_active[site_to_update + 1] = true;
		}
		return;
	}

	cerr << "NO TRANSITION IMPLEMENTED in apply_transition()";
	exit(EXIT_FAILURE);
}

bool PolymerSystemIO::is_monomer(int idx)
{
	// Returns true if a lone monomer is bound to the template at the site idx
		// Returns false if the item at idx is empty or part of a polymer

	if (idx >= length_limit || idx < 0)
	{
		cerr << "INVALID INDEX PASSED TO is_monomer";
		exit(EXIT_FAILURE);
	}
	if (idx == 0)
	{
		if (t_occupied_by[idx] == -1)
		{
			// If nothing bound to template
			return false;
		}
		if (t_c_bond[idx] == true &&
			bb_bond[t_occupied_by[idx]][idx] == true)
		{
			// If something is bound to the template
			// and it has a backbone bond to another polymer
			return false;
		}
		else {
			// If idx is not empty or part of a polymer then it is a monomer
			return true;
		}
	}

	if (t_occupied_by[idx] == -1)
	{
		return false;
	}
	if (t_c_bond[idx] == true &&
		(bb_bond[t_occupied_by[idx]][idx] == true || bb_bond[t_occupied_by[idx]][idx - 1] == true))
	{
		return false;
	}
	else {
		return true;
	}
}

void PolymerSystemIO::setup() {

	read_input();   // Read the input file
	init_length();  // Initialise the working length of the arrays
	init_output();  // Initialise streams input and output
	set_seed(seed); // Set the random number seed 
	set_rates();    // Calculate the rates based upon parameters
	if (t_activation_flag == false)
	{
		// If we are neglecting template activation, then set the initial condition for template activation to be all active
		fill(&t_active[0], &t_active[length_limit_max], true); // Initialise with fully active template.

	}
}

void PolymerSystemIO::read_input()
{
	ifstream fin;
	string line;

	string delimiter = ",";
	// Open an existing file
	fin.open(inputfilename);
	cout << inputfilename;
	if (fin.fail()) {
		// file could not be opened

		cerr << "Input file could not be opened. \n Exiting... \n\n";
		exit(EXIT_FAILURE);
	}
	int temp_counter = 0;
	while (!fin.eof()) {
		fin >> line;
		cout << line << "\n";

		if (temp_counter == 1)
		{
			line += ", "; // Add this to end of line to catch the last column
			int tt_contr = 0;
			size_t pos = 0;
			string token;
			while ((pos = line.find(delimiter)) != std::string::npos) {
				/*
				token = line.substr(0, pos);

				if (tt_contr == 0 ) { seed = stoi(token); }
				if (tt_contr == 1 ) { length_limit = stoi(token); }
				if (tt_contr == 2 ) { time_limit = stoi(token); }
				if (tt_contr == 3 ) { print_flag = (stoi(token) > 0); }
				if (tt_contr == 4 ) { print_valid_trans_flag = (stoi(token) > 0); }
				if (tt_contr == 5 ) { write_output_flag = (stoi(token) > 0); }
				if (tt_contr == 6 ) { write_traj_flag = (stoi(token) > 0); }
				if (tt_contr == 7 ) { traj_write_interval = stoi(token); }
				if (tt_contr == 8 ) { outputfilename = token; }
				if (tt_contr == 9 ) { outputtrajfilename = token; }
				if (tt_contr == 10) { breakable_backbone_flag = (stoi(token) > 0); }
				if (tt_contr == 11) { t_activation_flag = (stoi(token) > 0); }
				if (tt_contr == 12) { t_activation_flag = (stoi(token) > 0); }
				if (tt_contr == 13) { local_updates_flag = (stoi(token) > 0); }
				if (tt_contr == 14) { use_test_rates_flag = (stoi(token) > 0); }
				if (tt_contr == 15) { k0 = stold(token); }
				if (tt_contr == 16) { k = stold(token); }
				if (tt_contr == 17) { kact = stold(token); }
				if (tt_contr == 18) { ConcEff = stold(token); }
				if (tt_contr == 19) { Gact = stold(token); }
				if (tt_contr == 20) { Gbb = stold(token); }
				if (tt_contr == 21) { Conc[0] = stold(token); }
				if (tt_contr == 22) { Conc[1] = stold(token); }
				if (tt_contr == 23) { Gx[0] = stold(token); }
				if (tt_contr == 24) { Gx[1] = stold(token); }
				if (tt_contr == 25) { Ggen = stold(token); }
				if (tt_contr == 26) { Gend = stold(token);
				break;	*/

				token = line.substr(0, pos);

				if (tt_contr == 0) { seed = stoi(token); }
				if (tt_contr == 1) { length_limit = stoi(token); }
				if (tt_contr == 2) { time_limit = stoi(token); }
				if (tt_contr == 3) { max_pols = stoi(token); }
				if (tt_contr == 4) { k0 = stold(token); }
				if (tt_contr == 5) { k = stold(token); }
				if (tt_contr == 6) { kact = stold(token); }
				if (tt_contr == 7) { ConcEff = stold(token); }
				if (tt_contr == 8) { Gact = stold(token); }
				if (tt_contr == 9) { Gbb = stold(token); }
				if (tt_contr == 10) { Conc[0] = stold(token); }
				if (tt_contr == 11) { Conc[1] = stold(token); }
				if (tt_contr == 12) { Gx[0] = stold(token); }
				if (tt_contr == 13) { Gx[1] = stold(token); }
				if (tt_contr == 14) { Ggen = stold(token); }
				if (tt_contr == 15) { Gend = stold(token); }
				if (tt_contr == 16) { outputfilename = token; }
				if (tt_contr == 17) { outputtrajfilename = token; }
				if (tt_contr == 18) { print_flag = (stoi(token) > 0); }
				if (tt_contr == 19) { print_valid_trans_flag = (stoi(token) > 0); }
				if (tt_contr == 20) { write_output_flag = (stoi(token) > 0); }
				if (tt_contr == 21) { write_traj_flag = (stoi(token) > 0); }
				if (tt_contr == 22) { traj_write_interval = stoi(token); }
				if (tt_contr == 23) { breakable_backbone_flag = (stoi(token) > 0); }
				if (tt_contr == 24) { off_rate_discrim_flag = (stoi(token) > 0); }
				if (tt_contr == 25) { t_activation_flag = (stoi(token) > 0); }
				if (tt_contr == 26) { local_updates_flag = (stoi(token) > 0); }
				if (tt_contr == 27) {
					use_test_rates_flag = (stoi(token) > 0);
					break;
				}

				//std::cout << token << std::endl;
				line.erase(0, pos + delimiter.length());
				tt_contr++;

			}
			line = "";
		}

		++temp_counter;

	}



}

void PolymerSystemIO::init_length()
{

	if (length_limit > length_limit_max || length_limit < 0) {
		cerr << "\n\n Invalid template length. Template length must be an integer between 0 and " + to_string(length_limit_max) + ". \n Exiting...\n";
		exit(EXIT_FAILURE);
		return;
	}
	else
	{
		pols_limit = length_limit;
	}
}

void PolymerSystemIO::init_output()
{
	// Initialise the output streams to write to outputfile and trajectory file

	if (write_output_flag == true)
	{
		outputfile.open(outputfilename, ios::out | ios::app);

		if (outputfile.fail()) {
			// file could not be opened
			cerr << "Output file could not be opened. \n Exiting... \n\n";
			exit(EXIT_FAILURE);
		}
		if (outputfile.is_open())
		{
			for (size_t i = 1; i < length_limit; i++)
			{
				outputfile << i << ",";
			}
			outputfile << length_limit << endl;
		}
	}

	if (write_traj_flag == true)
	{
		outputtrajfile.open(outputtrajfilename, ios::out | ios::app);
		if (outputtrajfile.fail()) {
			// file could not be opened
			cerr << "Output trajectory file could not be opened. \n Exiting... \n\n";
			exit(EXIT_FAILURE);
		}
		if (outputtrajfile.is_open())
		{
			outputtrajfile << "Time (s),";
			outputtrajfile << "Template length,";
			outputtrajfile << "# TC bonds,";
			outputtrajfile << "# Active template sites (s),";
			outputtrajfile << "# total,";
			outputtrajfile << "# Polymers,";
			outputtrajfile << "# other way of counting Polymers,";
			outputtrajfile << "# Monomers,";
			outputtrajfile << "<pollength/siteidx> avg over polymers,";
			outputtrajfile << "<attached by/polymer length> avg over polymers,";
			outputtrajfile << "<spacing> avg over polymers \n";

		}
	}

}

void PolymerSystemIO::set_seed(int s)
{
	// Set seed to input value s (record of the seed is essential!)
	seed = s;
	// Update the RNG seed
	psRandNumGen.seed(seed);
	std::cout << std::setprecision(6);

	cout << "Resolution is " << 1 / (long double)psRandNumGen.max() << endl << endl;
}

void PolymerSystemIO::set_rates()
{
	// Initialises the values of the rate matrix
	

	// OFF RATE DISCRIMINATION
	// The rate of monomer unbinding is monomer dependent, and the rate of binding is independent
	if (off_rate_discrim_flag == true)
	{
		// Set the rates for all reactions
		for (size_t i = 0; i < c_alph_size; i++)
		{
			// Provisional rates can be set at the start if we don't care about the template sequence (Binary system only)
			rates[i] = k0 * Conc[i];   // Rate of adding copy monomer type i.
			rates[i + c_alph_size] = k0 * exp(Gx[i] + Ggen - Gact);   // Rate of removing copy monomer type i and deactivating ahead site.
			rates[i + 2 * c_alph_size] = k0 * exp(Gx[i] + Ggen);   // Rate of removing copy monomer type i without deactivating template.
			rates[i + 3 * c_alph_size] = k0 * exp(Gx[i] + Ggen + Gend);   // Rate of removing copy monomer type i at end of template.
			rates[i + 4 * c_alph_size] = k0 * exp(Gx[i]);   // Rate of removing monomer i upstream. 
			rates[i + 5 * c_alph_size] = k0 * ConcEff;   // Rate of adding monomer i upstream. 
			rates[i + 6 * c_alph_size + 4] = k0 * exp(Gx[i] + Ggen - Gact);   // Rate of terminating monomer i with deactivating forward site. 
			rates[i + 7 * c_alph_size + 4] = k0 * exp(Gx[i] + Ggen);   // Rate of terminating monomer i without deactivating +1 site. 
			rates[i + 8 * c_alph_size + 4] = k0 * exp(Gx[i] + Ggen + Gend);   // Rate of terminating monomer i without deactivating +1 site at end. 
			rates[i + 9 * c_alph_size + 4] = k0 * exp(Gx[i] - Gact);   // Terminating monomer i deact ahead in middle of polymer
			rates[i + 10 * c_alph_size + 4] = k0 * exp(Gx[i]);   // Terminating monomer i no deactivation ahead in mid of polymer.
			rates[i + 11 * c_alph_size + 4] = k0 * exp(Gx[i] - Gact);   // downstream unbind monomer i deact ahead in mid polymer
			rates[i + 12 * c_alph_size + 4] = k0 * exp(Gx[i]);  // downstream unbind monomer i no deact ahead in mid polymer
			rates[i + 13 * c_alph_size + 4] = k0 * exp(Gx[i] + Ggen - Gact);   //Leading downstream unbind monomer i deact ahead at leading edge
			rates[i + 14 * c_alph_size + 4] = k0 * exp(Gx[i] + Ggen);  // Leading downstream unbind monomer i no deact ahead at leading edge
			rates[i + 15 * c_alph_size + 4] = k0 * exp(Gx[i] + Ggen + Gend);  // Leading downstream unbind monomer i no deact ahead at leading edge at end
			rates[i + 16 * c_alph_size + 4] = k0 * ConcEff;  // Bind polymer downstream

		}

		rates[6 * c_alph_size] = k;                       // Rate of polymerisation
		rates[6 * c_alph_size + 1] = k * exp(Gbb - Ggen);   // Rate of depolymerisation

		rates[6 * c_alph_size + 2] = kact * exp(Gact);   // Rate of template activation
		rates[6 * c_alph_size + 3] = kact;   // Rate of template deactivation

	}

	// ON RATE DISCRIMINATION
	// The rate of monomer binding is monomer dependent, and the rate of UNbinding is independent
	else if (off_rate_discrim_flag == false)
	{
		// Set the rates for all reactions
		for (size_t i = 0; i < c_alph_size; i++)
		{
			// Provisional rates can be set at the start if we don't care about the template sequence (Binary system only)
			rates[i] = k0 * Conc[i] * exp(-Gx[i]);   // Rate of adding copy monomer type i.
			rates[i + c_alph_size] = k0 * exp(Ggen - Gact);   // Rate of removing copy monomer type i and deactivating ahead site.
			rates[i + 2 * c_alph_size] = k0 * exp(Ggen);   // Rate of removing copy monomer type i without deactivating template.
			rates[i + 3 * c_alph_size] = k0 * exp(Ggen + Gend);   // Rate of removing copy monomer type i at end of template.
			rates[i + 4 * c_alph_size] = k0;   // Rate of removing monomer i upstream. 
			rates[i + 5 * c_alph_size] = k0 * ConcEff * exp(-Gx[i]);   // Rate of adding monomer i upstream. 
			rates[i + 6 * c_alph_size + 4] = k0 * exp(Ggen - Gact);   // Rate of terminating monomer i with deactivating forward site. 
			rates[i + 7 * c_alph_size + 4] = k0 * exp(Ggen);   // Rate of terminating monomer i without deactivating +1 site. 
			rates[i + 8 * c_alph_size + 4] = k0 * exp(Ggen + Gend);   // Rate of terminating monomer i without deactivating +1 site at end. 
			rates[i + 9 * c_alph_size + 4] = k0 * exp(-Gact);   // Terminating monomer i deact ahead in middle of polymer
			rates[i + 10 * c_alph_size + 4] = k0;   // Terminating monomer i no deactivation ahead in mid of polymer.
			rates[i + 11 * c_alph_size + 4] = k0 * exp(-Gact);   // downstream unbind monomer i deact ahead in mid polymer
			rates[i + 12 * c_alph_size + 4] = k0;  // downstream unbind monomer i no deact ahead in mid polymer
			rates[i + 13 * c_alph_size + 4] = k0 * exp(Ggen - Gact);   //Leading downstream unbind monomer i deact ahead at leading edge
			rates[i + 14 * c_alph_size + 4] = k0 * exp(Ggen);  // Leading downstream unbind monomer i no deact ahead at leading edge
			rates[i + 15 * c_alph_size + 4] = k0 * exp(Ggen + Gend);  // Leading downstream unbind monomer i no deact ahead at leading edge at end
			rates[i + 16 * c_alph_size + 4] = k0 * ConcEff * exp(-Gx[i]);  // Bind polymer downstream

		}

		rates[6 * c_alph_size] = k;                       // Rate of polymerisation
		rates[6 * c_alph_size + 1] = k * exp(Gbb - Ggen);   // Rate of depolymerisation

		rates[6 * c_alph_size + 2] = kact * exp(Gact);   // Rate of template activation
		rates[6 * c_alph_size + 3] = kact;   // Rate of template deactivation

	}


	// If we are excluding template activation from the simulation then set these rates to 0
	if (t_activation_flag == false) {
		rates[6 * c_alph_size + 2] = 0;   // Rate of template activation
		rates[6 * c_alph_size + 3] = 0;   // Rate of template deactivation

		for (size_t i = 0; i < c_alph_size; i++)
		{
			rates[i + 6 * c_alph_size + 4] = 0;  // Rate of terminating i with deactivating template
			rates[i + c_alph_size] = 0;   // Rate of removing copy monomer type i with deactivating template.
			rates[i + 9 * c_alph_size + 4] = 0;
			rates[i + 11 * c_alph_size + 4] = 0;
			rates[i + 13 * c_alph_size + 4] = 0;

		}
	}

	// Test rates can be used for debugging and only allowing certain transitions
	if (use_test_rates_flag == true)
	{
		// Set the rates for all reactions
		for (size_t i = 0; i < c_alph_size; i++)
		{
			// Provisional rates can be set at the start if we don't care about the template sequence (Binary system only)
			rates[i] = 1;   // Rate of adding copy monomer type i.
			rates[i + c_alph_size] = 1;   // Rate of removing copy monomer type i and deactivating ahead site.
			rates[i + 2 * c_alph_size] = 1;   // Rate of removing copy monomer type i without deactivating template.
			rates[i + 3 * c_alph_size] = 1;   // Rate of removing copy monomer type i at end of template.
			rates[i + 4 * c_alph_size] = 1;   // Rate of removing monomer i upstream. 
			rates[i + 5 * c_alph_size] = 1;   // Rate of adding monomer i upstream. 
			rates[i + 6 * c_alph_size + 4] = 1;   // Rate of terminating monomer i with deactivating forward site. 
			rates[i + 7 * c_alph_size + 4] = 1;   // Rate of terminating monomer i without deactivating +1 site. 
			rates[i + 8 * c_alph_size + 4] = 1;   // Rate of terminating monomer i without deactivating +1 site at end. 
			rates[i + 9 * c_alph_size + 4] = 1;   // Terminating monomer i deact ahead in middle of polymer
			rates[i + 10 * c_alph_size + 4] = 1;   // Terminating monomer i no deactivation ahead in mid of polymer.
			rates[i + 11 * c_alph_size + 4] = 1;   // downstream unbind monomer i deact ahead in mid polymer
			rates[i + 12 * c_alph_size + 4] = 1;  // downstream unbind monomer i no deact ahead in mid polymer
			rates[i + 13 * c_alph_size + 4] = 1;   //Leading downstream unbind monomer i deact ahead at leading edge
			rates[i + 14 * c_alph_size + 4] = 1;  // Leading downstream unbind monomer i no deact ahead at leading edge
			rates[i + 15 * c_alph_size + 4] = 1;  // Leading downstream unbind monomer i no deact ahead at leading edge at end
			rates[i + 16 * c_alph_size + 4] = 1;  // Bind polymer downstream

		}

		rates[6 * c_alph_size] = 1;                       // Rate of polymerisation
		rates[6 * c_alph_size + 1] = 1;   // Rate of depolymerisation

		rates[6 * c_alph_size + 2] = 1;   // Rate of template activation
		rates[6 * c_alph_size + 3] = 1;   // Rate of template deactivation

		if (t_activation_flag == false) {
			rates[6 * c_alph_size + 2] = 0;   // Rate of template activation
			rates[6 * c_alph_size + 3] = 0;   // Rate of template deactivation

			for (size_t i = 0; i < c_alph_size; i++)
			{
				rates[i + 6 * c_alph_size + 4] = 0;  // Rate of terminating i with deactivating template
				rates[i + c_alph_size] = 0;   // Rate of removing copy monomer type i with deactivating template.
				rates[i + 9 * c_alph_size + 4] = 0;
				rates[i + 11 * c_alph_size + 4] = 0;
				rates[i + 13 * c_alph_size + 4] = 0;

			}
		}

	}
	cout << "Rate values: \n\n";
	for (size_t i = 0; i < num_rules; i++)
	{
		cout << i << ")";
		if (i < 10)
		{
			cout << " ";
		}
		cout << "  " << rule_labels[i] << " : " << rates[i] << "\n";

	}
	cout << endl;

}

void PolymerSystemIO::generate_initial_condition(int start_idx, int end_idx, int attached_by, int temp_row)
{
	// Generates a polymer that starts at start_idx, ends at end_idx, is connected by attached_by bonds, in row temp_row of cseq matrix.

	if (end_idx > length_limit || end_idx < 0 || start_idx > length_limit || start_idx < 0 || start_idx > end_idx || end_idx - start_idx < attached_by - 1)
	{
		cerr << "INVALID INITIAL CONDITIONS GIVEN";
		exit(EXIT_FAILURE);
		return;

	}
	// If we don't care about template activation then we need to set the template to open
	if (t_activation_flag == false)
	{
		fill(&t_active[0], &t_active[length_limit], true); // Initialise with fully active template.

	}
	t_active[0] = true;

	// Flag that a polymer info is in this row
	polymer_present[temp_row] = true;

	// From the start to the end
	for (int i = start_idx; i < end_idx; i++)
	{
		// Set the sequence
		c_seq[temp_row][i] = 0;
		// Polymerise
		bb_bond[temp_row][i] = true;
		// Attached by the last (attached_by) monomers should be set
		if (i > end_idx - attached_by)
		{
			t_c_bond[i] = true;
			t_occupied_by[i] = temp_row;
			t_active[i] = true;
		}
	}
	// Last few are tcbonded
	c_seq[temp_row][end_idx] = 0;
	t_c_bond[end_idx] = true;
	t_occupied_by[end_idx] = temp_row;
	t_active[end_idx] = true;

	if (end_idx < length_limit - 1 && t_activation_flag == true)
	{
		t_active[end_idx + 1] = true;
	}
}

void PolymerSystemIO::write_traj_output()
{
	print_state_vertically();

	// Calculates observables and then writes to traj file
	if (outputtrajfile.is_open())
	{
		double num_tcbonds = 0;
		double num_tactive = 0;
		double num_things = 0;
		double num_polymers = 0;
		double other_num_polymers = 0;
		double num_monomers = 0;

		double current_occupied_by = t_occupied_by[0];
		double idx_old = 0;
		double spacing_counter = 0;
		double spacing_total = 0;

		double polymer_completeness = 0;
		double polymer_bondedness = 0;
		for (size_t i = 0; i < length_limit; i++)
		{

			if (t_occupied_by[i] != current_occupied_by)
			{
				if (current_occupied_by == -1)
				{
					// Have crossed a space to a new polymer

					spacing_total += i - idx_old;
					spacing_counter += 1;
					current_occupied_by = t_occupied_by[i];
				}
				else if (t_occupied_by[i] == -1)
				{
					// Have reached the end of a polymer
					// Don't count the distance here
					current_occupied_by = t_occupied_by[i];
					idx_old = i;
				}
				else if (current_occupied_by != -1 && t_occupied_by[i] != -1)
				{
					// Two polymers are next to each other with distance = 0
					// Increase the counter
					spacing_counter += 1;
					idx_old = i;
					current_occupied_by = t_occupied_by[i];
				}

			}

			if (i == length_limit - 1 && idx_old != i && current_occupied_by == -1)
			{
				// Don't forget to count the space on the end!
				spacing_total += i - idx_old;
				spacing_counter += 1;
				current_occupied_by = t_occupied_by[i];
			}


			num_tcbonds += t_c_bond[i];
			num_tactive += t_active[i];
			num_things += polymer_present[i];
			num_monomers += is_monomer(i);
			if (polymer_present[i])
			{ // If there's a polymer present check if it's a polymer
				for (size_t j = 0; j < length_limit; j++)
				{
					if (t_occupied_by[j] == i)
					{
						// find the site that this polymer is attached to and check if it's NOT a monomer
						// Then it must be a polymer
						if (is_monomer(j) == false)
						{
							num_polymers += 1;

							// Now to calculate pol length/furthest site idx
							// and calculate the % of the polymer which is stuck down
							// We need to find polymer start and end idx, and tc_bond start and end idx
							/*  0		  s           b     e
										  0-0-0-0-0-0-0-0-0-0
													  | | | |
							*/
							double pol_start = -1; //s
							double pol_end = -1;    //e
							double bond_start = -1; //b
							double bond_end = -1;

							for (size_t idx = 0; idx < length_limit; idx++)
							{
								if (c_seq[i][idx] != -1 && pol_start == -1)
								{
									pol_start = idx;
								}
								if (c_seq[i][idx] == -1 && pol_start != -1 && pol_end == -1)
								{
									pol_end = idx - 1;
								}
								if (t_occupied_by[idx] == i && bond_start == -1)
								{
									bond_start = idx;
								}
								if (t_occupied_by[idx] != i && bond_start != -1 && bond_end == -1)
								{
									bond_end = idx - 1;
								}

								if (idx == length_limit - 1)
								{
									if (pol_end == -1)
									{
										pol_end = length_limit - 1;
									}
									if (bond_end == -1)
									{
										bond_end = length_limit - 1;
									}
									if (bond_start == -1)
									{
										bond_start = length_limit - 1;
									}
								}

								if (pol_start != -1 && pol_end != -1 && bond_start != -1 && bond_end != -1)
								{
									/*
									cout << "\n\n";
									cout << "Polymer        " << j << "\n";
									cout << "Pol starts at  " << pol_start << "\n";
									cout << "Pol ends at    " << pol_end << "\n";
									cout << "Bonds start at " << bond_start << "\n";
									cout << "Bonds end at   " << bond_end << "\n";
									cout << "Completeness   " << (pol_end - pol_start + 1) / (pol_end + 1) << " = " << (pol_end - pol_start + 1) << "/" << (pol_end + 1) << "\n";
									cout << "Bond fraction  " << (bond_end + 1 - bond_start) / (pol_end - pol_start + 1) << " = " << (bond_end + 1 - bond_start) << "/" << (pol_end - pol_start + 1) << "\n";
									*/
									break;
								}
							}
							polymer_completeness += (pol_end - pol_start + 1) / (pol_end + 1); // Polymer length as a fraction of the maximum length
							polymer_bondedness += (bond_end + 1 - bond_start) / (pol_end - pol_start + 1); // Polymer bonds as a fraction of total polymer length


							break;
						}
					}
				}

			}
		}
		other_num_polymers = num_things - num_monomers;


		outputtrajfile << time << ",";
		outputtrajfile << length_limit << ",";
		outputtrajfile << num_tcbonds << ",";
		outputtrajfile << num_tactive << ",";
		outputtrajfile << num_things << ",";
		outputtrajfile << num_polymers << ",";
		outputtrajfile << other_num_polymers << ",";
		outputtrajfile << num_monomers << ",";
		outputtrajfile << polymer_completeness / num_polymers << ",";
		outputtrajfile << polymer_bondedness / num_polymers << ",";
		outputtrajfile << (spacing_total / spacing_counter) << "\n";

		/*
		cout << "\n_____________________\nTrajectory: \n";
		cout << "time : " << time << "\n";
		cout << "length_limit : " << length_limit  << "\n";
		cout << "num_tcbonds : " << num_tcbonds  << "\n";
		cout << "num_tactive : " << num_tactive  << "\n";
		cout << "num_things : " << num_things  << "\n";
		cout << "num_polymers : " << num_polymers  << "\n";
		cout << "other_num_polymers : " << other_num_polymers  << "\n";
		cout << "num_monomers : " << num_monomers  << "\n";
		cout << "<polymer_completeness> : " << polymer_completeness / num_polymers << "\n";
		cout << "<polymer_bondedness> : " << polymer_bondedness / num_polymers  << "\n";
		cout << "<spacing> : " << (spacing_total / spacing_counter) << "\n";
		cout << "\n";
		*/



	}
}

void PolymerSystemIO::write_copy_output()
{
	// If a polymer is terminated
	// Then we write this to the output file

	if (outputfile.is_open())
	{
		int temp_idx = t_occupied_by[site_to_update];


		for (size_t i = 0; i < length_limit; i++)
		{

			if (c_seq[temp_idx][i] != -1)
			{
				outputfile << c_seq[temp_idx][i];
			}
			if (i == length_limit - 1)
			{
				outputfile << endl;
			}
			else
			{
				outputfile << ',';
			}

		}

	}

	pols_written++;
	/*
	cout << "__________________________";
	print_state_vertically();
	print_valid_transitions();
	print_rate_per_site();
	cout << "__________________________";
	*/
}

void PolymerSystemIO::print_state()
{
	// Ugly way of printing the state of the system
	cout << "      ";
	for (size_t idx = 0; idx < length_limit; idx++)
	{
		if (t_active[idx] == true)
		{
			cout << t_seq[idx] << " ";
		}
		else if (t_active[idx] == false)
		{
			cout << "X ";
		}


	}

	cout << "\n";
	cout << "      ";
	for (size_t idx = 0; idx < length_limit_max; idx++)
	{
		if (t_c_bond[idx] == true)
		{
			//cout << "| ";
			if (t_occupied_by[idx] > 9 || t_occupied_by[idx] == -1)
			{
				// formatting for double digits
				cout << t_occupied_by[idx];
			}
			else {
				cout << t_occupied_by[idx] << " ";
			}

		}
		else {
			cout << "  ";
		}
	}
	cout << "\n";



	for (size_t cop_i = 0; cop_i < pols_limit; cop_i++)
	{
		if (polymer_present[cop_i] == true)
		{
			if (cop_i < 10)
			{
				cout << " " << cop_i << " :  ";
			}
			else
			{
				cout << cop_i << " :  ";
			}

			size_t idx = 0;
			int trailing_flag = 0;
			while (idx < length_limit)
			{
				if (c_seq[cop_i][idx] == -1)
				{

					if (trailing_flag == 1)
					{
						trailing_flag = 2;
						break;
					}
					cout << " ";
				}
				if (c_seq[cop_i][idx] != -1) {

					cout << c_seq[cop_i][idx] << "";
					if (trailing_flag == 0)
					{
						trailing_flag = 1;
					}
				}

				if (bb_bond[cop_i][idx] == true)
				{
					cout << "-";
				}
				else
				{
					cout << " ";
				}
				idx++;
			}
			cout << "\n";
		}
	}

}

void PolymerSystemIO::print_state_vertically() {
	// Pretty way of printing state of system

	cout << "\n\n\n\n\n\n\n\n\n";
	cout << "\n State at T = " << time << "\n";
	bool layer_is_empty = false;

	for (int layer = length_limit; layer > 0; layer--)
	{
		layer_is_empty = true; // Flag to indicate whether the current layer needs be printed (avoids 200 lines of empty being printed)

		for (int idx = 0; idx < length_limit; idx++)
		{
			if (idx > 0 && idx - layer >= 0 && idx - layer < length_limit)
			{
				if (t_c_bond[idx] == true && t_occupied_by[idx - 1] != t_occupied_by[idx] && c_seq[t_occupied_by[idx]][idx - layer] != -1)
				{
					layer_is_empty = false;
					break;
				}
			}
			if (idx < length_limit - 1 && idx + layer >= 0 && idx + layer < length_limit)
			{
				// Print any leading polymer ends
				if (t_c_bond[idx] == true && t_occupied_by[idx + 1] != t_occupied_by[idx] && c_seq[t_occupied_by[idx]][idx + layer] != -1)
				{
					layer_is_empty = false;
					break;
				}
			}
		}

		if (layer_is_empty == false)
		{
			cout << "\n      ";
			for (int idx = 0; idx < length_limit; idx++)
			{
				if (idx > 0 && idx - layer >= 0 &&
					idx - layer < length_limit &&
					t_c_bond[idx] == true &&
					t_occupied_by[idx - 1] != t_occupied_by[idx] &&
					c_seq[t_occupied_by[idx]][idx - layer] != -1 &&
					t_c_bond[idx] == true &&
					t_occupied_by[idx - 1] != t_occupied_by[idx] &&
					c_seq[t_occupied_by[idx]][idx - layer] != -1)
				{
					// Print any lagging polymer ends
					cout << c_seq[t_occupied_by[idx]][idx - layer];
				}
				else {
					cout << " ";
				}


				if (idx < length_limit - 1 &&
					idx + layer >= 0 &&
					idx + layer < length_limit &&
					t_c_bond[idx] == true &&
					t_occupied_by[idx + 1] != t_occupied_by[idx] &&
					c_seq[t_occupied_by[idx]][idx + layer] != -1)
				{
					// Print any leading polymer ends
					cout << c_seq[t_occupied_by[idx]][idx + layer];
				}
				else {
					cout << " ";
				}


				cout << " ";
			}
		}
	}


	// DISPLAY CONNECTED POLYMERS
	cout << "\n      ";
	for (size_t idx = 0; idx < length_limit; idx++)
	{
		if (t_c_bond[idx] == true)
		{
			cout << c_seq[t_occupied_by[idx]][idx];

			if (bb_bond[t_occupied_by[idx]][idx] == true)
			{
				if (t_occupied_by[idx] != t_occupied_by[idx + 1]) {
					cout << "  ";
				}
				else {
					cout << "--";
				}
			}
			else {
				cout << "  ";
			}

		}
		else
		{
			cout << "   ";
		}

	}

	// DISPLAY TC BONDS
	cout << "\n      ";
	for (size_t idx = 0; idx < length_limit; idx++)
	{
		if (t_c_bond[idx] == true)
		{
			cout << "|  ";
		}
		else
		{
			cout << "   ";
		}

	}

	cout << "\n      ";


	for (size_t idx = 0; idx < length_limit; idx++)
	{
		if (t_active[idx] == true)
		{
			cout << t_seq[idx] << "  ";
		}
		else if (t_active[idx] == false)
		{
			cout << "X  ";
		}

	}
}

void PolymerSystemIO::print_valid_transitions()
{
	// Prints matrix which stores the valid transitions at each site
	cout << " \n\nValid transitions\n\n";
	cout << " T = " << time << "\n\n";
	for (size_t rule_i = 0; rule_i < num_rules; rule_i++)
	{
		cout << " " << rule_i;
		if (rule_i < 10)
		{
			cout << " ";
		}
		cout << " : ";

		for (size_t idx = 0; idx < length_limit; idx++)
		{
			if (valid_transitions[rule_i][idx] == true)
			{
				cout << "1 ";
			}
			else {
				cout << "0 ";
			}

		}
		cout << "   " << rule_labels[rule_i] << "\n";
	}
}

void PolymerSystemIO::print_valid_transitions_per_site(int site_idx) {
	// Prints valid transitions at a specific site
	cout << "\n Valid transitions at site " << site_idx << "\n";

	for (size_t rule_idx = 0; rule_idx < num_rules; rule_idx++)
	{

		if (valid_transitions[rule_idx][site_idx] == true) {
			cout << " Rule " << rule_idx;
			if (rule_idx < 10)
			{
				cout << " ";
			}
			cout << ": " << rates[rule_idx] << "    " << rule_labels[rule_idx];
			if (transition_to_update == rule_idx) { cout << " *** selected ***"; }

			cout << "\n";
		}
	}
	cout << endl;
}

void PolymerSystemIO::print_rate_per_site() {
	long double temp_sum = 0;
	for (size_t idx = 0; idx < length_limit; idx++)
	{
		temp_sum += total_rate_per_site[idx];
	}

	cout << "\nTotal rate per site\n";
	std::cout << std::left << std::setw(6) << "Site" << std::setw(14) << "Rate" << std::setw(14) << "Probability" << std::endl;



	for (size_t idx = 0; idx < length_limit; idx++)
	{

		cout << std::left << std::setw(6) << idx << std::setw(14) << total_rate_per_site[idx] << std::setw(14) << total_rate_per_site[idx] / temp_sum;
		if (site_to_update == idx) { cout << "  *** selected ***"; }

		cout << endl;
	}
}