/*
 * Experiments for monotone conjunctions
 * Dimitris Diochnos
 */

/****************************************************************************************/

#include <iostream>
#include <fstream>
#include <ctime>
#include <cmath>
#include <cstdlib>
#include <cassert>
#include <vector>
#include <algorithm>    // std::shuffle
#include <random>       // std::default_random_engine

/****************************************************************************************/

#define PROB_THRESHOLD_PER_THOUSAND 500
#define PROB_LAMBDA ((long double) PROB_THRESHOLD_PER_THOUSAND / 1000.0)
#define DIMENSION 100
#define EPSILON 0.05
#define DELTA 0.1
#define GAMMA 0.4
#define XI 0.9
#define NUM_EXPERIMENTS_PER_TARGET_SIZE 1000

/****************************************************************************************/

class Experiment {
public:
	std::vector<char> create_truth_assignment (const int dimension) {
		std::vector<char> ta;
		for (int i = 0; i < dimension; i++) {
			//char current_bit = static_cast<char> (rand () % 2);
			//assert ((current_bit == 0) || (current_bit == 1));
			int temp = rand () % 1000;
			char current_bit = (temp >= PROB_THRESHOLD_PER_THOUSAND) ? 0 : 1;
			//std::cout << "pushing back: " << static_cast<int> ( current_bit ) << std::endl;
			ta.push_back(current_bit);
		}
		return ta;
	}
	
	std::vector<std::vector<char> > create_sample (const int dimension, const int num_samples) {
		std::vector<std::vector<char> > the_sample;
		for (int i = 0; i < num_samples; i++) {
			the_sample.push_back( create_truth_assignment(dimension) );
		}
		return the_sample;
	}
	/*
	std::vector<std::vector<char> > create_EVO_sample (const int dimension) {
		create_sample (dimension, 0);
	}*/
	
	std::vector<char> create_target (const int dimension, const int size) {
		assert( (1 <= size) && (size <= dimension) );
		std::vector<int> my_shuffle;
		for (int i = 0; i < dimension; i++) {
			my_shuffle.push_back( i );
		}
		//std::shuffle ( my_shuffle.begin(), my_shuffle.end(),  std::default_random_engine(153));
		std::random_shuffle ( my_shuffle.begin(), my_shuffle.end());

		/*
		std::cout << "the shuffled vector is: ";
		for (int i = 0; i < my_shuffle.size(); i++) {
			std::cout << my_shuffle[i] << " ";
		}
		std::cout << std::endl;
		*/
		
		std::vector<char> the_target;
		for (int i = 0; i < dimension; i++) {
			the_target.push_back( static_cast<char> (0) );
		}
		for (int i = 0; i < size; i++) {
			the_target[my_shuffle[i]] = static_cast<char> (1);
		}
		return the_target;
	}
	
	bool evaluate (std::vector<char> const &target, std::vector<char> const &truth_assignment) {
		assert(target.size() == truth_assignment.size());
		
		for (int i = 0; i < target.size(); i++) {
			bool var_participates = (static_cast<int> (target[i]) == 1) ? true : false;
			if (var_participates == true) {
				if (static_cast<int>(truth_assignment[i]) == 0)
					return false;
			}
		}
		return true; // default value
	}
	
	int pac_sample_size (const long double epsilon, const long double delta, const int dimension) {
		long double md = (1.0/epsilon)*(dimension * log(2) + log(1.0/delta));
		return static_cast<int> (ceil(md));
	}
	
	int extended_pac_sample_size (const long double epsilon, const long double delta, const long double gamma, const long double xi, long double lower_bound_target_weight, const int dimension) {
		assert( gamma <= 0.5 );
		long double mymin       = epsilon;
		long double second_term = gamma * lower_bound_target_weight;
		long double third_term  = xi * lower_bound_target_weight / 2.0;
		if (second_term < mymin) {
			mymin = second_term;
		}
		if (third_term < mymin) {
			mymin = third_term;
		}
		long double md = (1.0/mymin)*(dimension * log(2) + log(1.0/delta));
		return static_cast<int> (ceil(md));
	}
	
	std::vector<char> create_initial_hypothesis (const int dimension) {
		std::vector<char> myguess;
		for (int i = 0; i < dimension; i++) {
			myguess.push_back( static_cast<char> (1) );
		}
		return myguess;
	}
	
	int get_hypothesis_length (std::vector<char> const &h) {
		int l = 0;
		for (int i = 0; i < h.size(); i++) {
			if ( static_cast<int>(h[i]) == 1 ) l++;
		}
		return l;
	}
	
	long double get_risk (std::vector<char> const &c, std::vector<char> const &h) {
		assert (c.size() == DIMENSION );
		assert ( c.size() == h.size() );
		int m = 0, u = 0, w = 0;
		for (int i = 0; i < c.size(); i++) {
			if ( static_cast<int> (c[i]) == 1 ) {
				if ( static_cast<int> (h[i]) == 1 ) m++;
				else u++;
			}
			else if ( static_cast<int> (h[i]) == 1 ) w++;
		}
        	// risk = lambda^{m}(lambda^{u} + lambda^{w} - 2lambda^{u+w})
		// See Proposition 4.1 from the paper
		//std::cout << "     m = " << m << "   u = " << u << "   w = " << w << "          " << std::flush;
		return pow(static_cast<long double>(PROB_LAMBDA), static_cast<long double>(m+u)) + pow(static_cast<long double>(PROB_LAMBDA), static_cast<long double>(m+w)) - static_cast<long double>(2.0) * pow(static_cast<long double>(PROB_LAMBDA), static_cast<long double>(m+u+w));
	}
	
	long double get_recall (std::vector<char> const &c, std::vector<char> const &h) {
		assert ( c.size() == h.size() );
		int m = 0, u = 0, w = 0;
		for (int i = 0; i < c.size(); i++) {
			if ( static_cast<int> (c[i]) == 1 ) {
				if ( static_cast<int> (h[i]) == 1 ) m++;
				else u++;
			}
			else if ( static_cast<int> (h[i]) == 1 ) w++;
		}
		if (w == 0) return static_cast<long double>(1.0);
		else {
			// recall = lambda^{w}      (due to Proposition 4.1)
			return pow(static_cast<long double>(PROB_LAMBDA), static_cast<long double>(w));
		}
	}
	
	long double get_precision (std::vector<char> const &c, std::vector<char> const &h) {
		assert ( c.size() == h.size() );
		int m = 0, u = 0, w = 0;
		for (int i = 0; i < c.size(); i++) {
			if ( static_cast<int> (c[i]) == 1 ) {
				if ( static_cast<int> (h[i]) == 1 ) m++;
				else u++;
			}
			else if ( static_cast<int> (h[i]) == 1 ) w++;
		}
		if (u == 0) return static_cast<long double>(1.0);
		else {
			// precision = lambda^{u}   (due to Proposition 4.1)
			return pow(static_cast<long double>(PROB_LAMBDA), static_cast<long double>(u));
		}
	}
	
	//---------------------------------------------------------------------
	
	std::vector<char> create_EVO_initial_hypothesis (const int dimension, const int frontier) {
		assert( (1 <= frontier) && (frontier <= dimension) );
		std::vector<char> myguess;
		for (int i = 0; i < dimension; i++) {
			myguess.push_back( static_cast<char> (0) );
		}
		return myguess;
	}
	
	long double get_correlation (std::vector<char> const &target, std::vector<char> const &hyp) {
		// Here we are exploiting again the fact that we have a closed formula for the correlation
		// This saves us computational resources as we do not have to simulate the whole process
		// and inspect directly the outcome in the case of success.
		return static_cast<long double> (1.0) - static_cast<long double> (2.0) * get_risk (target, hyp);
	}

	std::vector<std::vector<char> > get_neighborhood_of_interest (std::vector<char> const &target,
																  std::vector<char> const &hyp,
																  const int frontier,
																  const long double tolerance) {
		std::vector<std::vector<char> > all_neighbors;
		std::vector<std::vector<char> > bene, neut, dele;
		int hyp_length = get_hypothesis_length(hyp);
		
		// Create neighbors by adding one variable
		if (hyp_length < frontier) {
			// In this case, N+ exists
			for (int i = 0; i < hyp.size(); i++) {
				if (hyp[i] == static_cast<char>(0)) {
					// a new candidate is to flip this 0 to a 1
					std::vector<char> alien_friend = hyp;
					alien_friend[i] = static_cast<char>(1);
					all_neighbors.push_back(alien_friend);
				}
			}
		}
		
		// Create neighbors by deleting one variable
		if (hyp_length > 0) {
			// In this case, N- exists
			for (int i = 0; i < hyp.size(); i++) {
				if (hyp[i] == static_cast<char>(1)) {
					// a new candidate is to flip this 0 to a 1
					std::vector<char> alien_friend = hyp;
					alien_friend[i] = static_cast<char>(0);
					all_neighbors.push_back(alien_friend);
				}
			}
		}
		
		// Create neighbors by forming all possible swaps
		for (int i = 0; i < hyp.size(); i++) {
			if (hyp[i] == static_cast<char>(1)) {
				// This is a candidate for removal
				std::vector<char> alien_friend;
				for (int j = 0; j < hyp.size(); j++) {
					if (hyp[j] == static_cast<char>(0)) {
						// this can be flipped
						alien_friend = hyp;
						alien_friend[i] = 0;
						alien_friend[j] = 1;
						all_neighbors.push_back(alien_friend);
					}
				}
			}
		}
		
		/*
		 * All neighbors created (we have not added the hypothesis itself without any changes yet ...)
		 * Partition neighbors to bene, neut, deleterious
		 */
		long double base_corr = get_correlation(target, hyp);
		for (int i = 0; i < all_neighbors.size(); i++) {
			long double current_corr = get_correlation(target, all_neighbors[i]);
			if (current_corr > base_corr + tolerance) {
				// this is beneficial
				bene.push_back(all_neighbors[i]);
			} else if (current_corr >= base_corr - tolerance) {
				// this is neutral
				neut.push_back(all_neighbors[i]);
			} else {
				// this is deleterious
				dele.push_back(all_neighbors[i]);
			}
		}
		// Add the current hypothesis into the neighborhood
		neut.push_back(hyp);
		
		// The selection mechanism is the following:
		// If bene non-empty, this is the set of interest
		// otherwise, neut is non-empty (current hyp is always there) and this is the set of interest
		if (bene.size() > 0)
			return bene;
		assert(neut.size() > 0);
		return neut;
	}
	
	//---------------------------------------------------------------------
	
	std::vector<int> create_vector_with_target_sizes () {
		// {1, 5, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 50, 100, 500, 1000};
		std::vector<int> list;
		//list.push_back(1);
		
		list.push_back(2);
		list.push_back(3);
		list.push_back(4);
		list.push_back(5);
		list.push_back(6);
		list.push_back(7);
		list.push_back(8);
		list.push_back(9);
		list.push_back(10);
		list.push_back(11);
		list.push_back(12);

		return list;
	}
	
	//---------------------------------------------------------------------
	long double find_median_among_ints (std::vector<int> &h_array) {
		int array_length = h_array.size();
		std::sort (h_array.begin(), h_array.end());
		if (array_length % 2 == 1) {
			// odd length ==> median is right in the middle
			return static_cast<long double> ( h_array[static_cast<int> (floor(array_length / 2.0))] );
		} else {
			// even length ==> median is average of the two in the middle
			int index_1 = static_cast<int> (floor (array_length / 2.0));
			int index_2 = static_cast<int> (floor (array_length / 2.0)) - 1;
			long double avg = h_array[index_1] + h_array[index_2];
			avg = avg / 2.0;
			return avg;
		}
	}
	long double find_median_among_long_doubles (std::vector<long double> &error_array) {
		int array_length = error_array.size();
		std::sort (error_array.begin(), error_array.end());
		if (array_length % 2 == 1) {
			// odd length ==> median is right in the middle
			return static_cast<long double> ( error_array[static_cast<int> (floor(array_length / 2.0))] );
		} else {
			// even length ==> median is average of the two in the middle
			int index_1 = static_cast<int> (floor (array_length / 2.0));
			int index_2 = static_cast<int> (floor (array_length / 2.0)) - 1;
			long double avg = error_array[index_1] + error_array[index_2];
			avg = avg / 2.0;
			return avg;
		}
	}
	
	long double find_average_among_ints (std::vector<int> const &length_array) {
		int array_length = length_array.size();
		int sum = 0;
		for (int i = 0; i < array_length; i++) {
			sum += length_array[i];
		}
		return static_cast<long double> (sum) / static_cast<long double> (array_length);
	}
	long double find_average_among_long_doubles (std::vector<long double> const &error_array) {
		int array_length = error_array.size();
		long double sum = 0.0;
		for (int i = 0; i < array_length; i++) {
			sum += error_array[i];
		}
		return sum / static_cast<long double> (array_length);
	}
	
	long double find_min_among_long_doubles (std::vector<long double> const &floating_point_array) {
		int array_length = floating_point_array.size();
		assert(array_length >= 1);
		long double min = floating_point_array[0];
		for (int i = 1; i < array_length; i++) {
			if (floating_point_array[i] < min) {
				min = floating_point_array[i];
			}
		}
		return min;
	}
	long double find_max_among_long_doubles (std::vector<long double> const &floating_point_array) {
		int array_length = floating_point_array.size();
		assert(array_length >= 1);
		long double max = floating_point_array[0];
		for (int i = 1; i < array_length; i++) {
			if (floating_point_array[i] > max) {
				max = floating_point_array[i];
			}
		}
		return max;
	}
	long double find_max_among_ints (std::vector<int> const &int_array) {
		int array_length = int_array.size();
		assert(array_length >= 1);
		int max = int_array[0];
		for (int i = 1; i < array_length; i++) {
			if (int_array[i] > max) {
				max = int_array[i];
			}
		}
		return static_cast<long double> (max);
	}
	//---------------------------------------------------------------------

	int count_falsified_vars (std::vector<char> const &hypothesis, std::vector<char> const &assignment) {
		int num = 0;
		for (int i = 0; i < hypothesis.size(); i++) {
			if ( (static_cast<int> (hypothesis[i]) == 1) && (static_cast<int>(assignment[i]) == 0) )
				num++;
		}
		return num;
	}
	
	bool is_target_and_hypothesis_same (std::vector<char> const &target, std::vector<char> const &hypothesis) {
		assert ( target.size() == hypothesis.size() );
		for (int i = 0; i < target.size(); i++)
			if (target[i] != hypothesis[i])
				return false; // no need for typecasting ...
		return true;
	}
};

/****************************************************************************************/
/****************************************************************************************/

void driver () {
	Experiment e;
	
	std::cout << "sizeof(double) = " << sizeof(double) << std::endl;
	std::cout << "sizeof(long double) = " << sizeof(long double) << std::endl;
	
	//std::vector<char> test_char_vector = e.create_truth_assignment(10);
	//for (int i = 0; i < test_char_vector.size(); i++) {
	//	std::cout << static_cast<int> ( test_char_vector[i] ) << " ";
	//}
	//std::cout << std::endl;
	std::vector<std::vector<char> > sample = e.create_sample(10, 2);
	
	for (int i = 0; i < sample.size(); i++) {
		std::cout << "sample " << i+1 << ": ";
		for (int j = 0; j < sample[i].size(); j++) {
			std::cout << static_cast<int> (sample[i][j]) << " ";
		}
		std::cout << std::endl;
	}
	
	std::vector<char> target = e.create_target(10, 2);
	std::cout << "target in binary is   : ";
	for (int i = 0; i < target.size(); i++) {
		std::cout << static_cast<int> (target[i]) << " ";
	}
	std::cout << std::endl;
	
	std::cout << "Evaluation on s0: " << (e.evaluate(target, sample[0]) ? "true" : "false") << std::endl;
	std::cout << "Evaluation on s1: " << (e.evaluate(target, sample[1]) ? "true" : "false") << std::endl;
	
	std::cout << std::endl;
	int my_dimension = 100;
	int my_pac_sample_size = e.pac_sample_size(EPSILON, DELTA, my_dimension);
	std::cout << "PAC sample size for epsilon " << EPSILON << " and delta " << DELTA << " is (dimension = " << my_dimension << "): " << my_pac_sample_size << std::endl;
}

/****************************************************************************************/
/****************************************************************************************/

void run_PAC_experiment () {
	Experiment e;
	
	std::time_t t;
	std::tm* now;

	// Declarations
	long double current_median, current_average, current_min, current_max;
	std::vector<long double> hyp_sizes_avg, hyp_sizes_med, hyp_sizes_max;
	std::vector<long double> hyp_errors_avg, hyp_errors_med, hyp_errors_max;
	std::vector<long double> hyp_recalls_min, hyp_recalls_avg, hyp_recalls_med, hyp_recalls_max;
	std::vector<long double> hyp_precisions_min, hyp_precisions_avg, hyp_precisions_med;
	std::vector<std::vector<int> > hypotheses_sizes;
	std::vector<std::vector<long double> > hypotheses_errors;
	std::vector<std::vector<long double> > hypotheses_recalls;
	std::vector<std::vector<long double> > hypotheses_precisions;
	
	// Files
	std::cout << "open files now" << std::endl;
	std::ofstream fd_sizes_avg, fd_sizes_med, fd_sizes_max, fd_sizes; // store info about avg, median and all sizes
	std::ofstream fd_errors_avg, fd_errors_med, fd_errors_max, fd_errors; // store info about avg, median all errors
	std::ofstream fd_recalls_avg, fd_recalls_med, fd_recalls_min, fd_recalls_max, fd_recalls; // store info about recalls
	std::ofstream fd_precisions_avg, fd_precisions_med, fd_precisions_min, fd_precisions; // store info about precisions
	std::ofstream fd_log; // high level log information while conducting the experiment
	fd_sizes.open("PAC_stats_sizes_all.txt");
	fd_sizes_avg.open("PAC_stats_sizes_avg.txt");
	fd_sizes_med.open("PAC_stats_sizes_med.txt");
	fd_sizes_max.open("PAC_stats_sizes_max.txt");
	fd_errors.open("PAC_stats_risk_all.txt");
	fd_errors_avg.open("PAC_stats_risk_avg.txt");
	fd_errors_med.open("PAC_stats_risk_med.txt");
	fd_errors_max.open("PAC_stats_risk_max.txt");
	fd_recalls.open("PAC_stats_recalls_all.txt");
	fd_recalls_min.open("PAC_stats_recalls_min.txt");
	fd_recalls_avg.open("PAC_stats_recalls_avg.txt");
	fd_recalls_med.open("PAC_stats_recalls_med.txt");
	fd_recalls_max.open("PAC_stats_recalls_max.txt");
	fd_precisions.open("PAC_stats_precisions_all.txt");
	fd_precisions_min.open("PAC_stats_precisions_min.txt");
	fd_precisions_avg.open("PAC_stats_precisions_avg.txt");
	fd_precisions_med.open("PAC_stats_precisions_med.txt");
	fd_log.open("PAC_log.txt");
	std::cout << std::endl;
	std::vector<int> target_sizes = e.create_vector_with_target_sizes();
	
	// Book-keeping
	fd_log << "Dimension   = " << DIMENSION << std::endl;
	fd_log << "Epsilon     = " << EPSILON << std::endl;
	fd_log << "Delta       = " << DELTA << std::endl;
	fd_log << "Gamma       = " << GAMMA << std::endl;
	fd_log << "Xi          = " << XI << std::endl;
	fd_log << "Probability = " << PROB_LAMBDA << std::endl;
	//
	std::cout << "Dimension   = " << DIMENSION << std::endl;
	std::cout << "Epsilon     = " << EPSILON << std::endl;
	std::cout << "Delta       = " << DELTA << std::endl;
	std::cout << "Gamma       = " << GAMMA << std::endl;
	std::cout << "Xi          = " << XI << std::endl;
	std::cout << "Probability = " << PROB_LAMBDA << std::endl;
	
	int default_sample_size = e.pac_sample_size(EPSILON, DELTA, DIMENSION);
	fd_log << "PAC sample size            = " << default_sample_size << std::endl;
	fd_log << "Num experiments per target size = " << NUM_EXPERIMENTS_PER_TARGET_SIZE << std::endl;
	//
	t = std::time(0); now = std::localtime(&t);
	fd_log << "Experiments started: " << asctime(now);
		
	for (int j = 0; j < target_sizes.size(); j++) {
		// timestamp
		int current_target_size = target_sizes[j];
		t = std::time(0); now = std::localtime(&t);
		std::cout << "~~~~~~ " << "target size is " << current_target_size << " (case " << j+1 << " out of " << target_sizes.size() << ") ---- " << asctime(now); //<< " ~~~~~" << std::endl;
		fd_log << "- target size " << current_target_size << std::flush;
		
		//
		long double lower_bound_target_weight = pow(PROB_LAMBDA, current_target_size);
		std::cout << "       target size: " << current_target_size << " weight lower bound: " << lower_bound_target_weight << std::endl;
		int extended_sample_size = e.extended_pac_sample_size(EPSILON, DELTA, GAMMA, XI, lower_bound_target_weight, DIMENSION);
		fd_log << " with PAC sample size (extended) = " << extended_sample_size << std::endl;
		
		
		
		// Starting lines in files
		fd_sizes << current_target_size << " " << std::flush;
		fd_sizes_avg << current_target_size << " " << std::flush;
		fd_sizes_med << current_target_size << " " << std::flush;
		fd_sizes_max << current_target_size << " " << std::flush;
		fd_errors << current_target_size << " " << std::flush;
		fd_errors_avg << current_target_size << " " << std::flush;
		fd_errors_med << current_target_size << " " << std::flush;
		fd_errors_max << current_target_size << " " << std::flush;
		fd_recalls << current_target_size << " " << std::flush;
		fd_recalls_avg << current_target_size << " " << std::flush;
		fd_recalls_med << current_target_size << " " << std::flush;
		fd_recalls_min << current_target_size << " " << std::flush;
		fd_recalls_max << current_target_size << " " << std::flush;
		fd_precisions << current_target_size << " " << std::flush;
		fd_precisions_avg << current_target_size << " " << std::flush;
		fd_precisions_med << current_target_size << " " << std::flush;
		fd_precisions_min << current_target_size << " " << std::flush;
		
		
		
		// More Declarations -- specific to this particular target size
		std::vector<char> myhypothesis;
		std::vector<std::vector<char> > examples;
		std::vector<int> hyp_sizes_against_this_target_size;
		std::vector<long double> hyp_errors_against_this_target_size;
		std::vector<long double> hyp_recalls_against_this_target_size;
		std::vector<long double> hyp_precisions_against_this_target_size;
		
		
		for (int i = 0; i < NUM_EXPERIMENTS_PER_TARGET_SIZE; i++) {
			std::vector<char> mytarget = e.create_target(DIMENSION, current_target_size);
			assert(e.get_hypothesis_length(mytarget) == current_target_size);
			
			//std::cout << " drawing examples ..." << std::flush;
			std::cout << "\r       experiment " << i+1 << " out of " << NUM_EXPERIMENTS_PER_TARGET_SIZE << std::flush;
			//
			//examples = e.create_sample (DIMENSION, default_sample_size);
			examples = e.create_sample (DIMENSION, extended_sample_size);
			std::cout << " sample size = " << examples.size() << std::flush;
			//
			//std::cout << " done (" << examples.size() << " examples) --- experiment " << i+1 << " out of " << NUM_EXPERIMENTS_PER_TARGET_SIZE << std::endl;
			
			//std::cout << "   learning ..." << std::flush;
			//std::cout << "initializing the hypothesis to include all variables" << std::endl;
			myhypothesis = e.create_initial_hypothesis(DIMENSION);
			/*
			 * Now PAC learn based on the examples that you drew
			 */
			for (int j = 0; j < examples.size(); j++) {
				bool label = e.evaluate(mytarget, examples[j]);
				if (label == true) {
					// Delete all the variables that contradict the label
					for (int k = 0; k < myhypothesis.size(); k++) {
						if (( static_cast<int> (myhypothesis[k]) == 1) && ( static_cast<int> (examples[j][k]) == 0)) {
							myhypothesis[k] = 0;
						}
					}
				}
			}
			int current_hyp_size          = e.get_hypothesis_length(myhypothesis);
			long double current_risk      = e.get_risk      (mytarget, myhypothesis);
			long double current_recall    = e.get_recall    (mytarget, myhypothesis);
			long double current_precision = e.get_precision (mytarget, myhypothesis);
			hyp_sizes_against_this_target_size.push_back(current_hyp_size);
			hyp_errors_against_this_target_size.push_back(current_risk);
			hyp_recalls_against_this_target_size.push_back(current_recall);
			hyp_precisions_against_this_target_size.push_back(current_precision);

			
			// Write individual entries to files
			fd_sizes << current_hyp_size << " " << std::flush;
			fd_errors << current_risk << " " << std::flush;
			fd_recalls << current_recall << " " << std::flush;
			fd_precisions << current_precision << " " << std::flush;
			
			
			//std::cout << " done -- (hypothesis size: " << current_hyp_size << ", with risk: " << current_risk << ")" << std::endl;
			//std::cout << "I learned the following hypothesis (length " << e.get_hypothesis_length(myhypothesis) << "): " << std::flush;
			//for (int i = 0; i < myhypothesis.size(); i++) {
			//	std::cout << static_cast<int> (myhypothesis[i]) << " ";
			//}
			//std::cout << std::endl;
			mytarget.clear();
			myhypothesis.clear();
			examples.clear();
		}
		
		// First push them, then compute
		hypotheses_sizes.push_back(hyp_sizes_against_this_target_size);
		hypotheses_errors.push_back(hyp_errors_against_this_target_size);
		hypotheses_recalls.push_back(hyp_recalls_against_this_target_size);
		hypotheses_precisions.push_back(hyp_precisions_against_this_target_size);
		
		// Compute average and median hypothesis size against this target size
		current_median  = e.find_median_among_ints(hyp_sizes_against_this_target_size);
		current_average = e.find_average_among_ints(hyp_sizes_against_this_target_size);
		current_max     = e.find_max_among_ints(hyp_sizes_against_this_target_size);
		//
		hyp_sizes_med.push_back(current_median);
		hyp_sizes_avg.push_back(current_average);
		hyp_sizes_max.push_back(current_max);
		//
		fd_sizes_med << current_median << std::endl;
		fd_sizes_avg << current_average << std::endl;
		fd_sizes_max << current_max << std::endl;
		
		// Compute average, median and max hypothesis error (risk) against this target size
		current_median  = e.find_median_among_long_doubles(hyp_errors_against_this_target_size);
		current_average = e.find_average_among_long_doubles(hyp_errors_against_this_target_size);
		current_max     = e.find_max_among_long_doubles(hyp_errors_against_this_target_size);
		//
		hyp_errors_med.push_back(current_median);
		hyp_errors_avg.push_back(current_average);
		hyp_errors_max.push_back(current_max);
		//
		fd_errors_med << current_median << std::endl;
		fd_errors_avg << current_average << std::endl;
		fd_errors_max << current_max << std::endl;
		
		
		// Compute min, average and median hypothesis recall against this target size
		current_min     = e.find_min_among_long_doubles(hyp_recalls_against_this_target_size);
		current_median  = e.find_median_among_long_doubles(hyp_recalls_against_this_target_size);
		current_average = e.find_average_among_long_doubles(hyp_recalls_against_this_target_size);
		current_max     = e.find_max_among_long_doubles(hyp_recalls_against_this_target_size);
		//
		hyp_recalls_med.push_back(current_median);
		hyp_recalls_avg.push_back(current_average);
		hyp_recalls_min.push_back(current_min);
		hyp_recalls_max.push_back(current_max);
		//
		fd_recalls_med << current_median << std::endl;
		fd_recalls_avg << current_average << std::endl;
		fd_recalls_min << current_min << std::endl;
		fd_recalls_max << current_max << std::endl;
		
		// Compute min, average and median hypothesis precision against this target size
		current_min     = e.find_min_among_long_doubles(hyp_precisions_against_this_target_size);
		current_median  = e.find_median_among_long_doubles(hyp_precisions_against_this_target_size);
		current_average = e.find_average_among_long_doubles(hyp_precisions_against_this_target_size);
		//
		hyp_precisions_med.push_back(current_median);
		hyp_precisions_avg.push_back(current_average);
		hyp_precisions_min.push_back(current_min);
		//
		fd_precisions_med << current_median << std::endl;
		fd_precisions_avg << current_average << std::endl;
		fd_precisions_min << current_min << std::endl;
		
		//std::cout << std::endl << "       max risk = " << e.find_max_among_long_doubles(hyp_errors_against_this_target_size) << ",   min rc = " << e.find_min_among_long_doubles(hyp_recalls_against_this_target_size) << ", avg rc = " << e.find_average_among_long_doubles(hyp_recalls_against_this_target_size) << ", max rc = " << e.find_max_among_long_doubles(hyp_recalls_against_this_target_size) << ",   min prec = " << e.find_min_among_long_doubles(hyp_precisions_against_this_target_size) << std::endl;
		std::cout << std::endl;
		std::cout << std::fixed << std::setprecision(8);
		std::cout << "       RISK     :   min = " << e.find_min_among_long_doubles(hyp_errors_against_this_target_size);
		std::cout << ",   med = " << e.find_median_among_long_doubles(hyp_errors_against_this_target_size);
		std::cout << ",   avg = " << e.find_average_among_long_doubles(hyp_errors_against_this_target_size);
		std::cout << ",   max = " << e.find_max_among_long_doubles(hyp_errors_against_this_target_size);
		std::cout << std::endl;
		std::cout << "       RECALL   :   min = " << e.find_min_among_long_doubles(hyp_recalls_against_this_target_size);
		std::cout << ",   med = " << e.find_median_among_long_doubles(hyp_recalls_against_this_target_size);
		std::cout << ",   avg = " << e.find_average_among_long_doubles(hyp_recalls_against_this_target_size);
		std::cout << ",   max = " << e.find_max_among_long_doubles(hyp_recalls_against_this_target_size);
		std::cout << std::endl;
		std::cout << "       PRECISION:" << "   min = " << e.find_min_among_long_doubles(hyp_precisions_against_this_target_size);
		std::cout << ",   med = " << e.find_median_among_long_doubles(hyp_precisions_against_this_target_size);
		std::cout << ",   avg = " << e.find_average_among_long_doubles(hyp_precisions_against_this_target_size);
		std::cout << ",   max = " << e.find_max_among_long_doubles(hyp_precisions_against_this_target_size);
		std::cout << std::endl;
		
		
		// Cleaning up ....
		hyp_sizes_against_this_target_size.clear();
		hyp_errors_against_this_target_size.clear();
		hyp_recalls_against_this_target_size.clear();
		hyp_precisions_against_this_target_size.clear();
		
		fd_sizes << std::endl;
		fd_errors << std::endl;
		fd_recalls << std::endl;
		fd_precisions << std::endl;
		std::cout << std::endl;
	}
	
	fd_log << std::endl;
	t = std::time(0); now = std::localtime(&t);
	fd_log << "Experiments ended  : " << asctime(now) << std::endl;
	fd_log.close();
	
	hypotheses_sizes.clear();
	hypotheses_errors.clear();
	hypotheses_recalls.clear();
	hypotheses_precisions.clear();
	
	hyp_sizes_avg.clear();
	hyp_sizes_med.clear();
	hyp_sizes_max.clear();
	hyp_errors_avg.clear();
	hyp_errors_med.clear();
	hyp_errors_max.clear();
	
	hyp_recalls_med.clear();
	hyp_recalls_avg.clear();
	hyp_recalls_min.clear();
	hyp_recalls_max.clear();
	hyp_precisions_med.clear();
	hyp_precisions_avg.clear();
	hyp_precisions_min.clear();
	
	
	
	// Close the files
	std::cout << "Closing files now" << std::endl;
	fd_sizes_avg.close();
	fd_sizes_med.close();
	fd_sizes_max.close();
	fd_sizes.close();
	fd_errors_avg.close();
	fd_errors_med.close();
	fd_errors_max.close();
	fd_errors.close();
	fd_recalls_min.close();
	fd_recalls_med.close();
	fd_recalls_avg.close();
	fd_recalls_max.close();
	fd_recalls.close();
	fd_precisions_min.close();
	fd_precisions_med.close();
	fd_precisions_avg.close();
	fd_precisions.close();
	
	// final timestamp
	t = std::time(0); now = std::localtime(&t);
	fd_log << "Entire experiment ended: " << asctime(now);
	fd_log.close();
}


/****************************************************************************************/
/****************************************************************************************/


void run_EVO_experiment () {
	Experiment e;
	//driver();
	
	std::time_t t;
	std::tm* now;

	// Declarations
	long double current_median, current_average, current_min, current_max;
	std::vector<long double> hyp_sizes_avg, hyp_sizes_med, hyp_sizes_max;
	std::vector<long double> hyp_errors_avg, hyp_errors_med, hyp_errors_max;
	std::vector<long double> hyp_recalls_min, hyp_recalls_avg, hyp_recalls_med, hyp_recalls_max;
	std::vector<long double> hyp_precisions_min, hyp_precisions_avg, hyp_precisions_med;
	std::vector<std::vector<int> > hypotheses_sizes;
	std::vector<std::vector<long double> > hypotheses_errors;
	std::vector<std::vector<long double> > hypotheses_recalls;
	std::vector<std::vector<long double> > hypotheses_precisions;
	
	// Files
	std::cout << "open files now" << std::endl;
	std::ofstream fd_sizes_avg, fd_sizes_med, fd_sizes_max, fd_sizes; // store info about avg, median and all sizes
	std::ofstream fd_errors_avg, fd_errors_med, fd_errors_max, fd_errors; // store info about avg, median all errors
	std::ofstream fd_recalls_avg, fd_recalls_med, fd_recalls_min, fd_recalls_max, fd_recalls; // store info about recalls
	std::ofstream fd_precisions_avg, fd_precisions_med, fd_precisions_min, fd_precisions; // store info about precisions
	std::ofstream fd_log; // high level log information while conducting the experiment
	fd_sizes.open("EVO_stats_sizes_all.txt");
	fd_sizes_avg.open("EVO_stats_sizes_avg.txt");
	fd_sizes_med.open("EVO_stats_sizes_med.txt");
	fd_sizes_max.open("EVO_stats_sizes_max.txt");
	fd_errors.open("EVO_stats_risk_all.txt");
	fd_errors_avg.open("EVO_stats_risk_avg.txt");
	fd_errors_med.open("EVO_stats_risk_med.txt");
	fd_errors_max.open("EVO_stats_risk_max.txt");
	fd_recalls.open("EVO_stats_recalls_all.txt");
	fd_recalls_min.open("EVO_stats_recalls_min.txt");
	fd_recalls_avg.open("EVO_stats_recalls_avg.txt");
	fd_recalls_med.open("EVO_stats_recalls_med.txt");
	fd_recalls_max.open("EVO_stats_recalls_max.txt");
	fd_precisions.open("EVO_stats_precisions_all.txt");
	fd_precisions_min.open("EVO_stats_precisions_min.txt");
	fd_precisions_avg.open("EVO_stats_precisions_avg.txt");
	fd_precisions_med.open("EVO_stats_precisions_med.txt");
	fd_log.open("EVO_log.txt");
	std::cout << std::endl;
	std::vector<int> target_sizes = e.create_vector_with_target_sizes();

	
	// Book-keeping
	fd_log << "Dimension = " << DIMENSION << std::endl;
	fd_log << "Epsilon   = " << EPSILON << std::endl;
	fd_log << "Delta     = " << DELTA << std::endl;
	fd_log << "Gamma     = " << GAMMA << std::endl;
	fd_log << "Xi        = " << XI << std::endl;
	fd_log << "Probability = " << PROB_LAMBDA << std::endl;
	//
	std::cout << "Dimension   = " << DIMENSION << std::endl;
	std::cout << "Epsilon     = " << EPSILON << std::endl;
	std::cout << "Delta       = " << DELTA << std::endl;
	std::cout << "Gamma       = " << GAMMA << std::endl;
	std::cout << "Xi          = " << XI << std::endl;
	std::cout << "Probability = " << PROB_LAMBDA << std::endl;
	std::cout << std::endl;

	
	int default_sample_size = e.pac_sample_size(EPSILON, DELTA, DIMENSION);
	fd_log << "PAC sample size = " << default_sample_size << " (ignored b/c we are in the evolvability case)" << std::endl;
	fd_log << "Num experiments per target size = " << NUM_EXPERIMENTS_PER_TARGET_SIZE << std::endl;
	
	t = std::time(0); now = std::localtime(&t);
	fd_log << "Experiments started: " << asctime(now);
	
	for (int j = 0; j < target_sizes.size(); j++) {
		// timestamp
		int current_target_size = target_sizes[j];
		t = std::time(0); now = std::localtime(&t);
		std::cout << "~~~~~~ " << "target size is " << current_target_size << " (case " << j+1 << " out of " << target_sizes.size() << ") ---- " << asctime(now); //<< " ~~~~~" << std::endl;
		fd_log << "- target size " << current_target_size << std::endl;
		
		long double lower_bound_target_weight = pow(PROB_LAMBDA, current_target_size);
		std::cout << "       target size: " << current_target_size << " weight lower bound: " << lower_bound_target_weight << std::endl;
		
		
        	// This definition of frontier is for the uniform distribution only
		int frontier = static_cast<int> ( ceil(log(3.0/(2.0 * EPSILON))/log(2.0)) ); // EPSILON is prob of error (not how far off the correlation function is)
		/*
		long double critical_minimum = EPSILON;
		if (critical_minimum > GAMMA * lower_bound_target_weight) {
			critical_minimum = GAMMA * lower_bound_target_weight;
		}
		if (critical_minimum > XI * lower_bound_target_weight / 2.0) {
			critical_minimum = XI * lower_bound_target_weight / 2.0;
		}
		int frontier = static_cast<int> ( ceil(log(3.0/(2.0 * critical_minimum))/log(1.0/PROB_LAMBDA)) ); // Extended version for PAC learning with high recall and high precision
		 */
		fd_log << "Frontier = " << frontier << std::endl;
		std::cout << "       Frontier = " << frontier << std::endl;
		
		// This definition of tolerance is for the uniform distribution only;
        	// otherwise see the paper https://link.springer.com/chapter/10.1007/978-3-319-46379-7_7
		long double tolerance = static_cast<long double> ( pow( static_cast<long double> (2.0), static_cast<long double> (0.0-2.0*frontier)) );
		fd_log << "Tolerance ( 2^(-2q) ) = " << tolerance << std::endl;
		int num_generations = 2*frontier; // In general it is 3*q, but under the uniform distribution 2*q generations are enough!
		fd_log << "Num generations = " << num_generations << std::endl;



		// Starting lines in files
		fd_sizes << current_target_size << " " << std::flush;
		fd_sizes_avg << current_target_size << " " << std::flush;
		fd_sizes_med << current_target_size << " " << std::flush;
		fd_sizes_max << current_target_size << " " << std::flush;
		fd_errors << current_target_size << " " << std::flush;
		fd_errors_avg << current_target_size << " " << std::flush;
		fd_errors_med << current_target_size << " " << std::flush;
		fd_errors_max << current_target_size << " " << std::flush;
		fd_recalls << current_target_size << " " << std::flush;
		fd_recalls_avg << current_target_size << " " << std::flush;
		fd_recalls_med << current_target_size << " " << std::flush;
		fd_recalls_min << current_target_size << " " << std::flush;
		fd_recalls_max << current_target_size << " " << std::flush;
		fd_precisions << current_target_size << " " << std::flush;
		fd_precisions_avg << current_target_size << " " << std::flush;
		fd_precisions_med << current_target_size << " " << std::flush;
		fd_precisions_min << current_target_size << " " << std::flush;
		
		
		
		// More Declarations -- specific to this particular target size
		std::vector<char> myhypothesis;
		std::vector<std::vector<char> > examples;
		std::vector<int> hyp_sizes_against_this_target_size;
		std::vector<long double> hyp_errors_against_this_target_size;
		std::vector<long double> hyp_recalls_against_this_target_size;
		std::vector<long double> hyp_precisions_against_this_target_size;
		
		
		
		for (int i = 0; i < NUM_EXPERIMENTS_PER_TARGET_SIZE; i++) {
			std::vector<char> myhypothesis;
			std::vector<char> mytarget = e.create_target(DIMENSION, current_target_size);
			assert(e.get_hypothesis_length(mytarget) == current_target_size);
			
			//std::cout << " drawing examples ..." << std::flush;
			std::cout << "\r       experiment " << i+1 << " out of " << NUM_EXPERIMENTS_PER_TARGET_SIZE << std::flush;
			/*
			 * Evolve a hypothesis using the algorithm
			 *
			 * Below we exploit the structure theorems that were proved
			 * in the following two papers:
			 * -- On Evolvability: The Swapping Algorithm, Product Distributions, and Covariance
			 *    (by Diochnos and Turan)
			 * -- On the Evolution of Monotone Conjunctions: Drilling for Best Approximations
			 *    (by Diochnos)
			 *
			 * This approach speeds up the process in some cases and allows us to
			 * get the simulation results faster since we know how the final hypothesis
			 * will look like, when the target has size at most (frontier+1).
			 */
			//std::cout << "   learning ..." << std::flush;
			if (e.get_hypothesis_length (mytarget) <= frontier) {
				// The hypothesis is precisely the target (no mistakes)
				myhypothesis = mytarget;
			} else if (e.get_hypothesis_length(mytarget) == (frontier + 1)) {
				// Learn a best q-approximation
				myhypothesis = mytarget; // init to be the same thing
				// Now delete a variable at random
				int to_delete = rand() % (frontier+1);
				int temp_counter = 0;
				for (int j = 0; j < myhypothesis.size(); j++) {
					if (static_cast<int> (myhypothesis[j]) == 1) {
						if (temp_counter == to_delete) {
							// we delete this variable
							myhypothesis[j] = static_cast<char> (0);
						}
						temp_counter++; // increase the counter for the observations
					}
				}
				assert (e.get_hypothesis_length(myhypothesis) == frontier);
			} else {
				// Now it is better if we run the algorithm
				// so that we can simulate the randomness and approximate the consistency
				// of good and variables on the q-approximation
				//
				
				// Initialize the hypothesis to be the empty monotone conjunction
				myhypothesis = e.create_EVO_initial_hypothesis (DIMENSION, frontier);
				
				// Evolve by simulating the correlation queries (by the closed formula)!
				for (int gen = 0; gen < num_generations; gen++) {
					std::vector<std::vector<char> > critical_nhood = e.get_neighborhood_of_interest (mytarget, myhypothesis, frontier, tolerance);
					// Pick a solution from the critical neighborhood at random
					// Technically, there are also weights,
					// but in the face of the amortization argument
					// there is no real difference
					int survivor_index = rand() % critical_nhood.size();
					myhypothesis = critical_nhood[survivor_index]; // new hypothesis
				}
				
				assert( e.get_hypothesis_length(myhypothesis) == frontier ); // sanity check
				// evolution has finished
			}
			
			// Rest of the work
			int current_hyp_size          = e.get_hypothesis_length(myhypothesis);
			long double current_risk      = e.get_risk      (mytarget, myhypothesis);
			long double current_recall    = e.get_recall    (mytarget, myhypothesis);
			long double current_precision = e.get_precision (mytarget, myhypothesis);
			hyp_sizes_against_this_target_size.push_back(current_hyp_size);
			hyp_errors_against_this_target_size.push_back(current_risk);
			hyp_recalls_against_this_target_size.push_back(current_recall);
			hyp_precisions_against_this_target_size.push_back(current_precision);
			
			// Write individual entries to files
			fd_sizes << current_hyp_size << " " << std::flush;
			fd_errors << current_risk << " " << std::flush;
			fd_recalls << current_recall << " " << std::flush;
			fd_precisions << current_precision << " " << std::flush;


			//std::cout << " done -- (hypothesis size: " << current_hyp_size << ", with risk: " << current_risk << ")" << std::endl;
			//std::cout << "I learned the following hypothesis (length " << e.get_hypothesis_length(myhypothesis) << "): " << std::flush;
			//for (int i = 0; i < myhypothesis.size(); i++) {
			//	std::cout << static_cast<int> (myhypothesis[i]) << " ";
			//}
			//std::cout << std::endl;
			mytarget.clear();
			myhypothesis.clear();
			examples.clear();
		}
		
		// First push them, then compute
		hypotheses_sizes.push_back(hyp_sizes_against_this_target_size);
		hypotheses_errors.push_back(hyp_errors_against_this_target_size);
		hypotheses_recalls.push_back(hyp_recalls_against_this_target_size);
		hypotheses_precisions.push_back(hyp_precisions_against_this_target_size);
		
		// Compute average and median hypothesis size against this target size
		current_median  = e.find_median_among_ints(hyp_sizes_against_this_target_size);
		current_average = e.find_average_among_ints(hyp_sizes_against_this_target_size);
		current_max     = e.find_max_among_ints(hyp_sizes_against_this_target_size);
		//
		hyp_sizes_med.push_back(current_median);
		hyp_sizes_avg.push_back(current_average);
		hyp_sizes_max.push_back(current_max);
		//
		fd_sizes_med << current_median << std::endl;
		fd_sizes_avg << current_average << std::endl;
		fd_sizes_max << current_max << std::endl;
		
		// Compute average, median and max hypothesis error (risk) against this target size
		current_median  = e.find_median_among_long_doubles(hyp_errors_against_this_target_size);
		current_average = e.find_average_among_long_doubles(hyp_errors_against_this_target_size);
		current_max     = e.find_max_among_long_doubles(hyp_errors_against_this_target_size);
		//
		hyp_errors_med.push_back(current_median);
		hyp_errors_avg.push_back(current_average);
		hyp_errors_max.push_back(current_max);
		//
		fd_errors_med << current_median << std::endl;
		fd_errors_avg << current_average << std::endl;
		fd_errors_max << current_max << std::endl;
		
		
		// Compute min, average and median hypothesis recall against this target size
		current_min     = e.find_min_among_long_doubles(hyp_recalls_against_this_target_size);
		current_median  = e.find_median_among_long_doubles(hyp_recalls_against_this_target_size);
		current_average = e.find_average_among_long_doubles(hyp_recalls_against_this_target_size);
		current_max     = e.find_max_among_long_doubles(hyp_recalls_against_this_target_size);
		//
		hyp_recalls_med.push_back(current_median);
		hyp_recalls_avg.push_back(current_average);
		hyp_recalls_min.push_back(current_min);
		hyp_recalls_max.push_back(current_max);
		//
		fd_recalls_med << current_median << std::endl;
		fd_recalls_avg << current_average << std::endl;
		fd_recalls_min << current_min << std::endl;
		fd_recalls_max << current_max << std::endl;
		
		// Compute min, average and median hypothesis precision against this target size
		current_min     = e.find_min_among_long_doubles(hyp_precisions_against_this_target_size);
		current_median  = e.find_median_among_long_doubles(hyp_precisions_against_this_target_size);
		current_average = e.find_average_among_long_doubles(hyp_precisions_against_this_target_size);
		//
		hyp_precisions_med.push_back(current_median);
		hyp_precisions_avg.push_back(current_average);
		hyp_precisions_min.push_back(current_min);
		//
		fd_precisions_med << current_median << std::endl;
		fd_precisions_avg << current_average << std::endl;
		fd_precisions_min << current_min << std::endl;
		
		std::cout << std::endl;
		std::cout << std::fixed << std::setprecision(8);
		std::cout << "       RISK     :   min = " << e.find_min_among_long_doubles(hyp_errors_against_this_target_size);
		std::cout << ",   med = " << e.find_median_among_long_doubles(hyp_errors_against_this_target_size);
		std::cout << ",   avg = " << e.find_average_among_long_doubles(hyp_errors_against_this_target_size);
		std::cout << ",   max = " << e.find_max_among_long_doubles(hyp_errors_against_this_target_size);
		std::cout << std::endl;
		std::cout << "       RECALL   :   min = " << e.find_min_among_long_doubles(hyp_recalls_against_this_target_size);
		std::cout << ",   med = " << e.find_median_among_long_doubles(hyp_recalls_against_this_target_size);
		std::cout << ",   avg = " << e.find_average_among_long_doubles(hyp_recalls_against_this_target_size);
		std::cout << ",   max = " << e.find_max_among_long_doubles(hyp_recalls_against_this_target_size);
		std::cout << std::endl;
		std::cout << "       PRECISION:" << "   min = " << e.find_min_among_long_doubles(hyp_precisions_against_this_target_size);
		std::cout << ",   med = " << e.find_median_among_long_doubles(hyp_precisions_against_this_target_size);
		std::cout << ",   avg = " << e.find_average_among_long_doubles(hyp_precisions_against_this_target_size);
		std::cout << ",   max = " << e.find_max_among_long_doubles(hyp_precisions_against_this_target_size);
		std::cout << std::endl;
		
		
		// Cleaning up ....
		hyp_sizes_against_this_target_size.clear();
		hyp_errors_against_this_target_size.clear();
		hyp_recalls_against_this_target_size.clear();
		hyp_precisions_against_this_target_size.clear();
		
		fd_sizes << std::endl;
		fd_errors << std::endl;
		fd_recalls << std::endl;
		fd_precisions << std::endl;
		std::cout << std::endl;
	}

	// Finalizing things
	fd_log << std::endl;
	t = std::time(0); now = std::localtime(&t);
	fd_log << "Experiments ended  : " << asctime(now) << std::endl;
	
	/*
	std::cout << "Hypotheses sizes" << std::endl;
	for (int i = 0; i < hyp_sizes_med.size(); i++) {
		std::cout << target_sizes[i] << ":" << std::flush;
		std::cout << " median: " << hyp_sizes_med[i] << ", average: " << hyp_sizes_avg[i] << std::endl;
	}
	std::cout << "Hypotheses errors" << std::endl;
	for (int i = 0; i < hypotheses_errors.size(); i++) {
		std::cout << target_sizes[i] << ":" << std::flush;
		std::cout << " median: " << hyp_errors_med[i] << ", average: " << hyp_errors_avg[i] << std::endl;
	}*/
	
	hypotheses_sizes.clear();
	hypotheses_errors.clear();
	hypotheses_recalls.clear();
	hypotheses_precisions.clear();
	
	hyp_sizes_avg.clear();
	hyp_sizes_med.clear();
	hyp_sizes_max.clear();
	hyp_errors_avg.clear();
	hyp_errors_med.clear();
	hyp_errors_max.clear();
	
	hyp_recalls_med.clear();
	hyp_recalls_avg.clear();
	hyp_recalls_min.clear();
	hyp_recalls_max.clear();
	hyp_precisions_med.clear();
	hyp_precisions_avg.clear();
	hyp_precisions_min.clear();
	
	
	// Close the files
	std::cout << "Closing files now" << std::endl;
	fd_sizes_avg.close();
	fd_sizes_med.close();
	fd_sizes_max.close();
	fd_sizes.close();
	fd_errors_avg.close();
	fd_errors_med.close();
	fd_errors_max.close();
	fd_errors.close();
	fd_recalls_min.close();
	fd_recalls_med.close();
	fd_recalls_avg.close();
	fd_recalls_max.close();
	fd_recalls.close();
	fd_precisions_min.close();
	fd_precisions_med.close();
	fd_precisions_avg.close();
	fd_precisions.close();
	
	// final timestamp
	t = std::time(0); now = std::localtime(&t);
	fd_log << "Entire experiment ended: " << asctime(now);
	fd_log.close();
}


/****************************************************************************************/
/****************************************************************************************/

int main () {
	srand(time(NULL));
	//driver();
	
	run_PAC_experiment();
	
	//run_EVO_experiment();
	
	return 0;
}
