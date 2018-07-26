//created by John Stanco 7.1.18


#include "../include/ising_state.hpp"
#include "../include/markov_chain_mc.hpp"

template<class T1>
void print_ising_run(const size_t iter, 
				const size_t burn, 
				const size_t lag, 
				const double T, 
				const double J, 
				const double h,  
				const std::vector<uint> lat_spec,
				const double t,
				const T1 & init,
				const std::string method) {
	size_t dim = lat_spec.size();

	std::string mag;
	std::string size;
	if(J > 0){
		mag = "ferro";
	} else {
		mag = "antiferro";
	}

	std::cout << "\n--- ISING MCMC ---\n";
	std::cout << "Method:\t" << method << std::endl;
	std::cout << "Iter:\t" 	<< iter << std::endl;
	std::cout << "Burn:\t"	<< burn << std::endl;
	std::cout << "Lag:\t" 	<< lag 	<< std::endl;
	std::cout << "T:\t" 	<< T 	<< std::endl;
	std::cout << "J:\t" 	<< mag 	<< std::endl;
	std::cout << "h:\t" 	<< h 	<< std::endl;
	std::cout << "Size:\t";
	for (size_t i = 0; i < dim - 1; i++) {
		size += std::to_string(lat_spec[i]) + "x";
		std::cout << lat_spec[i] << "x";
	}
	size += std::to_string(lat_spec[dim-1]);
	std::cout << lat_spec[dim - 1]
				<< "\ntime:\t" 
				<< (float)t / CLOCKS_PER_SEC 
				<< " sec" << std::endl;
	printf("Mag:\t%lf\n", mean(init.mag()));
	printf("Energy:\t%lf\n", mean(init.energy()));
	std::cout << "-------------------\n";

	print(init.energy(),"../data/ising_energy_"+size+".dat");
	print(init.mag(),"../data/ising_mag_"+size+".dat");
	print(autocorr(init.mag()),"../data/ising_auto_"+size+".dat");
}


template<class T1>
void run_ising(const size_t iter, 
			const size_t burn, 
			const size_t lag, 
			const double T, 
			const double J, 
			const double h, 
			const std::vector<uint> lat_spec,
			const std::string method) {
	T1 init(lat_spec, T, J, h);
	markov_chain_mc<T1> mcmc;

	clock_t t = clock();
	T1 rslt = mcmc.run(init, iter);

	t = clock() - t;
	print_ising_run<T1>(iter, burn, lag, T, J, h, lat_spec, t, rslt, method);
}


static std::string switches[7] = { "-i", "-b", "-l", "-t", "-j", "-h", "-s" }; 


int display_usage() {
	std::cout << "usage: ising  [-b burn]  [-l lag]  [-j J]  [-h H]  -i <iter>  -t <temp>  -s <size_arr>\n";
	exit(0);
}


int display_no_arg(const std::string c) {
	std::cout << "Input Error: All switch identifiers must be followed by argument. ("<< c << ")\n";
	exit(0);
}


int display_incomplete() {
	std::cout << "Input Error: Must specify all arguments denoted with <...>\n";
	exit(0);
}


int is_switch(const std::string c) {
	for (size_t i = 0; i < 7; i++) {
		if( c == switches[ i ] ) { return 1; }
	}
	return 0;
}


void check_arg(int argc, char * argv[], const size_t i) {
	if (i + 1 >= argc || is_switch(argv[i + 1])){
		display_no_arg(argv[i]);
	}
}


int main(int argc, char *argv[]) {
	//mandatory arguments are iter, size, 
	if (argc == 1) {
		display_usage();
	} else {
		std::vector<uint> lat_spec;
		double J 	= 1;
		double h 	= 0;
		double T 	= 0;
		size_t iter = 0;
		size_t burn	= 0;
		size_t lag 	= 0;
		size_t n;
		/// for required arguments
		bool iter_flag = 0;
		bool temp_flag = 0;
		bool size_flag = 0;

		for (size_t i = 1; i < argc; i++) {
			std::string c = argv[ i ];
			if(c == "-i") {
				check_arg( argc, argv, i );
				sscanf( argv[ i + 1 ], "%zu", &iter );
				iter_flag = 1;
			} else if (c == "-b") {
				check_arg(argc, argv, i);
				sscanf(argv[i + 1], "%zu", &burn);
			} else if (c == "-l") {
				check_arg(argc, argv, i);
				sscanf(argv[i + 1], "%zu", &lag);
			} else if (c == "-t") {
				check_arg(argc, argv, i);
				sscanf(argv[ i + 1], "%lf", &T);
				temp_flag = 1;
			} else if (c == "-j") {
				check_arg(argc, argv, i);
				sscanf(argv[i + 1], "%lf", &J);
			} else if (c == "-h") {
				check_arg(argc, argv, i);
				sscanf(argv[i + 1], "%lf", &h);
			} else if (c == "-s") {
				check_arg( argc, argv, i );
				while (i + 1 < argc && !is_switch(argv[i + 1])) { //parse all non-switch arguments following this switch 
					sscanf( argv[ i + 1 ], "%zu", &n );
					lat_spec.push_back(n);
					i++;
					size_flag = 1;
				}
			} 
		}

		if (!(iter_flag && temp_flag && size_flag)) {
			display_incomplete();
		}

		run_ising<square_mh_ising>(iter, burn, lag, T, J, h, lat_spec, "MH");
	}
	return 1;
}