//created by John Stanco 8.2.18


/// User defined headers
#include "mc_methods.hpp"
#include "signal.hpp"
#include "help.hpp"
#include "graph.hpp"


namespace Ising{

	template<class U>
	class lattice_ising{
	protected:
		/// Stores connectivity
		U lat_graph;
		/// Stores spins 
		arma::Col<int> lat;
		double E;
		double M;
		double T;
		double B;
		double J;
		double h;

		double dM;
		double dE;

		size_t iter;
		size_t n_smpl;

		Signal::signal E_dat;
		Signal::signal M_dat;

		void fill_lat();
		virtual void init();
		void reset();
		double nn_energy(const size_t vtx);
		void calc_energy();
		void calc_mag();
		void update_lat(){ E += dE; M += dM; }
		void update_obs(){ E_dat(iter) = E; M_dat(iter) = std::abs(M); }
		void flip_spin(const size_t idx){ lat(idx) *= -1; }
	public:
		lattice_ising(double J, double h, double T, std::vector<size_t> latsize):
		J(J), T(T), h(h), lat_graph(U(latsize)), E(0), M(0) { init(); } 

		virtual void update(bool burn = false) = 0;
		void prepare(const size_t len_of_run);

		Signal::signal E_data() const { return E_dat; }
		Signal::signal M_data() const { return M_dat; }

		size_t n_site() const { return lat_graph.n_vtx(); }
		void print_data() const {}
	};


	template<class U>
	void lattice_ising<U>::fill_lat(){
		for(size_t i = 0; i < lat.n_elem; i++){
			lat(i) = rand_int(0,1) * 2 - 1;
		}
	}


	template<class U>
	void lattice_ising<U>::init(){
		B = pow(1 * T, -1);
		lat = arma::Col<int>(lat_graph.n_vtx());
		fill_lat();
		calc_energy();
		calc_mag();
	}


	template<class U>
	void lattice_ising<U>::reset(){
		fill_lat();
		calc_energy();
		calc_mag();
	}


	template<class U>
	double lattice_ising<U>::nn_energy(const size_t vtx){
		double eps = 0;
		std::vector<size_t> edges = lat_graph.get_edges(vtx);
		for(auto& e : edges){
			size_t idx = e;
			eps += lat(vtx) * lat(idx);
		}
		return -eps * J;
	}


	template<class U>
	void lattice_ising<U>::calc_energy(){
		E = 0;
		for(size_t i = 0; i < lat_graph.n_vtx(); i++){
			E += nn_energy(i) / 2;
			if(h){ E -= h * lat(i); }
		}
	}


	template<class U>
	void lattice_ising<U>::calc_mag(){
		M = 0;
		for(size_t i = 0; i < lat.n_elem; i++){
			M += lat(i);
		}
	}


	template<class U>
	void lattice_ising<U>::prepare(const size_t len_of_run){
		n_smpl = len_of_run;
		iter = 0;
		E_dat = Signal::signal(n_smpl);
		M_dat = Signal::signal(n_smpl);
	}


	template<class U>
	class local_ising : public lattice_ising<U>{
	protected:
		bool accept(const double p){ if(p>=1){return 1;} return flip(p); }
		virtual double acc_prob(const double dE) = 0;
	public:
		local_ising(double J, double h, double T, std::vector<size_t> latsize):
		lattice_ising<U>(J,h,T,latsize){}
		/// single_update
		void update(bool burn = false);		
	};


	template<class U>
	void local_ising<U>::update(bool burn){
		double flip_idx = rand_int(this->lat.n_elem);
		this->flip_spin(flip_idx);
		this->dE = 2 * this->nn_energy(flip_idx);
		this->dM = 2 * this->lat(flip_idx);

		if(this->accept(this->acc_prob(this->dE))){
			this->update_lat();
		} else { this->flip_spin(flip_idx); }
		if(!burn){
			this->update_obs();
			this->iter++;
		}
	}	


	template <class U>
	class mh_ising : public local_ising<U>{
	protected:
		double acc_prob(const double dE){ if(dE < 0){return 1;} return exp(-this->B*dE); }
	public:
		mh_ising(double J, double h, double T, std::vector<size_t> latsize) :
		local_ising<U>(J,h,T,latsize){}
	};


	template <class U>
	class glauber_ising : public local_ising<U>{
	protected:
		double acc_prob(const double dE){ double p=exp(-this->B*dE); return p/(1+p); }
	public:
		glauber_ising(double J, double h, double T, std::vector<size_t> latsize) :
		local_ising<U>(J,h,T,latsize){}
	};


	template <class U>
	class cluster_ising : public lattice_ising<U>{
	protected:
		size_t n_same;
		size_t n_diff;
		double flip_prob;
		std::vector<size_t> cluster;
		arma::Col<int> is_checked;

		void init();
		virtual void calc_flip_prob(const double T, const double J)=0;
		void reset_cluster();
		void add_to_cluster(const size_t spin_index);
		int is_added(const size_t spin_index){ return(is_checked(spin_index) == -1); }
		void expand_cluster(const size_t spin_index);
		int flip_cluster(const size_t spin_index);
	public:
		cluster_ising(double J, double h, double T, std::vector<size_t> latsize) :
		lattice_ising<U>(J,h,T,latsize) { init(); }
		void update(bool burn = false);
		void set_temp(const double T);
	};


	template <class U>
	void cluster_ising<U>::init(){
		n_same = 0;
		n_diff = 0;
		is_checked = arma::zeros<arma::Col<int>>(this->lat.n_elem);
	}


	template <class U>
	void cluster_ising<U>::reset_cluster(){
		cluster.clear();
		n_diff = 0;
		n_same = 0;
		is_checked.zeros();
	}


	template <class U>
	void cluster_ising<U>::add_to_cluster(const size_t spin_index){
		n_same -= is_checked(spin_index);
		is_checked(spin_index) = -1;
		cluster.push_back(spin_index);
	}


	template <class U>
	void cluster_ising<U>::expand_cluster(const size_t spin_index){
		add_to_cluster(spin_index);
		std::vector<size_t> nn = this->lat_graph.get_edges(spin_index);
		for(size_t idx : nn) {
			if (!is_added(idx)) {
				if(this->lat(idx) == this->lat(spin_index)) {
					if(flip(this->flip_prob)){
						expand_cluster(idx); 
					} else { n_same++; is_checked(idx)++; }
				} else {
					n_diff++;
				}
			}
		}
	}


	template <class U>
	int cluster_ising<U>::flip_cluster(const size_t spin_index){
		expand_cluster(spin_index);
		int n_flipped = cluster.size();
		int m_n = n_same - n_diff;
		for (size_t i : cluster) { this->flip_spin(i); }
		this->dM = 2 * n_flipped * this->lat(spin_index);
		this->dE = 2 * this->J * m_n;
		this->dE -= this->dM * this->h;
		reset_cluster();
		return m_n;
	}


	template <class U>
	void cluster_ising<U>::update(bool burn){
		int flip_index = rand_int(this->lat.n_elem);
		flip_cluster(flip_index);
		this->update_lat();
		if (!burn) {
			this->update_obs();
			this->iter++;
		}
	}


	template <class U>
	class wolff_ising : public cluster_ising<U>{
	protected:
		void calc_flip_prob(const double T, const double J){ this->flip_prob = 1-exp(-2*J/(1*T)); }
	public:
		wolff_ising(double J, double h, double T, std::vector<size_t> latsize) :
		cluster_ising<U>(J,h,T,latsize){ calc_flip_prob(T,J); }
	};

	/*
	template<class U>
	class worm_ising : public ising_state_mh<U> {
	protected:
		size_t ira;
		size_t masha;
	public:
		worm_ising(const std::vector<uint> r_max, 
							const double T, 
							const double J, 
							const double h ) : 
		flip_prob( calc_flip_prob(T, J)),
		ising_state<U>(r_max, T, J, h) {}

		void set_temp(const double _T){
			this->T = _T;
			this->B = pow(1 * _T, -1);
		}

		void update_lat(){
			size_t site;
			if(ira==masha){
				ira = rand_int(lat.n_site);
				masha = ira;
			}

			size_t nn = rand_int(2 * lat.dim);
		}
	};
	*/

	typedef mh_ising<Graph::hypercubic_lattice_graph> 	square_mh_ising;
	typedef glauber_ising<Graph::hypercubic_lattice_graph> square_glauber_ising;
	typedef wolff_ising<Graph::hypercubic_lattice_graph>	square_wolff_ising;

}



