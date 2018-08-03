//created by John Stanco 7.14.18

/// Std lib headers
#include <queue>
#include <functional>
#include <vector>
#include <limits>

/// User defined headers
#include "help.hpp"
#include "timer.hpp"
#include "signal.hpp"

#ifndef GRAPH_H
#define GRAPH_H

typedef arma::vec graphSignal;
typedef std::pair<size_t, double> my_pair_t;
typedef std::vector<my_pair_t> my_container_t;


namespace Graph{

	enum search_id{
		SEARCH,
		ADD_NN
	};

	class uw_ud_graph {
	public:
		struct node{
			bool is_searched;
			std::vector<size_t> edges;
			void add_edge(const size_t idx){ edges.push_back(idx); }
			node() : is_searched(false) {}
		};
	protected:
		size_t n_edge;
		size_t n_vert;
		int cc;
		/// Adjacency list
		std::vector<node> adj_list;
		std::vector<size_t> nn;
		std::queue<size_t> q;
		/// compute degree matrix from weights
		void compute_adj_list(const arma::mat& W);
		void dfs_init() { reset_search(); nn.clear(); }
		void bfs_init() { reset_search(); nn.clear(); q=std::queue<size_t>(); }
		void dijkstra_init() { reset_search(); nn.clear(); }
		void dfs_end() {}
		void bfs_end() {}
		void previsit(const size_t u) { search(adj_list[u]); nn.push_back(u); }
		void postvisit(const size_t u) {}
		/// explores connected component containing vertex u in depth-first fashion
		void explore_dfs(const size_t u);
		void explore_dfs(const size_t u, const size_t n);
		/// explores connected component containing vertex u in breadth-first fashion
		void explore_bfs(const size_t u);
		/// Computes connected components
		void compute_n_cc();
		void search(node& u) { u.is_searched = true; }
		void reset(node& u){ u.is_searched = false; }
		void reset_search() { for(auto& v : adj_list){ reset(v); } }
		bool is_searched(const size_t v) { return adj_list[v].is_searched; }
	public:
		uw_ud_graph() : cc(-1), n_edge(0), n_vert(0) {}
		uw_ud_graph(const size_t nvtx) : cc(0), n_edge(0), n_vert(nvtx) { adj_list=std::vector<node>(n_vert); }
		uw_ud_graph(const arma::mat W) : n_vert(W.n_cols), cc(-1) { compute_adj_list(W); }
		/// Depth-first search - returns connected component of 
		void dfs();
		std::vector<size_t> dfs(const size_t);
		/// Bread-first search
		void bfs();
		std::vector<size_t> bfs(const size_t);
		/// TODO - UNION find
		/// returns number of connected components
		const int n_cc() { if (cc < 0) { compute_n_cc(); } return cc; }
		void add_vtx(){ adj_list.push_back(node()); n_vert++; }
		void add_edge(const size_t u, const size_t v){ adj_list[u].add_edge(v); n_edge++; }
		std::vector<size_t> get_edges(const size_t v) const { return adj_list[v].edges; }
		/// returns set of all vertices within 'n' edges of vertex v
		std::vector<size_t> get_nn_upto(const size_t v, const size_t n);
		/// TODO - returns set of only n-nn of vertex v
		/*std::vector<size_t> get_nn_only(const size_t v, const size_t n);*/
		size_t n_vtx() const { return n_vert; }
	};

	


	class w_d_graph {
	public:
		struct node{
			bool is_searched;
			my_container_t edges;
			void add_edge(const size_t idx, const double wt){ edges.push_back(std::make_pair(idx, wt)); }
			node() : is_searched(false) {}
		};
	protected:
		size_t n_edge;
		size_t n_vert;
		int cc;
		std::vector<node> adj_list;
		std::vector<std::pair<size_t, double> > nn;
		std::queue<size_t> q;
		/// compute degree matrix from weights
		void compute_adj_list(const arma::mat& W);
		void dfs_init() { reset_search(); }
		void bfs_init() { reset_search(); q = std::queue<size_t>(); }
		void dijkstra_init() {	reset_search(); }
		void dfs_end() {}
		void bfs_end() {}
		void previsit(const size_t u) { search(adj_list[u]); }
		void postvisit(const size_t u) {}
		/// explores connected component containing vertex u in depth-first fashion
		void explore_dfs(const size_t u);
		void explore_dfs(const size_t u, const size_t n);
		/// explores connected component containing vertex u in breadth-first fashion
		void explore_bfs(const size_t n);
		/// Computes connected components
		void compute_n_cc();
		/// Mark node as searched
		void search(node& u) { u.is_searched = true; }
		void reset(node& u) { u.is_searched = false; }
		void reset_search() { for(auto& v : adj_list){ reset(v); } }
		bool is_searched(const size_t u){ return adj_list[u].is_searched; }
	public:
		w_d_graph() : cc(-1), n_edge(0), n_vert(0) {}
		w_d_graph(const size_t nvtx) : cc(0), n_edge(0), n_vert(nvtx) { adj_list=std::vector<node>(n_vert); }
		w_d_graph(const arma::mat W) : n_vert(W.n_cols), cc(-1) { compute_adj_list(W); }
		/// Depth-first search
		void dfs();

		/// Bread-first search
		void bfs();
		/// Dijkstra's shortest-path algorithm
		void dijkstra(const size_t& source);
		/// TODO - UNION find
		/// returns number of connected components
		const int n_cc() { if (cc < 0) { compute_n_cc(); } return cc; }
		void add_vtx(){ adj_list.push_back(node()); n_vert++; }
		void add_edge(const size_t u, const size_t v, const double wt){ adj_list[u].add_edge(v,wt); n_edge++; }
		my_container_t get_edges(const size_t v) const { return adj_list[v].edges; }
		/// returns set of all vertices within 'n' edges of vertex v
		my_container_t get_nn_upto(const size_t v, const size_t n);
		/// TODO - returns set of only n-nn of vertex v
		/*my_container_t get_nn_only(const size_t v, const size_t n);*/
		size_t n_vtx() const {return n_vert;}
	};


	/// simple Graph - undirected, no self-edges
	/// Implements Graph signal processing methods
	/// see https://arxiv.org/pdf/1211.0053.pdf

	class simple_graph : public w_d_graph {
	protected:
		/// Laplacian
		arma::mat L;						
		/// Laplacian eigenvectors
		arma::mat lpn_basis;		
		arma::mat lpn_basis_inv;
		/// Laplacian eigenvalues
		arma::vec lpn_eigvals;				

		/// returns degree of Graph vertex u by summing over all of its connections
		double row_sum( const arma::mat& wt_mtx, const size_t& u );
		arma::mat compute_degree_mat( const arma::mat& W );
		void compute_lpn_mtx( const arma::mat& W ) { L = compute_degree_mat(W) - W; }
		void compute_eigs();
		/// computes the Dirichlet 2-form -- functions as semi-norm for simple Graph
		double two_form(const graphSignal f) { return dot( f, L * f ); }
	public:
		simple_graph(const arma::mat W) : w_d_graph(W) { compute_lpn_mtx(W); compute_eigs(); }
		/// produces Graph signal in basis of Laplacian eigenvectors
		graphSignal FT(const graphSignal& f) { return lpn_basis_inv * f; }
		/// returns Graph signal in Cartesian basis
		graphSignal IFT(const graphSignal& f) { return lpn_basis * f; }
		/// computes the p-Dirichlet form - generalizes 2-form
		double p_dirichlet_form(const graphSignal& f, const size_t p );
		void print_mat(const arma::mat& mtx) { mtx.print(); std::cout << "\n"; }
		void print() {  print_mat(L); }
	};


	class bravais_lattice_graph : public uw_ud_graph{
		/// this will be used for the pimc simulation
		/// can define this in both 2d or 3d
	};


	class hypercubic_lattice_graph : public uw_ud_graph{
		/// this is used for the ising simulation
	protected:
		size_t dim;
		arma::uvec rmax;
		arma::uvec stride_len;

		arma::uvec to_vec(std::vector<size_t> const& x);
		void compute_stride();
		arma::uvec to_vec_indices(const size_t index);
		void compute_adj_list();
		void init();
	public:
		hypercubic_lattice_graph(std::vector<size_t> latsize):
		dim(latsize.size()), rmax(to_vec(latsize)) { init(); }
	};

	class triangular_lattice_graph : public uw_ud_graph{
	protected:
		arma::uvec rmax;
		arma::uvec stride_len;

		arma::uvec to_vec(std::vector<size_t> const& x);
		void compute_adj_list();
	public:
	};


	class kagome_lattice_graph : public uw_ud_graph{
	protected:
	public:
	};
}


#endif /* GRAPH_H */