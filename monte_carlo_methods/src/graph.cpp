//created by John Stanco 7.17.18

#include "../include/Graph.hpp"



/////// UD UW GRAPH ////////


void Graph::uw_ud_graph::compute_adj_list(const arma::mat& W){
	adj_list = std::vector<node>(n_vert);
	for(size_t u = 0; u < n_vert; u++) {
		for(size_t v = 0; v < n_vert; v++) {
			if(W(u, v)) { add_edge(u, v); }
		}
	}
}


/// explores connected component containing vertex u in depth-first fashion
void Graph::uw_ud_graph::explore_dfs(const size_t u) {
	previsit(u);
	for (size_t v : adj_list[u].edges){
		if (!is_searched(v)) { explore_dfs(v); }
	}
	postvisit(u);
}


void Graph::uw_ud_graph::explore_dfs(const size_t u, const size_t n){
	previsit(u);
	nn.push_back(u);
	if(n<1) { return; }
	for(size_t v : adj_list[u].edges){
		if (!is_searched(v)){ explore_dfs(v, n-1); }
	}
	postvisit(u);
}


/// explores connected component containing vertex u in breadth-first fashion
void Graph::uw_ud_graph::explore_bfs(const size_t n){
	q.push(n);
	previsit(n);
	while (!q.empty()) {
		size_t u = q.front();
		for (size_t v : adj_list[u].edges) {
			if (!is_searched(v)) { q.push(v); previsit(v); }
		}
		q.pop();
		postvisit(u);
	}
}


/// Depth-first search
void Graph::uw_ud_graph::dfs(){
	dfs_init();
	for(size_t v = 0; v < n_vert; v++){
		if(!is_searched(v)) { explore_dfs(v); }
	}
	dfs_end();
}


std::vector<size_t> Graph::uw_ud_graph::dfs(const size_t u){
	dfs_init();
	explore_dfs(u);
	dfs_end();
	return nn;
}


/// Bread-first search
void Graph::uw_ud_graph::bfs(){
	bfs_init();
	for(size_t v = 0; v < n_vert; v++){
		if (!is_searched(v)) { explore_bfs(v); }
	}
	bfs_end();
}


std::vector<size_t> Graph::uw_ud_graph::bfs(const size_t u){
	bfs_init();
	explore_bfs(u);
	bfs_end();
	return nn;
}


/// Computes connected components
void Graph::uw_ud_graph::compute_n_cc(){
	dfs_init();
	cc = 0;
	for(size_t v = 0; v < n_vert; v++){
		if(!is_searched(v)){ explore_dfs(v); cc++; }
	}
	dfs_end();
}


std::vector<size_t> Graph::uw_ud_graph::get_nn_upto(const size_t u, const size_t n){
	dfs_init();
	explore_dfs(u,n);
	dfs_end();
	return nn;
}


/////// W D GRAPH ///////


/// compute adjacency list from matrix
void Graph::w_d_graph::compute_adj_list(const arma::mat& W){
	adj_list = std::vector<node>(n_vert);
	for(size_t u = 0; u < n_vert; u++) {
		for(size_t v = 0; v < n_vert; v++) {
			if(W(u, v)) { add_edge(u, v, W(u,v)); }
		}
	}
}


/// explores connected component containing vertex u in depth-first fashion
void Graph::w_d_graph::explore_dfs( const size_t u ) {
	previsit(u);
	for (auto& v : adj_list[u].edges) {
		if (!is_searched(v.first)) { explore_dfs(v.first); }
	}
	postvisit(u);
}


/// explores connected component containing vertex u in breadth-first fashion
void Graph::w_d_graph::explore_bfs(const size_t n) {
	q.push(n);
	previsit(n);
	while (!q.empty()) {
		size_t u = q.front();
		for (auto& v : adj_list[u].edges) {
			if (!is_searched(v.first)) { q.push(v.first); previsit(v.first); }
		}
		q.pop();
		postvisit(u);
	}
}


/// Depth-first search
void Graph::w_d_graph::dfs() {
	dfs_init();
	for (size_t u = 0; u < n_vert; u++) {
		if (!is_searched(u)) { explore_dfs(u); }
	}
}


/// Bread-first search
void Graph::w_d_graph::bfs() {
	bfs_init();
	for (size_t u = 0; u < n_vert; u++) {
		if (!is_searched(u)) { explore_bfs(u); }
	}
	return bfs_end();
}


/// Computes connected components
void Graph::w_d_graph::compute_n_cc() {
	dfs_init();
	cc = 0;
	for (size_t u = 0; u < n_vert; u++) {
		if (!is_searched(u)){ explore_dfs(u); cc++; }
	}
	dfs_end();
}


/// implemements Dijkstra's shortest-path algorithm
void Graph::w_d_graph::dijkstra(const size_t& source) {
	dijkstra_init();
	arma::vec dist(n_vert);
	arma::vec prev(n_vert);
	auto cmp = [](my_pair_t l, my_pair_t r){ return l.second < r.second; };
	std::priority_queue<my_pair_t, my_container_t, decltype(cmp)> p_q(cmp);
	dist(source) = 0;
	prev(source) = source;
	p_q.push(std::make_pair(source, 0));

	for (size_t u = 1; u < n_vert; u++) {
		dist(u) = std::numeric_limits<double>::max();
	}

	while (!p_q.empty()) {
		size_t u = p_q.top().first;
		p_q.pop();
		for (auto& v : adj_list[u].edges) {
			double alt = dist(u) + v.second;
			if (alt < dist(v.first)) {
				dist(v.first) = alt;
				prev(v.first) = u;
				p_q.push(v);
			}
		}
	}
	//dist.print();
	//std::cout << std::endl;
	//prev.print();
	bfs_end();
}


/*
void Graph::w_d_graph::dfs_spec_depth(const size_t u, const size_t n){
	for(auto& v : adj_list[u].vs){
		//dfs_spec_depth
	}
}
*/



my_container_t Graph::w_d_graph::get_nn_upto(const size_t v, const size_t n){
	my_container_t nn;
	//dfs_spec_depth(v, n);
	return nn;
}



////// SIMPLE GRAPH //////


/// returns degree of Graph vertex u by summing over all of its connections
double Graph::simple_graph::row_sum(const arma::mat& wt_mtx, const size_t& u) {
	double sum = 0;
	for(size_t v = 0; v < wt_mtx.n_cols; v++) {
		sum += wt_mtx(u, v);
	}
	return sum;
}


/// compute degree matrix from weights
arma::mat Graph::simple_graph::compute_degree_mat(const arma::mat& W) {
	arma::mat D = arma::zeros<arma::mat>(n_vert, n_vert);
	for(size_t u = 0; u < n_vert; u++) {
		D(u, u) = row_sum(W, u);
	}
	return D;
}


void Graph::simple_graph::compute_eigs() {
	lpn_basis = arma::mat(n_vert, n_vert);
	lpn_eigvals = arma::vec(n_vert);
	arma::eig_sym(lpn_eigvals, lpn_basis, L);
	lpn_basis_inv = inv(lpn_basis);
}


/// computes the p-Dirichlet form - generalizes 2-form
double Graph::simple_graph::p_dirichlet_form(const graphSignal& f, const size_t p) {
	if (p == 2) { return two_form(f); }
	double p_form = 0;
	for (size_t u = 0; u < n_vert; u++) {
		double p_norm = 0;
		for (auto& v : adj_list[u].edges) {
			p_norm += pow(f(u) - f(v.first), 2) * v.second;
		}
		p_norm = pow(p_norm, p / 2);
		p_form += p_norm;
	}
	return p_form / (double)p;
}


///////////  LATTICES  ////////////


arma::uvec Graph::hypercubic_lattice_graph::to_vec(const std::vector<size_t>& x){
	arma::uvec y(x.size());
	for(size_t i = 0; i < x.size(); i++){
		y(i) = x[i];
	}
	return y;
}


void Graph::hypercubic_lattice_graph::compute_stride(){
	stride_len = arma::uvec(dim);
	stride_len(dim - 1) = 1;
	for(int i = dim - 2; i >= 0; i--){
		stride_len(i) = stride_len(i + 1) * rmax(i + 1);
	}
}


arma::uvec Graph::hypercubic_lattice_graph::to_vec_indices(const size_t index){
  size_t y_i = 0;
  arma::uvec xs(dim);
  for(int i = dim - 1; i > 0; i--){
   	xs(i) = ((index - y_i) / stride_len(i)) % rmax(i);
   	y_i += xs(i) * stride_len(i);
  }
  xs(0) = ((index - y_i) / stride_len(0));
  return xs;
}	


void Graph::hypercubic_lattice_graph::compute_adj_list(){
	adj_list = std::vector<node>(n_vert);
	size_t n_neighbors = 2*dim;
	size_t v;
	for(size_t u = 0; u < n_vert; u++){
		arma::uvec indices = to_vec_indices(u);
		for(size_t i = 0; i < n_neighbors; i += 2){
			size_t j = i / 2;
			size_t k = indices(j);
			v = u + ((k + 1) % rmax(j) - k) * stride_len(j);
			add_edge(u,v);
			v = u + ((k + rmax(j) - 1) % rmax(j) - k)*stride_len(j);
			add_edge(u,v);
		}
	}
}


void Graph::hypercubic_lattice_graph::init(){
	n_vert = prod(rmax);
	compute_stride();
	compute_adj_list();
}


void Graph::triangular_lattice_graph::compute_adj_list(){
	size_t dim = 2;
}