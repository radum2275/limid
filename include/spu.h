/*
 * spu.h
 *
 *  Created on: 5 Jul 2016
 *      Author: radu
 *
 * Copyright (c) 2016, International Business Machines Corporation
 * and University of California Irvine. All rights reserved.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 * DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
 * FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
 * DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
 * SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
 * CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
 * OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 * OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */

/// \file spu.h
/// \brief Single Policy Update (SPU) algorithm for LIMIDs
/// \author Radu Marinescu

#ifndef IBM_MERLIN_SPU_H_
#define IBM_MERLIN_SPU_H_

#include "limid.h"
#include "algorithm.h"

namespace merlin {

/**
 * Single Policy Update (SPU)
 *
 * Models supported: LIMID
 *
 *
 */
class spu : public limid, public algorithm {
public:
	typedef limid::findex findex;        ///< Factor index
	typedef limid::vindex vindex;        ///< Variable index
	typedef limid::flist flist;          ///< Collection of factor indices

	typedef std::pair<factor, factor> potential;	///< Potential (phi, psi)

	typedef struct {
		size_t from;			///< Id of the 'from' cluster
		size_t to;				///< Id of the 'to' cluster
		variable_set sep;		///< Separator
		potential forward;		///< Forward message from->to
		potential backward;		///< Backward message to->from
	} edge;

	///
	/// \brief Properties of the algorithm
	///
	MER_ENUM( Property , Order,Iter,Debug );

public:

	///
	/// \brief Default constructor.
	///
	spu() : limid() {
		set_properties();
	}

	///
	/// \brief Constructor with a graphical model.
	///
	spu(const limid& lm) : limid(lm), m_gmo(lm) {
		clear_factors();
		set_properties();
	}

	///
	/// \brief Clone the algorithm.
	/// \return the pointer to the new object containing the cloned algorithm.
	///
	virtual spu* clone() const {
		spu* lm = new spu(*this);
		return lm;
	}

	// Can be an optimization algorithm or a summation algorithm....
	double ub() const {
		throw std::runtime_error("Not implemented");
	}
	double lb() const {
		throw std::runtime_error("Not implemented");
	}
	std::vector<size_t> best_config() const {
		throw std::runtime_error("Not implemented");
	}

	double logZ() const {
		throw std::runtime_error("Not implemented");
	}
	double logZub() const {
		throw std::runtime_error("Not implemented");
	}
	double logZlb() const {
		throw std::runtime_error("Not implemented");
	}

	// No beliefs defined currently
	const factor& belief(size_t f) const {
		throw std::runtime_error("Not implemented");
	}
	const factor& belief(variable v) const {
		throw std::runtime_error("Not implemented");
	}
	const factor& belief(variable_set vs) const {
		throw std::runtime_error("Not implemented");
	}
	const vector<factor>& beliefs() const {
		throw std::runtime_error("Not implemented");
	}

	///
	/// \brief Set the properties of the algorithm.
	/// \param opt 	The string containing comma separated property value pairs
	///
	virtual void set_properties(std::string opt = std::string()) {
		if (opt.length() == 0) {
			set_properties("Order=MinFill,Iter=1,Debug=1");
			return;
		}
		std::vector<std::string> strs = merlin::split(opt, ',');
		for (size_t i = 0; i < strs.size(); ++i) {
			std::vector<std::string> asgn = merlin::split(strs[i], '=');
			switch (Property(asgn[0].c_str())) {
			case Property::Order:
				m_order.clear();
				m_order_method = OrderMethod(asgn[1].c_str());
				break;
			case Property::Iter:
				m_iterations = (atol(asgn[1].c_str()) == 0) ? false : true;
				break;
			case Property::Debug:
				m_debug = (atol(asgn[1].c_str()) == 0) ? false : true;
				break;
			default:
				break;
			}
		}
	}

	///
	/// \brief Eliminate a set of variables from a factor.
	/// \param F 	The reference of the factor to eliminate from
	///	\param vs 	The set of variables to be eliminated
	/// \param w 	The weight of the weighted elimination operator
	/// \return the factor resulted from eliminating the set of variables.
	///
	factor elim(const factor& F, const variable_set& vs) {
		return F.sum(vs);
	}

	///
	/// \brief Compute the marginal over a set of variables.
	/// \param F 	The reference of the factor to marginalize over
	/// \param vs 	The set of variables representing the scope of the marginal
	/// \param w 	The weight of the weighted elimination operator
	/// \return the factor representing the weighted marginal over the set of variables.
	///
	factor marg(const factor& F, const variable_set& vs) {
		return F.marginal(vs);
	}

	///
	/// \brief Initialize the bucket tree structure.
	///
	void init() {

		// Start the timer and store it
		m_start_time = timeSystem();

		// Check if LIMID
		bool is_limid = m_gmo.islimid();
		if (!is_limid) {
			throw std::runtime_error("SPU is only supported for LIMIDs.");
		}

		// Prologue
		std::cout << "Initialize solver ..." << std::endl;
		std::cout << " + models supported : LIMID" << std::endl;
		std::cout << " + algorithm        : SPU" << std::endl;

		// Construct the elimination ordering
		m_order = m_gmo.order(m_order_method);

		// Get the induced width of the order
		size_t wstar = m_gmo.induced_width(m_order);
		std::cout << " + elimination (u)  : ";
		std::copy(m_order.begin(), m_order.end(),
				std::ostream_iterator<size_t>(std::cout, " "));
		std::cout << std::endl;
		std::cout << " + induced width    : " << wstar << std::endl;
		assert(m_porder.empty() == true);
		std::cout << " + partial order    : N/A" << std::endl;
		std::cout << " + policy updates   :";
		m_decisions.clear();
		for (variable_order_t::reverse_iterator ri = m_order.rbegin();
				ri != m_order.rend(); ++ri) {
			if (m_vtypes[*ri] == 'd') {
				std::cout << " " << *ri;
				m_decisions.push_back(*ri); // build the ordered list of decisions
											// required for the iterative SPU updates
			}
		}
		std::cout << std::endl;

		// Initialize original variable and factor types
		m_vtypes = m_gmo.get_vtypes();
		m_ftypes = m_gmo.get_ftypes();

		// Get the original factors scopes
		vector<variable_set> fin;
		for (vector<factor>::const_iterator i = m_gmo.get_factors().begin();
				i != m_gmo.get_factors().end(); ++i) {
			fin.push_back((*i).vars());
		}

		// Mark factors depending on variable i
		vector<flist> vin;
		for (size_t i = 0; i < m_gmo.nvar(); ++i) {
			vin.push_back(m_gmo.with_variable(var(i)));
		}

		vector<flist> Orig(m_gmo.num_factors()); // origination info: which original factors are
		for (size_t i = 0; i < Orig.size(); ++i) {
			Orig[i] |= i;    					// included for the first time, and which newly
		}
		vector<flist> New(m_gmo.num_factors()); // created clusters feed into this cluster

		// First downward pass to initialize the bucket tree and backward messages
		std::cout << "Initializing bucket-tree ... " << std::endl;
		std::cout << "  - initial number of clique factors is: " << m_factors.size() << std::endl;
		m_clusters.resize(m_order.size()); // init clusters
		m_ctypes.resize(m_order.size(), false); // init cluster types
		for (variable_order_t::const_iterator x = m_order.begin(); x != m_order.end(); ++x) {

			std::cout << "  - create cluster for var "
				<< *x << " (" << m_vtypes[*x] << ")\n" ;

			// Get the current variable to process
			variable VX = var(*x);
			if (*x >= vin.size() || vin[*x].size() == 0)
				continue;  // check that we have some factors over this variable

			flist ids = vin[*x];  // list of factor IDs contained in this bucket

			if (m_debug) {
				std::cout << "  - factors in this bucket: " << ids.size() << std::endl;
				for (flist::const_iterator i = ids.begin(); i != ids.end(); ++i) {
					std::cout << "   factor id " << *i << " : " << fin[*i] << std::endl;
				}
			}

			// Merge all factor scopes in this bucket into a single scope
			vector<size_t> temp;
			typedef flist::const_iterator flistIt;
			size_t jj = *(ids.begin());
			for (flistIt i = ids.begin(); i != ids.end(); ++i) {
				size_t ii = *i;
				if (ii == jj) continue;
				fin[jj] |= fin[ii];     // combine into j
				erase(vin, ii, fin[ii]);
				fin[ii] = variable_set();  //   & remove i

				Orig[jj] |= Orig[ii];
				Orig[ii].clear(); // keep track of list of original factors in this cluster
				New[jj] |= New[ii];
				New[ii].clear(); //  list of new "message" clusters incoming to this cluster

				temp.push_back(ii);
			}

			for (size_t i = 0; i < temp.size(); ++i) {
				size_t ii = temp[i];
				ids /= ii;
			}

			// Sanity checks
			assert(ids.size() == 1);
			if (m_debug) {
				std::cout << "  After merging: " << ids.size() << std::endl;
				for (flist::const_iterator i = ids.begin(); i != ids.end(); ++i) {
					std::cout << "  Factor id " << *i << std::endl;
					std::cout << "  Scope: " << fin[*i] << std::endl;
				}
			}

			// Eliminate the bucket variable
			vector<findex> alphas;
			for (flistIt i = ids.begin(); i != ids.end(); ++i) {
				//
				// Create new cluster alpha over this set of variables; save function parameters also
				findex alpha = findex(-1);
				alpha = add_factor(factor(fin[*i])); // add the clique factor
				alphas.push_back(alpha);

				fin[*i] = fin[*i] - VX;

				// Add inter clusters edges
				for (flistIt j = New[*i].begin(); j != New[*i].end(); ++j) {
					add_edge(*j, alpha);
					m_schedule.push_back(std::make_pair(*j, alpha));
				}

				// Keep track of original factors
				m_originals.push_back(flist());
				m_originals[alpha] |= Orig[*i];
				m_clusters[*x] |= alpha; // map variable to bucket/cluster
				m_cluster2var[alpha] = *x; // map cluster to variable

				// Now incoming nodes to *i is just alpha
				Orig[*i].clear();
				New[*i].clear();
				New[*i] |= alpha;

				// Recompute and update adjacency
				insert(vin, *i, fin[*i]);
			}
		}
		// end for: variable elim order
		std::cout << "  - final number of clique factors is: " << m_factors.size() << std::endl;
		std::cout << "Done initializing the bucket-tree." << std::endl;

		// Set the separators and cluster scopes
		size_t max_sep_size = 0, max_clique_size = 0;
		size_t C = m_factors.size(); // number of clique factors
		m_separators.resize(C);
		for (size_t i = 0; i < C; ++i) m_separators[i].resize(C);
		m_scopes.resize(C); // clique scopes
		for (size_t i = 0; i < C; ++i) {
			m_scopes[i] = m_factors[i].vars();
			max_clique_size = std::max(max_clique_size, m_scopes[i].size());
			for (flist::const_iterator j = m_originals[i].begin();
					j != m_originals[i].end(); ++j) {
				if (m_gmo.get_factor(*j).get_type() == factor::FactorType::Decision) {
					m_ctypes[i] = true;
				}
			}
		}
		const std::vector<edge_id>& elist = edges();
		for (size_t i = 0; i < elist.size(); ++i) {
			findex a,b;
			a = elist[i].first;
			b = elist[i].second;
			if (a > b) continue;
			variable_set sep = m_factors[a].vars() & m_factors[b].vars();
			m_separators[a][b] = sep;
			m_separators[b][a] = sep;
			max_sep_size = std::max(max_sep_size, sep.size());
		}

		// Set the incoming and outgoing lists for each cluster/clique
		m_in.resize(C);
		m_out.resize(C);
		for (vector<std::pair<findex, findex> >::const_iterator i = m_schedule.begin();
				i != m_schedule.end(); ++i) {
			findex from = (*i).first;
			findex to = (*i).second;
			m_in[to] |= from;
			m_out[from] |= to;
		}

		// Initialize the root cluster(s)
		for (size_t i = 0; i < m_out.size(); ++i) {
			if ( m_out[i].empty() )
				m_roots |= i;
		}

		// Init forward and backward messages
		size_t N = m_schedule.size();
		m_forward.resize(N);
		m_backward.resize(N);
		m_edge_indeces.resize(C);
		for (size_t i = 0; i < C; ++i) m_edge_indeces[i].resize(C);
		for (size_t i = 0; i < N; ++i) {
			size_t from = m_schedule[i].first;
			size_t to = m_schedule[i].second;
			m_edge_indeces[from][to] = i;
		}

		// Init clique potentials as dummy potentials
		for (size_t c = 0; c < m_factors.size(); ++c) {
			m_factors[c] = factor(1.0); // init
		}

		// Initialize maximum expected utility
		m_meu = -infty(); // negative infinity (maximizing)

		// Output the bucket tree statistics
		std::cout << "Created bucket-tree with " << std::endl;
		std::cout << " - number of cliques:  " << C << std::endl;
		std::cout << " - number of edges:    " << elist.size() << std::endl;
		std::cout << " - max clique size:    " << max_clique_size << std::endl;
		std::cout << " - max separator size: " << max_sep_size << std::endl;
		std::cout << std::endl;

		// Now, we must redirect the edges in the bucket tree to the new root
		// which corresponds to the first decision variable in the order
		size_t d = m_decisions[0];
		m_root = m_clusters[d][0];

		// Updates the schedule so that messages are collected at the new root
		update_schedule(m_root);

		// Debug information
		if (m_debug) {
			std::cout << "[MERLIN DEBUG]\n";
			std::cout << "[DBG] Bucket-tree with " << m_factors.size() << " clusters and "
					<< elist.size() << " edges" << std::endl;
			for (size_t i = 0; i < elist.size(); ++i) {
				findex a,b;
				a = elist[i].first;
				b = elist[i].second;
				if (a > b) continue;
				std::cout << "  edge from "
						<< m_scopes[a] << " to "
						<< m_scopes[b] << " (a=" << a << ", b=" << b << ")"
						<< " sep: " << m_separators[a][b]
						<< std::endl;
			}

			std::cout << "[DBG] Forward propagation schedule:" << std::endl;
			for (size_t i = 0; i < m_schedule.size(); ++i) {
				std::cout << " msg " << m_schedule[i].first << " --> "
						<< m_schedule[i].second << std::endl;
			}
			std::cout << "[DBG] Backward propagation schedule:" << std::endl;
			vector<std::pair<findex, findex> >::reverse_iterator ri = m_schedule.rbegin();
			for (; ri != m_schedule.rend(); ++ri) {
				std::cout << " msg " << ri->second << " --> "
						<< ri->first << std::endl;
			}

			std::cout << "[DBG] Original factors per cluster:" << std::endl;
			for (size_t i = 0; i < m_originals.size(); ++i) {
				if (m_ctypes[i]) std::cout << "*cl " << i << " : ";
				else std::cout << " cl " << i << " : ";
				std::copy(m_originals[i].begin(), m_originals[i].end(),
						std::ostream_iterator<int>(std::cout, " "));
				std::cout << std::endl;
			}

			std::cout << "[DBG] Original factors:" << std::endl;
			for (size_t i = 0; i < m_gmo.num_factors(); ++i) {
				std::cout << " " << i << " (" << m_ftypes[i] << ") : "
					<< m_gmo.get_factor(i) << std::endl;
			}

			// _in and _out lists
			std::cout << "[DBG] IN list:" << std::endl;
			for (size_t i = 0; i < m_in.size(); ++i) {
				std::cout << "  in[" << i << "] = ";
				std::copy(m_in[i].begin(), m_in[i].end(),
						std::ostream_iterator<int>(std::cout, " "));
				std::cout << std::endl;
			}
			std::cout << "[DBG] OUT list:" << std::endl;
			for (size_t i = 0; i < m_out.size(); ++i) {
				std::cout << "  out[" << i << "] = ";
				std::copy(m_out[i].begin(), m_out[i].end(),
						std::ostream_iterator<int>(std::cout, " "));
				std::cout << std::endl;
			}
			std::cout << "[DBG] ROOTS: ";
			std::copy(m_roots.begin(), m_roots.end(),
					std::ostream_iterator<int>(std::cout, " "));
			std::cout << std::endl;

			// clique factors, forward and backward
			std::cout << "[DBG] clique_factors:" << std::endl;
			for (size_t i = 0; i < m_factors.size(); ++i) {
				std::cout << "[" << i << "]: " << m_factors[i] << std::endl;
			}
			std::cout << "[DBG] forward messages (top-down):" << std::endl;
			for (size_t i = 0; i < m_forward.size(); ++i) {
				std::cout << "(" << i << "): " << m_forward[i].first
					<< " | " << m_forward[i].second << std::endl;
			}
			std::cout << "[DBG] backward messages (bottop-up):" << std::endl;
			for (size_t i = 0; i < m_backward.size(); ++i) {
				std::cout << "(" << i << "): " << m_backward[i].first
					<< " | " << m_backward[i].second << std::endl;
			}
			std::cout << "[MERLIN DEBUG]\n";
		} // end if debug
	}

	///
	/// \brief Update the adjacencies with a new factor
	///
//	void insert(std::vector<flist>& adj, findex fid, const factor& f,
//			vector<vindex>::const_iterator x,
//			const variable_order_t& ord) {
//		for (vector<vindex>::const_iterator y = x+1;
//				y != ord.end(); ++y) {
//			variable VY = var(*y);
//			if (f.vars().contains(VY)) {
//				adj[*y] |= fid;
//				break;
//			}
//		}
//	}


	///
	/// \brief Run the SPU algorithm for LIMIDs.
	///
	virtual void run() {

		// Initialize the algorithm
		init();

		// SPU propagation
		propagate();

		// Output solution (UAI output format)
		std::cout << "Finished in " << (timeSystem() - m_start_time)
				<< " seconds" << std::endl;
	}

	///
	/// \brief Propagate the messages along the edges of the bucket tree.
	///
	/// The bucket tree is a particular kind of join trees in the sense that
	/// it doesn't enforce its cliques to be marginal and each variable is
	/// associated with exactly one clique/bucket/cluster. We also assume that
	/// the bucket tree is connected, namely there is exactly one root cluster.
	/// And the root is associated with the first decision variable in the
	/// ordering (again, we assume that the decision policies will be updated
	/// in a given ordering e.g., corresponding to the eliminatin order).
	///
	void propagate() {

		std::cout << "Begin message passing along the bucket tree ..." << std::endl;

		collect(m_root);
		std::cout << "  Finished initial collect to the root." << std::endl;
		std::cout << "  Initial random strategy has " << m_meu << " expected utility" << std::endl;
		for (size_t iter = 1; iter <= m_iterations; ++iter) {

			std::cout << " SPU iteration: " << iter << std::endl;
			size_t n = m_decisions.size();
			for (size_t i = 0; i < n; ++i) {
				size_t di = m_decisions[i];
				size_t j = (i == n-1 ? 0 : i+1);
				size_t dj = m_decisions[j];

				std::cout << "   update policy for " << di << " (next is " << dj << ")\n";

				update(di); // update the policy of decision 'di', and report EU
				propagate(di, dj); // update the messages between cliques di and dj
			}
		}
	}

	///
	/// \brief Propagate messages along the (unique) path between two clusters
	///	\param from 	The id of the source cluster
	/// \param to		The id of the target cluster
	///
	void propagate(size_t from, size_t to) {

	}

	///
	/// \brief Compute the incoming potential of a cluster
	/// \param a 	The index of the cluster
	/// \return the potential i.e., probability and expected utility of cluster *a*
	///
	potential incoming(findex a) {

		// Collect the probability, policy and utility factors
		factor phi = factor(1.0), psi = factor(0.0);
		for (flist::const_iterator j = m_originals[a].begin();
				j != m_originals[a].end(); ++j) {
			const factor& f = m_gmo.get_factor(*j);
			if (f.get_type() == factor::FactorType::Probability ||
				f.get_type() == factor::FactorType::Decision) {
				phi *= f;
			} else if (f.get_type() == factor::FactorType::Utility) {
				psi += f;
			}
		}

		// Collect forward messages to 'a'
		for (flist::const_iterator ci = m_in[a].begin();
				ci != m_in[a].end(); ++ci) {
			findex p = (*ci);
			size_t j = m_edge_indeces[p][a];
			phi *= m_forward[j].first;
			psi += m_forward[j].second;
		}

		return std::make_pair(phi, psi);
	}

	///
	/// \brief Redirect the edges in the bucket tree towards a new root
	///	\param cl 	The cluster id of the new root
	///
	///	Also updates the propagation schedule towards the new root
	///
	void update_schedule(size_t cl) {

		// Select the root and redirect edges outwards
		std::vector<int> visited;
		std::stack<int> dfs;
		dfs.push(cl);
		visited.resize(m_clusters.size(), false);

		while ( !dfs.empty() ) {

			int n = dfs.top();
			dfs.pop();

			// Get all unvisited neighbors of 'n' and direct them towards 'n'
			for (flist::const_iterator i = m_neighbors[n].begin();
					i != m_neighbors[n].end(); ++i) {

				size_t m = *i;
				if (visited[m] == false) {
					dfs.push(m);
					edge e;
					e.from = m;
					e.to = n;
					e.forward = std::make_pair(factor(1), factor(0));
					e.backward = std::make_pair(factor(1), factor(0));
					m_jointree.push_back(e);
					m_parents[m] = n;
					m_children[n] |= m;
				}

			}

			// Mark node 'n' as visited
			visited[n] = true;
		}

		// find the message order (schedule) by a bfs of the junction tree
		std::queue<size_t> bfs;
		bfs.push(m_root);
		while (!bfs.empty()) {
			size_t n = bfs.front();
			bfs.pop();

			size_t p = m_parents[n];
			if (p >= 0) {
				flist elist = p->edges();
				for (flist::const_iterator j = elist.begin();
						j != elist.end(); ++j) {
					const edge& e = m_jointree[*j];
					if (e.from == n && e.to == p) {
						m_schedule.push_back(*j); // the edge id
						break;
					}
				}
			}

			for (flist::const_iterator i = m_children[n].begin();
					i != m_children[n].end(); ++i) {
				bfs.push(*i);
			}
		}

		// reverse the order to start from the leaves
		std::reverse(m_schedule.begin(), m_schedule.end());
	}

	///
	/// \brief Collect all forward and backward messages toa cluster.
	/// \param a 	The index of the cluster
	/// \return the potential i.e., probability and utility factors of the cluster.
	///
	potential collect(findex a) {

		// Collect the original probability, policy and utility factors
		factor phi = factor(1.0), psi = factor(0.0);
		for (flist::const_iterator j = m_originals[a].begin();
				j != m_originals[a].end(); ++j) {
			const factor& f = m_gmo.get_factor(*j);
			if (f.get_type() == factor::FactorType::Probability ||
				f.get_type() == factor::FactorType::Decision) {
				phi *= f;
			} else if (f.get_type() == factor::FactorType::Utility) {
				psi += f;
			}
		}

		// Collect the forward messages to 'a'
		for (flist::const_iterator ci = m_in[a].begin();
				ci != m_in[a].end(); ++ci) {
			findex p = (*ci);
			size_t j = m_edge_indeces[p][a];
			phi *= m_forward[j].first;
			psi += m_forward[j].second;
		}

		// Collect the backward message to 'a'
		for (flist::const_iterator ci = m_out[a].begin();
				ci != m_out[a].end(); ++ci) {
			findex p = (*ci);
			size_t j = m_edge_indeces[a][p];
			phi *= m_backward[j].first;
			psi += m_backward[j].second;
		}

		return std::make_pair(phi, psi);
	}

	///
	/// \brief Retract the policy from its cluster
	///
	factor retract(findex a) {
		for (flist::const_iterator j = m_originals[a].begin();
				j != m_originals[a].end(); ++j) {
			const factor& f = m_gmo.get_factor(*j);
			if (f.get_type() == factor::FactorType::Decision) {
				return f;
			}
		}

		return factor(1.0); // dummy
	}

	///
	/// \brief Update the policy by overwriting the corresponding factor
	///
	void replace(findex a, const factor& d) {
		for (flist::const_iterator j = m_originals[a].begin();
				j != m_originals[a].end(); ++j) {
			const factor& f = m_gmo.get_factor(*j);
			if (f.get_type() == factor::FactorType::Decision) {
				m_gmo.set_factor(*j, d);
				break;
			}
		}
	}


	///
	/// \brief Forward (top-down) message passing along the edge of the bucket tree.
	///
	/// Updates all the forward messages (potentials) and computes the
	/// expected utility of the current strategy
	///
	void forward() {

		std::cout << "Begin forward (top-down) pass ..." << std::endl;

		double eu = -infty(); // init the expected utility
		for (variable_order_t::const_iterator x = m_order.begin(); x != m_order.end(); ++x) {

			std::cout << " - Eliminating var " << *x << " (" << m_vtypes[*x] << ")\n";

			// Generate forward messages
			variable VX = var(*x);
			findex a = m_clusters[*x][0]; // get source bucket of the variable
			if (m_out[a].size() > 0) {
				findex b = *(m_out[a].begin()); // destination bucket
				size_t ei = m_edge_indeces[a][b]; // edge index (for message)

				// Compute the probability factor
				potential pot = incoming(a);
				m_forward[ei].first = pot.first.sum(VX);
				// Compute the expected utility factor
				m_forward[ei].second = (pot.first * pot.second).sum(VX);
				m_forward[ei].second /= m_forward[ei].first;

				if (m_debug) {
					std::cout << "  forward msg (" << a << "," << b << "): elim = " << VX << "\n";
					std::cout << "    P: " << m_forward[ei].first << std::endl;
					std::cout << "    U: " << m_forward[ei].second << std::endl;
				}
			}

		} // done

		// Compute expected utility of the current strategy
		factor F(0.0);
		for (flist::const_iterator ci = m_roots.begin();
				ci != m_roots.end(); ++ci) {

			// Get the variable associated with the root cluster
			std::map<size_t, size_t>::iterator mi = m_cluster2var.find(*ci);
			assert(mi != m_cluster2var.end());
			size_t v = mi->second;
			variable VX = var(v);

			// Eliminate the root cluster variable (must be only one)
			potential pot = incoming(*ci);
			factor prob = pot.first.sum(VX);
			factor util = (pot.first*pot.second).sum(VX);
			util /= prob;

			F += (prob*util); // must be a constant
		}

		// Partition function or MAP/MMAP value
		eu = F.max();
		m_meu = std::max(eu, m_meu); // update the MEU

		std::cout << "Finished forward pass with EU: " << eu << std::endl;
	}

	///
	/// \brief Backward (bottom-up) message passing along the edges of the bucket tree.
	///
	void backward() {

		// Calculate 'backward' messages
		std::cout << "Begin backward (bottom-up) pass ..." << std::endl;

		vector<std::pair<findex, findex> >::reverse_iterator ri = m_schedule.rbegin();
		for (; ri != m_schedule.rend(); ++ri ) {

			// Select backward message m(b->a)
			findex a = (*ri).first;
			findex b = (*ri).second;
			size_t ei = m_edge_indeces[a][b]; // edge index

			// Update the policy if decision cluster; update MEU as well
			update(b);

			std::cout << " - Sending backward msg from " << b << " to " << a << std::endl;

			// Select variables to eliminate
			variable_set VX = m_scopes[b] - m_separators[a][b];

			// Collect the potential at "b"
			potential pot = collect(b);
			// Discard forward message on the edge
			pot.first /= m_forward[ei].first;
			pot.second -= m_forward[ei].second;

			// Update the backward message
			m_backward[ei].first = pot.first.sum(VX);
			// Compute the expected utility factor
			m_backward[ei].second = (pot.first * pot.second).sum(VX);
			m_backward[ei].second /= m_backward[ei].first;

			if (m_debug) {
				std::cout << "  - backward msg (" << b << "," << a << "): elim = " << VX << std::endl;
				std::cout << "      P: " << m_backward[ei].first << std::endl;
				std::cout << "      U: " << m_backward[ei].second << std::endl;
			}
		}

		std::cout << "Finished backward pass." << std::endl;
	}

	///
	/// \brief Update the policy of a decision variable.
	///
	/// \param a The cluster where the decision variable resides
	///
	void update(vindex a) {

		// Check if current cluster is a policy one
		if (m_ctypes[a] == false) return;
		std::cout << "Updating policy for cluster: " << a << std::endl;

		// Get the variable associated with the cluster/bucket
		std::map<size_t, size_t>::iterator mi = m_cluster2var.find(a);
		assert(mi != m_cluster2var.end());
		size_t v = mi->second;
		variable VX = var(v); // decision variable
		assert(m_vtypes[v] == true); // sanity checks

		// Retract the current policy, collect messages, compute contraction
		factor d = retract(a);
		assert(d.vars().contains(VX)); // sanity checks
		potential pot = collect(a);
		pot.first /= d; // do the actual retraction
		factor cont = (pot.first * pot.second).marginal(d.vars());

		// Compute the skewed distribution, ie. for each configuration of
		// the decision variable's parents, get the argmax in the corresponding
		// sliced/conditioned factor and set it as the degenerate distribution
		// in the target policy (for that particular configuration)
		factor nd = d; // to have the same scope
		size_t n = d.vars().size();
		variable_set scope = d.vars();
		for (size_t e = 0; e < d.numel(); ++e) {

			// Get the scope configuration correp. to the current index (big_endian)
			int i = e;
			size_t p = 0;
			std::vector<size_t> I(n); // configuration of source variables
			for (variable_set::const_iterator v = d.vars().begin();
					v != d.vars().end(); ++v) {
				I[p] = i % v->states();
				i -= I[p];
				i /= v->states();
				p++;
			}

			// SLice the factor on current parent configuration
			size_t old;
			factor tmp = d;
			p = 0;
			for (variable_set::const_iterator v = d.vars().begin();
					v != d.vars().end(); ++v) {
				if (*v != VX) {
					tmp = tmp.condition(*v, I[p]);
				} else {
					old = I[p];
				}

				p++;
			}

			// Find the optimal decision
			size_t opt = tmp.argmax();
			if (opt == old) {
				nd[e] = 1;
			} else {
				nd[e] = 0;
			}
		}

		// Update the policy
		replace(a, nd);

		// Compute the expected utility of the current strategy
		potential p = collect(a);
		factor F = (p.first*p.second).sum(m_scopes[a]);
		double eu = F.max();
		m_meu = std::max(eu, m_meu);

		std::cout << "Finished update step with EU: " << eu << std::endl;
	}

protected:
	// Members:

	limid m_gmo; 						///< Original influence diagram
	double m_meu;						///< Log maximum expected utility
	std::map<vindex, factor> m_policy;	///< Optimal decision policy
	OrderMethod m_order_method;			///< Variable ordering method
	variable_order_t m_order;			///< Variable order
	bool m_debug;						///< Internal debugging flag
	size_t m_iterations;				///< Number of SPU iterations (default 1)

private:
	// Bucket Tree local structures:

	vector<bool> m_ctypes;				///< The type of each cluster (true if contains policy factor, false otherwise)
	vector<flist> m_clusters;			///< Clusters corresponding to each of the variables
	vector<flist> m_originals;			///< Original factors (index) for each cluster
	vector<variable_set> m_scopes;		///< The scope (vars) for each cluster
	vector<flist> m_in;					///< Incoming messages (source clusters) to each cluster
	vector<flist> m_out; 				///< Outgoing messages (target clusters) from each cluster
	flist m_roots;						///< Root cluster(s)
	vector<potential> m_forward; 		///< Forward messages (indexed by edge id in the bucket tree)
	vector<potential> m_backward; 		///< Backward messages (indexed by edge id in the bucket tree)

	vector<std::pair<findex, findex> > m_schedule;	///< Propagation schedule
	vector<vector<size_t> > m_edge_indeces;			///< Edge indeces
	vector<vector<variable_set> > m_separators; 	///< Separators between clusters

	std::map<size_t, size_t> m_cluster2var;			///< Maps cluster id to a variable id

	std::vector<vindex> m_decisions;	///< Ordered list of decisions
	std::vector<int> m_parents;			///< Parent for each cluster in the bucket tree
	std::vector<flist> m_children;		///< Children of each cluster in the bucket tree
	size_t m_root;						///< Root of the bucket tree (co-incides with the first decision)
	std::vector<flist> m_neighbors;		///< For each cluster, list its neighbors in the bucket tree
	std::vector<edge> m_jointree;		///< Edges in the bucket tree / join tree
};


} // end namespace


#endif /* IBM_MERLIN_SPU_H_ */
