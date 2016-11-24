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
 * This is an implementation of the classic SPU algorithm for computing the
 * strategy that maximizes the expected utility for Limited Memory Influence
 * Diagrams (LIMIDs). However, it isn't an exact algorithm, namely it will only
 * reach a local maxima. There are however special cases for which SPU will find
 * the optimal strategy, but in general it doesn't guarantee optimality.
 *
 * TODO: We assume that the LIMID graph is connected, and the decision policies
 *       involve EXACTLY one decision variable (but several chance variables)!
 *
 * For computing optimal policies in LIMIDs, there is an extension called
 * Multiple Policy Update (MPU) which will be implemented soon.
 *
 *
 */
class spu : public limid, public algorithm {
public:
	typedef limid::findex findex;        ///< Factor index
	typedef limid::vindex vindex;        ///< Variable index
	typedef limid::flist flist;          ///< Collection of factor indices

	typedef std::pair<factor, factor> potential;					///< Potential (phi, psi)
	typedef vector<std::pair<vindex, vindex> > schedule;			///< Schedule
	typedef std::map<vindex, std::pair<size_t, findex> > policy;	///< Policy
	typedef vector<vector<size_t> > matrix;

	typedef struct jointree_edge_t {
		size_t from;			///< Id of the 'from' cluster
		size_t to;				///< Id of the 'to' cluster
		variable_set separator;	///< Separator
		potential message;		///< Message on the edge
		jointree_edge_t() : from(0), to(0) {};
		jointree_edge_t(const jointree_edge_t& e) : from(e.from), to(e.to),
			separator(e.separator), message(e.message) {};
		jointree_edge_t(size_t a, size_t b, const variable_set& sep,
			const potential& pot) : from(a), to(b), separator(sep), message(pot) {};
	} edge;

	///
	/// \brief Properties of the algorithm
	///
	MER_ENUM( Property , Order,Iter,Verbose,PureStrategy,InitPureStrategy );

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
			set_properties("Order=MinFill,Iter=5,Verbose=0,PureStrategy=1,InitPureStrategy=0");
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
				m_iterations = atol(asgn[1].c_str());
				break;
			case Property::Verbose:
				m_verbosity = atol(asgn[1].c_str());
				break;
			case Property::PureStrategy:
				m_pure_strategy = (atol(asgn[1].c_str()) == 0) ? false : true;
				break;
			case Property::InitPureStrategy:
				m_init_pure_strategy = (atol(asgn[1].c_str()) == 0) ? false : true;
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
		std::cout << " + verbosity        : " << m_verbosity << std::endl;

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
		std::cout << " + iterations       : " << m_iterations << std::endl;

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
		if (m_verbosity > 0) std::cout << "Initializing bucket-tree ... " << std::endl;
		m_clusters.resize(m_order.size()); // init clusters
		m_ctypes.resize(m_order.size(), false); // init cluster types
		for (variable_order_t::const_iterator x = m_order.begin(); x != m_order.end(); ++x) {

			if (m_verbosity >= 1) {
				std::cout << "  - create cluster for var "
					<< *x << " (" << m_vtypes[*x] << ")\n" ;
			}

			// Get the current variable to process
			variable VX = var(*x);
			if (*x >= vin.size() || vin[*x].size() == 0)
				continue;  // check that we have some factors over this variable

			flist ids = vin[*x];  // list of factor IDs contained in this bucket

			if (m_verbosity >= 2) {
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
			if (m_verbosity >= 2) {
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
				}

				// Keep track of original factors
				m_originals.push_back(flist());
				m_originals[alpha] |= Orig[*i];
				m_clusters[*x] |= alpha; // map variable to bucket/cluster

				// Now incoming nodes to *i is just alpha
				Orig[*i].clear();
				New[*i].clear();
				New[*i] |= alpha;

				// Recompute and update adjacency
				insert(vin, *i, fin[*i]);
			}
		}
		// end for: variable elim order
		if (m_verbosity > 0) std::cout << "Done initializing the bucket-tree." << std::endl;

		// Set the separators and cluster scopes
		size_t max_sep_size = 0, max_clique_size = 0;
		size_t C = m_factors.size(); // number of clique factors
		m_scopes.resize(C); // clique scopes
		m_num_clusters = C; // set the number of clusters
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

		// Map each of the decision variables to a decision policy (factor)
		for (size_t i = 0; i < m_decisions.size(); ++i) {
			vindex d = m_decisions[i];
			int cluster = -1, delta = -1;
			for (size_t c = 0; c < m_num_clusters; ++c) {
				if (m_ctypes[c] == true &&
					m_scopes[c].contains(var(d))) {
					cluster = c; // index of decision cluster
					for (flist::const_iterator j = m_originals[c].begin();
							j != m_originals[c].end(); ++j) {
						const factor& f = m_gmo.get_factor(*j);
						if (f.get_type() == factor::FactorType::Decision &&
							f.vars().contains(var(d))) {
							delta = *j; // index of decision policy factor
							break;
						}
					}
					break;
				}
			}

			assert(cluster >= 0 && delta >= 0);
			m_policy[d] = std::make_pair(cluster, delta);
		}

		// Init the join-tree edges and separators
		m_neighbors.resize(C);
		m_edge_indeces.resize(C);
		for (size_t i = 0; i < C; ++i) {
			m_edge_indeces[i].resize(C);
		}
		potential pot = std::make_pair(factor(1), factor(0));
		const std::vector<edge_id>& elist = edges();
		for (size_t i = 0; i < elist.size(); ++i) {
			findex a,b;
			a = elist[i].first;
			b = elist[i].second;
			variable_set sep = m_factors[a].vars() & m_factors[b].vars();
			max_sep_size = std::max(max_sep_size, sep.size());

			m_jointree.push_back(edge(a, b, sep, pot));
			m_edge_indeces[a][b] = i;
			m_neighbors[a] |= b;
			m_neighbors[b] |= a; // keep track of the adjacenies (undirected)
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

		// Select the first decision (ordered) as the new root of the join tree
		size_t d = m_decisions[0]; 	// first decision variable
		m_root = m_policy[d].first;	// and its corresp. cluster
		std::cout << "First decision in order is " << d << std::endl;
		std::cout << "Selected cluster " << m_root << " as root" << std::endl;
		std::cout << "Redirected join tree edges to new root" << std::endl;

		// Redirect edges to the new root, and update the message schedule
		redirect(m_root);

		// Initialize with either a pure or stochastic strategy
		if (m_init_pure_strategy) {
			for (size_t i = 0; i < m_decisions.size(); ++i) {
				vindex d = m_decisions[i];
				findex p = m_policy[d].second;
				factor delta = m_gmo.get_factor(p);
				factor ndelta = init_policy(d, 0, true, delta);
				m_gmo.set_factor(p, ndelta); // replace the factor
			}
		}

		// Debug information
		if (m_verbosity >= 2) {
			std::cout << "Bucket-tree with " << m_factors.size() << " clusters and "
					<< elist.size() << " edges" << std::endl;
			for (size_t ei = 0; ei < elist.size(); ++ei) {
				findex a,b;
				a = elist[ei].first;
				b = elist[ei].second;
				if (a > b) continue;
				std::cout << "  edge from "
						<< m_scopes[a] << " to "
						<< m_scopes[b] << " (a=" << a << ", b=" << b << ")"
						<< " sep: " << m_jointree[ei].separator
						<< std::endl;
			}

			std::cout << "Forward propagation schedule:" << std::endl;
			for (size_t i = 0; i < m_schedule.size(); ++i) {
				std::cout << " msg " << m_schedule[i].first << " --> "
						<< m_schedule[i].second << std::endl;
			}

			std::cout << "Original factors per cluster:" << std::endl;
			for (size_t i = 0; i < m_originals.size(); ++i) {
				if (m_ctypes[i]) std::cout << "*cl " << i << " : ";
				else std::cout << " cl " << i << " : ";
				std::copy(m_originals[i].begin(), m_originals[i].end(),
						std::ostream_iterator<int>(std::cout, " "));
				std::cout << std::endl;
			}

			std::cout << "Original factors:" << std::endl;
			for (size_t i = 0; i < m_gmo.num_factors(); ++i) {
				std::cout << " " << i << " (" << m_ftypes[i] << ") : "
					<< m_gmo.get_factor(i) << std::endl;
			}

			// Display the *parents* and *children* lists
			std::cout << "Join tree PARENTS:" << std::endl;
			for (size_t i = 0; i < m_parents.size(); ++i) {
				std::cout << "  parents[" << i << "] = " << m_parents[i] << std::endl;
			}
			std::cout << "Join tree CHILDREN:" << std::endl;
			for (size_t i = 0; i < m_children.size(); ++i) {
				std::cout << "  children[" << i << "] = ";
				std::copy(m_children[i].begin(), m_children[i].end(),
						std::ostream_iterator<int>(std::cout, " "));
				std::cout << std::endl;
			}

			std::cout << "Join tree ROOT: " << m_root << " (assume connected graph)" << std::endl;
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

		std::cout << "Begin SPU updates along the join tree ..." << std::endl;

		// Evaluate initial policy
		m_meu = collect();
		double prev_eu = -infty();

		if (m_verbosity == 0) {
			std::cout << "  [SPU] init " << m_meu << std::endl;
		}

		for (size_t iter = 1; iter <= m_iterations; ++iter) {

			if (m_verbosity >= 1) {
				std::cout << "  Begin SPU iteration " << iter << " ..." << std::endl;
			}

			size_t n = m_decisions.size();
			for (size_t i = 0; i < n; ++i) {
				size_t di = m_decisions[i];
				size_t j = (i == n-1 ? 0 : i+1);
				size_t dj = m_decisions[j];

				update_policy(di); // update the policy of decision 'di', and report EU

				if (m_verbosity == 0) {
					std::cout << "  [SPU] " << iter << " D" << di << " : " << m_meu << std::endl;
				}


				propagate(di, dj); // update the messages between cliques di and dj
			}

			if (m_verbosity >= 1) {
				std::cout << "  Finished SPU iteration in "
					<< (timeSystem() - m_start_time) << " seconds" << std::endl;
			}

			// Check if converged
			if (m_meu <= prev_eu) {
				std::cout << "Converged after " << iter << " iterations" << std::endl;
				break; // converged
			}

			prev_eu = m_meu;
		}

		// Output solution (UAI output format)
		std::cout << "Finished SPU updates in " << (timeSystem() - m_start_time)
				<< " seconds" << std::endl;
	}

	schedule find_path(size_t from, size_t to) {

		// Find the LCA of the two clusters in the join tree
		std::list<size_t> path1, path2;
		path1.push_front(from);
		int n = m_parents[from];
		while (n != -1) {
			path1.push_front(n);
			n = m_parents[n];
		}

		path2.push_front(to);
		n = m_parents[to];
		while (n != -1) {
			path2.push_front(n);
			n = m_parents[n];
		}

		// Find the last element that's common to the two lists
		int lca = -1;
		std::list<size_t>::iterator it1 = path1.begin();
		std::list<size_t>::iterator it2 = path2.begin();
		while (it1 != path1.end() && it2 != path2.end()) {
			if (*it1 == *it2) {
				lca = *it1;
				++it1;
				++it2;
			} else {
				break;
			}
		}

		// Assemble the path (the actual message schedule)
		assert(lca != -1);
		schedule path;
		size_t c = from;
		while ( c != lca ) {
			int p = m_parents[c];
			path.push_back(std::make_pair(c, p));
			c = p;
		}

		c = to;
		schedule temp;
		while ( c != lca ) {
			int p = m_parents[c];
			temp.push_back(std::make_pair(p, c));
			c = p;
		}
		std::reverse(temp.begin(), temp.end());

		// Update the path
		std::copy(temp.begin(), temp.end(), std::back_inserter(path));

		return path;
	}

	///
	/// \brief Propagate messages between the corresponding clusters of two decisions
	///	\param di 	The source variable index
	/// \param dj	The target variable index
	///
	void propagate(vindex di, vindex dj) {

		// Get the corresponding clusters in the join tree
		findex from = m_policy[di].first;
		findex to = m_policy[dj].first;

		// Check if same clusters
		if (from == to) {
			return; // nothing to propagate
		}

		if (m_verbosity >= 1) {
			std::cout << "   [PROPAGATE] messages between clusters: "
				<< from << " and " << to << std::endl;
		}

		// Find unique path between the two clusters
		schedule path = find_path(from, to);
		if (m_verbosity >= 2) {
			std::cout << "      Path found:";
			for (schedule::const_iterator i = path.begin(); i != path.end(); ++i) {
				std::cout << " (" << i->first << " --> " << i->second << ")";
			}
			std::cout << std::endl;
		}

		// Re-compute the messages along the path
		for (schedule::const_iterator i = path.begin(); i != path.end(); ++i) {

			// Select the message m(a->b)
			findex a = i->first;
			findex b = i->second;
			size_t ei = m_edge_indeces[a][b]; // edge index
			m_jointree[ei].message = message(a, b);

			if (m_verbosity >= 2) {
				std::cout << "    - Sending message from " << a << " to " << b
					<< " m(" << a << "," << b << "): elim = " << (m_scopes[a] - m_jointree[ei].separator) << std::endl;
				std::cout << "        P: " << m_jointree[ei].message.first << std::endl;
				std::cout << "        U: " << m_jointree[ei].message.second << std::endl;
			}
		}

		if (m_verbosity >= 1) {
			std::cout << "    Finished propagation in "
				<< (timeSystem() - m_start_time) << " seconds" << std::endl;
		}
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

		// Collect the messages from the neighbors of 'a'
		for (flist::const_iterator i = m_neighbors[a].begin();
				i != m_neighbors[a].end(); ++i) {
			findex b = (*i);
			size_t ej = m_edge_indeces[b][a];
			phi *= m_jointree[ej].message.first;
			psi += m_jointree[ej].message.second;
		}

		return std::make_pair(phi, psi);
	}

	///
	/// \brief Redirect the edges in the bucket tree towards a new root
	///	\param cl 	The cluster id of the new root
	///
	///	Also updates the propagation schedule towards the new root
	///
	void redirect(size_t root) {

		// Clear in/out lists, and schedule
		m_schedule.clear();
		m_parents.clear();
		m_children.clear();

		// Select the root and redirect edges outwards
		std::vector<int> visited;
		std::stack<int> dfs;
		dfs.push(root);
		visited.resize(m_num_clusters, false);
		m_parents.resize(m_num_clusters, -1);	// reset parents
		m_children.resize(m_num_clusters, flist()); // reset children lists
		while ( !dfs.empty() ) {

			int n = dfs.top();
			dfs.pop();

			// Get all unvisited neighbors of 'n' and direct them towards 'n'
			for (flist::const_iterator i = m_neighbors[n].begin();
					i != m_neighbors[n].end(); ++i) {

				size_t m = *i;
				if (visited[m] == false) {
					dfs.push(m);
					m_parents[m] = n;
					m_children[n] |= m;
				}
			}

			// Mark node 'n' as visited
			visited[n] = true;
		}

		// Find the message propagation order (schedule) by a BFS over the join tree
		std::queue<size_t> bfs;
		bfs.push(m_root);
		while (!bfs.empty()) {
			size_t n = bfs.front();
			bfs.pop();

			for (flist::const_iterator i = m_children[n].begin();
					i != m_children[n].end(); ++i) {
				size_t c = *i;
				size_t ei = m_edge_indeces[c][n];
				assert(m_jointree[ei].from == c && m_jointree[ei].to == n);
				m_schedule.push_back(std::make_pair(c, n));
				bfs.push(c);
			}
		}

		// reverse the order to start from the leaves
		std::reverse(m_schedule.begin(), m_schedule.end());
	}

	///
	/// \brief Send message from cluter *a* to adjacent cluster *b*
	/// \param a 	The source cluster index
	///	\param b	The target cluster index
	///
	potential message(size_t a, size_t b) {

		// Sanity checks
		assert(a != b);

		// Select the message m(a->b)
		size_t ei = m_edge_indeces[a][b]; // edge index

		// Get the incoming messages to 'a' (except for 'b')
		potential pot = incoming(a);

		// Discard message on the edge (b->a)
		size_t bi = m_edge_indeces[b][a];
		pot.first /= m_jointree[bi].message.first;
		pot.second -= m_jointree[bi].message.second;

		// Select variables to eliminate
		variable_set VX = m_scopes[a] - m_jointree[ei].separator;

		// Update the message
		potential result;
		result.first = pot.first.sum(VX);
		result.second = (pot.first * pot.second).sum(VX);
		result.second /= result.first;

		return result;
	}

	///
	/// \brief Propagate messages towards the root of the join tree
	///
	/// Full message propagation from leaves to the root of the join tree. Each
	/// message is a potential having a probability (phi) and a utility (psi)
	/// component. The join tree edges are assumed to be directed towards the
	/// root, and only those messages are updated. Once the root receives all
	/// messages, we also compute the expected utility of the current strategy,
	/// which is then returned as a result.
	///
	double collect() {

		if (m_verbosity >= 1) {
			std::cout << "  [COLLECT] Begin message passing towards the root ..." << std::endl;
		}

		for (schedule::const_iterator i = m_schedule.begin();
				i != m_schedule.end(); ++i ) {

			// Select the message m(a->b)
			findex a = i->first;
			findex b = i->second;
			size_t ei = m_edge_indeces[a][b]; // edge index
			m_jointree[ei].message = message(a, b);

			if (m_verbosity >= 2) {
				std::cout << "  - Sending message from " << a << " to " << b
					<< " m(" << a << "," << b << "): elim = " << (m_scopes[a] - m_jointree[ei].separator) << std::endl;
				std::cout << "      P: " << m_jointree[ei].message.first << std::endl;
				std::cout << "      U: " << m_jointree[ei].message.second << std::endl;
			}
		}

		// Collect the original probability, policy and utility factors
		potential pot = incoming(m_root);
		double eu = (pot.first * pot.second).sum(m_scopes[m_root]).sum();

		if (m_verbosity >= 1) {
			std::cout << "  Finished collect to the root in "
					<< (timeSystem() - m_start_time) << " seconds" << std::endl;
		}

		if (m_verbosity >= 2) {
			std::cout << "  Current strategy (expected utility): " << eu << std::endl;
			std::cout << "  Current strategy (decision policies):" << std::endl;
			for (size_t i = 0; i < m_decisions.size(); ++i) {
				vindex d = m_decisions[i];
				findex p = m_policy[d].second;
				std::cout << "    " << d << " (d) : " << m_gmo.get_factor(p) << std::endl;
			}
		}

		return eu; // expected utility of the current policy
	}

	///
	/// \brief Retract the policy from a cluster
	/// \param a 	The index of the cluster
	/// \param d	The index of the decision variable
	/// \return The decision policy factor associated with the decision variable
	///
	factor retract(findex a, vindex d) {
		assert(a == m_policy[d].first); // make sure we have the same cluster
		findex i = m_policy[d].second;
		return m_gmo.get_factor(i); // get the decision policy factor

//		for (flist::const_iterator j = m_originals[a].begin();
//				j != m_originals[a].end(); ++j) {
//			const factor& f = m_gmo.get_factor(*j);
//			if (f.get_type() == factor::FactorType::Decision &&
//					f.vars().contains(var(d))) {
//				return f;
//			}
//		}
//
//		return factor(1.0); // dummy
	}

	///
	/// \brief Update the policy by overwriting the corresponding factor
	///
	/// \param a 	The index of the cluster
	/// \param d 	The index of the decision variable
	/// \param p	The policy factor
	///
	void replace(findex a, vindex d, const factor& p) {
		assert(a == m_policy[d].first); // make sure we have the same cluster
		findex i = m_policy[d].second;
		m_gmo.set_factor(i, p);

//		for (flist::const_iterator j = m_originals[a].begin();
//				j != m_originals[a].end(); ++j) {
//			const factor& f = m_gmo.get_factor(*j);
//			if (f.get_type() == factor::FactorType::Decision &&
//				f.vars().contains(var(d))) {
//				m_gmo.set_factor(*j, p);
//				break;
//			}
//		}
	}

	///
	/// \brief Create a pure (random) decision policy
	///
	///	\param d 		The index of the decision variable
	///	\param val		The value of the decision variable (action)
	///	\param random	Flag indicating a uniform random pure policy
	///	\param delta	The policy factor
	///	\return the policy factor representing either a fixed pure policy or
	/// a uniform random pure policy.
	///
	factor init_policy(vindex d, size_t val, const bool random, const factor& delta) {

		variable VX = var(d); // decision variable
		size_t opt = val;
		if (random) opt = randi(VX.states()); // select randomly a value

		// Compute the pure policy factor, ie. for each configuration of
		// the decision variable's parents, get the corresponding
		// sliced/conditioned factor and set the table entry to 1 if the value
		// corresponding to the decision variable is the same as the one sampled
		// otherwise set the table entry to 0
		factor ndelta = delta; // to have the same scope
		size_t n = delta.vars().size();
		variable_set scope = delta.vars();
		size_t pd; // index of the decision variable in policy scope
		for (size_t e = 0; e < delta.numel(); ++e) {

			// Get the scope configuration correp. to the current index (big_endian)
			int i = e;
			size_t p = 0;
			std::vector<size_t> I(n); // configuration of source variables
			for (variable_set::const_iterator v = delta.vars().begin();
					v != delta.vars().end(); ++v) {
				if (*v == VX) pd = p;
				I[p] = i % v->states();
				i -= I[p];
				i /= v->states();
				p++;
			}

			// Slice the factor on current parent configuration
			p = 0; // reset index
			factor tmp = delta;
			for (variable_set::const_iterator v = delta.vars().begin();
					v != delta.vars().end(); ++v) {
				if (*v != VX) {
					tmp = tmp.condition(*v, I[p]);
				}

				p++;
			}

			// Find the optimal decision
			size_t old = I[pd];
			if (opt == old) {
				ndelta[e] = 1;
			} else {
				ndelta[e] = 0;
			}
		}

		return ndelta;
	}

	///
	/// \brief Update the policy of a decision variable.
	///
	/// \param d 	The decision variable index
	///
	void update_policy(vindex d) {

		// Get the cluster containing the variable's decision policy
		assert(m_vtypes[d] == 'd'); // sanity checks
		variable VX = var(d); // decision variable
		findex a = m_policy[d].first; // cluster where decision resides

		if (m_verbosity >= 1) {
			if (m_pure_strategy) {
				std::cout << "   [UPDATE] pure policy for decision " << d << std::endl;
			} else {
				std::cout << "   [UPDATE] stochastic policy for decision " << d << std::endl;
			}
		}

		// Retract the current policy, collect messages, compute contraction
		factor delta = retract(a, d);
		assert(delta.vars().contains(VX)); // sanity checks
		factor dummy = delta;
		dummy.fill(1.0);
		replace(a, d, dummy);
		potential pot = incoming(a);
		factor cont = (pot.first * pot.second).marginal(delta.vars());

		if (m_verbosity >= 1) {
			std::cout << "    Target cluster: " << a << std::endl;
			std::cout << "    Current policy: " << delta << std::endl;
			std::cout << "    Contraction: " << cont << std::endl;
		}

		factor ndelta = delta; // to have the same scope
		if (m_pure_strategy) {
			// Compute the skewed distribution, ie. for each configuration of
			// the decision variable's parents, get the argmax in the corresponding
			// sliced/conditioned factor and set it as the degenerate distribution
			// in the target policy (for that particular configuration)
			size_t n = delta.vars().size();
			variable_set scope = delta.vars();
			size_t pd; // index of the decision variable in policy scope
			for (size_t e = 0; e < delta.numel(); ++e) {

				// Get the scope configuration correp. to the current index (big_endian)
				int i = e;
				size_t p = 0;
				std::vector<size_t> I(n); // configuration of source variables
				for (variable_set::const_iterator v = delta.vars().begin();
						v != delta.vars().end(); ++v) {
					if (*v == VX) pd = p;
					I[p] = i % v->states();
					i -= I[p];
					i /= v->states();
					p++;
				}

				// Slice the factor on current parent configuration
				p = 0; // reset index
				factor tmp = cont;
				for (variable_set::const_iterator v = delta.vars().begin();
						v != delta.vars().end(); ++v) {
					if (*v != VX) {
						tmp = tmp.condition(*v, I[p]);
					}

					p++;
				}

				// Find the optimal decision
				size_t old = I[pd];
				size_t opt = tmp.argmax();
				if (opt == old) {
					ndelta[e] = 1;
				} else {
					ndelta[e] = 0;
				}
			}
		} else {
			// Compute the skewed distribution, ie. for each configuration of
			// the decision variable's parents, get the argmax in the corresponding
			// sliced/conditioned factor and set it as the degenerate distribution
			// in the target policy (for that particular configuration)
			factor ndelta = delta; // to have the same scope
			size_t n = delta.vars().size();
			variable_set scope = delta.vars();
			size_t pd; // index of the decision variable in policy scope
			for (size_t e = 0; e < delta.numel(); ++e) {

				// Get the scope configuration correp. to the current index (big_endian)
				int i = e;
				size_t p = 0;
				std::vector<size_t> I(n); // configuration of source variables
				for (variable_set::const_iterator v = delta.vars().begin();
						v != delta.vars().end(); ++v) {
					if (*v == VX) pd = p;
					I[p] = i % v->states();
					i -= I[p];
					i /= v->states();
					p++;
				}

				// Slice the factor on current parent configuration
				p = 0; // reset index
				factor tmp = cont;
				for (variable_set::const_iterator v = delta.vars().begin();
						v != delta.vars().end(); ++v) {
					if (*v != VX) {
						tmp = tmp.condition(*v, I[p]);
					}

					p++;
				}

				// Now, we have the factor defined over the decision variable
				tmp = tmp.normalize(); // normalize it

				// Find the optimal decision
				size_t old = I[pd];
				double dd = tmp[old];
				ndelta[e] = tmp[old];
			}
		}

		if (m_verbosity >= 1) {
			std::cout << "    Updated policy: " << ndelta << std::endl;
			std::cout << "    Finished update step." << std::endl;
		}

		// Update the policy
		replace(a, d, ndelta);

		// Compute the expected utility of the current strategy
		potential p = incoming(a);
		double eu = (p.first*p.second).sum(m_scopes[a]).max();
		m_meu = std::max(eu, m_meu);

		if (m_verbosity >= 1) {
			std::cout << "    Current strategy (expected utility): " << eu << std::endl;
			std::cout << "    Current strategy (decision policies):" << std::endl;
			for (size_t i = 0; i < m_decisions.size(); ++i) {
				vindex d = m_decisions[i];
				findex p = m_policy[d].second;
				std::cout << "     " << d << " (d) : " << m_gmo.get_factor(p) << std::endl;
			}
		}
	}

	///
	/// \bried Brute-force enumeration of all possible strategies.
	///
	/// For simplicity, assume for now that the decision variables are parentless.
	///
	void brute_force() {

		// Initialize the algorithm
		init();

		m_meu = -infty();
		std::cout << "Finding optimal policy by brute-force enumeration ..." << std::endl;

		// Safety checks.
		bool ok = true;
		for (policy::const_iterator pi = m_policy.begin();
				pi != m_policy.end(); ++pi) {
			size_t fid = pi->second.second;
			const factor& f = m_gmo.get_factor(fid);
			if (f.vars().size() > 1) ok = false;
			if (m_verbosity > 0) {
				std::cout << " d " << pi->first << " : " << f << std::endl;
			}
		}

		if (!ok) {
			throw std::runtime_error("Brute force supported only for parentless decisions.");
		}

		// Initialize the decision variables and their values
		vector<int> opt;
		vector<vindex> scope = m_decisions;
		vector<int> values(scope.size(), 0);
		values[values.size()-1] = -1;

		// Enumerate all possible assignments
		int i, num_policies = 0;
		while (true) {

			// Enumerate the scope variables
			for (i = scope.size() - 1; i >= 0; --i) {
				vindex v = scope[i];
				int last = var(v).states() - 1;
				if (values[i] < last) break;
				values[i] = 0;
			}

			if (i < 0) break; // done;
			++values[i];

			// Now all variables in scope have a specific value combination
			for (size_t j = 0; j < values.size(); ++j) {
				variable_set vs(var(scope[j]));
				factor f(vs, 0.0);
				f[values[j]] = 1.0;

				size_t a = m_policy[scope[j]].first; // get the cluster
				replace(a, scope[j], f); // replace the policy
			}

			// Evaluate the current policy
			double eu = collect();
			std::cout << "  [ ";
			std::copy(values.begin(), values.end(),
					std::ostream_iterator<int>(std::cout, " "));
			std::cout << "] : " << eu << std::endl;

			if (eu > m_meu) {
				m_meu = eu;
				opt = values;
			}

			++num_policies;
		}

		std::cout << "Finished evaluating " << num_policies << " policies in "
			<< (timeSystem() - m_start_time) << " seconds" << std::endl;
		std::cout << "MEU is " << m_meu << std::endl;
		std::cout << "OPT policy is [ ";
		std::copy(opt.begin(), opt.end(), std::ostream_iterator<int>(std::cout, " "));
		std::cout << "] : " << m_meu << std::endl;
	}

protected:
	// Members:

	limid m_gmo; 						///< Original influence diagram
	double m_meu;						///< Log maximum expected utility
	OrderMethod m_order_method;			///< Variable ordering method
	variable_order_t m_order;			///< Variable order
	size_t m_verbosity;					///< Verbosity level (0=default, 1,2, 3)
	size_t m_iterations;				///< Number of SPU iterations (default 1)
	bool m_pure_strategy;				///< Compute pure or stochastic strategies (default 1=pure)
	bool m_init_pure_strategy;			///< Initialize with pure, or stochastic policies

private:
	// Bucket/Join Tree local structures:

	size_t m_num_clusters;				///< Number of clusters

	vector<bool> m_ctypes;				///< The type of each cluster (true if contains policy factor, false otherwise)
	vector<flist> m_clusters;			///< Clusters corresponding to each of the variables
	vector<flist> m_originals;			///< Original factors (index) for each cluster
	vector<variable_set> m_scopes;		///< The scope (vars) for each cluster
	schedule m_schedule;				///< Initial message propagation schedule towards the root
	policy m_policy;					///< Map decision variable to pair <cluster, policy factor>

	matrix m_edge_indeces;				///< Stores the index of the edge between adjacent clusters

	vector<edge> m_jointree;			///< All edges in the join(bucket) tree (i.e., a->b and b->a)
	vector<vindex> m_decisions;			///< Ordered list of decisions
	vector<int> m_parents;				///< Parent for each cluster in the bucket tree
	vector<flist> m_children;			///< Children of each cluster in the bucket tree
	size_t m_root;						///< Root of the bucket tree (ie, the first decision)
	vector<flist> m_neighbors;			///< For each cluster, list its neighbors in the bucket tree
};


} // end namespace


#endif /* IBM_MERLIN_SPU_H_ */
