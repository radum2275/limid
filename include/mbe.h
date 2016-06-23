/*
 * be.h
 *
 *  Created on: 21 Jul 2015
 *      Author: radu
 *
 * Copyright (c) 2015, International Business Machines Corporation
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

/// \file be.h
/// \brief Mini-Bucket Elimination algorithm for IDs
/// \author Radu Marinescu

#ifndef IBM_MERLIN_MBE_H_
#define IBM_MERLIN_MBE_H_

#include "limid.h"
#include "algorithm.h"

namespace merlin {

/**
 * Mini-Bucket Elimination (MBE)
 *
 * Models supported: ID, LIMID
 *
 * For standard IDs, mini-bucket elimination (MBE) assumes a constrained elimination
 * order that respects the partial order induced by the temporal order of the
 * decisions (ie, porder is given as input). The input probability and utility
 * factors are partitioned into buckets, each corresponding to a variable. Then,
 * each bucket is processed by a variable elimination procedure which eliminates
 * the bucket variable from the combination of factors from that bucket. Chance
 * buckets typically generate two messages, a probability one and an expected
 * utility ones (note that the utility messages are divided by the compiled
 * probability of that chance bucket). Decision buckets typically generated
 * maximum expected utility messages (any probability components residing in
 * these buckets are actually constants when viewed as functions of the decision
 * variables). The optimal decision policy is recovered by performing a backward
 * step that argmax'es the decision buckets. Notice that parents set of each
 * decision variables is automatically computed during this backward step.
 *
 * For LIMIDs, the mini-bucket elimination works on an unconstrained elimination
 * order. Moreover, the parents sets of each decision variable are fixed and
 * given as input. The valuation algebra doesn't assume division.
 *
 */
class mbe : public limid, public algorithm {
public:
	typedef limid::findex findex;        ///< Factor index
	typedef limid::vindex vindex;        ///< Variable index
	typedef limid::flist flist;          ///< Collection of factor indices

	///
	/// \brief Properties of the algorithm
	///
	MER_ENUM( Property , Order,iBound,Debug );

	MER_ENUM( Operator , Sum,Max,Min );

public:

	///
	/// \brief Default constructor.
	///
	mbe() : limid() {
		set_properties();
	}

	///
	/// \brief Constructor with a graphical model.
	///
	mbe(const limid& lm) : limid(lm), m_gmo(lm) {
		clear_factors();
		set_properties();
	}

	///
	/// \brief Clone the algorithm.
	/// \return the pointer to the new object containing the cloned algorithm.
	///
	virtual mbe* clone() const {
		mbe* lm = new mbe(*this);
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
	/// \brief Set the mini-bucket i-bound parameter.
	///
	void set_ibound(size_t i) {
		m_ibound = i ? i : std::numeric_limits<size_t>::max();
	}

	///
	/// \brief Get the mini-bucket i-bound parameter.
	///
	size_t get_ibound() const {
		return m_ibound;
	}

	///
	/// \brief Set the properties of the algorithm.
	/// \param opt 	The string containing comma separated property value pairs
	///
	virtual void set_properties(std::string opt = std::string()) {
		if (opt.length() == 0) {
			set_properties("Order=MinFill,iBound=2,Debug=1");
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
			case Property::iBound:
				set_ibound(atol(asgn[1].c_str()));
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
	/// \brief Eliminate a set of variables either by summation or maximization.
	/// \param f 	The input factor (const)
	/// \param vs 	The set of variables to be eliminated (eliminator)
	/// \param op	The elimination operator (sum or max)
	/// \return a new factor that does not contain the eliminator.
	///
	factor elim(const factor& f, const variable_set& vs, const Operator op) {
		if (op == Operator::Sum) {
			return f.sum(vs);
		} else if (op == Operator::Max) {
			return f.max(vs);
		} else if (op == Operator::Min) {
			return f.min(vs);
		} else {
			throw std::runtime_error("Unkown elimination operator.");
		}
	}

	///
	/// \brief Initialize the weighted mini-buckets algorithm.
	///
	void init() {

		// Start the timer and store it
		m_start_time = timeSystem();

		bool islimid = m_gmo.islimid();
		if (islimid) {
			throw std::runtime_error("MBE is only supported for standard IDs.");
		}

		// Prologue
		std::cout << "Initialize solver ..." << std::endl;
		std::cout << " + models supported : ID" << std::endl;
		std::cout << " + algorithm        : MBE" << std::endl;
		std::cout << " + i-bound          : " << m_ibound << std::endl;

		if (m_order.size() == 0) { // if we need to construct an elimination ordering
			m_order = m_gmo.order(m_order_method);
		}

		// Get the induced width of the order
		size_t wstar = m_gmo.induced_width(m_order);
		std::cout << " + elimination      : ";
		std::copy(m_order.begin(), m_order.end(),
				std::ostream_iterator<size_t>(std::cout, " "));
		std::cout << std::endl;
		std::cout << " + induced width    : " << wstar << std::endl;
		if (m_porder.empty() == false) {
			std::cout << " + partial order    : ";
			std::copy(m_porder.begin(), m_porder.end(),
					std::ostream_iterator<size_t>(std::cout, " "));
			std::cout << std::endl;
			std::cout << " + decisions (ord)  : ";
			for (size_t i = 0; i < m_porder.size(); ++i) {
				if (m_vtypes[m_porder[i]] == 'd')
					std::cout << m_porder[i] << " ";
			}
			std::cout << std::endl;
		}

		std::cout << "Initialization complete in "
			<< (timeSystem() - m_start_time) << " seconds." << std::endl;
	}

	///
	/// \brief Update the adjacencies with a new factor
	///
	void insert(std::vector<flist>& adj, findex fid, const factor& f,
			vector<vindex>::const_iterator x,
			const variable_order_t& ord) {
		for (vector<vindex>::const_iterator y = x+1;
				y != ord.end(); ++y) {
			variable VY = var(*y);
			if (f.vars().contains(VY)) {
				adj[*y] |= fid;
				break;
			}
		}
	}

	///
	/// \brief Create a mini-bucket partitioning of a set of factors
	/// \param ids 	Ordered list of factor indeces
	/// \return the mini-bucket partitioning such that each mini-bucket contains
	/// at most i-bound distinct variables.
	///
	std::vector<flist> partition(const flist& ids, const std::vector<factor>& factors) {

		if (m_debug) {
			std::cout << "    Begin MB partitioning ..." << std::endl;
			std::cout << "     initial factor indeces: ";
			std::copy(ids.begin(), ids.end(), std::ostream_iterator<int>(std::cout, " "));
			std::cout << std::endl;
			std::cout << "     initial factors: " << std::endl;
			for (flist::const_iterator i = ids.begin(); i != ids.end(); ++i) {
				std::cout << "      " << *i << " : " << factors[*i] << std::endl;
			}
		}

		// Mini-bucket partition
		std::vector<flist> par;

		// Greedy partitioning
		std::multimap<size_t, findex> scores;
		for (flist::const_iterator i = ids.begin(); i < ids.end(); ++i) {
			size_t key = factors[*i].nvar(); // scope size
			scores.insert(std::make_pair(key, *i)); // (scope size, findex)
		}

		size_t pos = 0;
		flist mb; // initialize current mini-bucket
		std::multimap<size_t, findex>::iterator top = scores.begin();
		if (top != scores.end()) {
			mb |= top->second;
			scores.erase(top);
			par.push_back(mb);
		}

		while (!scores.empty()) {
			top = scores.begin(); // smallest scope first

			// Get the actual scopes in the current mini-bucket
			variable_set vs = factors[top->second].vars();
			for (flist::const_iterator j = par[pos].begin(); j != par[pos].end(); ++j)
				vs |= factors[*j].vars();

			// Check if new factor fits in the current mini-bucket
			if (vs.size() <= m_ibound) {
				par[pos] |= top->second; // extend current mini-bucket
			} else {
				par.push_back(flist());
				par[++pos] |= top->second;
			}

			scores.erase(top);
		}

		if (m_debug) {
			std::cout << "     number of mini-buckets: " << par.size() << std::endl;
			for (size_t j = 0; j < par.size(); ++j) {
				std::cout << "      mini-bucket " << j << ": ";
				std::copy(par[j].begin(), par[j].end(), std::ostream_iterator<int>(std::cout, " "));
				std::cout << std::endl;
			}
			std::cout << "    End MB partitioning." << std::endl;
		}

		return par;
	}

	///
	/// \brief Run bucket elimination for IDs.
	///
	virtual void run() {

		// Initialize the algorithm
		init();

		// Get the input factors
		std::vector<factor> fin(m_gmo.get_factors());
		findex fid = fin.size();
		flist roots; // constant factors
		std::map<findex, int> ftypes; // factor types

		if (m_debug) {
			std::cout << "Partition factors into buckets ..." << std::endl;
		}

		// Mark factors depending on variable i
		vector<flist> vin(m_gmo.nvar());
		vector<bool> used(fin.size(), false);
		for (vector<vindex>::const_iterator x = m_order.begin();
				x != m_order.end(); ++x) {
			// Mark all remaining factors depending on variable *x
			variable VX = var(*x);
			for (size_t i = 0; i < fin.size(); ++i) {
				if (used[i] == false && fin[i].vars().contains(VX)) {
					vin[*x] |= i;
					used[i] = true;
				}
			}

			if (m_debug) {
				std::cout << " Bucket " << *x << ":   ";
				std::copy(vin[*x].begin(), vin[*x].end(),
						std::ostream_iterator<size_t>(std::cout, " "));
				std::cout << std::endl;
				for (size_t j = 0; j < vin[*x].size(); ++j) {
					std::cout << "   " << vin[*x][j] << " " << fin[vin[*x][j]] << std::endl;
				}
			}
		}

		if (m_debug) {
			std::cout << "Finished initializing the buckets." << std::endl;
		}

		// Forward pass: eliminate variables one at a time
		size_t max_phi_scope = 0, max_psi_scope = 0;
		std::cout << "Begin variable elimination ..." << std::endl;
		for (vector<vindex>::const_iterator x = m_order.begin();
				x != m_order.end(); ++x) {

			// Get the corresponding variable object
			variable VX = var(*x);
			if (vin[*x].size() == 0)
				continue;  // check that we have some factors over this variable

			// Partition the factors into probabilities (phi's) and utilities (psi's)
			flist ids = vin[*x];  // list of all factor IDs contained in this bucket
			flist phi, psi;
			for (flist::const_iterator i = ids.begin(); i != ids.end(); ++i) {
				findex id = *i;
				if (fin[id].get_type() == factor::FactorType::Probability) {
					phi |= id;
				} else if (fin[id].get_type() == factor::FactorType::Utility) {
					psi |= id;
				}
			}

			// Process the bucket of the current variable
			if (m_vtypes[*x] == 'c') { // chance variable

				std::cout << "  Eliminating (C) variable " << *x << std::endl;

				if (m_debug) {
					for (flist::const_iterator i = phi.begin();
							i != phi.end(); ++i) {
						std::cout << "    PHI: " << *i << " " << fin[*i] << std::endl;
					}
					for (flist::const_iterator i = psi.begin();
							i != psi.end(); ++i) {
						std::cout << "    PSI: " << *i << " " << fin[*i] << std::endl;
					}
				}

				// Create mini-bucket partitioning of the probability factors only (phi)
				std::vector<flist> mini_buckets = partition(phi, fin);

				// Combine the factors in each mini-bucket (used latter)
				vector<factor> comb(mini_buckets.size());
				vector<findex> temp(mini_buckets.size());
				for (size_t i = 0; i < mini_buckets.size(); ++i) {
					comb[i] = factor(1.0);
					for (flist::const_iterator j = mini_buckets[i].begin();
							j != mini_buckets[i].end(); ++j) {
						comb[i] *= fin[*j];
					}

					// Eliminate the chance variable by summation
					Operator op = (i == 0 ? Operator::Sum : Operator::Max);
					factor f = elim(comb[i], VX, op);
					temp[i] = fid; // store the id of the new factor
					f.set_type(factor::FactorType::Probability); // probability factor
					fin.push_back(f); // store the new factor
					insert(vin, fid, f, x, m_order); // insert the factor in a lower bucket
					if (fin[fid].nvar() == 0) roots |= fid; // keep track of constants separately
					fid++;

					if (m_debug) {
						std::cout << "    Prob: " << f << std::endl;
					}

					max_phi_scope = std::max(max_phi_scope, f.nvar());
				}

				// Process each utility factor separately
				for (flist::const_iterator j = psi.begin(); j != psi.end(); ++j) {

					size_t i = 0;
					size_t max_sz = 0;

					// Find the mini-bucket with which it shares the most variables
					for (size_t l = 0; l < comb.size(); ++l) {
						variable_set inter = fin[*j].vars() & comb[l].vars();
						if (inter.size() > max_sz) {
							max_sz = inter.size();
							i = l;
						}
					}

					factor g = comb[i] * fin[*j];
					g = elim(g, VX, Operator::Sum);	// eliminate by summation
					g = g / fin[temp[i]];			// divide by
					g.set_type(factor::FactorType::Utility); // utility factor
					fin.push_back(g);				// store the new factor
					insert(vin, fid, g, x, m_order); 	// recompute and update adjacency
					if (fin[fid].nvar() == 0) roots |= fid; // keep track of constants separately
					fid++;

					if (m_debug) {
						std::cout << "    Util: " << g << std::endl;
					}

					max_psi_scope = std::max(max_psi_scope, g.nvar());
				}

			} else if (m_vtypes[*x] == 'd') { // decision variable

				std::cout << "  Eliminating (D) variable " << *x << std::endl;

				// Process each probability factor separately
				for (flist::const_iterator i = phi.begin(); i != phi.end(); ++i) {
					factor f = fin[*i].slice(VX, 0); // condition on any value
					f.set_type(factor::FactorType::Probability); // probability factor
					fin.push_back(f);				 // store the new factor
					insert(vin, fid, f, x, m_order);  // recompute and update adjacency
					if (fin[fid].nvar() == 0) roots |= fid; // keep track of constants separately
					fid++;

					if (m_debug) {
						std::cout << "    Prob: " << f << std::endl;
					}

					max_phi_scope = std::max(max_phi_scope, f.nvar());
				}

				// Create mini-bucket partitioning of the utility factors only (psi)
				std::vector<flist> mini_buckets = partition(psi, fin);

				// Process the utility mini-buckets
				for (size_t i = 0; i < mini_buckets.size(); ++i) {
					flist mb = mini_buckets[i];
					factor comb(0.0);
					for (flist::const_iterator j = mb.begin();
							j != mb.end(); ++j) {
						comb += fin[*j];
					}

					factor g = comb.max(VX); 	// eliminate by maximization
					g.set_type(factor::FactorType::Utility); // utility factor
					fin.push_back(g);			// store the new factor
					insert(vin, fid, g, x, m_order); // recompute and update adjacency
					if (fin[fid].nvar() == 0) roots |= fid; // keep track of constants separately
					fid++;

					if (m_debug) {
						std::cout << "    Util: " << g << std::endl;
					}

					max_psi_scope = std::max(max_psi_scope, g.nvar());
				}
			}
		} // end for

		// Collect all probability and utility factors (constants)
		factor P(1.0), U(0.0);
		for (size_t i = 0; i < roots.size(); ++i) {
			findex id = roots[i];
			if (fin[id].get_type() == factor::FactorType::Probability) {
				P *= fin[id];
			} else if (fin[id].get_type() == factor::FactorType::Utility) {
				U += fin[id];
			}
		}
		factor F = P*U;
		m_meu = F.max();

		std::cout << "End variable elimination." << std::endl;
		std::cout << "Max phi and psi scopes: " << max_phi_scope << " and " << max_psi_scope << std::endl;
		std::cout << "Upper Bound on MEU value is " << m_meu << "\n";
		std::cout << "CPU time is " << (timeSystem() - m_start_time) << " seconds" << std::endl;

		// Assemble the decision policy by going backward.
		std::cout << "Begin building policy ..." << std::endl;

		// Backward pass: create optimal decision policy
		for (vector<vindex>::const_reverse_iterator x = m_order.rbegin();
				x != m_order.rend(); ++x) {

			// Get the corresponding variable object
			variable VX = var(*x);
			if (m_vtypes[*x] != 'd')
				continue;  // skip over chance variables

			assert(m_vtypes[*x] == 'd');
			flist ids = vin[*x];  // list of all factor IDs contained in this bucket
			factor P(1.0), U(0.0);
			for (flist::const_iterator i = ids.begin(); i != ids.end(); ++i) {
				findex id = *i;
				if (fin[id].get_type() == factor::FactorType::Probability) {
					P *= fin[id];
				} else if (fin[id].get_type() == factor::FactorType::Utility) {
					U += fin[id];
				}
			}

			factor F = P*U;
			std::cout << "  Policy for decision " << *x << " is: " << F << std::endl;
			m_policy[*x] = F;
		}

		std::cout << "End building policy." << std::endl;
		std::cout << "Done." << std::endl << std::endl;
	}


	///
	/// \brief Run bucket elimination for IDs.
	///
	virtual void run2() {

		// Initialize the algorithm
		init();

		// Get the input factors
		std::vector<factor> fin(m_gmo.get_factors());
		findex fid = fin.size();
		flist roots; // constant factors
		std::map<findex, int> ftypes; // factor types

		if (m_debug) {
			std::cout << "Partition factors into buckets ..." << std::endl;
		}

		// Mark factors depending on variable i
		vector<flist> vin(m_gmo.nvar());
		vector<bool> used(fin.size(), false);
		for (vector<vindex>::const_iterator x = m_order.begin();
				x != m_order.end(); ++x) {
			// Mark all remaining factors depending on variable *x
			variable VX = var(*x);
			for (size_t i = 0; i < fin.size(); ++i) {
				if (used[i] == false && fin[i].vars().contains(VX)) {
					vin[*x] |= i;
					used[i] = true;
				}
			}

			if (m_debug) {
				std::cout << " Bucket " << *x << ":   ";
				std::copy(vin[*x].begin(), vin[*x].end(),
						std::ostream_iterator<size_t>(std::cout, " "));
				std::cout << std::endl;
				for (size_t j = 0; j < vin[*x].size(); ++j) {
					std::cout << "   " << vin[*x][j] << " " << fin[vin[*x][j]] << std::endl;
				}
			}
		}

		if (m_debug) {
			std::cout << "Finished initializing the buckets." << std::endl;
		}

		// Forward pass: eliminate variables one at a time
		std::cout << "Begin variable elimination ..." << std::endl;
		for (vector<vindex>::const_iterator x = m_order.begin();
				x != m_order.end(); ++x) {

			// Get the corresponding variable object
			variable VX = var(*x);
			if (vin[*x].size() == 0)
				continue;  // check that we have some factors over this variable

			// Partition the factors into probabilities (phi's) and utilities (psi's)
			flist ids = vin[*x];  // list of all factor IDs contained in this bucket

			// Process the bucket of the current variable
			if (m_vtypes[*x] == 'c') { // chance variable


				std::cout << "  Eliminating (C) variable " << *x << std::endl;

				if (m_debug) {
					for (flist::const_iterator i = ids.begin();
							i != ids.end(); ++i) {
						if (fin[*i].get_type() == factor::FactorType::Probability) {
							std::cout << "    PHI: " << *i << " " << fin[*i] << std::endl;
						} else {
							std::cout << "    PSI: " << *i << " " << fin[*i] << std::endl;
						}
					}
				}

				// Create mini-bucket partitioning of the factors in this bucket
				std::vector<flist> mini_buckets = partition(ids, fin);

				// Process each mini-bucket
				for (size_t i = 0; i < mini_buckets.size(); ++i) {

					flist phi, psi;
					for (flist::const_iterator j = mini_buckets[i].begin();
							j != mini_buckets[i].end(); ++j) {
						findex id = *j;
						if (fin[id].get_type() == factor::FactorType::Probability) {
							phi |= id;
						} else if (fin[id].get_type() == factor::FactorType::Utility) {
							psi |= id;
						}
					}

					// Combine all probability factors in the current mini-bucket
					factor comb(1.0);
					for (flist::const_iterator j = phi.begin();
							j != phi.end(); ++j) {
						comb *= fin[*j];
					}

					// Eliminate the chance variable by summation
					factor f = elim(comb, VX, Operator::Sum);
					f.set_type(factor::FactorType::Probability); // probability factor
					fin.push_back(f); // store the new factor
					insert(vin, fid, f, x, m_order); // insert the factor in a lower bucket
					if (fin[fid].nvar() == 0) roots |= fid; // keep track of constants separately
					fid++;

					if (m_debug) {
						std::cout << "    Prob: " << f << std::endl;
					}

					// Process each utility factor separately
					for (flist::const_iterator j = psi.begin(); j != psi.end(); ++j) {
						factor g = comb * fin[*j];
						g = elim(g, VX, Operator::Sum);	// eliminate by summation
						g = g / f;						// divide by
						g.set_type(factor::FactorType::Utility); // utility factor
						fin.push_back(g);				// store the new factor
						insert(vin, fid, g, x, m_order); 	// recompute and update adjacency
						if (fin[fid].nvar() == 0) roots |= fid; // keep track of constants separately
						fid++;

						if (m_debug) {
							std::cout << "    Util: " << g << std::endl;
						}
					}
				}

			} else if (m_vtypes[*x] == 'd') { // decision variable

				std::cout << "  Eliminating D variable " << *x << std::endl;

				flist phi, psi;
				for (flist::const_iterator i = ids.begin(); i != ids.end(); ++i) {
					findex id = *i;
					if (fin[id].get_type() == factor::FactorType::Probability) {
						phi |= id;
					} else if (fin[id].get_type() == factor::FactorType::Utility) {
						psi |= id;
					}
				}

				// Process each probability factor separately
				for (flist::const_iterator i = phi.begin(); i != phi.end(); ++i) {
					factor f = fin[*i].slice(VX, 0); // condition on any value
					f.set_type(factor::FactorType::Probability); // probability factor
					fin.push_back(f);				 // store the new factor
					insert(vin, fid, f, x, m_order);  // recompute and update adjacency
					if (fin[fid].nvar() == 0) roots |= fid; // keep track of constants separately
					fid++;

					if (m_debug) {
						std::cout << "    Prob: " << f << std::endl;
					}
				}

				// Create mini-bucket partitioning of the utility factors only
				std::vector<flist> mini_buckets = partition(psi, fin);

				// Process the utility factors
				for (size_t i = 0; i < mini_buckets.size(); ++i) {
					flist mb = mini_buckets[i];
					factor comb(0.0);
					for (flist::const_iterator j = mb.begin(); j != mb.end(); ++j) {
						comb += fin[*j];
					}
					factor g = comb.max(VX); 	// eliminate by maximization
					g.set_type(factor::FactorType::Utility); // utility factor
					fin.push_back(g);			// store the new factor
					insert(vin, fid, g, x, m_order); // recompute and update adjacency
					if (fin[fid].nvar() == 0) roots |= fid; // keep track of constants separately
					fid++;

					if (m_debug) {
						std::cout << "    Util: " << g << std::endl;
					}
				}
			}
		} // end for

		// Collect all probability and utility factors (constants)
		factor P(1.0), U(0.0);
		for (size_t i = 0; i < roots.size(); ++i) {
			findex id = roots[i];
			if (fin[id].get_type() == factor::FactorType::Probability) {
				P *= fin[id];
			} else if (fin[id].get_type() == factor::FactorType::Utility) {
				U += fin[id];
			}
		}

		factor F = P*U;
		m_meu = F.max();

		std::cout << "End variable elimination." << std::endl;
		std::cout << "Upper Bound on MEU value is " << m_meu << "\n";
		std::cout << "CPU time is " << (timeSystem() - m_start_time) << " seconds" << std::endl;

		// Assemble the decision policy by going backward.
		std::cout << "Begin building policy ..." << std::endl;

		// Backward pass: create optimal decision policy
		for (vector<vindex>::const_reverse_iterator x = m_order.rbegin();
				x != m_order.rend(); ++x) {

			// Get the corresponding variable object
			variable VX = var(*x);
			if (m_vtypes[*x] != 'd')
				continue;  // skip over chance variables

			assert(m_vtypes[*x] == 'd');
			flist ids = vin[*x];  // list of all factor IDs contained in this bucket
			factor P(1.0), U(0.0);
			for (flist::const_iterator i = ids.begin(); i != ids.end(); ++i) {
				findex id = *i;
				if (fin[id].get_type() == factor::FactorType::Probability) {
					P *= fin[id];
				} else if (fin[id].get_type() == factor::FactorType::Utility) {
					U += fin[id];
				}
			}

			factor F = P*U;
			std::cout << "  Policy for decision " << *x << " is: " << F << std::endl;
			m_policy[*x] = F;
		}

		std::cout << "End building policy." << std::endl;
		std::cout << "Done." << std::endl << std::endl;
	}


protected:
	// Members:

	size_t m_ibound;					///< Mini-bucket i-bound
	limid m_gmo; 						///< Original influence diagram
	double m_meu;						///< Maximum expected utility (upper bound)
	std::map<vindex, factor> m_policy;	///< Optimal decision policy
	OrderMethod m_order_method;			///< Variable ordering method
	variable_order_t m_order;			///< Variable order
	bool m_debug;						///< Internal debugging flag

};


} // end namespace


#endif /* IBM_MERLIN_BE_H_ */
