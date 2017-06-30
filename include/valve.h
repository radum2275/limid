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
/// \brief Valuation algebra based variable elimination for IDs
/// \author Radu Marinescu

#ifndef IBM_MERLIN_VALVE_H_
#define IBM_MERLIN_VALVE_H_

#include "limid.h"
#include "algorithm.h"
#include "valuation.h"

namespace merlin {

/**
 * Valuation Algebra based Variable Elimination (VE)
 *
 * Models supported: ID, LIMID
 *
 * For standard IDs, variable elimination (VE) assumes a constrained elimination
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
 * For LIMIDs, the variable elimination works on an unconstrained elimination
 * order. Moreover, the parents sets of each decision variable are fixed and
 * given as input. The valuation algebra doesn't assume division.
 *
 */
class valve : public limid, public algorithm {
public:
	typedef limid::findex findex;        ///< Factor index
	typedef limid::vindex vindex;        ///< Variable index
	typedef limid::flist flist;          ///< Collection of factor indices

	///
	/// \brief Properties of the algorithm
	///
	MER_ENUM( Property , Order,Debug );
	MER_ENUM( Operator , Sum,Max,Min );

public:

	///
	/// \brief Default constructor.
	///
	valve() : limid() {
		set_properties();
	}

	///
	/// \brief Constructor with a graphical model.
	///
	valve(const limid& lm) : limid(lm), m_gmo(lm) {
		clear_factors();
		set_properties();
	}

	///
	/// \brief Clone the algorithm.
	/// \return the pointer to the new object containing the cloned algorithm.
	///
	virtual valve* clone() const {
		valve* lm = new valve(*this);
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
			set_properties("Order=MinFill,Debug=1");
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
	/// \param v
	///			the input valuation (const)
	/// \param vs
	///			the set of variables to be eliminated (eliminator)
	/// \param op
	///			the elimination operator (sum or max)
	///
	/// \return a new valuation that does not contain the eliminator.
	///
	valuation elim(const valuation& v, const variable_set& vs, const Operator op) {
		if (op == Operator::Sum) {
			return v.sum(vs);
		} else if (op == Operator::Max) {
			return v.max(vs);
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

		// Check if LIMID
		bool is_limid = m_gmo.islimid();
		if (is_limid) {
			throw std::runtime_error("BE is only supported for standard IDs.");
		}

		// Prologue
		std::cout << "Initialize solver ..." << std::endl;
		std::cout << " + models supported : ID" << std::endl;
		std::cout << " + algorithm        : VALVE (with valuations)" << std::endl;

		// Construct the elimination ordering
		m_order = m_gmo.order(m_order_method);

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
	void insert(std::vector<flist>& adj, findex fid, const valuation& f,
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
	/// \brief Run variable elimination for IDs.
	///
	virtual void run() {

		// Initialize the algorithm
		init();

		// Get the input factors and convert them into valuations
		std::vector<valuation> fin;
		const std::vector<factor>& temp = m_gmo.get_factors();
		for (std::vector<factor>::const_iterator it = temp.begin();
				it != temp.end(); ++it) {
			const factor& f = (*it);
			if (f.get_type() == factor::FactorType::Probability) {
				fin.push_back(valuation(f, factor(0.0)));
			} else if (f.get_type() == factor::FactorType::Utility) {
				fin.push_back(valuation(factor(1.0), f));
			}
		}

		findex fid = fin.size();
		flist roots; // constant valuations
		std::map<findex, int> ftypes; // factor types

		// Brute force to ensure correctness of results
		if (m_debug) {
			std::cout << "*** Start brute force computation ***" << std::endl;
			valuation Comb;
			for (size_t i = 0; i < fin.size(); ++i) {
				Comb *= fin[i];
			}
			std::cout << "Combined valuation: " << Comb << std::endl;
			for (vector<vindex>::const_iterator x = m_order.begin();
					x != m_order.end(); ++x) {
				if (m_vtypes[*x] == 'c') {
					Comb = elim(Comb, var(*x), Operator::Sum);
				} else if (m_vtypes[*x] == 'd') {
					Comb = elim(Comb, var(*x), Operator::Max);
				}
			}
			std::cout << "Brute force MEU is " << Comb.util().max() << std::endl;
			std::cout << "Final valuation: " << Comb << std::endl;
			std::cout << "*** End brute force computation ***" << std::endl;
		}

		if (m_debug) {
			std::cout << "Partition valuations into buckets ..." << std::endl;
		}

		// Partition into buckets: mark factors depending on variable i
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
				continue;  // check that we have some valuations over this variable

			// Get the list of all valuation IDs contained in this bucket
			flist ids = vin[*x];

			// Process the bucket of the current variable
			Operator op = (m_vtypes[*x] == 'c') ? Operator::Sum : Operator::Max;
			if (m_vtypes[*x] == 'c') { // chance variable
				std::cout << "  Eliminating (C) variable " << *x << std::endl;
			} else if (m_vtypes[*x] == 'd') {
				std::cout << "  Eliminating (D) variable " << *x << std::endl;
			}


			// Combine all valuations
			valuation comb;
			for (flist::const_iterator i = ids.begin();
					i != ids.end(); ++i) {
				comb *= fin[*i];
			}

			// Eliminate chance variables by summation, decisions by maximization
			valuation f = elim(comb, VX, op); // eliminate by summation/maximization
			fin.push_back(f); // store the new factor
			insert(vin, fid, f, x, m_order); // recompute and update adjacency
			if (fin[fid].nvar() == 0)
				roots |= fid; // keep track of constants separately
			fid++;

			if (m_debug) {
				std::cout << "    Old valuation: " << comb << std::endl;
				std::cout << "    New valuation: " << f << std::endl;
			}
		} // end for

		// Combine all remaining valuations
		valuation res;
		for (size_t i = 0; i < roots.size(); ++i) {
			findex id = roots[i];
			res *= fin[id];
		}

		// MEU is the util component of the combined valuation
		m_meu = res.util().max();

		std::cout << "End variable elimination." << std::endl;
		std::cout << "MEU value is " << m_meu << "\n";
		std::cout << "CPU time is " << timeSystem() - m_start_time << " seconds" << std::endl;

		// Memory usage
		double mem_usage = 0;
		for (size_t i = 0; i < fin.size(); ++i) {
			mem_usage += ((double)fin[i].numel() * sizeof(double) / (1024 * 1024)); // MB
		}
		std::cout << "Memory usage is " << mem_usage << " MBytes" << std::endl;

		// Assemble the decision policy by going backward.
		std::cout << "Begin building optimal policy ..." << std::endl;

		// Backward pass: create optimal decision policy
		mem_usage = 0;
		for (vector<vindex>::const_reverse_iterator x = m_order.rbegin();
				x != m_order.rend(); ++x) {

			// Get the corresponding variable object
			variable VX = var(*x);
			if (m_vtypes[*x] != 'd')
				continue;  // skip over chance variables

			assert(m_vtypes[*x] == 'd');
			flist ids = vin[*x];  // list of all valuation IDs contained in this bucket
			valuation V;
			for (flist::const_iterator i = ids.begin(); i != ids.end(); ++i) {
				findex id = *i;
				V *= fin[id];
			}

			std::cout << "  Policy for decision " << *x << " is: " << V.util() << std::endl;
			m_policy[*x] = V.util();
			mem_usage += ((double)V.numel() * sizeof(double) / (1024*1024)); // MBytes
		}

		std::cout << "End building optimal policy." << std::endl;
		std::cout << "Estimated memory usage is " << mem_usage << " MBytes" << std::endl;
		std::cout << "Done." << std::endl << std::endl;
	}


protected:
	// Members:

	limid m_gmo; 						///< Original influence diagram
	double m_meu;						///< Log maximum expected utility
	std::map<vindex, factor> m_policy;	///< Optimal decision policy
	OrderMethod m_order_method;			///< Variable ordering method
	variable_order_t m_order;			///< Variable order
	bool m_debug;						///< Internal debugging flag

};


} // end namespace


#endif /* IBM_MERLIN_BE_H_ */
