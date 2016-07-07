/*
 * limid.h
 *
 *  Created on: 20 Jul 2015
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

/// \file limid.h
/// \brief A limited memory influence diagram
/// \author Radu Marinescu


#ifndef IBM_MERLIN_LIMID_H_
#define IBM_MERLIN_LIMID_H_

#include "graphical_model.h"

namespace merlin {

///
/// \brief Limited memory influence diagram base model.
///
/// Simplest form of an *influence diagram* model, namely just a collection of
/// chance and decision variables, and probability, utility and decision factors
//  (defined on subsets of variables). For standard influence diagrams, the
//  scope of the decision factors is not known a priori (it is actually
//  computed automatically during the maximum expected utility computation). For
//  limited memory influence diagrams (LIMIDs) the scope of the decision factors
//  is known a priori.
/// Internally, factors and variables are mainly referenced by integer indices, with
///  0 <= f < nFactors() and 0 <= v < nvar().  However, many interface functions are called
///  using a variable object, Var(label,dim).  To convert, var(v) gives the vth Var object
///  and the internal function _vindex(V) gives the index corresponding to variable object V.
///
class limid: public graphical_model {
public:

	///
	/// \brief Constructs an empty LIMID.
	///
	limid() : graphical_model(), m_forgetful(true) {};

	///
	/// \brief Copy constructor.
	///
	limid(const limid& lm) :
		graphical_model((graphical_model&) lm), m_vtypes(lm.m_vtypes),
		m_ftypes(lm.m_ftypes), m_forgetful(lm.m_forgetful), m_porder(lm.m_porder) {};

	///
	/// \brief Assignment operator (deep copy).
	/// \param lm 	The LIMID model to be copied from
	/// \return a reference to the modified object containing the copied content.
	///
	limid& operator=(const limid& lm) {
		graphical_model::operator=((graphical_model&) lm); // copy graphical model elements over
		m_vtypes = lm.m_vtypes;				// copy the variable types
		m_ftypes = lm.m_ftypes;				// copy the factor types
		m_forgetful = lm.m_forgetful;			// copy the LIMID flag
		m_porder = lm.m_porder;				// copy the partial order (ID only)
		return *this;
	};

	///
	/// \brief Clone the LIMID model.
	/// \return the pointer to the new object representing the copy of
	/// the current model.
	///
	virtual limid* clone() {
		limid* lm = new limid(*this);
		return lm;
	};

	///
	/// \brief Destructor.
	///
	virtual ~limid() {};

public:

	///
	/// \brief Read the LIMID model from a file in the UAI format.
	/// \param file_name 	The full path to the file
	///
	void read(const char* file_name) {
		std::ifstream is(file_name);
		if (is.fail()) {
			std::cout << "Error while reading input file: " << file_name << std::endl;
			throw std::runtime_error("Input file error");
		}

		size_t nvar, ncliques, csize, v, nval, ndec = 0, max_csize = 0;
		char st[20];
		is >> st;
		if ( strcasecmp(st, "LIMID") == 0 ) {
			m_forgetful = true;
		} else if (strcasecmp(st, "ID") == 0) {
			m_forgetful = false;
		} else {
			throw std::runtime_error("Only UAI Limid-format files are supported currently");
		}

		// Prologue
		std::cout << "Reading " << file_name << std::endl;

		// Read number of variables
		is >> nvar;
		std::vector<size_t> dims(nvar);
		for (size_t i = 0; i < nvar; i++)
			is >> dims[i];	// read domain size
		m_vtypes.resize(nvar);
		for (size_t i = 0; i < nvar; i++) {
			is >> m_vtypes[i];
			if (m_vtypes[i] == 'd') ++ndec;
		}

		// Read the temporal order of the decisions (ID only)
		if (m_forgetful == false) {
			m_porder.resize(nvar);
			for (size_t i = 0; i < nvar; i++)
				is >> m_porder[i];
		}

		// Read the factor scopes (probability, utility, decision)
		is >> ncliques;
		std::vector<std::vector<variable> > cliques(ncliques);
		std::vector<variable_set> sets(ncliques);
		vector<char> ftypes(ncliques);
		for (size_t i = 0; i < ncliques; i++) {
			is >> ftypes[i];
			is >> csize;
			cliques[i].reserve(csize);
			for (size_t j = 0; j < csize; j++) {
				is >> v;
				variable V(v, dims[v]);
				cliques[i].push_back(V);
				sets[i] |= V;
			}
			max_csize = std::max(csize, max_csize);
		}

		// Read the factor tables (ignore decisions as they have size 0)
		vector<char> temp;
		std::vector<factor> tables;
		for (size_t i = 0; i < ncliques; i++) {
			is >> nval;
			if (ftypes[i] == 'd' && m_forgetful == false) continue;
			if (ftypes[i] != 'd') assert(nval == sets[i].num_states());

			factor f(sets[i], 0.0); // preallocate memory and convert from given order, bigEndian
			if (ftypes[i] != 'd') {
				convert_index ci(cliques[i], false, true);
				for (size_t j = 0; j < nval; j++) {
					size_t k = ci.convert(j);
					is >> f[k];
				}

				// Re-code 0 probability entries to 1e-06
//				if (m_ftypes[i] == 'p') {
//					for (size_t j = 0; j < nval; ++j)
//						if (tables[i][j] == 0.0)
//							tables[i][j] = 1e-06;
//				}

				if (ftypes[i] == 'p') {
					f.set_type(factor::FactorType::Probability);
				} else if (ftypes[i] == 'u') {
					f.set_type(factor::FactorType::Utility);
				}

				tables.push_back(f);
				temp.push_back(ftypes[i]);
			} else {
				// initialize with a uniform random policy
				vindex d = cliques[i].back();
				double v = 1.0/(double)dims[d];
				f.fill(v);
				f.set_type(factor::FactorType::Decision);
				tables.push_back(f);
				temp.push_back(ftypes[i]);
			}
		}

		m_factors = tables;
		m_ftypes = temp;
		fixup();

		// Log statistics
		size_t num_prob = 0, num_util = 0, num_dec = 0;
		for (size_t i = 0; i < m_ftypes.size(); ++i) {
			if (m_ftypes[i] == 'p') num_prob++;
			else if (m_ftypes[i] == 'u') num_util++;
			else if (m_ftypes[i] == 'd') num_dec++;
		}

		std::cout << " + model type          : " << (m_forgetful ? "LIMID" : "ID") << std::endl;
		std::cout << " + number of variables : " << nvar << std::endl;
		std::cout << " + decision variables  : " << ndec << std::endl;
		std::cout << " + chance variables    : " << (nvar - ndec) << std::endl;
		std::cout << " + max domain size     : " << *(std::max_element(dims.begin(), dims.end())) << std::endl;
		std::cout << " + max scope size      : " << max_csize << std::endl;
		std::cout << " + number of factors   : " << m_factors.size() << std::endl;
		std::cout << " + probability factors : " << num_prob << std::endl;
		std::cout << " + utility factors     : " << num_util << std::endl;
		std::cout << " + policy factors (dec): " << num_dec << std::endl;
	}

	///
	/// \brief Output operator.
	/// \param out 	The reference of an output stream
	/// \param gm 	The reference of a graphical model
	/// \return the reference of the modified output stream containing the
	/// 	formatted content of the graphical model.
	///
	friend std::ostream& operator<<(std::ostream& out, const limid& lm) {
		out << "Dumping decision network content with " << lm.nvar() << " variables and "
				<< lm.num_factors() << " factors: " << std::endl;
		out << (lm.m_forgetful ? "LIMID" : "ID") << std::endl;
		for (size_t j = 0; j < lm.get_factors().size(); j++)
			out << " " << j << " " << lm.m_ftypes[j]
				<< " " << lm.get_factors()[j] << std::endl;
		return out;
	}

	///
	/// \brief Check if is ID or LIMID.
	///
	bool islimid() {
		return m_forgetful;
	}

    ///
    /// \brief Find a variable elimination order.
	///
	/// For standard IDs, the ordering must respect the partial order induced
	/// by the temporal order of the decisions. For LIMIDs there is no such
	/// restriction of the temporal order of the decisions.
    /// \param ord_type 	The ordering method
    /// \return the variable ordering corresponding to the method, such that
    ///		the first variable in the ordering is eliminated first.
    ///
	virtual variable_order_t order(OrderMethod ord_type) const {
		variable_order_t order;
		order.resize(nvar());

		if (m_forgetful) { // Limited memory influence diagrams (LIMID)
			if (ord_type == OrderMethod::Random) {	// random orders are treated here
				for (size_t i = 0; i < nvar(); i++)
					order[i] = var(i).label();	//   build a list of all the variables
				std::random_shuffle(order.begin(), order.end());//   and randomly permute them
				return order;											//   then return
			}

			std::vector<variable_set> adj = mrf();
			//vector<set<Var> > adj(adj1.size());
			//for (size_t i=0;i<adj1.size();++i) adj[i] = set<Var>(adj1[i].begin(),adj1[i].end());

			typedef std::pair<double, size_t> NN;
			typedef std::multimap<double, size_t> sMap;
			sMap scores;
			std::vector<sMap::iterator> reverse(nvar());

			for (size_t v = 0; v < nvar(); v++) 		// get initial scores
				reverse[v] = scores.insert(NN(order_score(adj, v, ord_type), v));

			for (size_t ii = 0; ii < nvar(); ++ii) {// Iterate through, selecting variables
				sMap::iterator first = scores.begin();// Choose a random entry from among the smallest
				sMap::iterator last = scores.upper_bound(first->first);
				std::advance(first, randi(std::distance(first, last)));
				size_t i = first->second;

				order[ii] = var(i).label();  	       // save its label in the ordering
				scores.erase(reverse[i]);					// remove it from our list
				variable_set vi = adj[i]; // go through adjacent variables (copy: adj may change)
				variable_set fix;					//   and keep track of which need updating
				for (variable_set::const_iterator j = vi.begin(); j != vi.end(); ++j) {
					size_t v = _vindex(*j);
					adj[v] |= vi;             // and update their adjacency structures
					adj[v] /= var(i);
					if (fix.size() < scores.size()) {
						if (ord_type == OrderMethod::MinWidth
								|| ord_type == OrderMethod::WtMinWidth)
							fix |= adj[v]; //var(v);	// (width methods only need v, not nbrs)
						else
							fix |= adj[v];	// come back and recalculate their scores
					}
				}
				for (variable_set::const_iterator j = fix.begin(); j != fix.end(); ++j) {
					size_t jj = j->label();
					scores.erase(reverse[jj]);	// remove and update (score,index) pairs
					reverse[jj] = scores.insert(NN(order_score(adj, jj, ord_type), jj));
				}
			}
			return order;
		} else { // Standard influence diagrams (ID)

//			std::cout << "Building constrained elimination order:" << std::endl;

			// Safety checks
			assert(m_porder.size() == order.size());

			// Select the bundles of chance variables preceeding each decision
			vector<vindex> decisions; // ordered decision variables
			vector<vector<vindex> > bundles;
			vector<vindex> temp;
			for (vector<vindex>::const_reverse_iterator ri = m_porder.rbegin();
					ri != m_porder.rend(); ++ri) {
				vindex v = *ri;
				if (m_vtypes[v] == 'd') { // decision variable
					decisions.push_back(v);
					bundles.push_back(temp); // some bundles can be empty
					temp.clear();
				} else { // chance variable
					temp.push_back(v);
				}
			}
			bundles.push_back(temp); // add the last bundle of chance variables
			std::reverse(decisions.begin(), decisions.end());

			// Dump bundles
//			for (size_t b = 0; b < bundles.size(); ++b) {
//				std::cout << "  Bundle " << b << ": ";
//				std::copy(bundles[b].begin(), bundles[b].end(),
//						std::ostream_iterator<int>(std::cout, " "));
//				std::cout << std::endl;
//			}
//			// And decisions
//			std::cout << "  Decisions: ";
//			std::copy(decisions.begin(), decisions.end(),
//					std::ostream_iterator<int>(std::cout, " "));
//			std::cout << std::endl;

			size_t m = decisions.size(); // number of decisions
			assert( bundles.size() == m+1);
			order.clear();
			if (ord_type == OrderMethod::Random) {	// random orders are treated here
				for (size_t i = 0; i < bundles.size(); i++) {
					std::random_shuffle(bundles[i].begin(), bundles[i].end());
					std::copy(bundles[i].begin(), bundles[i].end(),
							std::back_inserter(order));
					if (decisions.empty() == false) {
						vindex d = decisions.back();
						order.push_back(d);
						decisions.pop_back();
					}
				}
				return order;		//   then return
			}

			std::vector<variable_set> adj = mrf();
			//vector<set<Var> > adj(adj1.size());
			//for (size_t i=0;i<adj1.size();++i) adj[i] = set<Var>(adj1[i].begin(),adj1[i].end());

			typedef std::pair<double, size_t> NN;
			typedef std::multimap<double, size_t> sMap;

			// Iterate through bundles of chance variables
			for (size_t b = 0; b < bundles.size(); ++b) {
//				std::cout << "  Ordering bundle " << b << std::endl;

				vector<vindex> chance_vars = bundles[b];
				sMap scores;
				std::vector<sMap::iterator> reverse(nvar(), scores.end());

				for (size_t j = 0; j < chance_vars.size(); j++) { 		// get initial scores
					vindex v = chance_vars[j];
					reverse[v] = scores.insert(NN(order_score(adj, v, ord_type), v));
				}

				for (size_t ii = 0; ii < chance_vars.size(); ++ii) {// Iterate through, selecting variables
					sMap::iterator first = scores.begin();// Choose a random entry from among the smallest
					sMap::iterator last = scores.upper_bound(first->first);
					std::advance(first, randi(std::distance(first, last)));
					size_t i = first->second;

					order.push_back(i);       // save its label in the ordering
					scores.erase(reverse[i]); // remove it from our list
					variable_set vi = adj[i]; // go through adjacent variables (copy: adj may change)
					variable_set fix;         //   and keep track of which need updating
					for (variable_set::const_iterator j = vi.begin(); j != vi.end(); ++j) {
						size_t v = _vindex(*j);
						adj[v] |= vi;             // and update their adjacency structures
						adj[v] /= var(i);
						if (fix.size() < scores.size()) {
							if (ord_type == OrderMethod::MinWidth
									|| ord_type == OrderMethod::WtMinWidth)
								fix |= adj[v]; //var(v);	// (width methods only need v, not nbrs)
							else
								fix |= adj[v];	// come back and recalculate their scores
						}
					}
					for (variable_set::const_iterator j = fix.begin(); j != fix.end(); ++j) {
						size_t jj = j->label();
						if (reverse[jj] != scores.end()) {
							scores.erase(reverse[jj]);	// remove and update (score,index) pairs
							reverse[jj] = scores.insert(NN(order_score(adj, jj, ord_type), jj));
						}
					}
				}

				// add the current decision variable
				if (decisions.empty() == false) {
					vindex d = decisions.back();
					decisions.pop_back();
					order.push_back(d);
				}
			}

			return order;
		}
	}


	void test() {
		variable_order_t temp;
		temp = order(OrderMethod::MinFill);
		std::cout << "Elimination order is: ";
		for (size_t i = 0; i < temp.size(); ++i) {
			std::cout << " " << temp[i];
		}
		std::cout << std::endl;
	}

	///
	/// \brief Return the variable types.
	///
	const vector<char>& get_vtypes() const {
		return m_vtypes;
	}

	///
	/// \brief Return the factor types.
	///
	const vector<char>& get_ftypes() const {
		return m_ftypes;
	}

protected:

	// Members:

	vector<char> m_vtypes;				///< Variable types ('c' = chance, 'd' = decision)
	vector<char> m_ftypes;				///< Factor types ('p' = probability, 'u' = utility, 'd' = decision)
	bool m_forgetful;					///< LIMID flag (ie, forgetful)
	vector<vindex> m_porder;			///< Partial order induced by the temporal order of decisions (ID only)

};

}// end namespace

#endif /* IBM_MERLIN_LIMID_H_ */
