/*
 * valuation.h
 *
 *  Created on: 28 Jun 2017
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

/// \file valuation.h
/// \brief A valuation algebra for influence diagrams
/// \author Radu Marinescu

#ifndef IBM_MERLIN_VALUATION_H_
#define IBM_MERLIN_VALUATION_H_

#include <float.h>

#include "factor.h"

namespace merlin {

///
/// \brief A valuation algebra for influence diagrams.
///
/// Factor based valuation algebra for influence diagrams and limids.
///
///

class valuation {

protected:
	factor m_prob;		// probability component
	factor m_util;		// utility component

public:

	///
	/// \brief Default constructor
	///
	valuation() : m_prob(factor(1)), m_util(factor(0)){

	}

	///
	/// \brief Constructor
	///
	valuation(const factor& p, const factor& u)
		: m_prob(p), m_util(u) {

	}

	///
	/// \brief Copy constructor
	///
	valuation(const valuation& v)
		: m_prob(v.m_prob), m_util(v.m_util) {

	}

	///
	/// \brief Assignment operator
	///
	valuation& operator=(const valuation& v) {
		m_prob = v.m_prob;
		m_util = v.m_util;
		return (*this);
	}

	///
	/// \brief Destructor
	///
	virtual ~valuation() {};

	///
	/// \brief Clone from pointer.
	///
	virtual valuation* clone() {
		valuation *v = new valuation(*this);
		return v;
	};

	///
	/// \brief Get the probability component
	///
	factor prob() {
		return m_prob;
	}

	///
	/// \brief Get the utility component
	///
	factor util() {
		return m_util;
	}

	///
	/// \brief Get the probability component (const reference)
	///
	const factor& prob() const {
		return m_prob;
	}

	///
	/// \brief Get the utility component (const reference)
	///
	const factor& util() const {
		return m_util;
	}

	///
	/// \brief Combination operator
	///
	valuation& operator*(const valuation& a) {
		factor p = m_prob*a.m_prob;
		factor u = m_prob*a.m_util + a.m_prob*m_util;
		m_prob = p;
		m_util = u;

		return *this;
	}

	///
	/// \brief Combination operator (const)
	///
	valuation operator*(const valuation& a) const {
		valuation v;
		v.m_prob = m_prob*a.m_prob;
		v.m_util = m_prob*a.m_util + a.m_prob*m_util;

		return v;
	}

	///
	/// \brief Combination operator
	///
	valuation& operator*=(const valuation& a) {
		factor p = m_prob*a.m_prob;
		factor u = m_prob*a.m_util + a.m_prob*m_util;
		m_prob = p;
		m_util = u;

		return *this;
	}

	///
	/// Elimination by summation.
	///
	/// Sum out a subset of variables from the factor's scope.
	/// \param sum_out 	The set of variables to be summed-out
	/// \return a copy of the factor resulting from the variable elimination operation.
	///
	valuation sum(variable_set const& sum_out) const {
		variable_set t1 = m_prob.vars() - sum_out;
		variable_set t2 = m_util.vars() - sum_out;
		return valuation(m_prob.marginal(t1), m_util.marginal(t2));
	};

	///
	/// Elimination by weighted sum.
	///
	/// Sum out a subset of variables from the factor's scope.
	/// \param sum_out 	The set of variables to be summed-out
	/// \return a copy of the factor resulting from the variable elimination operation.
	///
	valuation sum_power(variable_set const& sum_out, double w) const {
		variable_set t1 = m_prob.vars() - sum_out;
		variable_set t2 = m_util.vars() - sum_out;
		return valuation(m_prob.marginal(t1, w), m_util.marginal(t2, w));
	};

	///
	/// \brief Elimination by maximization
	///
	valuation max(variable_set const& max_out) const {
		variable_set t1 = m_prob.vars() - max_out;
		variable_set t2 = m_util.vars() - max_out;
		return valuation(m_prob.maxmarginal(t1), m_util.maxmarginal(t2));
	}

	///
	/// \brief Scope of the valuation (union of probability and utility scopes)
	///
	variable_set vars() const {
		return (m_prob.vars() + m_util.vars());
	}

	///
	/// \brief Scope of the valuation (union of probability and utility scopes)
	///
	variable_set variables() const {
		return (m_prob.vars() + m_util.vars());
	}

	///
	/// \brief Size of the factor's scope.
	///
	/// \return the number of variables in the factor's scope.
	///
	const size_t nvar() const {
		return this->vars().size();
	};

	///
	/// \brief Size of the valuation's tables.
 	///
 	/// \return the size of the tables storing the valuation values.
 	///
	size_t numel() const {
		return (m_prob.numel() + m_util.numel());
	};

	///
	/// \brief Output operator (friend).
	///
	/// Write the (formatted) content of the valuation to an output stream.
	/// \param out
	///			the output stream
	/// \param f
	///			the valuation to be written out
	///
	/// \return a reference to the modified output stream containing the
	///		content of the valuation received as input.
	///
	friend std::ostream& operator<<(std::ostream& out, const valuation& f) {
		out << "Valuation over " << f.variables() << ": "
				<< f.m_prob << " | " << f.m_util;
		return out;
	};

};

} // end namespace

#endif /* IBM_MERLIN_VALUATION_H_ */
