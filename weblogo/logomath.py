# -------------------------------- WebLogo --------------------------------

#  Copyright (c) 2003-2004 The Regents of the University of California.
#  Copyright (c) 2005 Gavin E. Crooks
#  Copyright (c) 2006-2011, The Regents of the University of California, through
#  Lawrence Berkeley National Laboratory (subject to receipt of any required
#  approvals from the U.S. Dept. of Energy).  All rights reserved.

#  This software is distributed under the new BSD Open Source License.
#  <http://www.opensource.org/licenses/bsd-license.html>
#
#  Redistribution and use in source and binary forms, with or without
#  modification, are permitted provided that the following conditions are met:
#
#  (1) Redistributions of source code must retain the above copyright notice,
#  this list of conditions and the following disclaimer.
#
#  (2) Redistributions in binary form must reproduce the above copyright
#  notice, this list of conditions and the following disclaimer in the
#  documentation and or other materials provided with the distribution.
#
#  (3) Neither the name of the University of California, Lawrence Berkeley
#  National Laboratory, U.S. Dept. of Energy nor the names of its contributors
#  may be used to endorse or promote products derived from this software
#  without specific prior written permission.
#
#  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
#  AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
#  IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
#  ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
#  LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
#  CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
#  SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
#  INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
#  CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
#  ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
#  POSSIBILITY OF SUCH DAMAGE.

import random
from math import exp, log, sqrt
from typing import Tuple

import numpy as np
import scipy.optimize
from numpy import asarray, float64, shape, zeros
from scipy.special import digamma, gamma, gammaincc, polygamma


class Dirichlet(object):
    """The Dirichlet probability distribution. The Dirichlet is a continuous
    multivariate probability distribution across non-negative unit length
    vectors. In other words, the Dirichlet is a probability distribution of
    probability distributions. It is conjugate to the multinomial
    distribution and is widely used in Bayesian statistics.

    The Dirichlet probability distribution of order K-1 is

     p(theta_1,...,theta_K) d theta_1 ... d theta_K =
        (1/Z) prod_i=1,K theta_i^{alpha_i - 1} delta(1 -sum_i=1,K theta_i)

    The normalization factor Z can be expressed in terms of gamma functions:

      Z = {prod_i=1,K Gamma(alpha_i)} / {Gamma( sum_i=1,K alpha_i)}

    The K constants, alpha_1,...,alpha_K, must be positive. The K parameters,
    theta_1,...,theta_K are nonnegative and sum to 1.

    Status:
        Alpha
    """

    __slots__ = (
        "alpha",
        "_total",
        "_mean",
    )

    def __init__(self, alpha: np.ndarray) -> None:
        """
        Args:
            - alpha  -- The parameters of the Dirichlet prior distribution.
                        A vector of non-negative real numbers.
        """
        # TODO: Check that alphas are positive
        # TODO: what if alpha's not one dimensional?
        self.alpha = asarray(alpha, float64)

        self._total = sum(alpha)
        self._mean = self.alpha / self._total

    def sample(self) -> np.ndarray:
        """Return a randomly generated probability vector.

        Random samples are generated by sampling K values from gamma
        distributions with parameters a=\alpha_i, b=1, and renormalizing.

        Ref:
            A.M. Law, W.D. Kelton, Simulation Modeling and Analysis (1991).
        Authors:
            Gavin E. Crooks <gec@compbio.berkeley.edu> (2002)
        """
        alpha = self.alpha
        K = len(alpha)
        theta = zeros((K,), float64)

        for k in range(K):
            theta[k] = random.gammavariate(alpha[k], 1.0)
        theta /= sum(theta)

        return theta

    def mean(self) -> np.ndarray:
        return self._mean

    def covariance(self):
        alpha = self.alpha
        A = sum(alpha)
        # A2 = A * A
        K = len(alpha)
        cv = zeros((K, K), float64)

        for i in range(K):
            cv[i, i] = alpha[i] * (1.0 - alpha[i] / A) / (A * (A + 1.0))

        for i in range(K):
            for j in range(i + 1, K):
                v = -alpha[i] * alpha[j] / (A * A * (A + 1.0))
                cv[i, j] = v
                cv[j, i] = v
        return cv

    def mean_x(self, x: np.ndarray) -> float:
        x = asarray(x, float64)
        if shape(x) != shape(self.alpha):
            raise ValueError("Argument must be same dimension as Dirichlet")
        return sum(x * self.mean())

    def variance_x(self, x: np.ndarray) -> float:
        x = asarray(x, float64)
        if shape(x) != shape(self.alpha):
            raise ValueError("Argument must be same dimension as Dirichlet")

        cv = self.covariance()
        var = np.dot(np.dot(np.transpose(x), cv), x)
        return var

    def mean_entropy(self) -> float:
        """Calculate the average entropy of probabilities sampled
        from this Dirichlet distribution.

        Returns:
            The average entropy.

        Ref:
            Wolpert & Wolf, PRE 53:6841-6854 (1996) Theorem 7
            (Warning: this paper contains typos.)
        Status:
            Alpha
        Authors:
            GEC 2005

        """
        # TODO: Optimize
        alpha = self.alpha
        A = float(sum(alpha))
        ent = 0.0
        for a in alpha:
            if a > 0:
                ent += -1.0 * a * digamma(1.0 + a)  # FIXME: Check
        ent /= A
        ent += digamma(A + 1.0)
        return ent

    def variance_entropy(self) -> float:
        """Calculate the variance of the Dirichlet entropy.

        Ref:
            Wolpert & Wolf, PRE 53:6841-6854 (1996) Theorem 8
            (Warning: this paper contains typos.)
        """
        alpha = self.alpha
        A = float(sum(alpha))
        A2 = A * (A + 1)
        L = len(alpha)

        dg1 = zeros((L), float64)
        dg2 = zeros((L), float64)
        tg2 = zeros((L), float64)

        for i in range(L):
            dg1[i] = digamma(alpha[i] + 1.0)
            dg2[i] = digamma(alpha[i] + 2.0)
            tg2[i] = polygamma(1, alpha[i] + 2.0)

        dg_Ap2 = digamma(A + 2.0)
        tg_Ap2 = polygamma(1, A + 2.0)

        mean = self.mean_entropy()
        var = 0.0

        for i in range(L):
            for j in range(L):
                if i != j:
                    var += (
                        ((dg1[i] - dg_Ap2) * (dg1[j] - dg_Ap2) - tg_Ap2)
                        * (alpha[i] * alpha[j])
                        / A2
                    )
                else:
                    var += (
                        ((dg2[i] - dg_Ap2) ** 2 + (tg2[i] - tg_Ap2))
                        * (alpha[i] * (alpha[i] + 1.0))
                        / A2
                    )

        var -= mean ** 2
        return var

    def mean_relative_entropy(self, pvec: np.ndarray) -> float:
        ln_p = np.log(pvec)
        return -self.mean_x(ln_p) - self.mean_entropy()

    def variance_relative_entropy(self, pvec: np.ndarray) -> float:
        ln_p = np.log(pvec)
        return self.variance_x(ln_p) + self.variance_entropy()

    def interval_relative_entropy(
        self, pvec: np.ndarray, frac: float
    ) -> Tuple[float, float]:
        mean = self.mean_relative_entropy(pvec)
        variance = self.variance_relative_entropy(pvec)
        sd = sqrt(variance)

        # If the variance is small, use the standard 95%
        # confidence interval: mean +/- 1.96 * sd
        if mean / sd > 3.0:
            return max(0.0, mean - sd * 1.96), mean + sd * 1.96

        g = Gamma.from_mean_variance(mean, variance)
        low_limit = g.inverse_cdf((1.0 - frac) / 2.0)
        high_limit = g.inverse_cdf(1.0 - (1.0 - frac) / 2.0)

        return low_limit, high_limit


class Gamma(object):
    """The gamma probability distribution. (Not to be confused with the
    gamma function.)


    """

    __slots__ = "alpha", "beta"

    def __init__(self, alpha: float, beta: float) -> None:
        if alpha <= 0.0:
            raise ValueError("alpha must be positive")
        if beta <= 0.0:
            raise ValueError("beta must be positive")
        self.alpha = alpha
        self.beta = beta

    @classmethod
    def from_shape_scale(cls, shape: float, scale: float) -> "Gamma":
        return cls(shape, 1.0 / scale)

    @classmethod
    def from_mean_variance(cls, mean: float, variance: float) -> "Gamma":
        alpha = mean ** 2 / variance
        beta = alpha / mean
        return cls(alpha, beta)

    def mean(self) -> float:
        return self.alpha / self.beta

    def variance(self) -> float:
        return self.alpha / (self.beta ** 2)

    def sample(self) -> float:
        return random.gammavariate(self.alpha, 1.0 / self.beta)

    def pdf(self, x: float) -> float:
        if x == 0.0:
            return 0.0
        a = self.alpha
        b = self.beta
        return (x ** (a - 1.0)) * exp(-b * x) * (b ** a) / gamma(a)

    def cdf(self, x: float) -> float:
        return 1.0 - gammaincc(self.alpha, self.beta * x)

    def inverse_cdf(self, p: float) -> float:
        def rootof(x: float) -> float:
            return self.cdf(exp(x)) - p

        root = scipy.optimize.newton(rootof, log(self.mean()))
        return exp(root)
