#!/usr/bin/env python

#  Copyright (c) 2005 Gavin E. Crooks <gec@threeplusone.com>
#
#  This software is distributed under the MIT Open Source License.
#  <http://www.opensource.org/licenses/mit-license.html>
#
#  Permission is hereby granted, free of charge, to any person obtaining a
#  copy of this software and associated documentation files (the "Software"),
#  to deal in the Software without restriction, including without limitation
#  the rights to use, copy, modify, merge, publish, distribute, sublicense,
#  and/or sell copies of the Software, and to permit persons to whom the
#  Software is furnished to do so, subject to the following conditions:
#
#  The above copyright notice and this permission notice shall be included
#  in all copies or substantial portions of the Software.
#
#  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
#  IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
#  FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
#  AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
#  LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
#  OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
#  THE SOFTWARE.
#


""" Various bits of useful math not in the standard python library.

Special Functions :

- entropy()                     -- The entropy of a probability vector
- log2()                        -- Base 2 logarithms.

"""

__all__ = ('entropy', 'log2')

from math import log, exp


def log2(x):
    """ Return the base 2 logarithm of x """
    return log(x, 2)


def entropy(pvec, base=exp(1)):
    """ The entropy S = -Sum_i p_i ln p_i
        pvec is a frequency vector, not necessarily normalized.
    """
    # TODO: Optimize
    if len(pvec) == 0:
        raise ValueError("Zero length vector")

    total = 0.0
    ent = 0.0
    for p in pvec:
        if p > 0:  # 0 log(0) =0
            total += p
            ent += - log(float(p)) * p
        elif p < 0:
            raise ValueError("Negative probability")

    ent = (ent / total) + log(total)
    ent /= log(base)

    return ent
