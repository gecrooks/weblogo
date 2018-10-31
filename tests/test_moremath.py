#!/usr/bin/env python

#  Copyright (c) 2006, The Regents of the University of California, through
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

import unittest
from math import log

from weblogo.moremath import entropy


class test_misc_math(unittest.TestCase):
    def test_entropy(self):
        ent = entropy((1., 1.))
        self.assertAlmostEqual(ent, log(2))

    def test_entropy_with_flat_distribution(self):
        for n in range(1, 100):
            pvec = [1. / n for i in range(0, n)]
            ent = entropy(pvec)
            self.assertAlmostEqual(ent, log(n))

    def test_entropy_unnormalized(self):
        for n in range(1, 100):
            pvec = [1. for i in range(0, n)]
            ent = entropy(pvec)
            self.assertAlmostEqual(ent, log(n))

    def test_entropy_with_integer(self):
        ent = entropy((1, 1, 1, 0))
        self.assertAlmostEqual(ent, log(3))

    def test_entropy_with_short_pvec(self):
        ent = entropy((1,))
        self.assertEqual(ent, 0)
        ent = entropy((1, 0, 0, 0, 0))
        self.assertEqual(ent, 0)

    def test_entropy_invalid_pvec(self):
        self.assertRaises(ValueError, entropy, ())
        self.assertRaises(ValueError, entropy, (1, -1))

    def test_entropy_base(self):
        ent = entropy((2, 2, 2, 2, 0), 2)
        self.assertAlmostEqual(ent, 2)


if __name__ == '__main__':
    unittest.main()
