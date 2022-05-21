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

from weblogo.utils import (
    ArgumentError,
    Token,
    crc32,
    crc64,
    group_count,
    invert_dict,
    isblank,
    isfloat,
    isint,
    remove_whitespace,
    resource_filename,
    resource_stream,
    resource_string,
)


class test_utils(unittest.TestCase):
    def test_isfloat(self) -> None:
        self.assertTrue(isfloat("0.5"))
        self.assertTrue(isfloat(" 0"))
        self.assertTrue(isfloat("+1000000000  "))
        self.assertTrue(isfloat("2"))
        self.assertTrue(isfloat("0000.2323"))
        self.assertTrue(isfloat("0.1e-23"))
        self.assertTrue(isfloat(" -0.5e+23"))
        self.assertFalse(isfloat(None))
        self.assertFalse(isfloat(""))
        self.assertFalse(isfloat("asdad"))
        self.assertFalse(isfloat("q34sd"))
        self.assertFalse(isfloat("92384.kjdfghiksw"))
        self.assertFalse(isfloat("adf!@#nn"))

    def test_isint(self) -> None:
        self.assertTrue(isint("0"))
        self.assertTrue(isint("-1"))
        self.assertTrue(isint("10"))
        self.assertTrue(isint("100101012234"))
        self.assertTrue(isint("000"))
        self.assertFalse(isint(None))
        self.assertFalse(isint(""))
        self.assertFalse(isint("asdad"))
        self.assertFalse(isint("q34sd"))
        self.assertFalse(isint("0.23"))
        self.assertFalse(isint("adf!@#nn"))

    def test_remove_whitespace(self) -> None:
        self.assertEqual(
            remove_whitespace("  kjashd askjdh askjdh\tasdf"), "kjashdaskjdhaskjdhasdf"
        )

    def test_isblank(self) -> None:
        blank = ("", " ", "\n", "\t \n\n")
        not_blank = (" a",)
        for s in blank:
            self.assertTrue(isblank(s))
        for s in not_blank:
            self.assertFalse(isblank(s))

        self.assertFalse(isblank(123))

    def test_group_count(self) -> None:
        test = "aaabbbbcccddea"
        out = group_count(test)
        self.assertEqual(
            tuple(out), (("a", 3), ("b", 4), ("c", 3), ("d", 2), ("e", 1), ("a", 1))
        )

    def test_token(self) -> None:
        t = Token("kind", "some data", 4, 3)
        str(t)
        r = repr(t)
        t2 = eval(r)
        self.assertEqual(t2.typeof, "kind")

    def test_invert_dict(self) -> None:
        d = dict(a=3, b=4)
        invd = invert_dict(d)
        self.assertTrue(3 in invd)
        self.assertEqual(invd[3], "a")

    def test_crc64(self) -> None:
        self.assertEqual(crc64("IHATEMATH"), "E3DCADD69B01ADD1")

    def test_crc32(self) -> None:
        self.assertEqual(crc32("Test the CRC-32 of this string."), "%08X" % 1571220330)

    def test_ArgumentValueError(self) -> None:
        message = "Some message"
        component = "whatsit"
        try:
            raise ArgumentError(message, component)
        except ArgumentError as err:
            self.assertEqual(err.msg, message)
            self.assertEqual(err.key, component)
        try:
            raise ArgumentError(message, component, 10)
        except ArgumentError as err:
            self.assertEqual(err.msg, message)
            self.assertEqual(err.key, component)
            self.assertEqual(err.value, 10)

    def test_resource(self) -> None:
        fn = resource_filename(__name__, "data/cap.fa", __file__)
        self.assertTrue(fn.endswith("data/cap.fa"))
        f = resource_stream(__name__, "data/cap.fa", __file__)
        f.close()
        s = resource_string(__name__, "data/cap.fa", __file__).decode()
        self.assertTrue(s.startswith(">aldB"))


tfile = """line 0
line 1
Blah
line 3
line 4
"""

if __name__ == "__main__":
    unittest.main()
