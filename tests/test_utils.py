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

from weblogo.utils import (
    ArgumentError,
    Token,
    group_count,
    invert_dict,
    isblank,
    isfloat,
    isint,
    remove_whitespace,
)


def test_isfloat() -> None:
    assert isfloat("0.5")
    assert isfloat(" 0")
    assert isfloat("+1000000000  ")
    assert isfloat("2")
    assert isfloat("0000.2323")
    assert isfloat("0.1e-23")
    assert isfloat(" -0.5e+23")
    assert not isfloat(None)
    assert not isfloat("")
    assert not isfloat("asdad")
    assert not isfloat("q34sd")
    assert not isfloat("92384.kjdfghiksw")
    assert not isfloat("adf!@#nn")


def test_isint() -> None:
    assert isint("0")
    assert isint("-1")
    assert isint("10")
    assert isint("100101012234")
    assert isint("000")
    assert not isint(None)
    assert not isint("")
    assert not isint("asdad")
    assert not isint("q34sd")
    assert not isint("0.23")
    assert not isint("adf!@#nn")


def test_remove_whitespace() -> None:
    assert remove_whitespace("  kjashd askjdh askjdh\tasdf") == "kjashdaskjdhaskjdhasdf"


def test_isblank() -> None:
    blank = ("", " ", "\n", "\t \n\n")
    not_blank = (" a",)
    for s in blank:
        assert isblank(s)
    for s in not_blank:
        assert not isblank(s)

    assert not isblank(123)


def test_group_count() -> None:
    test = "aaabbbbcccddea"
    out = group_count(test)
    assert tuple(out) == (("a", 3), ("b", 4), ("c", 3), ("d", 2), ("e", 1), ("a", 1))


def test_token() -> None:
    t = Token("kind", "some data", 4, 3)
    str(t)
    r = repr(t)
    t2 = eval(r)
    assert t2.typeof == "kind"


def test_token_default_offset() -> None:
    t = Token("kind", "some data", 4)
    s = str(t)
    assert ":" not in s.split()[0]  # no offset in coord


def test_stdrepr_explicit_name() -> None:
    from weblogo.utils import stdrepr

    t = Token("kind", "some data", 4, 3)
    r = stdrepr(t, name="CustomName")
    assert r.startswith("CustomName(")


def test_invert_dict() -> None:
    d = dict(a=3, b=4)
    invd = invert_dict(d)
    assert 3 in invd
    assert invd[3] == "a"


def test_ArgumentValueError() -> None:
    message = "Some message"
    component = "whatsit"
    try:
        raise ArgumentError(message, component)
    except ArgumentError as err:
        assert err.msg == message
        assert err.key == component
    try:
        raise ArgumentError(message, component, 10)
    except ArgumentError as err:
        assert err.msg == message
        assert err.key == component
        assert err.value == 10


tfile = """line 0
line 1
Blah
line 3
line 4
"""
