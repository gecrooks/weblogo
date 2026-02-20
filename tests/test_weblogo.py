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

from math import log, sqrt
from typing import Tuple

import numpy as np
import pytest
from scipy.stats import entropy

from weblogo import (
    LogoData,
    LogoFormat,
    LogoOptions,
    equiprobable_distribution,
    parse_prior,
)
from weblogo.color import Color
from weblogo.colorscheme import ColorScheme, IndexColor, RefSeqColor, SymbolColor
from weblogo.logomath import Dirichlet, Gamma
from weblogo.seq import (
    Alphabet,
    unambiguous_dna_alphabet,
    unambiguous_protein_alphabet,
    unambiguous_rna_alphabet,
)
from weblogo.utils import ArgumentError


def test_logoformat_options() -> None:
    LogoOptions()


def test_logoformat_errors() -> None:
    logodata = LogoData()
    logodata.length = 100

    # Negactive logo_margin
    with pytest.raises(ArgumentError):
        logooptions = LogoOptions()
        logooptions.logo_margin = -1
        LogoFormat(logodata, logooptions)

    # logo_start before start of sequecne
    with pytest.raises(ArgumentError):
        logooptions = LogoOptions()
        logooptions.first_index = 10
        logooptions.logo_start = -10
        LogoFormat(logodata, logooptions)

    # logo_end before logo_staRT
    with pytest.raises(ArgumentError):
        logooptions = LogoOptions()
        logooptions.first_index = 1
        logooptions.logo_start = 10
        logooptions.logo_end = -10
        LogoFormat(logodata, logooptions)

    # logo_end past lenght of sequence
    with pytest.raises(ArgumentError):
        logooptions = LogoOptions()
        logooptions.first_index = 1
        logooptions.logo_start = 10
        logooptions.logo_end = 200
        LogoFormat(logodata, logooptions)

    # No alphabet
    with pytest.raises(ArgumentError):
        logooptions = LogoOptions()
        logooptions.first_index = 1
        logooptions.logo_start = 10
        logooptions.logo_end = 20
        LogoFormat(logodata, logooptions)

    with pytest.raises(ArgumentError):
        logooptions = LogoOptions()
        logooptions.yaxis_scale = -1
        LogoFormat(logodata, logooptions)


def test_logoformats() -> None:
    # Make sure all different logo option code gets run
    logodata = LogoData()
    logodata.alphabet = unambiguous_rna_alphabet
    logodata.length = 100

    logooptions = LogoOptions()
    logooptions.fineprint = ""
    logooptions.xaxis_label = ""
    logooptions.yaxis_label = "Label"
    LogoFormat(logodata, logooptions)

    logooptions.yaxis_label = ""
    logooptions.unit_name = "probability"
    LogoFormat(logodata, logooptions)

    logooptions.show_yaxis = False
    LogoFormat(logodata, logooptions)

    logooptions.yaxis_label = "Label"
    logooptions.show_ends = True
    logooptions.show_xaxis = True
    LogoFormat(logodata, logooptions)

    logooptions.rotate_numbers = True
    LogoFormat(logodata, logooptions)

    logooptions.show_xaxis = False
    LogoFormat(logodata, logooptions)

    logodata.alphabet = Alphabet("ABCD")
    LogoFormat(logodata, logooptions)


def test_parse_prior_none() -> None:
    assert parse_prior(None, unambiguous_protein_alphabet) is None
    assert parse_prior("none", unambiguous_protein_alphabet) is None


def test_parse_prior_equiprobable() -> None:
    assert np.all(
        20.0 * equiprobable_distribution(20)
        == parse_prior(
            "equiprobable", unambiguous_protein_alphabet, weight=20.0
        )
    )

    assert np.all(
        1.2 * equiprobable_distribution(3)
        == parse_prior(" equiprobablE  ", Alphabet("123"), 1.2)
    )


def test_parse_prior_percentage() -> None:
    # print(parse_prior('50%', unambiguous_dna_alphabet, 1.))
    assert np.all(
        equiprobable_distribution(4)
        == parse_prior("50%", unambiguous_dna_alphabet, 1.0)
    )

    assert np.all(
        equiprobable_distribution(4)
        == parse_prior(" 50.0 % ", unambiguous_dna_alphabet, 1.0)
    )

    assert np.all(
        np.array((0.3, 0.2, 0.2, 0.3), np.float64)
        == parse_prior(" 40.0 % ", unambiguous_dna_alphabet, 1.0)
    )


def test_parse_prior_float() -> None:
    assert np.all(
        equiprobable_distribution(4)
        == parse_prior("0.5", unambiguous_dna_alphabet, 1.0)
    )

    assert np.all(
        equiprobable_distribution(4)
        == parse_prior(" 0.500 ", unambiguous_dna_alphabet, 1.0)
    )

    assert np.all(
        np.array((0.3, 0.2, 0.2, 0.3), np.float64)
        == parse_prior(" 0.40 ", unambiguous_dna_alphabet, 1.0)
    )


def test_parse_prior_auto() -> None:
    assert np.all(
        2.0 * equiprobable_distribution(4)
        == parse_prior("auto", unambiguous_dna_alphabet)
    )
    assert np.all(
        2.0 * equiprobable_distribution(4)
        == parse_prior("automatic", unambiguous_dna_alphabet)
    )

    parse_prior("automatic", unambiguous_protein_alphabet)
    parse_prior("E. coli", unambiguous_dna_alphabet)


def test_parse_prior_weight() -> None:
    assert np.all(
        2.0 * equiprobable_distribution(4)
        == parse_prior("automatic", unambiguous_dna_alphabet)
    )
    assert np.all(
        123.123 * equiprobable_distribution(4)
        == parse_prior("auto", unambiguous_dna_alphabet, 123.123)
    )


def test_parse_prior_explicit() -> None:
    s = "{'A':10, 'C':40, 'G':40, 'T':10}"
    p = np.array((10, 40, 40, 10), np.float64) * 2.0 / 100.0
    assert np.all(p == parse_prior(s, unambiguous_dna_alphabet))


def test_parse_prior_error() -> None:
    with pytest.raises(ValueError):
        parse_prior("0.5", unambiguous_protein_alphabet, weight=-10000.0)

    with pytest.raises(ValueError):
        s = "{'A':10, 'C':40, 'G':40, 'T':10}"
        parse_prior(s, unambiguous_protein_alphabet)

    with pytest.raises(ValueError):
        s = "{'A':'ljkasd', 'C':40, 'G':40, 'T':10}"
        parse_prior(s, unambiguous_dna_alphabet)

    with pytest.raises(ValueError):
        s = "asjnd"
        parse_prior(s, unambiguous_dna_alphabet)


def test_logooptions_create() -> None:
    opt = LogoOptions()
    opt.small_fontsize = 10
    repr(opt)

    opt = LogoOptions(logo_title="sometitle")  # type: ignore
    assert opt.logo_title == "sometitle"


def test_symbol_color() -> None:
    sc = SymbolColor("ABC", "black", "Because")
    assert sc.description == "Because"
    assert sc.symbol_color(0, "A", 0) == Color.by_name("black")
    assert sc.symbol_color(1, "D", 0) is None


def test_index_color() -> None:
    ic = IndexColor([1, 3], "black", "Because")
    assert ic.description == "Because"
    assert ic.symbol_color(0, "A", 0) is None
    assert ic.symbol_color(1, "A", 0) == Color.by_name("black")


def test_ref_seq_color() -> None:
    rc = RefSeqColor("ABC", "black", "Because")
    assert rc.description == "Because"

    assert rc.symbol_color(0, "A", 0) == Color.by_name("black")
    assert rc.symbol_color(1, "A", 0) is None
    assert rc.symbol_color(2, "A", 0) is None

    assert rc.symbol_color(0, "B", 0) is None
    assert rc.symbol_color(1, "B", 0) == Color.by_name("black")
    assert rc.symbol_color(2, "B", 0) is None

    assert rc.symbol_color(0, "C", 0) is None
    assert rc.symbol_color(1, "C", 0) is None
    assert rc.symbol_color(2, "C", 0) == Color.by_name("black")


def test_colorscheme() -> None:
    cs = ColorScheme(
        [
            SymbolColor("G", "orange"),
            SymbolColor("TU", "red"),
            SymbolColor("C", "blue"),
            SymbolColor("A", "green"),
        ],
        title="title",
        description="description",
    )

    assert cs.symbol_color(1, "G", 1) == Color.by_name("orange")
    assert cs.symbol_color(1, "T", 1) == Color.by_name("red")
    assert cs.symbol_color(1, "C", 1) == Color.by_name("blue")
    assert cs.symbol_color(1, "A", 1) == Color.by_name("green")
    assert cs.symbol_color(1, "X", 1) == cs.default_color

    cs = ColorScheme(
        [
            SymbolColor("G", "orange"),
            SymbolColor("TU", "red"),
            SymbolColor("C", "blue"),
            SymbolColor("A", "green"),
        ],
        title="title",
        description="description",
        alphabet=Alphabet("GTUCA"),
    )
    with pytest.raises(KeyError):
        cs.symbol_color(1, "X", 1)


def test_colorscheme_string_rejected() -> None:
    logodata = LogoData()
    logodata.alphabet = unambiguous_dna_alphabet
    logodata.length = 100
    options = LogoOptions()
    options.color_scheme = "monochrome"  # type: ignore[assignment]
    with pytest.raises(TypeError):
        LogoFormat(logodata, options)


def test_color_names() -> None:
    names = Color.names()
    assert len(names) == 147
    for n in names:
        c = Color.by_name(n)
        assert c is not None


def test_color_components() -> None:
    white = Color.by_name("white")
    assert 1.0 == white.red
    assert 1.0 == white.green
    assert 1.0 == white.blue

    c = Color(0.3, 0.4, 0.2)
    assert 0.3 == c.red
    assert 0.4 == c.green
    assert 0.2 == c.blue

    c = Color(0, 128, 0)
    assert 0.0 == c.red
    assert 128.0 / 255.0 == c.green
    assert 0.0 == c.blue


def test_color_from_rgb() -> None:
    white = Color.by_name("white")
    assert white == Color(1.0, 1.0, 1.0)
    assert white == Color(255, 255, 255)
    assert white == Color.from_rgb(1.0, 1.0, 1.0)
    assert white == Color.from_rgb(255, 255, 255)


def test_color_from_hsl() -> None:
    red = Color.by_name("red")
    lime = Color.by_name("lime")
    saddlebrown = Color.by_name("saddlebrown")
    darkgreen = Color.by_name("darkgreen")
    blue = Color.by_name("blue")
    Color.by_name("green")
    assert red == Color.from_hsl(0, 1.0, 0.5)
    assert lime == Color.from_hsl(120, 1.0, 0.5)
    assert blue == Color.from_hsl(240, 1.0, 0.5)
    assert Color.by_name("gray") == Color.from_hsl(0, 0, 0.5)
    assert saddlebrown == Color.from_hsl(25, 0.76, 0.31)
    assert darkgreen == Color.from_hsl(120, 1.0, 0.197)


def test_color_by_name() -> None:
    white = Color.by_name("white")
    assert white == Color.by_name("white")
    assert white == Color.by_name("WHITE")
    assert white == Color.by_name(" wHiTe \t\n\t")
    assert Color(255, 255, 240) == Color.by_name("ivory")
    assert Color(70, 130, 180) == Color.by_name("steelblue")
    assert Color(0, 128, 0) == Color.by_name("green")


def test_color_from_invalid_name() -> None:
    with pytest.raises(ValueError):
        Color.by_name("not_a_color")


def test_color_clipping() -> None:
    red = Color.by_name("red")
    assert red == Color(255, 0, 0)
    assert red == Color(260, -10, 0)
    assert red == Color(1.1, -0.0, -1.0)
    assert Color(1.0001, 213.0, 1.2).red == 1.0
    assert Color(-0.001, -2183.0, -1.0).red == 0.0
    assert Color(1.0001, 213.0, 1.2).green == 1.0
    assert Color(-0.001, -2183.0, -1.0).green == 0.0
    assert Color(1.0001, 213.0, 1.2).blue == 1.0
    assert Color(-0.001, -2183.0, -1.0).blue == 0.0


def test_color_fail_on_mixed_type() -> None:
    with pytest.raises(TypeError):
        Color.from_rgb(1, 1, 1.0)
    with pytest.raises(TypeError):
        Color.from_rgb(1.0, 1, 1.0)


def test_color_red() -> None:
    # Check Usage comment in Color
    red = Color.by_name("red")
    assert red == Color(255, 0, 0)
    assert red == Color(1.0, 0.0, 0.0)
    assert red == Color.from_rgb(1.0, 0.0, 0.0)
    assert red == Color.from_rgb(255, 0, 0)
    assert red == Color.from_hsl(0.0, 1.0, 0.5)
    assert red == Color.from_string("red")
    assert red == Color.from_string("RED")
    assert red == Color.from_string("#F00")
    assert red == Color.from_string("#FF0000")
    assert red == Color.from_string("rgb(255, 0, 0)")
    assert red == Color.from_string("rgb(100%, 0%, 0%)")
    assert red == Color.from_string("hsl(0, 100%, 50%)")


def test_color_from_string() -> None:
    Color(128, 0, 128)  # purple
    red = Color(255, 0, 0)
    skyblue = Color(135, 206, 235)

    red_strings = (
        "red",
        "ReD",
        "RED",
        "   Red \t",
        "#F00",
        "#FF0000",
        "rgb(255, 0, 0)",
        "rgb(100%, 0%, 0%)",
        "hsl(0, 100%, 50%)",
    )
    for s in red_strings:
        assert red == Color.from_string(s)

    skyblue_strings = (
        "skyblue",
        "SKYBLUE",
        "  \t\n SkyBlue  \t",
        "#87ceeb",
        "rgb(135,206,235)",
    )
    for s in skyblue_strings:
        assert skyblue == Color.from_string(s)

    with pytest.raises(ValueError):
        Color.from_string("#not_a_color")
    with pytest.raises(ValueError):
        Color.from_string("rgb(not_a_color)")
    with pytest.raises(ValueError):
        Color.from_string("hsl(not_a_color)")
    with pytest.raises(ValueError):
        Color.from_string("not_a_color")


def test_color_equality() -> None:
    c1 = Color(123, 99, 12)
    c2 = Color(123, 99, 12)
    assert c1 == c2
    assert c1 != "not_a_color"


def test_gamma_create() -> None:
    a = 1.213
    b = 3.210
    g = Gamma(a, b)
    assert g.alpha == a
    assert g.beta == b


def test_gamma_mean_variance() -> None:
    g = Gamma.from_mean_variance(2.0, 3.0)
    assert g.mean() == 2.0
    assert g.variance() == 3.0

    g = Gamma.from_mean_variance(2.0123, 3.01283)
    assert g.mean() == 2.0123
    assert g.variance() == 3.01283


def test_gamma_from_shape_scale() -> None:
    g = Gamma.from_shape_scale(1.0, 8.0)
    assert g.alpha == 1.0
    assert g.beta == 1.0 / 8.0


def test_gamma_invalid_args() -> None:
    with pytest.raises(ValueError):
        Gamma(1.0, -1.0)
    with pytest.raises(ValueError):
        Gamma(0.0, 1.0)
    with pytest.raises(ValueError):
        Gamma(1.0, 0.0)


def test_gamma_sample() -> None:
    m = 10.0
    v = 2.0
    g = Gamma.from_mean_variance(m, v)
    # print(g.alpha, g.beta)
    S = 1000
    total = 0.0
    for s in range(S):
        total += g.sample()
    mean = total / S

    # The estimated mean will differ from true mean by a small amount

    error = 4.0 * sqrt(g.variance() / S)
    # print(mean, m, error)
    assert abs(mean - m) < error


def test_gamma_pdf() -> None:
    m = 3.0
    v = 2.0
    g = Gamma.from_mean_variance(m, v)
    upper = 30.0

    norm = integrate(g.pdf, 0, upper)
    assert norm == pytest.approx(1.0)

    def fx(x: float) -> float:
        return x * g.pdf(x)

    mean = integrate(fx, 0, upper)
    assert mean == pytest.approx(m)

    def fx2(x: float) -> float:
        return x * x * g.pdf(x)

    x2 = integrate(fx2, 0, upper)
    var = x2 - mean**2
    assert var == pytest.approx(v)


def test_gamma_cdf() -> None:
    m = 3.0
    v = 2.0
    g = Gamma.from_mean_variance(m, v)
    # Numerical integration
    S = 1000
    M = 10.0
    total_p = 0.0
    epsilon = 1e-4
    last = 0.0
    for s in range(S):
        x = s * M / S
        p = g.pdf(x) * M / S
        total_p += (last - p) / 2.0
        last = p
        # print(x, total_p, g.cdf(x))

        assert (total_p - g.cdf(x)) < epsilon


def test_gamma_inverse_cdf() -> None:
    g = Gamma.from_mean_variance(2.34, 4)
    assert 3.9 == pytest.approx(g.inverse_cdf(g.cdf(3.9)))
    assert 1.92 == pytest.approx(g.inverse_cdf(g.cdf(1.92)))

    g = Gamma.from_mean_variance(10.34, 2)
    assert 3.9 == pytest.approx(g.inverse_cdf(g.cdf(3.9)))
    assert 10.92 == pytest.approx(g.inverse_cdf(g.cdf(10.92)))

    g = Gamma.from_mean_variance(10.34, 2)
    assert 0.975 == pytest.approx(g.cdf(g.inverse_cdf(0.975)))
    assert 0.025 == pytest.approx(g.cdf(g.inverse_cdf(0.025)))

    g = Gamma.from_mean_variance(1.34, 4)
    assert 0.975 == pytest.approx(g.cdf(g.inverse_cdf(0.975)))
    assert 0.025 == pytest.approx(g.cdf(g.inverse_cdf(0.025)))


def test_dirichlet_init() -> None:
    Dirichlet(
        (
            1,
            1,
            1,
            1,
        )
    )


def test_dirichlet_random() -> None:
    def do_test(alpha: Tuple[float, ...], samples: int = 1000) -> None:
        ent = np.zeros((samples,), np.float64)
        # alpha = ones( ( K,), Float64 ) * A/K

        # pt = zeros( (len(alpha) ,), Float64)
        d = Dirichlet(alpha)
        for s in range(samples):
            p = d.sample()
            # print(p)
            # pt +=p
            ent[s] = entropy(p)

        # print(pt/samples)

        m = mean(ent)
        v = var(ent)

        dm = d.mean_entropy()
        dv = d.variance_entropy()

        # print(alpha, ':', m, v, dm, dv)
        error = 4.0 * sqrt(v / samples)
        assert abs(m - dm) < error
        assert abs(v - dv) < error  # dodgy error estimate

    do_test((1.0, 1.0))
    do_test((2.0, 1.0))
    do_test((3.0, 1.0))
    do_test((4.0, 1.0))
    do_test((5.0, 1.0))
    do_test((6.0, 1.0))

    do_test((1.0, 1.0))
    do_test((20.0, 20.0))
    do_test((1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0))
    do_test((0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1))
    do_test((0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01))
    do_test((2.0, 6.0, 1.0, 1.0))


def test_dirichlet_mean() -> None:
    alpha = np.ones((10,), np.float64) * 23.0
    d = Dirichlet(alpha)
    m = d.mean()
    assert m[2] == pytest.approx(1.0 / 10)
    assert sum(m) == pytest.approx(1.0)


def test_dirichlet_covariance() -> None:
    alpha = np.ones((4,), np.float64)
    d = Dirichlet(alpha)
    cv = d.covariance()
    assert cv.shape == (4, 4)
    assert cv[0, 0] == pytest.approx(1.0 * (1.0 - 1.0 / 4.0) / (4.0 * 5.0))
    assert cv[0, 1] == pytest.approx(-1 / (4.0 * 4.0 * 5.0))


def test_dirichlet_mean_x() -> None:
    alpha = (1.0, 2.0, 3.0, 4.0)
    xx = (2.0, 2.0, 2.0, 2.0)
    m = Dirichlet(alpha).mean_x(xx)
    assert m == 2.0

    xx2 = (2.0, 2.0, 2.0, 2.0, 2.0)
    with pytest.raises(ValueError):
        Dirichlet(alpha).mean_x(xx2)

    alpha = (1.0, 1.0, 1.0, 1.0)
    xx = (2.0, 3.0, 4.0, 3.0)
    m = Dirichlet(alpha).mean_x(xx)
    assert m == 3.0


def test_dirichlet_variance_x() -> None:
    alpha = (1.0, 1.0, 1.0, 1.0)
    xx = (2.0, 2.0, 2.0, 2.0)
    v = Dirichlet(alpha).variance_x(xx)
    assert v == pytest.approx(0.0)

    alpha = (1.0, 2.0, 3.0, 4.0)
    xx = (2.0, 0.0, 1.0, 10.0)
    v = Dirichlet(alpha).variance_x(xx)
    # print(v)
    # TODO: Don't actually know if this is correct

    xx2 = (2.0, 2.0, 2.0, 2.0, 2.0)
    with pytest.raises(ValueError):
        Dirichlet(alpha).variance_x(xx2)


def test_dirichlet_relative_entropy() -> None:
    alpha = (2.0, 10.0, 1.0, 1.0)
    d = Dirichlet(alpha)
    pvec = (0.1, 0.2, 0.3, 0.4)

    rent = d.mean_relative_entropy(pvec)
    vrent = d.variance_relative_entropy(pvec)
    low, high = d.interval_relative_entropy(pvec, 0.95)

    # print()
    # print('> ', rent, vrent, low, high)

    # This test can fail randomly, but the precision from a few
    # thousand samples is low. Increasing samples, 1000->2000
    samples = 2000
    sent = np.zeros((samples,), np.float64)

    for s in range(samples):
        post = d.sample()
        e = -entropy(post)
        for k in range(4):
            e += -post[k] * log(pvec[k])
        sent[s] = e
    sent.sort()
    assert abs(sent.mean() - rent) < 4.0 * sqrt(vrent)
    assert sent.std() == pytest.approx(sqrt(vrent), abs=0.05)
    assert abs(low - sent[int(samples * 0.025)]) < 0.2
    assert abs(high - sent[int(samples * 0.975)]) < 0.2


def test_from_URL_fileopen_URLscheme() -> None:
    """test for http, https, or ftp scheme"""
    from weblogo.logo import _from_URL_fileopen

    broken_url = "file://foo.txt"
    with pytest.raises(ValueError):
        _from_URL_fileopen(broken_url)


def mean(a) -> float:  # type: ignore
    return sum(a) / len(a)


def var(a) -> float:  # type: ignore
    return (sum(a * a) / len(a)) - mean(a) ** 2


def integrate(f, a, b, n=1000):  # type: ignore
    """
    Numerically integrate the function 'f' from 'a' to 'b' using a discretization with 'n' points.

    Args:
    - f -- A function that eats a float and returns a float.
    - a -- Lower integration bound (float)
    - b -- Upper integration bound (float)
    - n -- number of sample points (int)

    Status :
        Alpha (very primitive.)
    """
    h = (b - a) / (n - 1.0)
    total = 0.0
    for i in range(n):
        total += f(a + (i) * h)
    result = h * (total - 0.5 * f(a) - 0.5 * f(b))
    return result


def test_pdf_formatter() -> None:
    """Test that the PDF formatter produces valid PDF output."""
    from weblogo.logo_formatter import pdf_formatter
    from weblogo.seq import Seq, SeqList

    seqs = SeqList(
        [Seq("AACGTAG"), Seq("AAGGTAC"), Seq("AACGTAG"), Seq("GAAGTAC")],
        unambiguous_dna_alphabet,
    )
    logodata = LogoData.from_seqs(seqs)
    logooptions = LogoOptions()
    logooptions.logo_title = "Test"
    logoformat = LogoFormat(logodata, logooptions)

    pdf = pdf_formatter(logodata, logoformat)
    assert isinstance(pdf, bytes)
    assert len(pdf) > 0
    assert pdf[:5] == b"%PDF-"


def test_txt_formatter() -> None:
    """Test that the text formatter produces output."""
    from weblogo.logo_formatter import txt_formatter
    from weblogo.seq import Seq, SeqList

    seqs = SeqList(
        [Seq("AACGTAG"), Seq("AAGGTAC"), Seq("AACGTAG"), Seq("GAAGTAC")],
        unambiguous_dna_alphabet,
    )
    logodata = LogoData.from_seqs(seqs)
    logoformat = LogoFormat(logodata, LogoOptions())

    txt = txt_formatter(logodata, logoformat)
    assert isinstance(txt, bytes)
    assert len(txt) > 0


def test_csv_formatter() -> None:
    """Test that the CSV formatter produces output."""
    from weblogo.logo_formatter import csv_formatter
    from weblogo.seq import Seq, SeqList

    seqs = SeqList(
        [Seq("AACGTAG"), Seq("AAGGTAC"), Seq("AACGTAG"), Seq("GAAGTAC")],
        unambiguous_dna_alphabet,
    )
    logodata = LogoData.from_seqs(seqs)
    logoformat = LogoFormat(logodata, LogoOptions())

    csv = csv_formatter(logodata, logoformat)
    assert isinstance(csv, bytes)
    assert b"," in csv
