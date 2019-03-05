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
import pytest

from math import log, sqrt
from pkg_resources import resource_stream

from numpy import array, float64, ones, zeros, all

from weblogo import LogoOptions, equiprobable_distribution, LogoData, LogoFormat
from weblogo import parse_prior, GhostscriptAPI
from weblogo.color import Color
from weblogo.colorscheme import ColorScheme, RefSeqColor, SymbolColor, IndexColor
from weblogo.logomath import Dirichlet, Gamma
from weblogo.seq import (Alphabet, unambiguous_protein_alphabet, unambiguous_dna_alphabet)
from scipy.stats import entropy
from weblogo.utils import ArgumentError
from weblogo.seq import unambiguous_rna_alphabet


def data_stream(name):
    return resource_stream(__name__, 'data/' + name)


class test_logoformat(unittest.TestCase):
    def test_options(self):
        LogoOptions()


def test_logoformat_errors():
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


def test_logoformats():
    # Make sure all different logo option code gets run
    logodata = LogoData()
    logodata.alphabet = unambiguous_rna_alphabet
    logodata.length = 100

    logooptions = LogoOptions()
    logooptions.fineprint = None
    logooptions.xaxis_label = True
    logooptions.yaxis_label = "Label"
    LogoFormat(logodata, logooptions)

    logooptions.yaxis_label = ''
    logooptions.unit_name = 'probability'
    LogoFormat(logodata, logooptions)

    logooptions.show_yaxis = False
    LogoFormat(logodata, logooptions)

    logooptions.yaxis_label = 'Label'
    logooptions.show_ends = True
    logooptions.show_xaxis = True
    LogoFormat(logodata, logooptions)

    logooptions.rotate_numbers = True
    LogoFormat(logodata, logooptions)

    logooptions.show_xaxis = False
    LogoFormat(logodata, logooptions)

    logodata.alphabet = "ABCD"
    LogoFormat(logodata, logooptions)


class test_ghostscript(unittest.TestCase):
    def test_version(self):
        GhostscriptAPI().version()


class test_parse_prior(unittest.TestCase):
    def test_parse_prior_none(self):
        self.assertEqual(None,
                         parse_prior(None, unambiguous_protein_alphabet))
        self.assertEqual(None,
                         parse_prior('none', unambiguous_protein_alphabet))
        self.assertEqual(None,
                         parse_prior('noNe', None))

    def test_parse_prior_equiprobable(self):
        self.assertTrue(all(20. * equiprobable_distribution(20) ==
                            parse_prior('equiprobable', unambiguous_protein_alphabet, weight=20.)))

        self.assertTrue(
                all(1.2 * equiprobable_distribution(3)
                    == parse_prior(' equiprobablE  ', Alphabet('123'), 1.2)))

    def test_parse_prior_percentage(self):
        # print(parse_prior('50%', unambiguous_dna_alphabet, 1.))
        self.assertTrue(all(equiprobable_distribution(4)
                            == parse_prior('50%', unambiguous_dna_alphabet, 1.)))

        self.assertTrue(all(equiprobable_distribution(4)
                            == parse_prior(' 50.0 % ', unambiguous_dna_alphabet, 1.)))

        self.assertTrue(all(array((0.3, 0.2, 0.2, 0.3), float64)
                            == parse_prior(' 40.0 % ', unambiguous_dna_alphabet, 1.)))

    def test_parse_prior_float(self):
        self.assertTrue(all(equiprobable_distribution(4)
                            == parse_prior('0.5', unambiguous_dna_alphabet, 1.)))

        self.assertTrue(all(equiprobable_distribution(4)
                            == parse_prior(' 0.500 ', unambiguous_dna_alphabet, 1.)))

        self.assertTrue(all(array((0.3, 0.2, 0.2, 0.3), float64)
                            == parse_prior(' 0.40 ', unambiguous_dna_alphabet, 1.)))

    def test_auto(self):
        self.assertTrue(all(2. * equiprobable_distribution(4) ==
                            parse_prior('auto', unambiguous_dna_alphabet)))
        self.assertTrue(all(2. * equiprobable_distribution(4) ==
                            parse_prior('automatic', unambiguous_dna_alphabet)))

        parse_prior('automatic', unambiguous_protein_alphabet)
        parse_prior('E. coli', unambiguous_dna_alphabet)

    def test_weight(self):
        self.assertTrue(all(2. * equiprobable_distribution(4) ==
                            parse_prior('automatic', unambiguous_dna_alphabet)))
        self.assertTrue(all(123.123 * equiprobable_distribution(4) ==
                            parse_prior('auto', unambiguous_dna_alphabet, 123.123)))

    def test_explicit(self):
        s = "{'A':10, 'C':40, 'G':40, 'T':10}"
        p = array((10, 40, 40, 10), float64) * 2. / 100.
        self.assertTrue(all(
                p == parse_prior(s, unambiguous_dna_alphabet)))


def test_parse_prior_error():
    with pytest.raises(ValueError):
        parse_prior('0.5', unambiguous_protein_alphabet, weight=-10000.0)

    with pytest.raises(ValueError):
        s = "{'A':10, 'C':40, 'G':40, 'T':10}"
        parse_prior(s, unambiguous_protein_alphabet)

    with pytest.raises(ValueError):
        s = "{'A':'ljkasd', 'C':40, 'G':40, 'T':10}"
        parse_prior(s, unambiguous_dna_alphabet)

    with pytest.raises(ValueError):
        s = "asjnd"
        parse_prior(s, unambiguous_dna_alphabet)


class test_logooptions(unittest.TestCase):
    def test_create(self):
        opt = LogoOptions()
        opt.small_fontsize = 10
        repr(opt)

        opt = LogoOptions(title="sometitle")
        self.assertEqual(opt.title, "sometitle")


class test_colorscheme(unittest.TestCase):
    def test_symbol_color(self):
        sc = SymbolColor("abc", "black", "Because")
        self.assertEqual(sc.description, "Because")
        self.assertEqual(sc.symbol_color(0, "A", 0), Color.by_name("black"))
        self.assertEqual(sc.symbol_color(1, "D", 0), None)

    def test_index_color(self):
        ic = IndexColor([1, 3], "black", "Because")
        self.assertEqual(ic.description, "Because")
        self.assertEqual(ic.symbol_color(0, "A", 0), None)
        self.assertEqual(ic.symbol_color(1, "A", 0), Color.by_name("black"))

    def test_ref_seq_color(self):
        rc = RefSeqColor("abc", "black", "Because")
        self.assertEqual(rc.description, "Because")

        self.assertEqual(rc.symbol_color(0, "A", 0), Color.by_name("black"))
        self.assertEqual(rc.symbol_color(1, "A", 0), None)
        self.assertEqual(rc.symbol_color(2, "A", 0), None)

        self.assertEqual(rc.symbol_color(0, "B", 0), None)
        self.assertEqual(rc.symbol_color(1, "B", 0), Color.by_name("black"))
        self.assertEqual(rc.symbol_color(2, "B", 0), None)

        self.assertEqual(rc.symbol_color(0, "C", 0), None)
        self.assertEqual(rc.symbol_color(1, "C", 0), None)
        self.assertEqual(rc.symbol_color(2, "C", 0), Color.by_name("black"))

    def test_colorscheme(self):
        cs = ColorScheme([
            SymbolColor("G", "orange"),
            SymbolColor("TU", "red"),
            SymbolColor("C", "blue"),
            SymbolColor("A", "green")
        ],
                title="title",
                description="description",
        )

        self.assertEqual(cs.symbol_color(1, 'G', 1), Color.by_name("orange"))
        self.assertEqual(cs.symbol_color(1, 'T', 1), Color.by_name("red"))
        self.assertEqual(cs.symbol_color(1, 'C', 1), Color.by_name("blue"))
        self.assertEqual(cs.symbol_color(1, 'A', 1), Color.by_name("green"))
        self.assertEqual(cs.symbol_color(1, 'X', 1), cs.default_color)

        cs = ColorScheme([
            SymbolColor("G", "orange"),
            SymbolColor("TU", "red"),
            SymbolColor("C", "blue"),
            SymbolColor("A", "green")
        ],
                title="title",
                description="description",
                alphabet="GTUCA"
        )
        self.assertRaises(KeyError, cs.symbol_color, 1, 'X', 1)


class test_color(unittest.TestCase):
    def test_color_names(self):
        names = Color.names()
        self.assertEqual(len(names), 147)
        for n in names:
            c = Color.by_name(n)
            self.assertTrue(c is not None)

    def test_color_components(self):
        white = Color.by_name("white")
        self.assertEqual(1.0, white.red)
        self.assertEqual(1.0, white.green)
        self.assertEqual(1.0, white.blue)

        c = Color(0.3, 0.4, 0.2)
        self.assertEqual(0.3, c.red)
        self.assertEqual(0.4, c.green)
        self.assertEqual(0.2, c.blue)

        c = Color(0, 128, 0)
        self.assertEqual(0.0, c.red)
        self.assertEqual(128. / 255., c.green)
        self.assertEqual(0.0, c.blue)

    def test_color_from_rgb(self):
        white = Color.by_name("white")
        self.assertEqual(white, Color(1., 1., 1.))
        self.assertEqual(white, Color(255, 255, 255))
        self.assertEqual(white, Color.from_rgb(1., 1., 1.))
        self.assertEqual(white, Color.from_rgb(255, 255, 255))

    def test_color_from_hsl(self):
        red = Color.by_name("red")
        lime = Color.by_name("lime")
        saddlebrown = Color.by_name("saddlebrown")
        darkgreen = Color.by_name("darkgreen")
        blue = Color.by_name("blue")
        Color.by_name("green")
        self.assertEqual(red, Color.from_hsl(0, 1.0, 0.5))
        self.assertEqual(lime, Color.from_hsl(120, 1.0, 0.5))
        self.assertEqual(blue, Color.from_hsl(240, 1.0, 0.5))
        self.assertEqual(Color.by_name("gray"), Color.from_hsl(0, 0, 0.5))
        self.assertEqual(saddlebrown, Color.from_hsl(25, 0.76, 0.31))
        self.assertEqual(darkgreen, Color.from_hsl(120, 1.0, 0.197))

    def test_color_by_name(self):
        white = Color.by_name("white")
        self.assertEqual(white, Color.by_name("white"))
        self.assertEqual(white, Color.by_name("WHITE"))
        self.assertEqual(white, Color.by_name(" wHiTe \t\n\t"))
        self.assertEqual(Color(255, 255, 240), Color.by_name("ivory"))
        self.assertEqual(Color(70, 130, 180), Color.by_name("steelblue"))
        self.assertEqual(Color(0, 128, 0), Color.by_name("green"))

    def test_color_from_invalid_name(self):
        self.assertRaises(ValueError, Color.by_name, "not_a_color")

    def test_color_clipping(self):
        red = Color.by_name("red")
        self.assertEqual(red, Color(255, 0, 0))
        self.assertEqual(red, Color(260, -10, 0))
        self.assertEqual(red, Color(1.1, -0., -1.))
        self.assertEqual(Color(1.0001, 213.0, 1.2).red, 1.0)
        self.assertEqual(Color(-0.001, -2183.0, -1.0).red, 0.0)
        self.assertEqual(Color(1.0001, 213.0, 1.2).green, 1.0)
        self.assertEqual(Color(-0.001, -2183.0, -1.0).green, 0.0)
        self.assertEqual(Color(1.0001, 213.0, 1.2).blue, 1.0)
        self.assertEqual(Color(-0.001, -2183.0, -1.0).blue, 0.0)

    def test_color_fail_on_mixed_type(self):
        self.assertRaises(TypeError, Color.from_rgb, 1, 1, 1.0)
        self.assertRaises(TypeError, Color.from_rgb, 1.0, 1, 1.0)

    def test_color_red(self):
        # Check Usage comment in Color
        red = Color.by_name("red")
        self.assertEqual(red, Color(255, 0, 0))
        self.assertEqual(red, Color(1., 0., 0.))
        self.assertEqual(red, Color.from_rgb(1., 0., 0.))
        self.assertEqual(red, Color.from_rgb(255, 0, 0))
        self.assertEqual(red, Color.from_hsl(0., 1., 0.5))
        self.assertEqual(red, Color.from_string("red"))
        self.assertEqual(red, Color.from_string("RED"))
        self.assertEqual(red, Color.from_string("#F00"))
        self.assertEqual(red, Color.from_string("#FF0000"))
        self.assertEqual(red, Color.from_string("rgb(255, 0, 0)"))
        self.assertEqual(red, Color.from_string("rgb(100%, 0%, 0%)"))
        self.assertEqual(red, Color.from_string("hsl(0, 100%, 50%)"))

    def test_color_from_string(self):
        Color(128, 0, 128)  # purple
        red = Color(255, 0, 0)
        skyblue = Color(135, 206, 235)

        red_strings = ("red",
                       "ReD",
                       "RED",
                       "   Red \t",
                       "#F00",
                       "#FF0000",
                       "rgb(255, 0, 0)",
                       "rgb(100%, 0%, 0%)",
                       "hsl(0, 100%, 50%)")
        for s in red_strings:
            self.assertEqual(red, Color.from_string(s))

        skyblue_strings = ("skyblue",
                           "SKYBLUE",
                           "  \t\n SkyBlue  \t",
                           "#87ceeb",
                           "rgb(135,206,235)"
                           )
        for s in skyblue_strings:
            self.assertEqual(skyblue, Color.from_string(s))

        self.assertRaises(ValueError, Color.from_string, '#not_a_color')
        self.assertRaises(ValueError, Color.from_string, 'rgb(not_a_color)')
        self.assertRaises(ValueError, Color.from_string, 'hsl(not_a_color)')
        self.assertRaises(ValueError, Color.from_string, 'not_a_color')

    def test_color_equality(self):
        c1 = Color(123, 99, 12)
        c2 = Color(123, 99, 12)
        self.assertEqual(c1, c2)
        self.assertNotEqual(c1, "not_a_color")


class test_gamma(unittest.TestCase):
    def test_create(self):
        a = 1.213
        b = 3.210
        g = Gamma(a, b)
        self.assertEqual(g.alpha, a)
        self.assertEqual(g.beta, b)

    def test_mean_variance(self):
        g = Gamma.from_mean_variance(2.0, 3.0)
        self.assertEqual(g.mean(), 2.0)
        self.assertEqual(g.variance(), 3.0)

        g = Gamma.from_mean_variance(2.0123, 3.01283)
        self.assertEqual(g.mean(), 2.0123)
        self.assertEqual(g.variance(), 3.01283)

    def test_from_shape_scale(self):
        g = Gamma.from_shape_scale(1.0, 8.0)
        self.assertEqual(g.alpha, 1.0)
        self.assertEqual(g.beta, 1.0 / 8.0)

    def test_invalid_args(self):
        self.assertRaises(ValueError, Gamma, 1.0, -1.0)
        self.assertRaises(ValueError, Gamma, 0.0, 1.0)
        self.assertRaises(ValueError, Gamma, 1.0, 0.0)

    def test_sample(self):
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

        error = 4. * sqrt(g.variance() / S)
        # print(mean, m, error)
        self.assertTrue(abs(mean - m) < error)

    def test_pdf(self):
        m = 3.0
        v = 2.0
        g = Gamma.from_mean_variance(m, v)
        upper = 30.

        norm = integrate(g.pdf, 0, upper)
        self.assertAlmostEqual(norm, 1.0)

        def fx(x): return x * g.pdf(x)

        mean = integrate(fx, 0, upper)
        self.assertAlmostEqual(mean, m)

        def fx2(x): return x * x * g.pdf(x)

        x2 = integrate(fx2, 0, upper)
        var = x2 - mean ** 2
        self.assertAlmostEqual(var, v)

    def test_cdf(self):
        m = 3.0
        v = 2.0
        g = Gamma.from_mean_variance(m, v)
        # Numerical integration
        S = 1000
        M = 10.
        total_p = 0.0
        epsilon = 1e-4
        last = 0.0
        for s in range(S):
            x = s * M / S
            p = g.pdf(x) * M / S
            total_p += (last - p) / 2.0
            last = p
            # print(x, total_p, g.cdf(x))

            self.assertTrue((total_p - g.cdf(x)) < epsilon)

    def test_inverse_cdf(self):
        g = Gamma.from_mean_variance(2.34, 4)
        self.assertAlmostEqual(3.9, g.inverse_cdf(g.cdf(3.9)))
        self.assertAlmostEqual(1.92, g.inverse_cdf(g.cdf(1.92)))

        g = Gamma.from_mean_variance(10.34, 2)
        self.assertAlmostEqual(3.9, g.inverse_cdf(g.cdf(3.9)))
        self.assertAlmostEqual(10.92, g.inverse_cdf(g.cdf(10.92)))

        g = Gamma.from_mean_variance(10.34, 2)
        self.assertAlmostEqual(0.975, g.cdf(g.inverse_cdf(0.975)))
        self.assertAlmostEqual(0.025, g.cdf(g.inverse_cdf(0.025)))

        g = Gamma.from_mean_variance(1.34, 4)
        self.assertAlmostEqual(0.975, g.cdf(g.inverse_cdf(0.975)))
        self.assertAlmostEqual(0.025, g.cdf(g.inverse_cdf(0.025)))


class test_Dirichlet(unittest.TestCase):
    def test_init(self):
        Dirichlet((1, 1, 1, 1,))

    def test_random(self):

        def do_test(alpha, samples=1000):
            ent = zeros((samples,), float64)
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
            error = 4. * sqrt(v / samples)
            self.assertTrue(abs(m - dm) < error)
            self.assertTrue(abs(v - dv) < error)  # dodgy error estimate

        do_test((1., 1.))
        do_test((2., 1.))
        do_test((3., 1.))
        do_test((4., 1.))
        do_test((5., 1.))
        do_test((6., 1.))

        do_test((1., 1.))
        do_test((20., 20.))
        do_test((1., 1., 1., 1., 1., 1., 1., 1., 1., 1.))
        do_test((.1, .1, .1, .1, .1, .1, .1, .1, .1, .1))
        do_test((.01, .01, .01, .01, .01, .01, .01, .01, .01, .01))
        do_test((2.0, 6.0, 1.0, 1.0))

    def test_mean(self):
        alpha = ones((10,), float64) * 23.
        d = Dirichlet(alpha)
        m = d.mean()
        self.assertAlmostEqual(m[2], 1. / 10)
        self.assertAlmostEqual(sum(m), 1.0)

    def test_covariance(self):
        alpha = ones((4,), float64)
        d = Dirichlet(alpha)
        cv = d.covariance()
        self.assertEqual(cv.shape, (4, 4))
        self.assertAlmostEqual(cv[0, 0], 1.0 * (1.0 - 1. / 4.0) / (4.0 * 5.0))
        self.assertAlmostEqual(cv[0, 1], - 1 / (4. * 4. * 5.))

    def test_mean_x(self):
        alpha = (1.0, 2.0, 3.0, 4.0)
        xx = (2.0, 2.0, 2.0, 2.0)
        m = Dirichlet(alpha).mean_x(xx)
        self.assertEqual(m, 2.0)

        xx2 = (2.0, 2.0, 2.0, 2.0, 2.0)
        self.assertRaises(ValueError, Dirichlet(alpha).mean_x, xx2)

        alpha = (1.0, 1.0, 1.0, 1.0)
        xx = (2.0, 3.0, 4.0, 3.0)
        m = Dirichlet(alpha).mean_x(xx)
        self.assertEqual(m, 3.0)

    def test_variance_x(self):
        alpha = (1.0, 1.0, 1.0, 1.0)
        xx = (2.0, 2.0, 2.0, 2.0)
        v = Dirichlet(alpha).variance_x(xx)
        self.assertAlmostEqual(v, 0.0)

        alpha = (1.0, 2.0, 3.0, 4.0)
        xx = (2.0, 0.0, 1.0, 10.0)
        v = Dirichlet(alpha).variance_x(xx)
        # print(v)
        # TODO: Don't actually know if this is correct

        xx2 = (2.0, 2.0, 2.0, 2.0, 2.0)
        self.assertRaises(ValueError, Dirichlet(alpha).variance_x, xx2)

    def test_relative_entropy(self):
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
        sent = zeros((samples,), float64)

        for s in range(samples):
            post = d.sample()
            e = -entropy(post)
            for k in range(4):
                e += - post[k] * log(pvec[k])
            sent[s] = e
        sent.sort()
        self.assertTrue(abs(sent.mean() - rent) < 4. * sqrt(vrent))
        self.assertAlmostEqual(sent.std(), sqrt(vrent), 1)
        self.assertTrue(abs(low - sent[int(samples * 0.025)]) < 0.2)
        self.assertTrue(abs(high - sent[int(samples * 0.975)]) < 0.2)


class _from_URL_fileopen_Tests(unittest.TestCase):
    def test_URLscheme(self):
        """test for http, https, or ftp scheme"""
        from weblogo.logo import _from_URL_fileopen
        broken_url = "file://foo.txt"
        self.assertRaises(ValueError, _from_URL_fileopen, (broken_url))


def mean(a):
    return sum(a) / len(a)


def var(a):
    return (sum(a * a) / len(a)) - mean(a) ** 2


def integrate(f, a, b, n=1000):
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


if __name__ == '__main__':
    unittest.main()
