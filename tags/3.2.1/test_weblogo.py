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

import weblogolib
from weblogolib import  *
from weblogolib import  parse_prior, GhostscriptAPI
from weblogolib.color import *
from weblogolib.colorscheme import *
from StringIO import StringIO
import sys

from numpy import array, asarray, float64, ones, zeros, int32,all,any, shape
import numpy as na

from corebio import seq_io
from corebio.seq import *

from subprocess import *
from pkg_resources import resource_stream

from corebio.moremath import entropy
from math import log, sqrt

from weblogolib.logomath import Dirichlet, Gamma

def testdata_stream( name ): 
    return resource_stream(__name__, 'test_weblogo/data/'+name)    

class test_logoformat(unittest.TestCase) :

    def test_options(self) :
        options = LogoOptions()

           
class test_ghostscript(unittest.TestCase) :
    def test_version(self) :
        version = GhostscriptAPI().version
    

    
class test_parse_prior(unittest.TestCase) :
    def assertTrue(self, bool) :
        self.assertEquals( bool, True)
    
    def test_parse_prior_none(self) :
        self.assertEquals( None, 
            parse_prior(None, unambiguous_protein_alphabet ) )
        self.assertEquals( None, 
            parse_prior( 'none', unambiguous_protein_alphabet ) )        
        self.assertEquals( None, 
            parse_prior( 'noNe', None) )                

    def test_parse_prior_equiprobable(self) :            
        self.assertTrue( all(20.*equiprobable_distribution(20)  ==
            parse_prior( 'equiprobable',  unambiguous_protein_alphabet ) ) )
                        
        self.assertTrue( 
            all( 1.2* equiprobable_distribution(3) 
            == parse_prior( ' equiprobablE  ',  Alphabet('123'), 1.2 )  ) )

    def test_parse_prior_percentage(self) :  
        #print       parse_prior( '50%', unambiguous_dna_alphabet, 1. )          
        self.assertTrue( all( equiprobable_distribution(4) 
            == parse_prior( '50%', unambiguous_dna_alphabet, 1. )  ) )
 
        self.assertTrue( all( equiprobable_distribution(4) 
            == parse_prior( ' 50.0 % ', unambiguous_dna_alphabet, 1. )  ) )

        self.assertTrue( all( array( (0.3,0.2,0.2,0.3), float64) 
            == parse_prior( ' 40.0 % ', unambiguous_dna_alphabet, 1. )  ) )

    def test_parse_prior_float(self) :  
        self.assertTrue( all( equiprobable_distribution(4) 
            == parse_prior( '0.5', unambiguous_dna_alphabet, 1. )  ) )
 
        self.assertTrue( all( equiprobable_distribution(4) 
            == parse_prior( ' 0.500 ', unambiguous_dna_alphabet, 1. )  ) )

        self.assertTrue( all( array( (0.3,0.2,0.2,0.3), float64) 
            == parse_prior( ' 0.40 ', unambiguous_dna_alphabet, 1. )  ) )

    def test_auto(self) :
        self.assertTrue( all(4.*equiprobable_distribution(4)  ==
            parse_prior( 'auto',  unambiguous_dna_alphabet ) ) )
        self.assertTrue( all(4.*equiprobable_distribution(4)  ==
            parse_prior( 'automatic',  unambiguous_dna_alphabet ) ) )
       
    def test_weight(self) :
        self.assertTrue( all(4.*equiprobable_distribution(4)  ==
            parse_prior( 'automatic',  unambiguous_dna_alphabet ) ) )
        self.assertTrue( all(123.123*equiprobable_distribution(4)  ==
            parse_prior( 'auto',  unambiguous_dna_alphabet , 123.123) ) )
 
    def test_explicit(self) :
        s = "{'A':10, 'C':40, 'G':40, 'T':10}"
        p = array( (10, 40, 40,10), float64)*4./100.
        self.assertTrue( all(
            p == parse_prior( s,  unambiguous_dna_alphabet ) ) )
        
        
class test_logooptions(unittest.TestCase) :
    def test_create(self) :
        opt = LogoOptions()
        opt.small_fontsize =10
        options = repr(opt)
        
        opt = LogoOptions(title="sometitle")
        assert opt.title == "sometitle"
        
        
class test_seqlogo(unittest.TestCase) :
    # FIXME: The version of python used by Popen may not be the
    # same as that used to run this test.
    def _exec(self,  args, outputtext, returncode =0, stdin=None) :
        if not stdin : 
            stdin = testdata_stream("cap.fa")
        args = ["./weblogo"] + args
        p = Popen(args,stdin=stdin,stdout=PIPE, stderr=PIPE)    
        (out, err) = p.communicate()
        if returncode ==0 and p.returncode >0 :
            print err
        self.assertEquals(returncode, p.returncode)
        if returncode == 0 : self.assertEquals( len(err), 0)
    
        for item in outputtext :
            self.failUnless(item in out)  


                 
    def test_malformed_options(self) :
        self._exec( ["--notarealoption"], [], 2)
        self._exec( ["extrajunk"], [], 2)
        self._exec( ["-I"], [], 2)
    
    def test_help_option(self) :
        self._exec( ["-h"], ["options"])
        self._exec( ["--help"], ["options"])

    def test_version_option(self) :
        self._exec( ['--version'], weblogolib.__version__)

    
    def test_default_build(self) :
        self._exec( [], ["%%Title:        Sequence Logo:"] )

    
    # Format options
    def test_width(self) :
        self._exec( ['-W','1234'], ["/stack_width         1234"] )
        self._exec( ['--stack-width','1234'], ["/stack_width         1234"] )

    def test_height(self) :
        self._exec( ['-W','1000'], ["/stack_height        5000"] )
        self._exec( ['-W','1000', '--aspect-ratio', '2'], ["/stack_height        2000"] )
        
        
    def test_stacks_per_line(self) :
        self._exec( ['-n','7'], ["/stacks_per_line     7 def"] )    
        self._exec( ['--stacks-per-line','7'], ["/stacks_per_line     7 def"] )    
    
        
    def test_title(self) :            
        self._exec( ['-t', '3456'], ['/logo_title         (3456) def',
            '/show_title         True def'])
        self._exec( ['-t', ''], ['/logo_title         () def',
            '/show_title         False def'])
        self._exec( ['--title', '3456'], ['/logo_title         (3456) def',
            '/show_title         True def'])
        
    def test_annotate(self) :
        self._exec( ["--annotate","1,2,3,4"], [], 2)
        self._exec( ["--annotate","1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,,"],[])

        

        


class test_colorscheme(unittest.TestCase) :
    
    def test_colorgroup(self) :
        cr = ColorGroup( "ABC", "black", "Because")
        self.assertEquals( cr.description, "Because")       

    def test_colorscheme(self) :
        cs = ColorScheme([
                ColorGroup("G", "orange"),
                ColorGroup("TU", "red"),
                ColorGroup("C",  "blue"),
                ColorGroup("A",  "green") 
            ],
            title = "title",
            description = "description",
            ) 

        self.assertEquals( cs.color('A'), Color.by_name("green"))
        self.assertEquals( cs.color('X'), cs.default_color)



class test_color(unittest.TestCase) :

    def test_color_names(self) :
        names = Color.names()
        self.failUnlessEqual( len(names), 147)
        
        for n in names:
            c = Color.by_name(n)
            self.assertTrue( c != None )

       
    def test_color_components(self) :
        white = Color.by_name("white")
        self.failUnlessEqual( 1.0, white.red)
        self.failUnlessEqual( 1.0, white.green)
        self.failUnlessEqual( 1.0, white.blue)
    
    
        c = Color(0.3, 0.4, 0.2)
        self.failUnlessEqual( 0.3, c.red)
        self.failUnlessEqual( 0.4, c.green)
        self.failUnlessEqual( 0.2, c.blue)

        c = Color(0,128,0)
        self.failUnlessEqual( 0.0, c.red)
        self.failUnlessEqual( 128./255., c.green)
        self.failUnlessEqual( 0.0, c.blue)


    def test_color_from_rgb(self) :
        white = Color.by_name("white")
        
        self.failUnlessEqual(white, Color(1.,1.,1.) )        
        self.failUnlessEqual(white, Color(255,255,255) )        
        self.failUnlessEqual(white, Color.from_rgb(1.,1.,1.) )        
        self.failUnlessEqual(white, Color.from_rgb(255,255,255) )        

    
    def test_color_from_hsl(self) :
        red = Color.by_name("red")
        lime = Color.by_name("lime")
        saddlebrown = Color.by_name("saddlebrown")
        darkgreen = Color.by_name("darkgreen")
        blue = Color.by_name("blue")
        green = Color.by_name("green")
            
        self.failUnlessEqual(red, Color.from_hsl(0, 1.0,0.5) )          
        self.failUnlessEqual(lime, Color.from_hsl(120, 1.0, 0.5) ) 
        self.failUnlessEqual(blue, Color.from_hsl(240, 1.0, 0.5) )
        self.failUnlessEqual(Color.by_name("gray"), Color.from_hsl(0,0,0.5) )

        self.failUnlessEqual(saddlebrown, Color.from_hsl(25, 0.76, 0.31) ) 

        self.failUnlessEqual(darkgreen, Color.from_hsl(120, 1.0, 0.197) )           
                    
    
    def test_color_by_name(self):
        white = Color.by_name("white")
        self.failUnlessEqual(white, Color.by_name("white"))
        self.failUnlessEqual(white, Color.by_name("WHITE"))
        self.failUnlessEqual(white, Color.by_name(" wHiTe \t\n\t"))
        
        
        self.failUnlessEqual(Color(255,255,240), Color.by_name("ivory"))
        self.failUnlessEqual(Color(70,130,180), Color.by_name("steelblue"))        
        
        self.failUnlessEqual(Color(0,128,0), Color.by_name("green"))        
        

    def test_color_from_invalid_name(self):
        self.failUnlessRaises( ValueError, Color.by_name, "not_a_color")            


    def test_color_clipping(self):
        red = Color.by_name("red")
        self.failUnlessEqual(red, Color(255,0,0) ) 
        self.failUnlessEqual(red, Color(260,-10,0) ) 
        self.failUnlessEqual(red, Color(1.1,-0.,-1.) )        

        self.failUnlessEqual( Color(1.0001,   213.0,  1.2).red, 1.0 )
        self.failUnlessEqual( Color(-0.001, -2183.0, -1.0).red, 0.0 )
        self.failUnlessEqual( Color(1.0001,   213.0,  1.2).green, 1.0 )
        self.failUnlessEqual( Color(-0.001, -2183.0, -1.0).green, 0.0 )        
        self.failUnlessEqual( Color(1.0001,   213.0,  1.2).blue, 1.0 )
        self.failUnlessEqual( Color(-0.001, -2183.0, -1.0).blue, 0.0 )
    
    
    def test_color_fail_on_mixed_type(self):
        self.failUnlessRaises( TypeError, Color.from_rgb, 1,1,1.0 )            
        self.failUnlessRaises( TypeError, Color.from_rgb, 1.0,1,1.0 )            
    
    def test_color_red(self) :
        # Check Usage comment in Color
        red = Color.by_name("red")
        self.failUnlessEqual( red , Color(255,0,0) )
        self.failUnlessEqual( red,  Color(1., 0., 0.) )
    
        self.failUnlessEqual( red , Color.from_rgb(1.,0.,0.) )
        self.failUnlessEqual( red , Color.from_rgb(255,0,0) )
        self.failUnlessEqual( red , Color.from_hsl(0.,1., 0.5) )
    
        self.failUnlessEqual( red , Color.from_string("red") )
        self.failUnlessEqual( red , Color.from_string("RED") )
        self.failUnlessEqual( red , Color.from_string("#F00") )
        self.failUnlessEqual( red , Color.from_string("#FF0000") ) 
        self.failUnlessEqual( red , Color.from_string("rgb(255, 0, 0)") )
        self.failUnlessEqual( red , Color.from_string("rgb(100%, 0%, 0%)") )
        self.failUnlessEqual( red , Color.from_string("hsl(0, 100%, 50%)") )

    
    def test_color_from_string(self) :
        purple = Color(128,0,128)
        red    = Color(255,0,0)
        skyblue = Color(135,206,235)
        
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
            self.failUnlessEqual( red, Color.from_string(s) ) 
            
        skyblue_strings = ("skyblue", 
                        "SKYBLUE", 
                        "  \t\n SkyBlue  \t", 
                        "#87ceeb", 
                        "rgb(135,206,235)"
                        )
                      
        for s in skyblue_strings:        
            self.failUnlessEqual( skyblue, Color.from_string(s) ) 
    
    
    
    def test_color_equality(self):
        c1 = Color(123,99,12)
        c2 = Color(123,99,12)

        self.failUnlessEqual(c1,c2)        
        
        
#

class test_gamma(unittest.TestCase) :
    def test_create(self) :
        a = 1.213
        b = 3.210
        g = Gamma(a, b)
        self.assertEquals( g.alpha, a)
        self.assertEquals( g.beta, b)
        
    def test_mean_variance(self) :
        g = Gamma.from_mean_variance(2.0, 3.0)
        self.assertEquals( g.mean(), 2.0)
        self.assertEquals( g.variance(), 3.0)
                
        g = Gamma.from_mean_variance(2.0123, 3.01283)
        self.assertEquals( g.mean(), 2.0123)
        self.assertEquals( g.variance(), 3.01283)
        
    def test_from_shape_scale(self):
        g = Gamma.from_shape_scale(1.0, 8.0)
        self.assertEquals( g.alpha, 1.0)
        self.assertEquals( g.beta, 1.0/8.0)

    def test_invalid_args(self):
        self.failUnlessRaises( ValueError, Gamma, 1.0, -1.0   )
        self.failUnlessRaises( ValueError, Gamma, 0.0,  1.0   )
        self.failUnlessRaises( ValueError, Gamma, 1.0, 0.0   )

        
    def test_sample(self) :
        m = 10.0
        v = 2.0
        g = Gamma.from_mean_variance(m, v)
        #print g.alpha, g.beta
        S = 1000
        total = 0.0
        for s in range(S):
            total += g.sample()
        mean = total/S
        
        # The estimated mean will differ from true mean by a small amount

        error = 4. * sqrt( g.variance()/S)
        # print mean, m, error
        self.assertTrue( abs(mean - m) < error)

        
    def test_pdf(self) :
        m = 3.0
        v = 2.0
        g = Gamma.from_mean_variance(m, v)
        upper = 30.
                    
        norm = integrate( g.pdf, 0, upper)
        self.assertAlmostEqual(norm , 1.0)
        
        def fx(x) : return x * g.pdf(x)
        mean = integrate( fx, 0, upper) 
        self.assertAlmostEqual( mean, m)

        def fx2(x) : return x*x * g.pdf(x)
        x2 = integrate( fx2, 0, upper)
        var = x2 - mean**2
        self.assertAlmostEqual( var, v)
         
    def test_cdf(self) :
        m = 3.0
        v = 2.0
        g = Gamma.from_mean_variance(m, v)    
        # Numerical integration
        S = 1000
        M = 10.
        total_p = 0.0
        epsilon = 1e-4
        last = 0.0
        for s in range(S) :
            x = s*M/S
            p = g.pdf(x) *M/S
            total_p += (last-p) /2.0
            last = p
            #print x, total_p, g.cdf(x)
            
            self.assertTrue( (total_p - g.cdf(x)) <epsilon )

    def test_inverse_cdf(self) :
        g = Gamma.from_mean_variance( 2.34, 4)   
        self.assertAlmostEquals( 3.9, g.inverse_cdf(g.cdf(3.9) ) )
        self.assertAlmostEquals( 1.92, g.inverse_cdf(g.cdf(1.92) ) )

        g = Gamma.from_mean_variance( 10.34, 2)   
        self.assertAlmostEquals( 3.9, g.inverse_cdf(g.cdf(3.9) ) )
        self.assertAlmostEquals( 10.92, g.inverse_cdf(g.cdf(10.92) ) )

        g = Gamma.from_mean_variance( 10.34, 2)  
        self.assertAlmostEquals( 0.975, g.cdf(g.inverse_cdf(0.975) ) )
        self.assertAlmostEquals( 0.025, g.cdf(g.inverse_cdf(0.025) ) )

        g = Gamma.from_mean_variance( 1.34, 4)  
        self.assertAlmostEquals( 0.975, g.cdf(g.inverse_cdf(0.975) ) )
        self.assertAlmostEquals( 0.025, g.cdf(g.inverse_cdf(0.025) ) )
       
       
       
   
class test_Dirichlet(unittest.TestCase) :
    
    def test_init(self) :
        d = Dirichlet( ( 1,1,1,1,) )
        
    
    def test_random(self) :


        def do_test( alpha, samples = 1000) :
            ent = zeros( (samples,), float64)
            #alpha = ones( ( K,), Float64 ) * A/K
            
            #pt = zeros( (len(alpha) ,), Float64)
            d = Dirichlet(alpha)
            for s in range(samples) :
                p = d.sample()
                #print p
                #pt +=p
                ent[s] = entropy(p)
            
            #print pt/samples
            
            m = mean(ent)
            v = var(ent)
            
            dm = d.mean_entropy()
            dv = d.variance_entropy()

           #print alpha, ':', m, v, dm, dv 
            error = 4. * sqrt(v/samples)
            self.assertTrue( abs(m-dm) < error)
            self.assertTrue( abs(v-dv) < error) # dodgy error estimate 


        do_test( (1., 1.) )
        do_test( (2., 1.) )                                        
        do_test( (3., 1.) )                                        
        do_test( (4., 1.) )                                        
        do_test( (5., 1.) )                                        
        do_test( (6., 1.) )                                        

        do_test( (1., 1.) )
        do_test( (20., 20.) )
        do_test( (1., 1., 1., 1., 1., 1., 1., 1., 1., 1.) )
        do_test( (.1, .1, .1, .1, .1, .1, .1, .1, .1, .1) )
        do_test( (.01, .01, .01, .01, .01, .01, .01, .01, .01, .01) )
        do_test( (2.0, 6.0, 1.0, 1.0) )
    
    
    def test_mean(self) :
        alpha = ones( ( 10,), float64 ) * 23.
        d = Dirichlet(alpha)
        m = d.mean()
        self.assertAlmostEqual( m[2], 1./10)
        self.assertAlmostEqual( sum(m), 1.0)
    
    def test_covariance(self) :
        alpha = ones( ( 4,), float64 ) 
        d = Dirichlet(alpha)
        cv = d.covariance()
        self.assertEqual( cv.shape, (4,4)  )
        self.assertAlmostEqual( cv[0,0], 1.0 * (1.0 - 1./4.0)/ (4.0 * 5.0)  )
        self.assertAlmostEqual( cv[0,1],  - 1 / ( 4. * 4. * 5.) )

    def test_mean_x(self) :
        alpha = (1.0, 2.0, 3.0, 4.0)
        xx = (2.0, 2.0, 2.0, 2.0)
        m = Dirichlet(alpha).mean_x(xx)
        self.assertEquals( m, 2.0)
        
        alpha = (1.0, 1.0, 1.0, 1.0)
        xx = (2.0, 3.0, 4.0, 3.0)
        m = Dirichlet(alpha).mean_x(xx)
        self.assertEquals( m, 3.0)        

    def test_variance_x(self) :
        alpha = (1.0, 1.0, 1.0, 1.0)
        xx = (2.0, 2.0, 2.0, 2.0)
        v = Dirichlet(alpha).variance_x(xx)
        self.assertAlmostEquals( v, 0.0)        

        alpha = (1.0, 2.0, 3.0, 4.0)
        xx = (2.0, 0.0, 1.0, 10.0)
        v = Dirichlet(alpha).variance_x(xx)
        #print v
        # TODO: Don't actually know if this is correct
        
    def test_relative_entropy(self):
        alpha = (2.0, 10.0, 1.0, 1.0)
        d = Dirichlet(alpha)
        pvec = (0.1, 0.2, 0.3, 0.4)
        
        rent = d.mean_relative_entropy(pvec)
        vrent = d.variance_relative_entropy(pvec)
        low, high = d.interval_relative_entropy(pvec, 0.95)
        
        #print
        #print '> ', rent, vrent, low, high
     
        # This test can fail randomly, but the precision from a few
        # thousand samples is low. Increasing samples, 1000->2000
        samples = 2000
        sent = zeros( (samples,), float64) 
        
        for s in range(samples) :
            post = d.sample()
            e = -entropy(post)
            for k in range(4) :
                e += - post[k] * log(pvec[k])
            sent[s] = e
        sent.sort()
        self.assertTrue( abs(sent.mean() - rent) < 4.*sqrt(vrent) )
        self.assertAlmostEqual( sent.std(), sqrt(vrent), 1 )
        self.assertTrue( abs(low-sent[ int( samples *0.025)])<0.2 )
        self.assertTrue( abs(high-sent[ int( samples *0.975)])<0.2 )
        
        #print '>>', mean(sent), var(sent), sent[ int( samples *0.025)] ,sent[ int( samples *0.975)] 
       
       
   
def mean( a) :
    return sum(a)/ len(a)
    
def var(a) :
    return (sum(a*a) /len(a) ) - mean(a)**2
    
    
    
#
def integrate(f, a, b, n=1000): 
    """Numerically integrate the function 'f' from 'a' to 'b' using
    a discretization with 'n' points.

    Args:    
    - f -- A function that eats a float and returns a float.
    - a -- Lower integration bound (float)
    - b -- Upper integration bound (float)
    - n -- number of sample points (int)

    Status :
        Alpha (very primitive.)
    """
    h = (b-a)/(n-1.0); 
    total = 0.0; 
    for i in range(n): 
        total += f(a+(i)*h); 
    result  = h*(total - 0.5*f(a) -0.5*f(b)); 
    return result  

        
if __name__ == '__main__':
    unittest.main()
