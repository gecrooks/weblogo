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

Constants :

- euler_gamma  = 0.577215...
- catalan      = 0.915965...
- golden_ratio = 1.618033...
- bits_per_nat = log2(e) = 1/log(2) 
- sqrt_2pi     = 2.50662...
    
Special Functions :



- gamma()                       -- Gamma function.
- lngamma()                     -- Logarithm of the gamma function
- factorial()                   -- Factorial function.
- digamma()                     -- Digamma function (logarithmic derivative of gamma).
- trigamma()                    -- Trigamma function (derivative of digamma).

- cgamma()                       -- Complex math version of counterparts above.
- clngamma()                                   
- cdigamma()                     
- ctrigamma()                    

- entropy()                     -- The entropy of a probability vector
- incomplete_gamma()            -- The 'upper' incomplete gamma function.
- normalized_incomplete_gamma() -- 
- log2()                          -- Base 2 logarithms.
- argmin()
- argmax()


"""



__all__ = ('euler_gamma', 'catalan', 'golden_ratio', 'bits_per_nat', 'sqrt_2pi',
            'gamma', 'lngamma', 'factorial', 'digamma', 'trigamma',
            'cgamma', 'clngamma', 'cdigamma', 'ctrigamma',            
            'entropy', 'log2',            
            'incomplete_gamma', 'normalized_incomplete_gamma', 
            'argmax', 'argmin'
            )            

from math import *
import cmath as cm

import random
from itertools import izip, count

# Some mathematical constants
euler_gamma  = 0.57721566490153286060651
catalan      = 0.91596559417721901505460
golden_ratio = 1.6180339887498948482046
bits_per_nat = 1.44269504088896340735992468100 # = log_2(e) = 1/log(2) 
sqrt_2pi     = 2.5066282746310005024157652848110





# The Lanczos approximation for the gamma function is
#                          
#                -(z + g + 1/2)               (z + 1/2)                   
# Gamma(z+1) =  e               * (z + g + 1/2)        * Sqrt(2Pi) * C 
#                                                                        
#
#                      c[1]    c[2]    c[3]
#          C = [c[0] + ----- + ----- + ----- + ...   ]
#                      z + 1   z + 2   z + 3 
#
#
#  To calculate digamma and trigamma functions we take an analytic derivative
#  of the Lanczos approximation.
#
#  Gamma(z)  = Gamma(z+1)/z
#  Digamma(z) = D ln Gamma(z)
#  Trigamma(z) = D Digamma(z)         

# These Lanczos constants are from
# "A note on the computation of the convergent
# Lanczos complex Gamma approximation." Paul Godfrey (2001)
# http://my.fit.edu/~gabdo/gamma.txt


__lanczos_gamma = 607./128.
__lanczos_coefficients = ( 
       0.99999999999999709182,
      57.156235665862923517,
     -59.597960355475491248,
      14.136097974741747174,
      -0.49191381609762019978,
        .33994649984811888699e-4,
        .46523628927048575665e-4,
       -.98374475304879564677e-4,
        .15808870322491248884e-3,
       -.21026444172410488319e-3,
        .21743961811521264320e-3,
       -.16431810653676389022e-3,
        .84418223983852743293e-4,
       -.26190838401581408670e-4,
        .36899182659531622704e-5)

__factorial =(
    1.,
    1.,
    2.,
    6.,
    24.,
    120.,
    720.,
    5040.,
    40320.,
    362880.,
    3628800.,
    39916800.,
    479001600.,
    6227020800.,
    87178291200.,
    1307674368000.,
    20922789888000.,
    355687428096000.,
    6402373705728000.,
    121645100408832000.,
    2432902008176640000.,
    51090942171709440000.,
    1124000727777607680000.,
    25852016738884976640000.,
    620448401733239439360000.,
    15511210043330985984000000.,
    403291461126605635584000000.,
    10888869450418352160768000000.,
    304888344611713860501504000000.,
    8841761993739701954543616000000.,
    265252859812191058636308480000000.,
    8222838654177922817725562880000000.,
    263130836933693530167218012160000000. )
    
def gamma(x) :
    """The gamma function. Returns exact results for small integers. Will
    overflow for modest sized arguments. Use lngamma(z) instead. 
    
    See: Eric W. Weisstein. "Gamma Function." From MathWorld, A Wolfram Web Resource.
         http://mathworld.wolfram.com/GammaFunction.html

    """
    return cgamma(x).real


def cgamma(z) :
    """Gamma function with complex arguments and results."""
    
    z = complex(z)
    n = floor(z.real)   
    if  n == z  :
        if n <= 0 :
            return complex(1.0/0.0) # Infinity
        elif n <= len(__factorial) :
            return complex(__factorial[int(n)-1])
        
    zz = z
    if z.real < 0.5 :   # Reflection formula
        zz = 1-z
    
        
    g = __lanczos_gamma
    c = __lanczos_coefficients 
    
    zz = zz - 1.    
    zh = zz + 0.5
    zgh = zh + g
    zp = zgh** (zh*0.5) # trick for avoiding FP overflow above z=141

    ss = 0.0
    for k in range(len(c)-1,0,-1):
        ss += c[k]/(zz+k)
    
    f = (sqrt_2pi*(c[0]+ss)) * (( zp*cm.exp(-zgh)) *zp)

    if z.real<0.5 :
        f  = pi /( cm.sin(pi*z) *f)

    return complex(f)

def lngamma(x) :
    """The logarithm of the gamma function. 
    """
    return clngamma(x).real

   
def clngamma(z) :
    """The logarithm of the gamma function for a complex argument
    """
    z = complex(z)
    
    # common case optimization
    n = floor(z.real)   
    if  n == z  :
        if n <= 0 :
            return complex(1.0/0.0) # Infinity
        elif n <= len(__factorial) :
            return complex(__factorial[int(n)-1])
        
    zz = z
    if z.real < 0.5 :
        zz = 1-z
    
        
    g = __lanczos_gamma
    c = __lanczos_coefficients 
    
    zz = zz - 1.    
    zh = zz + 0.5
    zgh = zh + g
    zp = zgh** (zh*0.5) # trick for avoiding FP overflow above z=141

    ss = 0.0
    for k in range(len(c)-1,0,-1):
        ss += c[k]/(zz+k)
    
    f = (sqrt_2pi*(c[0]+ss)) * (( zp*cm.exp(-zgh)) *zp)

    if z.real<0.5 :
        f  = pi /( cm.sin(pi*z) *f)

    return cm.log(f)
          
                    
def factorial(z) :
    """ The factorial function. 
    factorial(z) == gamma(z+1)
    """
    return gamma(z+1)
              
def digamma(x) :   
    """The digamma function, the logarithmic derivative of the gamma function.
             digamma(z) = d/dz ln( gamma(z) )


     See: Eric W. Weisstein. "Digamma Function." From MathWorld--
     A Wolfram Web Resource. http://mathworld.wolfram.com/DigammaFunction.html
    """
    return cdigamma(x).real
          
          
def cdigamma(z) :
    """Digamma function with complex arguments"""
    z = complex(z)

    g = __lanczos_gamma
    c = __lanczos_coefficients 
        
    zz = z
    if z.real < 0.5 :
        zz = 1 -z
        
    n=0.
    d=0.
    for k in range(len(c)-1,0,-1):
        dz =1./(zz+(k+1)-2);
        dd =c[k] * dz
        d = d + dd 
        n = n - dd * dz

    d = d + c[0]
    gg = zz + g - 0.5
    f = cm.log(gg) + (n/d - g/gg)

    if z.real<0.5 :
        f -= pi / cm.tan( pi * z)
        
    return f


def trigamma(x) :
    """The trigamma function, the derivative of the digamma function.
            trigamma(z) = d/dz digamma(z) = d/dz d/dz ln( gamma(z) )
    
    See: Eric W. Weisstein. "Trigamma Function." From MathWorld--
    A Wolfram Web Resource. http://mathworld.wolfram.com/TrigammaFunction.html
    """
    return ctrigamma(x).real    

def ctrigamma(z):
    "The trigamma function with complex arguments and return value."
    z = complex(z)
    g = __lanczos_gamma
    c = __lanczos_coefficients 
        
    t1=0.
    t2=0.
    t3=0.
    for k in range(len(c)-1,0,-1):
        dz =1./(z+k);       
        dd1 = c[k]* dz
        t1 += dd1
        dd2 = dd1 * dz
        t2 += dd2
        t3 += dd2 * dz

    t1 += c[0]
    c =  - (t2*t2)/(t1*t1)  +2*t3/t1

    result = 1./(z*z)
    gg = z + g + 0.5
    result += - (z+0.5)/ (gg*gg)
    result += 2./gg

    result += c

    return result

def incomplete_gamma(a,x) :
    """The 'upper' incomplete gamma function:

                            oo
                             -
                            |    -t  a-1
    incomplete_gamma(a,x) = |   e   t   dt.
                            |
                           -                         
                            x

    In Mathematica, Gamma[a,x].

    Note that, very confusingly, the phrase 'incomplete gamma fucntion'
    can also refer to the same integral between 0 and x, (the 'lower'
    incomplete gamma function) or to the normalized versions,
    normalized_incomplete_gamma() )


    See: Eric W. Weisstein. "Gamma Function." From MathWorld, A Wolfram Web Resource.
         http://mathworld.wolfram.com/IncompleteGammaFunction.html

    Bugs :
        This implementation is not very accurate for some arguments. 
    """   
    return  normalized_incomplete_gamma(a,x) * gamma(a)

    
def normalized_incomplete_gamma(a,x) :
    """The upper, incomplete gamma function normalized so that the limiting
    values are zero and one.
    
     Q(a,x) = incomplete_gamma(a,x) / gamma(a) 

    See: 
        incomplete_gamma()
    Bugs :
        This implementation is not very accurate for some arguments. 
    """
    maxiter = 100
    epsilon = 1.48e-8
    small = 1e-30
    
        
    if a<=0 or x<0 : 
        raise ValueError("Invalid arguments")
    if x == 0.0 : return 1.0
    
    if x<= a+1 :
        # Use the series representation
        term = 1./a
        total = term
        for n in range(1,maxiter) :
            term *= x/(a+n)
            total += term
            if abs(term/total) < epsilon : 
                return 1. - total * exp(-x+a*log(x) - lngamma(a) )
        raise RuntimeError(
            "Failed to converge after %d iterations." % (maxiter) )
    else :
        # Use the continued fraction representation
        total = 1.0
        b = x + 1. -a
        c = 1./small
        d = 1./b
        h = d
        for i in range(1, maxiter) :
            an = -i * (i-a)
            b = b+2.
            d = an * d + b
            if abs(d) < small : d = small
            c = b + an /c
            if abs(c) < small : c= small
            d = 1./d
            term = d * c
            h = h * term
            if abs( term-1.) < epsilon :
                return h * exp(-x+a*log(x) - lngamma(a) )
        raise RuntimeError(
            "Failed to converge after %d iterations." % (maxiter) )


 
def log2( x) :
    """ Return the base 2 logarithm of x """
    return log(x,2)


def entropy( pvec, base= exp(1) ) :
    """ The entropy S = -Sum_i p_i ln p_i
        pvec is a frequency vector, not necessarily normalized. 
    """
    # TODO: Optimize
    if len(pvec) ==0 : 
        raise ValueError("Zero length vector")
        
        
    total = 0.0
    ent = 0.0
    for p in pvec:
        if p>0 : # 0 log(0) =0 
            total += p
            ent += - log(float(p)) *p
        elif p<0:
            raise ValueError("Negative probability")
    
    
    ent = (ent/total) + log(total)
    ent /= log(base)

    return ent
   
    


    
def argmax( alist) :
    """Return the index of the last occurrence of the maximum value in the list."""
    return max(izip(alist, count() ))[1]

def argmin( alist) :
    """Return the index of the first occurrence of the minimum value in the list."""
    return min(izip(alist, count() ))[1]

     
    
           


