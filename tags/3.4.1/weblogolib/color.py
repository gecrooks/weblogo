

#  Copyright (c) 2005 Gavin E. Crooks
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

""" Color specifications using CSS2 (Cascading Style Sheet) syntax."""
from __future__ import division


class Color(object):
    """ Color specifications using CSS2 (Cascading Style Sheet) syntax.
    
    http://www.w3.org/TR/REC-CSS2/syndata.html#color-units

    Usage:
    
    red = Color(255,0,0)
    red = Color(1., 0., 0.)
    red = Color.by_name("red")
    red = Color.from_rgb(1.,0.,0.)
    red = Color.from_rgb(255,0,0)
    red = Color.from_hsl(0.,1., 0.5)
    
    red = Color.from_string("red") 
    red = Color.from_string("RED") 
    red = Color.from_string("#F00") 
    red = Color.from_string("#FF0000") 
    red = Color.from_string("rgb(255, 0, 0)")
    red = Color.from_string("rgb(100%, 0%, 0%)")
    red = Color.from_string("hsl(0, 100%, 50%)")

    """
    def __init__(self, red, green, blue) :
        if not (type(red) == type(green) == type(blue)):
            raise TypeError("Mixed floats and integers?")
        # Convert integer RBG values in [0, 255] to floats in [0, 1]
        if isinstance(red, int):
            red /= 255.
        if isinstance(green, int):
            green /= 255.
        if isinstance(blue, int):
            blue /= 255.
        # Clip RBG values to [0, 1]
        self.red = max(0., min(red, 1.0))
        self.green = max(0., min(green, 1.0))
        self.blue = max(0., min(blue, 1.0))

    @staticmethod
    def names():
        "Return a list of standard color names."
        return _std_colors.keys()

    @classmethod
    def from_rgb(cls, r, g, b):
        return cls(r,g,b)

    @classmethod
    def from_hsl(cls, hue_angle, saturation, lightness ):
        def hue_to_rgb( v1, v2, vH) :
            if vH < 0.0 : vH += 1.0
            if vH > 1.0 : vH -= 1.0
            if vH*6.0 < 1.0 : return (v1 + (v2 - v1) * 6.0 * vH)
            if vH*2.0 < 1.0 : return v2
            if vH*3.0 < 2.0 : return (v1 + (v2 - v1) * ((2.0/3.0) - vH) * 6.0)
            return v1      
    
        hue =  (((hue_angle % 360.) + 360.) % 360.)/360.
        
        if not (saturation >= 0.0 and saturation <=1.0) :
            raise ValueError("Out-of-range saturation %f"% saturation)
        if not (lightness >= 0.0 and lightness <=1.0) :
            raise ValueError("Out-of-range lightness %f"% lightness)
                
        if saturation == 0 :
            # greyscale
            return cls.from_rgb( lightness, lightness, lightness)
    
        if lightness < 0.5 :
            v2 = lightness * (1.0+ saturation)
        else :
            v2 = (lightness + saturation) - (saturation* lightness)       

        v1 = 2.0 * lightness - v2
        r = hue_to_rgb( v1, v2, hue + (1./3.) )
        g = hue_to_rgb( v1, v2, hue )
        b = hue_to_rgb( v1, v2, hue - (1./3.) )

        return cls(r,g,b)


    @staticmethod
    def by_name(string):
        s = string.strip().lower().replace(' ', '')
        try:
            return _std_colors[s]
        except KeyError:
            raise ValueError("Unknown color name: %s"%  s) 

    @classmethod
    def from_string(cls, string):
        def to_frac(string) :
            # string can be "255" or "100%"
            if string[-1]=='%': 
                return float(string[0:-1])/100.
            else:
                return float(string)/255.
                
        s = string.strip().lower().replace(' ', '').replace('_', '')
        
        if s in _std_colors : # "red"
            return _std_colors[s]

        if s[0] == "#" :    # "#fef"
            if len(s) == 4 :
                r = int(s[1]+s[1],16)
                g = int(s[2]+s[2],16)
                b = int(s[3]+s[3],16)
                return cls(r,g,b)
            elif len(s) ==7 :   # "#ff00aa"
                r = int(s[1:3],16)
                g = int(s[3:5],16)
                b = int(s[5:7],16)
                return cls(r,g,b)
            else :
                raise ValueError("Cannot parse string: %s" % s)

        if s[0:4] == 'rgb(' and s[-1] == ')' :
            rgb = s[4:-1].split(",")
            if len(rgb) != 3 :
                raise ValueError("Cannot parse string a: %s" % s)
            return cls( to_frac(rgb[0]), to_frac(rgb[1]), to_frac(rgb[2]))

        if s[0:4] == 'hsl(' and s[-1] == ')' :
            hsl = s[4:-1].split(",")
            if len(hsl) != 3 :
                raise ValueError("Cannot parse string a: %s" % s)
            return cls.from_hsl( int(hsl[0]), to_frac(hsl[1]), to_frac(hsl[2]))
        
        raise ValueError("Cannot parse string: %s" % s)

    def __eq__(self, other) :
        if not isinstance(other, self.__class__):
            return False
        req = int(0.5+255.*self.red) == int(0.5+255.*other.red)
        beq = int(0.5+255.*self.blue) == int(0.5+255.*other.blue)
        geq = int(0.5+255.*self.green) == int(0.5+255.*other.green)
        return req and beq and geq

    def __repr__(self):
        return "Color(%f,%f,%f)" % (self.red, self.green, self.blue)


_std_colors = dict(
    aliceblue = Color(240,248,255),     #f0f8ff
    antiquewhite = Color(250,235,215),     #faebd7
    aqua = Color(0,255,255),     #00ffff
    aquamarine = Color(127,255,212),     #7fffd4
    azure = Color(240,255,255),     #f0ffff
    beige = Color(245,245,220),     #f5f5dc
    bisque = Color(255,228,196),     #ffe4c4
    black = Color(0,0,0),     #000000
    blanchedalmond = Color(255,235,205),     #ffebcd
    blue = Color(0,0,255),     #0000ff
    blueviolet = Color(138,43,226),     #8a2be2
    brown = Color(165,42,42),     #a52a2a
    burlywood = Color(222,184,135),     #deb887
    cadetblue = Color(95,158,160),     #5f9ea0
    chartreuse = Color(127,255,0),     #7fff00
    chocolate = Color(210,105,30),     #d2691e
    coral = Color(255,127,80),     #ff7f50
    cornflowerblue = Color(100,149,237),     #6495ed
    cornsilk = Color(255,248,220),     #fff8dc
    crimson = Color(220,20,60),     #dc143c
    cyan = Color(0,255,255),     #00ffff
    darkblue = Color(0,0,139),     #00008b
    darkcyan = Color(0,139,139),     #008b8b
    darkgoldenrod = Color(184,134,11),     #b8860b
    darkgray = Color(169,169,169),     #a9a9a9
    darkgreen = Color(0,100,0),     #006400
    darkgrey = Color(169,169,169),     #a9a9a9
    darkkhaki = Color(189,183,107),     #bdb76b
    darkmagenta = Color(139,0,139),     #8b008b
    darkolivegreen = Color(85,107,47),     #556b2f
    darkorange = Color(255,140,0),     #ff8c00
    darkorchid = Color(153,50,204),     #9932cc
    darkred = Color(139,0,0),     #8b0000
    darksalmon = Color(233,150,122),     #e9967a
    darkseagreen = Color(143,188,143),     #8fbc8f
    darkslateblue = Color(72,61,139),     #483d8b
    darkslategray = Color(47,79,79),     #2f4f4f
    darkslategrey = Color(47,79,79),     #2f4f4f
    darkturquoise = Color(0,206,209),     #00ced1
    darkviolet = Color(148,0,211),     #9400d3
    deeppink = Color(255,20,147),     #ff1493
    deepskyblue = Color(0,191,255),     #00bfff
    dimgray = Color(105,105,105),     #696969
    dimgrey = Color(105,105,105),     #696969    
    dodgerblue = Color(30,144,255),     #1e90ff
    firebrick = Color(178,34,34),     #b22222
    floralwhite = Color(255,250,240),     #fffaf0
    forestgreen = Color(34,139,34),     #228b22
    fuchsia = Color(255,0,255),     #ff00ff
    gainsboro = Color(220,220,220),     #dcdcdc
    ghostwhite = Color(248,248,255),     #f8f8ff
    gold = Color(255,215,0),     #ffd700
    goldenrod = Color(218,165,32),     #daa520
    gray = Color(128,128,128),     #808080
    green = Color(0,128,0),     #008000
    greenyellow = Color(173,255,47),     #adff2f
    grey = Color(128,128,128),     #808080
    honeydew = Color(240,255,240),     #f0fff0
    hotpink = Color(255,105,180),     #ff69b4
    indianred = Color(205,92,92),     #cd5c5c
    indigo = Color(75,0,130),     #4b0082
    ivory = Color(255,255,240),     #fffff0
    khaki = Color(240,230,140),     #f0e68c
    lavender = Color(230,230,250),     #e6e6fa
    lavenderblush = Color(255,240,245),     #fff0f5
    lawngreen = Color(124,252,0),     #7cfc00
    lemonchiffon = Color(255,250,205),     #fffacd
    lightblue = Color(173,216,230),     #add8e6
    lightcoral = Color(240,128,128),     #f08080
    lightcyan = Color(224,255,255),     #e0ffff
    lightgoldenrodyellow = Color(250,250,210),     #fafad2
    lightgray = Color(211,211,211),     #d3d3d3
    lightgreen = Color(144,238,144),     #90ee90
    lightgrey = Color(211,211,211),     #d3d3d3
    lightpink = Color(255,182,193),     #ffb6c1
    lightsalmon = Color(255,160,122),     #ffa07a
    lightseagreen = Color(32,178,170),     #20b2aa
    lightskyblue = Color(135,206,250),     #87cefa
    lightslategray = Color(119,136,153),     #778899
    lightslategrey = Color(119,136,153),     #778899
    lightsteelblue = Color(176,196,222),     #b0c4de
    lightyellow = Color(255,255,224),     #ffffe0
    lime = Color(0,255,0),     #00ff00
    limegreen = Color(50,205,50),     #32cd32
    linen = Color(250,240,230),     #faf0e6
    magenta = Color(255,0,255),     #ff00ff
    maroon = Color(128,0,0),     #800000
    mediumaquamarine = Color(102,205,170),     #66cdaa
    mediumblue = Color(0,0,205),     #0000cd
    mediumorchid = Color(186,85,211),     #ba55d3
    mediumpurple = Color(147,112,219),     #9370db
    mediumseagreen = Color(60,179,113),     #3cb371
    mediumslateblue = Color(123,104,238),     #7b68ee
    mediumspringgreen = Color(0,250,154),     #00fa9a
    mediumturquoise = Color(72,209,204),     #48d1cc
    mediumvioletred = Color(199,21,133),     #c71585
    midnightblue = Color(25,25,112),     #191970
    mintcream = Color(245,255,250),     #f5fffa
    mistyrose = Color(255,228,225),     #ffe4e1
    moccasin = Color(255,228,181),     #ffe4b5
    navajowhite = Color(255,222,173),     #ffdead
    navy = Color(0,0,128),     #000080
    oldlace = Color(253,245,230),     #fdf5e6
    olive = Color(128,128,0),     #808000
    olivedrab = Color(107,142,35),     #6b8e23
    orange = Color(255,165,0),     #ffa500
    orangered = Color(255,69,0),     #ff4500
    orchid = Color(218,112,214),     #da70d6
    palegoldenrod = Color(238,232,170),     #eee8aa
    palegreen = Color(152,251,152),     #98fb98
    paleturquoise = Color(175,238,238),     #afeeee
    palevioletred = Color(219,112,147),     #db7093
    papayawhip = Color(255,239,213),     #ffefd5
    peachpuff = Color(255,218,185),     #ffdab9
    peru = Color(205,133,63),     #cd853f
    pink = Color(255,192,203),     #ffc0cb
    plum = Color(221,160,221),     #dda0dd
    powderblue = Color(176,224,230),     #b0e0e6
    purple = Color(128,0,128),     #800080
    red = Color(255,0,0),     #ff0000
    rosybrown = Color(188,143,143),     #bc8f8f
    royalblue = Color(65,105,225),     #4169e1
    saddlebrown = Color(139,69,19),     #8b4513
    salmon = Color(250,128,114),     #fa8072
    sandybrown = Color(244,164,96),     #f4a460
    seagreen = Color(46,139,87),     #2e8b57
    seashell = Color(255,245,238),     #fff5ee
    sienna = Color(160,82,45),     #a0522d
    silver = Color(192,192,192),     #c0c0c0
    skyblue = Color(135,206,235),     #87ceeb
    slateblue = Color(106,90,205),     #6a5acd
    slategray = Color(112,128,144),     #708090
    slategrey = Color(112,128,144),     #708090
    snow = Color(255,250,250),     #fffafa
    springgreen = Color(0,255,127),     #00ff7f
    steelblue = Color(70,130,180),     #4682b4
    tan = Color(210,180,140),     #d2b48c
    teal = Color(0,128,128),     #008080
    thistle = Color(216,191,216),     #d8bfd8
    tomato = Color(255,99,71),     #ff6347
    turquoise = Color(64,224,208),     #40e0d0
    violet = Color(238,130,238),     #ee82ee
    wheat = Color(245,222,179),     #f5deb3
    white = Color(255,255,255),     #ffffff
    whitesmoke = Color(245,245,245),     #f5f5f5
    yellow = Color(255,255,0),     #ffff00
    yellowgreen = Color(154,205,50)     #9acd32
    )
