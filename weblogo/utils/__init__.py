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
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
#  OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
#  THE SOFTWARE.
#


"""Extra utilities and core classes not in standard python.
"""


from typing import Any, Iterable, List, TextIO

import pkg_resources

__all__ = (
    "isblank",
    "isfloat",
    "isint",
    "ischar",
    "remove_whitespace",
    "invert_dict",
    "stdrepr",
    "Token",
    "deoptparse",
    "crc32",
    "crc64",
    "ArgumentError",
    "group_count",
    "resource_string",
    "resource_stream",
    "resource_filename",
)


def isblank(s: Any) -> bool:
    """Is this whitespace or an empty string?"""
    if isinstance(s, str):
        if not s:
            return True
        else:
            return s.isspace()
    else:
        return False


def isfloat(s: Any) -> bool:
    """Does this object represent a floating point number?"""
    try:
        float(s)
        return True
    except (ValueError, TypeError):
        return False


def isint(s: Any) -> bool:
    """Does this object represent an integer?"""
    try:
        int(s)
        return True
    except (ValueError, TypeError):
        return False


def ischar(s: Any) -> bool:
    """Does this object represent a character?"""
    # Adapted from: https://stackoverflow.com/a/14321721 and
    # https://mail.python.org/pipermail/python-list/2007-March/425058.html
    return isinstance(s, str) and bool(s) and s == len(s) * s[0]


def remove_whitespace(astring: str) -> str:
    """Remove all whitespace from a string."""
    # TODO: Is this horrible slow?
    return "".join(astring.split())


def invert_dict(dictionary: dict) -> dict:
    """Constructs a new dictionary with inverted mappings so that keys
    become values and vice versa. If the values of the original dictionary
    are not unique then only one of the original keys will be included
    in the new dictionary.
    """
    return dict((value, key) for key, value in dictionary.items())


def stdrepr(obj: Any, attributes: List[str] = None, name: str = None) -> str:
    """Create a standard representation of an object."""
    if name is None:
        name = obj.__class__.__name__
    if attributes is None:
        attributes = obj.__class__.__slots__
    args = []
    for a in attributes:
        if a[0] == "_":
            continue  # pragma: no cover
        args.append("%s=%s" % (a, repr(getattr(obj, a))))
    arg_str = ",\n".join(args).replace("\n", "\n    ")

    return "%s(\n    %s\n)" % (name, arg_str)


def group_count(i: Iterable) -> list:
    """An iteration that returns tuples of items and the number of consecutive
    occurrences. Thus group_count('aabbbc') yields ('a',2), ('b',3), ('c',1)
    """
    from itertools import groupby

    return [(item, sum(1 for n in group)) for item, group in groupby(i)]


class Token:
    """Represents the items returned by a file scanner, normally processed
    by a parser.

    Attributes :
    o typeof    -- a string describing the kind of token
    o data      -- the value of the token
    o lineno    -- the line of the file on which the data was found (if known)
    o offset    -- the offset of the data within the line (if known)
    """

    __slots__ = ["typeof", "data", "lineno", "offset"]

    def __init__(
        self, typeof: str, data: str = None, lineno: int = -1, offset: int = -1
    ) -> None:
        self.typeof = typeof
        self.data = data
        self.lineno = lineno
        self.offset = offset

    def __repr__(self) -> str:
        return stdrepr(self)

    def __str__(self) -> str:
        coord = str(self.lineno)
        if self.offset != -1:
            coord += ":" + str(self.offset)
        coord = coord.ljust(7)
        return (coord + "  " + self.typeof + " : ").ljust(32) + str(self.data or "")


def crc32(string: str) -> str:
    """Return the standard CRC32 checksum as a hexidecimal string."""
    import binascii

    return "%08X" % binascii.crc32(string.encode())


_crc64_table = None


def crc64(string: str) -> str:
    """Calculate ISO 3309 standard cyclic redundancy checksum.
    Used, for example, by SWISS-PROT.

    Returns : The CRC as a hexadecimal string.

    Reference:
    o W. H. Press, S. A. Teukolsky, W. T. Vetterling, and B. P. Flannery,
     "Numerical recipes in C", 2nd ed., Cambridge University Press.
     Pages 896ff.
    """
    # Adapted from biopython, which was adapted from bioperl
    global _crc64_table
    if _crc64_table is None:
        # Initialisation of CRC64 table
        table = []
        for i in range(256):
            k = i
            part_h = 0
            for j in range(8):
                rflag = k & 1
                k >>= 1
                if part_h & 1:
                    k |= 1 << 31  # pragma: no cover
                part_h >>= 1
                if rflag:
                    part_h ^= 0xD8000000
            table.append(part_h)
        _crc64_table = tuple(table)

    crcl = 0
    crch = 0
    for c in string:
        shr = (crch & 0xFF) << 24
        temp1h = crch >> 8
        temp1l = (crcl >> 8) | shr
        idx = (crcl ^ ord(c)) & 0xFF
        crch = temp1h ^ _crc64_table[idx]
        crcl = temp1l

    return "%08X%08X" % (crch, crcl)


# End crc64


class ArgumentError(ValueError):
    """A subclass of ValueError raised when a function receives an argument
    that has the right type but an inappropriate value, and the situation is
    not described by a more precise exception such as IndexError. The name of
    the argument or component at fault and (optionally) the value
    are also stored.
    """

    def __init__(self, message: str, key: str, value: Any = None) -> None:
        """Args:
        - msg -- An error message.
        - key -- The name of the argument or component at fault.
        - value -- Optional value of the argument.
        """
        ValueError.__init__(self, key, message)
        # Changed .message to .msg because of deprecation warning in python 2.6
        self.msg = message
        self.key = key
        self.value = value


# end class ArgumentError


# TODO: Replace with direct calls to pkg_resources


def resource_string(modulename: str, resource: str, basefilename: str = None) -> bytes:
    """Locate and return a resource as a string.
    >>> f = resource_string( __name__, 'somedatafile', __file__)
    """
    return pkg_resources.resource_string(modulename, resource)


def resource_stream(modulename: str, resource: str, basefilename: str = None) -> TextIO:
    """Locate and return a resource as a stream.
    >>> f = resource_stream( __name__, 'somedatafile', __file__)
    """
    return open(resource_filename(modulename, resource, basefilename))


def resource_filename(modulename: str, resource: str, basefilename: str = None) -> str:
    """Locate and return a resource filename.
    >>> f = resource_filename( __name__, 'somedatafile', __file__)
    """
    return pkg_resources.resource_filename(modulename, resource)
