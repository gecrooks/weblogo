#  Copyright (c) 2005 Gavin E. Crooks
#  Copyright (c) 2006 John Gilman

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
#  LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
#  FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS
#  IN THE SOFTWARE.


"""
Arrays indexed by alphabetic strings.
"""

from array import array
from typing import TYPE_CHECKING, Any, List, TextIO, Tuple, Union

import numpy as np

from .seq import (
    Alphabet,
    Seq,
    unambiguous_dna_alphabet,
    unambiguous_protein_alphabet,
    unambiguous_rna_alphabet,
)
from .utils import ischar, isint

if TYPE_CHECKING:
    from numpy.typing import ArrayLike, DTypeLike

__all__ = "AlphabeticArray", "Motif"


class AlphabeticArray(object):
    """An alphabetic array. Wraps a numpy array so that each dimension
    can be associated with an alphabet and indexed with characters or strings.

    Attributes :
    - alphabets -- A sequence of alphabets used to index the array
    - array     -- The underlying array object that is indexed.

    Examples :

    >>> from weblogo.seq import *
    >>> from weblogo.matrix import AlphabeticArray
    >>>
    >>> str(protein_alphabet)
    'ACDEFGHIKLMNOPQRSTUVWYBJZX*-'
    >>> matrix = AlphabeticArray( (protein_alphabet, protein_alphabet) )
    >>>
    >>> # Index by character or integer:
    >>> matrix['A', 'C'] = 10
    >>> matrix[0,1]
    10
    >>>
    >>> # Different alphabets on each dimension:
    >>> import numpy as na
    >>> a234 = zeros( shape = (2,3,4) )
    >>> alpha = ( "AB", "ABC", "ABCD")
    >>> aa = AlphabeticArray(alpha,a234)
    >>> aa['A', 'B', 'C'] = 22
    >>>
    >>> # String indices are converted to integer index arrays:
    ...
    >>> aa['A', 'B', 'ABCD']
    array([ 0,  0, 22,  0])


    Authors:
    o GEC 2005, JXG 2006

    """

    # Design note: Subclassing numpy arrays is hard, so instead we
    # build this proxy wrapper.

    __slots__ = ["alphabets", "array"]

    def __init__(
        self,
        alphabets: Union[Alphabet, Tuple[Union[Alphabet, None], ...], str],
        values: "ArrayLike" = None,
        dtype: "DTypeLike" = None,
    ):
        """
        Args:
        - alphabets -- a list of alphabets (as string or Alphabet objects) to
                    be used to convert strings into indices. The lengths of
                    the alphabets match the shape of the indexed array.
                    Alternatively, an integer or None in the list indicate a
                    non-alphabetic dimension. If None the dimension length is
                    taken from values argument.
        - values -- An array of values to be indexed. If None a new
                 array is created. If this argument is not a numpy array
                 then the alphabet list must be explicit (cannot contain
                 None.)
        - dtype -- An optional numpy type code.
        """

        # A dummy object to be used in place of None in the alphabets list
        # so that we get meaningful error messages if we try to index a
        # nonalphabetic dimension with a string.
        class NullAlphabet(Alphabet):
            def ord(self, key: str) -> int:
                raise IndexError(
                    "This dimension does not have an alphabet"
                )  # pragma: no cover

            def ords(self, string: Union["Seq", str]) -> array:
                raise IndexError(
                    "This dimension does not have an alphabet"
                )  # pragma: no cover

        alpha: List[Union[Alphabet, None]] = []
        shape: List[Union[int, None]] = []
        for a in alphabets:
            if isinstance(a, str):
                a = Alphabet(a)

            if a is None:
                shape.append(None)
                alpha.append(NullAlphabet(letters=""))
            elif isinstance(a, Alphabet):
                shape.append(len(a))
                alpha.append(a)
            else:
                shape.append(int(a))  # pragma: no cover
                alpha.append(None)  # pragma: no cover

        # shape = tuple(shape) # CHECKME, line not needed?
        if values is None:
            values = np.zeros(shape=shape, dtype=dtype)
        else:
            values = np.asarray(values, dtype=dtype)
            vshape = values.shape
            if len(shape) != len(vshape):
                raise ValueError(
                    "The values array is the wrong shape."
                )  # pragma: no cover
            for s1, s2 in zip(shape, vshape):
                if s1 is not None and s1 != s2:
                    raise ValueError(
                        "The values array is the wrong shape."
                    )  # pragma: no cover
        self.array = values
        self.alphabets = tuple(alpha)

    def __getitem__(self, key: Any) -> Any:
        return self.array.__getitem__(self._ordkey(key))

    def __setitem__(self, key: Any, value: Any) -> None:
        self.array.__setitem__(self._ordkey(key), value)

    def _ordkey(self, key: Any) -> Union[int, np.ndarray, slice, Tuple, Alphabet]:
        """Convert string indices into integers. Handles characters, strings
        slices with strings, and tuples of the same. Anything else is
        unchanged.
        """

        def norm(
            key: Any, alpha: Alphabet
        ) -> Union[int, np.ndarray, slice, Tuple, Alphabet, None]:
            if key is None:
                return None
            elif isinstance(key, str) or isinstance(key, Alphabet):
                assert alpha is not None
                key = str(key)
                if len(key) == 1:
                    return alpha.ord(key)
                if len(key) == 0:
                    return None  # pragma: no cover
                return np.asarray(alpha.ords(key))
            elif isinstance(key, slice):
                start = norm(key.start, alpha)  # pragma: no cover
                stop = norm(key.stop, alpha)  # pragma: no cover
                step = key.step  # pragma: no cover
                return slice(start, stop, step)  # pragma: no cover
            else:
                return key

        if isinstance(key, tuple):
            return tuple([norm(k, a) for k, a in zip(key, self.alphabets)])  # type: ignore
        else:
            return norm(key, self.alphabets[0])  # type: ignore # FIXME

    def index(self, keys: Any) -> np.ndarray:
        """Return an array of shape (len(key1), len(key2), ...) whose values
        are indexed by keys.

        a.outerindex( (I,J,K) )[i,j,k] == a.array[I[i],J[j],K[k]]

        """
        # TODO: Above docstring is not very clear.
        # Deep voodoo using numpy indexing
        keys = self._ordkey(keys)
        outerkeys = []
        for i, k in enumerate(keys):
            if k is None:
                k = range(0, self.array.shape[i])
            k = np.asarray(k)
            for j in range(len(keys) - i - 1):
                k = k[..., np.newaxis]
            outerkeys.append(k)
        return self.array.__getitem__(tuple(outerkeys))

    def reindex(self, new_alphabet: Alphabet) -> "AlphabeticArray":
        """Create a new AlphabeticArray with the given alphabet. The new
        alphabet must be a subset of the current alphabet. Useful for
        extracting a submatrix or for permuting the alphabet.
        """
        new_array = self.index(new_alphabet)
        return AlphabeticArray(new_alphabet, new_array)

    # The following code is designed to proxy all attributes
    # of the wrapped array. But I'm not entirely sure that this will work as
    # intended.
    def __getattr__(self, name: str) -> str:
        try:
            return self.array.__getattr__(self, name)
        except AttributeError:
            return getattr(self.array, name)

    def __setattr__(self, name: str, value: str) -> None:
        try:
            return object.__setattr__(self, name, value)
        except AttributeError:  # pragma: no cover
            return setattr(self.array, name, value)  # pragma: no cover


# End class AlphabeticArray


class Motif(AlphabeticArray):
    """A two dimensional array where the second dimension is indexed by an
    Alphabet. Used to represent sequence motifs and similar information.


    Attr:
    - alphabet     -- An Alphabet
    - array        -- A numpy array
    - name         -- The name of this motif (if any) as a string.
    - description  -- The description, if any.

    """

    _TRANSFAC_DELIM_LINES = ["XX", "//"]

    def __init__(
        self,
        alphabet: Alphabet,
        array: "ArrayLike" = None,
        dtype: "DTypeLike" = None,
        name: str = None,
        description: str = None,
        scale: float = None,
    ):
        AlphabeticArray.__init__(self, (None, alphabet), array, dtype)
        self.name = name
        self.description = description
        self.scale = scale

    @property
    def alphabet(self) -> Alphabet:
        assert self.alphabets[1] is not None
        return self.alphabets[1]

    def reindex(self, alphabet: Union[Alphabet, Tuple[Alphabet, ...], str]) -> "Motif":
        return Motif(alphabet, AlphabeticArray.reindex(self, (None, alphabet)))  # type: ignore

    # These methods alter self, and therefore do not return a value.
    # (Compare to Seq objects, where the data is immutable and
    #  therefore methods return a new Seq.)
    # TODO: Should reindex (above) also act on self?

    # Deprecate?
    def reverse(self) -> None:
        """Reverse sequence data"""
        self.array = self.array[::-1]  # view into the original numpy array

    # Deprecate?
    def complement(self) -> None:
        """Complement nucleic acid sequence."""
        from weblogo.seq import Alphabet, Seq

        alphabet = self.alphabet
        complement_alphabet = Alphabet(Seq(str(alphabet), alphabet).complement())
        self.alphabets = (None, complement_alphabet)

        assert alphabet is not None
        m = self.reindex(alphabet)
        self.alphabets = (None, alphabet)
        self.array = m.array

    # Deprecate?
    def reverse_complement(self) -> None:
        """Complements and reverses nucleic acid
        sequence (i.e. the other strand of a DNA sequence.)
        """
        self.reverse()
        self.complement()

    @classmethod
    def read_transfac(
        cls, fin: TextIO, alphabet: Union[Alphabet, str] = None
    ) -> "Motif":
        """Parse a TRANSFAC-format PWM from a file.
        Returns a Motif object, representing the provided
        PWM along with an inferred or provided alphabet.
        """

        items = []

        start = False
        for line in fin:
            if line.isspace() or line[0] == "#":
                continue  # pragma: no cover

            stuff = line.split()

            if stuff[0] == "PO" or stuff[0] == "P0":
                start = True

            # 'XX' delimiters may precede the first motif
            if start:
                if stuff[0] in cls._TRANSFAC_DELIM_LINES:
                    break
                else:
                    items.append(stuff)

        if len(items) < 2:
            raise ValueError("Vacuous file.")

        # Is the first line a header line?
        header = items.pop(0)
        hcols = len(header)
        rows = len(items)
        cols = len(items[0])
        if not (
            header[0] == "PO"
            or header[0] == "P0"
            or hcols == cols - 1
            or hcols == cols - 2
        ):
            raise ValueError("Missing header line!")  # pragma: no cover

        # Do all lines (except the first) contain the same number of items?
        cols = len(items[0])
        for i in range(1, len(items)):
            if cols != len(items[i]):
                raise ValueError(
                    "Inconsistant length, row: {}".format(i)
                )  # pragma: no cover

        # Vertical or horizontal arrangement?
        if header[0] == "PO" or header[0] == "P0":
            header.pop(0)

        position_header = True

        for h in header:
            if not ischar(h):
                raise ValueError(
                    "Expected a single character per header "
                    'item, but got "{}" as one item'.format(h)
                )  # pragma: no cover
            if not isint(h):
                position_header = False

        alphabet_header = False if position_header else True

        # Check row headers
        if alphabet_header:
            for i, r in enumerate(items):
                if not isint(r[0]) and r[0][0] != "P":
                    raise ValueError(
                        "Expected position " "as first item on line {}".format(i)
                    )  # pragma: no cover
                r.pop(0)
                defacto_alphabet_str = "".join(header)
        else:
            a = []  # pragma: no cover
            for i, r in enumerate(items):  # pragma: no cover
                if not ischar(r[0]) and r[0][0] != "P":  # pragma: no cover
                    raise ValueError(
                        "Expected position "  # pragma: no cover
                        "as first item on line {}".format(i)
                    )  # pragma: no cover
                a.append(r.pop(0))  # pragma: no cover
            defacto_alphabet_str = "".join(a)  # pragma: no cover

        # Check defacto_alphabet
        defacto_alphabet = Alphabet(defacto_alphabet_str)

        if alphabet:
            alphabet = Alphabet(str(alphabet))
            if not defacto_alphabet.alphabetic(str(alphabet)):
                # Allow alphabet to be a superset of defacto_alphabet
                alphabet = defacto_alphabet

        else:
            alphabets = (
                unambiguous_rna_alphabet,
                unambiguous_dna_alphabet,
                unambiguous_protein_alphabet,
            )
            for aa in alphabets:
                if defacto_alphabet.alphabetic(str(aa)):
                    alphabet = aa
                    break
            if not alphabet:
                alphabet = defacto_alphabet  # pragma: no cover

        # The last item of each row may be extra cruft. Remove
        if len(items[0]) == len(header) + 1:
            for r in items:
                r.pop()

        # items should now be a list of lists of numbers (as strings)
        rows = len(items)
        cols = len(items[0])
        matrix = np.zeros((rows, cols), dtype=np.float64)
        for rr in range(rows):
            for cc in range(cols):
                matrix[rr, cc] = float(items[rr][cc])

        if position_header:
            matrix.transpose()  # pragma: no cover

        return Motif(defacto_alphabet, matrix).reindex(alphabet)
