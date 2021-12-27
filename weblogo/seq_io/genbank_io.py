#!/usr/bin/env python


"""Read GenBank flat files.

Currently only reads sequence data and not annotations.

"""

from typing import Iterator, TextIO

from ..seq import Alphabet, Seq, SeqList
from ..utils import isblank

names = ("genbank",)
extensions = ("gb", "genbank", "gbk")


def read(fin: TextIO, alphabet: Alphabet = None) -> SeqList:
    """Read and parse a file of genbank records.

    Args:
    fin -- A stream or file to read
    alphabet -- The expected alphabet of the data, if given

    Returns:
    SeqList -- A list of sequences

    Raises:
    ValueError -- If the file is unparsable
    """
    seqs = [s for s in iterseq(fin, alphabet)]
    return SeqList(seqs)


def iterseq(fin: TextIO, alphabet: Alphabet = None) -> Iterator[Seq]:
    """Iterate over genbank records

    Args:
    fin -- A stream or file to read
    alphabet -- The expected alphabet of the data, if given

    Yields:
    Seq -- One alphabetic sequence at a time.

    Raises:
    ValueError -- If the file is unparsable
    """
    alphabet = Alphabet(alphabet)

    header, block, data = range(3)
    state = header
    seq: list = []
    for L, line in enumerate(fin):
        if isblank(line):
            continue
        if state == header:
            if not line.startswith("LOCUS"):
                raise ValueError("Cannot find start of record at line %d" % L)
            state = block
        elif state == block:
            if line.startswith("ORIGIN") or line.startswith("//"):
                state = data
        elif state == data:
            if line.startswith("//"):
                yield Seq("".join(seq), alphabet)
                seq = []
                state = block
            else:
                seq.extend(line.split()[1:])
