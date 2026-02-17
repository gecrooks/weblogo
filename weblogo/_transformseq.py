#!/usr/bin/env python

#  Copyright (c) 2006, The Regents of the University of California, through
#  Lawrence Berkeley National Laboratory (subject to receipt of any required
#  approvals from the U.S. Dept. of Energy).  All rights reserved.

#  This software is distributed under the new BSD Open Source License.
#  <http://www.opensource.org/licenses/bsd-license.html>
#  Please see the LICENSE.txt file that should have been included
#  as part of this package.


"""A tool for converting multiple sequence alignments from
one format to another.
"""

import argparse
import sys

from weblogo import seq_io
from weblogo._cli import _lookup
from weblogo.seq import Seq, SeqList, nucleic_alphabet
from weblogo.transform import mask_low_complexity

__version__ = "1.0.0"
description = """ A tool for converting multiple sequence alignments from
one format to another. """


def main() -> None:
    # ------ Parse Command line ------
    parser = _build_argument_parser()
    opts = parser.parse_args(sys.argv[1:])

    seqs = opts.reader.read(opts.fin)

    if opts.trans_seg:
        seqs = SeqList([mask_low_complexity(s) for s in seqs])

    if opts.subsample is not None:
        from random import random

        frac = opts.subsample
        ss = []
        for s in seqs:
            if random() < frac:
                ss.append(s)
        seqs = SeqList(ss)

    if opts.reverse:
        seqs = SeqList([s.reverse() for s in seqs])

    if opts.complement:
        seqs = SeqList([Seq(s, alphabet=nucleic_alphabet) for s in seqs])
        seqs = SeqList([s.complement() for s in seqs])

    opts.writer.write(opts.fout, seqs)


def _build_argument_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        usage="%(prog)s [options]  < sequence_data.fa > output.fa",
        description=description,
    )
    parser.add_argument(
        "--version", action="version", version=__version__
    )

    io_grp = parser.add_argument_group("Input/Output Options")

    io_grp.add_argument(
        "-f",
        "--fin",
        dest="fin",
        type=argparse.FileType("r"),
        default=sys.stdin,
        help="Sequence input file (default: stdin)",
        metavar="FILENAME",
    )

    io_grp.add_argument(
        "--format-in",
        dest="reader",
        type=_lookup(seq_io.format_names(), "format"),
        default=seq_io,
        help=f"Multiple sequence alignment format: ({', '.join([f.names[0] for f in seq_io.formats])})",  # type: ignore
        metavar="FORMAT",
    )

    io_grp.add_argument(
        "-o",
        "--fout",
        dest="fout",
        type=argparse.FileType("w"),
        default=sys.stdout,
        help="Output file (default: stdout)",
        metavar="FILENAME",
    )

    trans_grp = parser.add_argument_group("Transformations")

    trans_grp.add_argument(
        "--seg",
        dest="trans_seg",
        action="store_true",
        default=False,
        help="Mask low complexity regions in protein sequences.",
    )

    trans_grp.add_argument(
        "--subsample",
        dest="subsample",
        type=float,
        default=None,
        help="Return a random subsample of the sequences.",
        metavar="FRACTION",
    )

    trans_grp.add_argument(
        "--reverse",
        dest="reverse",
        action="store_true",
        default=False,
        help="reverse sequences",
    )

    trans_grp.add_argument(
        "--complement",
        dest="complement",
        action="store_true",
        default=False,
        help="complement DNA sequences",
    )

    # Writers
    out_formats = []
    for f in seq_io.formats:
        if hasattr(f, "write"):
            out_formats.append(f)
    out_choices = {}
    for f in out_formats:
        out_choices[f.names[0]] = f  # type: ignore
    out_names = [f.names[0] for f in out_formats]  # type: ignore

    io_grp.add_argument(
        "-F",
        "--format-out",
        dest="writer",
        type=_lookup(out_choices, "output format"),
        default=seq_io.fasta_io,
        help=f"Multiple sequence alignment output format: ({', '.join(out_names)}) (Default: fasta)",
        metavar="FORMAT",
    )

    return parser
