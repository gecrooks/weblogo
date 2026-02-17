#!/usr/bin/env python

# -------------------------------- WebLogo --------------------------------

#  Copyright (c) 2003-2004 The Regents of the University of California.
#  Copyright (c) 2005 Gavin E. Crooks
#  Copyright (c) 2006-2011, The Regents of the University of California, through
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

# WebLogo Command Line Interfaceg


import argparse
import atexit
import os
import sys
from contextlib import ExitStack
from io import StringIO
from os import PathLike
from typing import Any, Callable

import importlib_resources

from . import (
    LogoData,
    LogoFormat,
    LogoOptions,
    default_formatter,
    description,
    formatters,
    parse_prior,
    read_seq_data,
    release_description,
    seq_io,
    std_alphabets,
    std_color_schemes,
    std_sizes,
    std_units,
)
from .colorscheme import ColorScheme, SymbolColor
from .logo import _seq_formats, _seq_names
from .seq import Seq, SeqList, nucleic_alphabet


# ====================== Main: Parse Command line =============================
def main() -> None:
    """WebLogo command line interface"""

    # ------ Parse Command line ------
    parser = _build_argument_parser()
    opts = parser.parse_args(sys.argv[1:])

    if opts.serve:
        httpd_serve_forever(opts.port)  # Never returns?    # pragma: no cover
        sys.exit(0)  # pragma: no cover

    # ------ Create Logo ------
    try:
        data = _build_logodata(opts)
        format = _build_logoformat(data, opts)

        formatter = opts.formatter
        logo = formatter(data, format)

        opts.fout.buffer.write(logo)

    except ValueError as err:
        print("Error:", err, file=sys.stderr)
        sys.exit(2)
    except KeyboardInterrupt:  # pragma: no cover
        sys.exit(0)


def httpd_serve_forever(port: int = 8080) -> None:
    """Start a webserver on a local port."""

    import http.server as server
    import http.server as cgiserver

    class __HTTPRequestHandler(cgiserver.CGIHTTPRequestHandler):
        # Modify CGIHTTPRequestHandler so that it will run the cgi script directly,
        # instead of exec'ing
        # This bypasses the need for the cgi script to have execute permissions set,
        # since distutils install does not preserve permissions.
        def is_cgi(self) -> bool:
            self.have_fork = False  # Prevent CGIHTTPRequestHandler from using fork
            if self.path == "/create.cgi":
                self.cgi_info = "", "create.cgi"
                return True
            return False

        def is_python(
            self, path: str | PathLike
        ) -> bool:  # Let CGIHTTPRequestHandler know that cgi script is python
            return True

    # Add current directory to PYTHONPATH. This is
    # so that we can run the standalone server
    # without having to run the install script.
    pythonpath = os.getenv("PYTHONPATH", "")
    pythonpath += os.pathsep + os.path.abspath(sys.path[0])  # .split()[0]
    os.environ["PYTHONPATH"] = pythonpath

    file_manager = ExitStack()
    atexit.register(file_manager.close)
    ref = importlib_resources.files("weblogo") / "htdocs"
    path = file_manager.enter_context(importlib_resources.as_file(ref))

    os.chdir(path)

    HandlerClass = __HTTPRequestHandler
    ServerClass = server.HTTPServer
    httpd = ServerClass(("", port), HandlerClass)
    print(f"WebLogo server running at http://localhost:{port}/")

    try:
        httpd.serve_forever()
    except KeyboardInterrupt:
        sys.exit(0)


def _build_logodata(options: Any) -> LogoData:
    motif_flag = False

    fin = options.fin

    if options.upload is None:
        if fin is None:
            fin = StringIO(sys.stdin.read())
    else:
        if fin is None:
            from .logo import _from_URL_fileopen

            fin = _from_URL_fileopen(options.upload)
        else:
            raise ValueError("error: options --fin and --upload are incompatible")

    try:
        # Try reading data in transfac format first.
        from .matrix import Motif

        motif = Motif.read_transfac(fin, alphabet=options.alphabet)
        motif_flag = True
    except ValueError as motif_err:
        # Failed reading Motif, try reading as multiple sequence data.
        if options.input_parser == "transfac":
            raise motif_err  # Adding transfac as str insted of parser is a bit of a ugly kludge
        seqs = read_seq_data(
            fin,
            options.input_parser.read,
            alphabet=options.alphabet,
            ignore_lower_case=options.ignore_lower_case,
        )

    if motif_flag:
        if options.ignore_lower_case:
            raise ValueError(
                "error: option --ignore-lower-case incompatible with matrix input"
            )
        if options.reverse or options.revcomp:
            motif.reverse()
        if options.complement or options.revcomp:
            motif.complement()

        prior = parse_prior(options.composition, motif.alphabet, options.weight)
        data = LogoData.from_counts(motif.alphabet, motif.array, prior)
    else:
        if options.reverse or options.revcomp:
            seqs = SeqList([s.reverse() for s in seqs], seqs.alphabet)

        if options.complement or options.revcomp:
            if not nucleic_alphabet.alphabetic(str(seqs.alphabet)):
                raise ValueError(
                    "non-nucleic sequence cannot be complemented"
                )  # pragam: no cover
            aaa = seqs.alphabet
            seqs.alphabet = nucleic_alphabet
            seqs = SeqList(
                [Seq(s, seqs.alphabet).complement() for s in seqs], seqs.alphabet
            )
            seqs.alphabet = aaa

        a = seqs.alphabet
        assert a is not None
        prior = parse_prior(options.composition, a, options.weight)
        data = LogoData.from_seqs(seqs, prior)

    return data


def _build_logoformat(logodata: LogoData, opts: Any) -> LogoFormat:
    """Extract and process relevant option values and return a
    LogoFormat object."""

    args = {}
    direct_from_opts = [
        # Logo Data Options.
        "alphabet",
        "unit_name",
        "first_index",
        "logo_start",
        "logo_end",
        # Logo Format Options.
        "stack_width",
        "stacks_per_line",
        "logo_title",
        "logo_label",
        "show_xaxis",
        "xaxis_label",
        "annotate",
        "rotate_numbers",
        "number_interval",
        "yaxis_scale",
        "show_yaxis",
        "yaxis_label",
        "show_ends",
        "fineprint",
        "yaxis_tic_interval",
        "show_errorbars",
        "reverse_stacks",
        # Color Options.
        "color_scheme",
        "default_color",
        # Font Format Options.
        "fontsize",
        "title_fontsize",
        "small_fontsize",
        "number_fontsize",
        "text_font",
        "logo_font",
        "title_font",
        # Advanced Format Options.
        "stack_aspect_ratio",
        "show_boxes",
        "resolution",
        "scale_width",
        "debug",
        "errorbar_fraction",
        "errorbar_width_fraction",
        "errorbar_gray",
    ]

    for k in direct_from_opts:
        args[k] = opts.__dict__[k]

    if opts.colors:
        color_scheme = ColorScheme()
        for color, symbols, desc in opts.colors:
            try:
                color_scheme.rules.append(SymbolColor(symbols, color, desc))
            except ValueError:
                raise ValueError(f"error: option --color: invalid value: '{color}'")

        args["color_scheme"] = color_scheme

    if opts.annotate:
        args["annotate"] = opts.annotate.split(",")

    logooptions = LogoOptions()
    for a, v in args.items():
        setattr(logooptions, a, v)

    theformat = LogoFormat(logodata, logooptions)
    return theformat


# ========================== Helpers ==========================


def _parse_bool(value: str) -> bool:
    """Parse a boolean string value for argparse."""
    v = value.lower()
    choices = {
        "no": False,
        "false": False,
        "0": False,
        "yes": True,
        "true": True,
        "1": True,
    }
    if v not in choices:
        raise argparse.ArgumentTypeError(
            f"invalid choice: '{value}' (choose from 'yes' or 'no', 'true' or 'false')"
        )
    return choices[v]


def _lookup(choices_dict: dict, label: str) -> Callable:
    """Create an argparse type function that looks up values in a dict."""

    def parse(value: str) -> Any:
        v = value.lower()
        if v not in choices_dict:
            raise argparse.ArgumentTypeError(
                f"invalid choice: '{value}' (choose from '{"', '".join(choices_dict)}')"
            )
        return choices_dict[v]

    parse.__name__ = label
    return parse


# ========================== OPTIONS ==========================
def _build_argument_parser() -> argparse.ArgumentParser:
    defaults = LogoOptions()
    parser = argparse.ArgumentParser(
        usage="%(prog)s [options]  < sequence_data.fa > sequence_logo.pdf",
        description=description,
    )
    parser.add_argument(
        "--version", action="version", version=release_description
    )

    io_grp = parser.add_argument_group("Input/Output Options")
    data_grp = parser.add_argument_group("Logo Data Options")
    trans_grp = parser.add_argument_group(
        "Transformations", "Optional transformations of the sequence data."
    )
    format_grp = parser.add_argument_group(
        "Logo Format Options",
        "These options control the format and display of the logo.",
    )
    color_grp = parser.add_argument_group(
        "Color Options",
        "Colors can be specified using CSS2 syntax. e.g. 'red', '#FF0000', etc",
    )
    font_grp = parser.add_argument_group(
        "Font Format Options",
        "These options provide control over the font sizes and types.",
    )
    advanced_grp = parser.add_argument_group(
        "Advanced Format Options",
        "These options provide fine control over the display of the logo.",
    )
    server_grp = parser.add_argument_group(
        "WebLogo Server", "Run a standalone webserver on a local port."
    )

    # ========================== IO OPTIONS ==========================

    io_grp.add_argument(
        "-f",
        "--fin",
        dest="fin",
        type=argparse.FileType("r"),
        default=None,
        help="Sequence input file (default: stdin)",
        metavar="FILENAME",
    )

    io_grp.add_argument(
        "--upload",
        dest="upload",
        default=None,
        help="Upload input file from URL",
        metavar="URL",
    )

    io_grp.add_argument(
        "-D",
        "--datatype",
        dest="input_parser",
        type=_lookup(_seq_formats(), "datatype"),
        default=seq_io,
        help="Type of multiple sequence alignment or position"
        f" weight matrix file: ({', '.join(_seq_names())})",
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

    io_grp.add_argument(
        "-F",
        "--format",
        dest="formatter",
        type=_lookup(formatters, "format"),
        metavar="FORMAT",
        help="Format of output: pdf (default), png, jpeg, svg, "
        "logodata, csv",
        default=default_formatter,
    )

    # ========================== Data OPTIONS ==========================

    data_grp.add_argument(
        "-A",
        "--sequence-type",
        dest="alphabet",
        type=_lookup(std_alphabets, "sequence type"),
        help="The type of sequence data: 'protein', 'rna' or 'dna'.",
        metavar="TYPE",
    )

    data_grp.add_argument(
        "-a",
        "--alphabet",
        dest="alphabet",
        help="The set of symbols to count, e.g. 'AGTC'. "
        "All characters not in the alphabet are ignored. "
        "If neither the alphabet nor sequence-type are specified then weblogo "
        "will examine the input data and make an educated guess. "
        "See also --sequence-type, --ignore-lower-case",
    )

    data_grp.add_argument(
        "-U",
        "--units",
        dest="unit_name",
        choices=list(std_units.keys()),
        default=defaults.unit_name,
        help="A unit of entropy ('bits' (default), 'nats', 'digits'), or a unit of"
        "free energy ('kT', 'kJ/mol', 'kcal/mol'), or 'probability' for"
        " probabilities",
        metavar="UNIT",
    )

    data_grp.add_argument(
        "--composition",
        dest="composition",
        default="auto",
        help="The expected composition of the sequences: 'auto' (default), "
        "'equiprobable', 'none' (do not perform any compositional "
        "adjustment), a CG percentage, a species name (e.g. 'E. coli', "
        "'H. sapiens'), or an explicit distribution (e.g. \"{'A':10, 'C':40,"
        " 'G':40, 'T':10}\"). The automatic option uses a typical "
        "distribution for proteins and equiprobable distribution for "
        "everything else. ",
        metavar="COMP.",
    )

    data_grp.add_argument(
        "--weight",
        dest="weight",
        type=float,
        default=None,
        help="The weight of prior data.  Default depends on alphabet length",
        metavar="NUMBER",
    )

    data_grp.add_argument(
        "-i",
        "--first-index",
        dest="first_index",
        type=int,
        default=1,
        help="Index of first position in sequence data (default: 1)",
        metavar="INDEX",
    )

    data_grp.add_argument(
        "-l",
        "--lower",
        dest="logo_start",
        type=int,
        help="Lower bound of sequence to display",
        metavar="INDEX",
    )

    data_grp.add_argument(
        "-u",
        "--upper",
        dest="logo_end",
        type=int,
        help="Upper bound of sequence to display",
        metavar="INDEX",
    )

    # ========================== Transformation OPTIONS ==========================

    # FIXME Add test?
    trans_grp.add_argument(
        "--ignore-lower-case",
        dest="ignore_lower_case",
        action="store_true",
        default=False,
        help="Disregard lower case letters and only count upper case letters"
        " in sequences.",
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
        help="complement nucleic sequences",
    )

    trans_grp.add_argument(
        "--revcomp",
        dest="revcomp",
        action="store_true",
        default=False,
        help="reverse complement nucleic sequences",
    )

    # ========================== FORMAT OPTIONS ==========================

    format_grp.add_argument(
        "-s",
        "--size",
        dest="stack_width",
        type=_lookup(std_sizes, "logo size"),
        metavar="LOGOSIZE",
        default=defaults.stack_width,
        help="Specify a standard logo size (small, medium (default), large)",
    )

    format_grp.add_argument(
        "-n",
        "--stacks-per-line",
        dest="stacks_per_line",
        type=int,
        help="Maximum number of logo stacks per logo line. (default: %(default)s)",
        default=defaults.stacks_per_line,
        metavar="COUNT",
    )

    format_grp.add_argument(
        "-t",
        "--title",
        dest="logo_title",
        help="Logo title text.",
        default=defaults.logo_title,
        metavar="TEXT",
    )

    format_grp.add_argument(
        "--label",
        dest="logo_label",
        help="A figure label, e.g. '2a'",
        default=defaults.logo_label,
        metavar="TEXT",
    )

    format_grp.add_argument(
        "-X",
        "--show-xaxis",
        type=_parse_bool,
        default=defaults.show_xaxis,
        metavar="YES/NO",
        help="Display sequence numbers along x-axis? (default: %(default)s)",
    )

    format_grp.add_argument(
        "-x",
        "--xlabel",
        dest="xaxis_label",
        default=defaults.xaxis_label,
        help="X-axis label",
        metavar="TEXT",
    )

    format_grp.add_argument(
        "--annotate",
        dest="annotate",
        default=None,
        help="A comma separated list of custom stack annotations, "
        "e.g. '1,3,4,5,6,7'.  Annotation list must be same length as "
        "sequences.",
        metavar="TEXT",
    )

    format_grp.add_argument(
        "--rotate-numbers",
        dest="rotate_numbers",
        type=_parse_bool,
        default=defaults.rotate_numbers,
        help="Draw X-axis numbers with vertical orientation (default: %(default)s).",
        metavar="YES/NO",
    )

    format_grp.add_argument(
        "--number-interval",
        dest="number_interval",
        type=float,
        default=defaults.number_interval,
        help=f"Distance between numbers on X-axis (default: {defaults.number_interval})",
        metavar="NUMBER",
    )

    format_grp.add_argument(
        "-S",
        "--yaxis",
        dest="yaxis_scale",
        type=float,
        help="Height of yaxis in units. (Default: Maximum value with "
        "uninformative prior.)",
        metavar="NUMBER",
    )

    format_grp.add_argument(
        "-Y",
        "--show-yaxis",
        type=_parse_bool,
        dest="show_yaxis",
        default=defaults.show_yaxis,
        metavar="YES/NO",
        help="Display entropy scale along y-axis? (default: %(default)s)",
    )

    format_grp.add_argument(
        "-y",
        "--ylabel",
        dest="yaxis_label",
        help="Y-axis label (default depends on plot type and units)",
        metavar="TEXT",
    )

    format_grp.add_argument(
        "-E",
        "--show-ends",
        type=_parse_bool,
        default=defaults.show_ends,
        metavar="YES/NO",
        help="Label the ends of the sequence? (default: %(default)s)",
    )

    format_grp.add_argument(
        "-P",
        "--fineprint",
        dest="fineprint",
        default=defaults.fineprint,
        help="The fine print (default: weblogo version)",
        metavar="TEXT",
    )

    format_grp.add_argument(
        "--ticmarks",
        dest="yaxis_tic_interval",
        type=float,
        default=defaults.yaxis_tic_interval,
        help="Distance between ticmarks (default: %(default)s)",
        metavar="NUMBER",
    )

    format_grp.add_argument(
        "--errorbars",
        dest="show_errorbars",
        type=_parse_bool,
        default=defaults.show_errorbars,
        metavar="YES/NO",
        help="Display error bars? (default: %(default)s)",
    )

    format_grp.add_argument(
        "--reverse-stacks",
        dest="reverse_stacks",
        type=_parse_bool,
        default=defaults.show_errorbars,
        metavar="YES/NO",
        help="Draw stacks with largest letters on top? (default: %(default)s)",
    )

    # ========================== Color OPTIONS ==========================

    color_scheme_choices = list(std_color_schemes.keys())
    color_scheme_choices.sort()
    color_grp.add_argument(
        "-c",
        "--color-scheme",
        dest="color_scheme",
        type=_lookup(std_color_schemes, "color scheme"),
        metavar="SCHEME",
        default=None,  # Auto
        help=f"Specify a standard color scheme ({', '.join(color_scheme_choices)})",
    )

    color_grp.add_argument(
        "-C",
        "--color",
        dest="colors",
        action="append",
        metavar="COLOR SYMBOLS DESCRIPTION ",
        nargs=3,
        default=[],
        help="Specify symbol colors, e.g. --color black AG 'Purine' "
        "--color red TC 'Pyrimidine' ",
    )

    color_grp.add_argument(
        "--default-color",
        dest="default_color",
        metavar="COLOR",
        default=defaults.default_color,
        help="Symbol color if not otherwise specified.",
    )

    # ========================== Font options =========================

    font_grp.add_argument(
        "--fontsize",
        dest="fontsize",
        type=float,
        default=defaults.fontsize,
        help=f"Regular text font size in points (default: {defaults.fontsize})",
        metavar="POINTS",
    )

    font_grp.add_argument(
        "--title-fontsize",
        dest="title_fontsize",
        type=float,
        default=defaults.title_fontsize,
        help=f"Title text font size in points (default: {defaults.title_fontsize})",
        metavar="POINTS",
    )

    font_grp.add_argument(
        "--small-fontsize",
        dest="small_fontsize",
        type=float,
        default=defaults.small_fontsize,
        help=f"Small text font size in points (default: {defaults.small_fontsize})",
        metavar="POINTS",
    )

    font_grp.add_argument(
        "--number-fontsize",
        dest="number_fontsize",
        type=float,
        default=defaults.number_fontsize,
        help=f"Axis numbers font size in points (default: {defaults.number_fontsize})",
        metavar="POINTS",
    )

    font_grp.add_argument(
        "--text-font",
        dest="text_font",
        default=defaults.text_font,
        help=f"Specify font for labels (default: {defaults.text_font})",
        metavar="FONT",
    )

    font_grp.add_argument(
        "--logo-font",
        dest="logo_font",
        default=defaults.text_font,
        help=f"Specify font for logo (default: {defaults.logo_font})",
        metavar="FONT",
    )

    font_grp.add_argument(
        "--title-font",
        dest="title_font",
        default=defaults.title_font,
        help=f"Specify font for title (default: {defaults.title_font})",
        metavar="FONT",
    )

    # ========================== Advanced options =========================

    advanced_grp.add_argument(
        "-W",
        "--stack-width",
        dest="stack_width",
        type=float,
        default=defaults.stack_width,
        help=f"Width of a logo stack (default: {defaults.stack_width})",
        metavar="POINTS",
    )

    advanced_grp.add_argument(
        "--aspect-ratio",
        dest="stack_aspect_ratio",
        type=float,
        default=defaults.stack_aspect_ratio,
        help=f"Ratio of stack height to width (default: {defaults.stack_aspect_ratio})",
        metavar="POINTS",
    )

    advanced_grp.add_argument(
        "--box",
        dest="show_boxes",
        type=_parse_bool,
        default=False,
        metavar="YES/NO",
        help="Draw boxes around symbols? (default: no)",
    )

    advanced_grp.add_argument(
        "--resolution",
        dest="resolution",
        type=float,
        default=600,
        help="Bitmap resolution in dots per inch (DPI).  (Default: 600 DPI)"
        " Low resolution bitmaps (DPI<300) are antialiased.",
        metavar="DPI",
    )

    advanced_grp.add_argument(
        "--scale-width",
        dest="scale_width",
        type=_parse_bool,
        default=True,
        metavar="YES/NO",
        help="Scale the visible stack width by the fraction of symbols in the"
        " column?  (I.e. columns with many gaps of unknowns are narrow.)  "
        "(Default: yes)",
    )

    advanced_grp.add_argument(
        "--debug",
        type=_parse_bool,
        default=defaults.debug,
        metavar="YES/NO",
        help="Output additional diagnostic information. (Default: %(default)s)",
    )

    advanced_grp.add_argument(
        "--errorbar-fraction",
        dest="errorbar_fraction",
        type=float,
        default=defaults.errorbar_fraction,
        help=f"Sets error bars display proportion (default: {defaults.errorbar_fraction})",
        metavar="NUMBER",
    )

    advanced_grp.add_argument(
        "--errorbar-width-fraction",
        dest="errorbar_width_fraction",
        type=float,
        default=defaults.errorbar_width_fraction,
        help=f"Sets error bars width display proportion (default: {defaults.errorbar_width_fraction})",
        metavar="NUMBER",
    )

    advanced_grp.add_argument(
        "--errorbar-gray",
        dest="errorbar_gray",
        type=float,
        default=defaults.errorbar_gray,
        help=f"Sets error bars' gray scale percentage (default: {defaults.errorbar_gray})",
        metavar="NUMBER",
    )

    # ========================== Server options =========================
    server_grp.add_argument(
        "--serve",
        dest="serve",
        action="store_true",
        default=False,
        help="Start a standalone WebLogo server for creating sequence logos.",
    )

    server_grp.add_argument(
        "--port",
        dest="port",
        type=int,
        default=8080,
        help="Listen to this local port. (Default: %(default)s)",
        metavar="PORT",
    )

    return parser
