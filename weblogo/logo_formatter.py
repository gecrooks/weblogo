"""
Logo formatting. Each formatter is a function f(data, format) that draws
a representation of the logo. The main graphical formatter is eps_formatter. A mapping
'formatters' containing all available formatters . Each formatter returns binary data. The eps
and data formats can decoded to strings, e.g. eps_as_string = eps_data.decode()
"""

import shutil
from subprocess import Popen, PIPE
from math import log
from string import Template
import os

from .utils import resource_string
from .logo import LogoData, LogoFormat
from .color import Color

__all__ = ['pdf_formatter', 'jpeg_formatter', 'svg_formatter', 'png_formatter',
           'png_print_formatter', 'txt_formatter', 'eps_formatter', 'formatters',
           'default_formatter', 'GhostscriptAPI']

std_units = {
    "bits": 1. / log(2),
    "nats": 1.,
    "digits": 1. / log(10),
    "kT": 1.,
    "kJ/mol": 8.314472 * 298.15 / 1000.,
    "kcal/mol": 1.987 * 298.15 / 1000.,
    "probability": None,
}
"""Some text"""


def pdf_formatter(logodata: LogoData, logoformat: LogoFormat) -> bytes:
    """ Generate a logo in PDF format."""
    eps = eps_formatter(logodata, logoformat).decode()
    gs = GhostscriptAPI()
    return gs.convert('pdf', eps, logoformat.logo_width, logoformat.logo_height)


def _bitmap_formatter(logodata: LogoData, logoformat: LogoFormat, device: str) -> bytes:
    eps = eps_formatter(logodata, logoformat).decode()
    gs = GhostscriptAPI()
    return gs.convert(device, eps,
                      logoformat.logo_width, logoformat.logo_height, logoformat.resolution)


def jpeg_formatter(logodata: LogoData, logoformat: LogoFormat) -> bytes:
    """ Generate a logo in JPEG format."""
    return _bitmap_formatter(logodata, logoformat, device="jpeg")


def svg_formatter(logodata: LogoData, logoformat: LogoFormat) -> bytes:
    """ Generate a logo in Scalable Vector Graphics (SVG) format.
    Requires the program 'pdf2svg' be installed.
    """
    pdf = pdf_formatter(logodata, logoformat)

    command = shutil.which('pdf2svg')
    if command is None:
        raise EnvironmentError("Scalable Vector Graphics (SVG) format requires the program"
                               "'pdf2svg'.")  # pragma: no cover

    import tempfile
    fpdfi, fname_pdf = tempfile.mkstemp(suffix=".pdf")
    fsvgi, fname_svg = tempfile.mkstemp(suffix=".svg")
    try:

        fpdf2 = open(fname_pdf, 'w')
        fpdf2.buffer.write(pdf)
        fpdf2.seek(0)

        args = [command, fname_pdf, fname_svg]
        p = Popen(args)
        (out, err) = p.communicate()

        fsvg = open(fname_svg)
        return fsvg.read().encode()
    finally:
        os.remove(fname_svg)
        os.remove(fname_pdf)


def png_formatter(logodata: LogoData, logoformat: LogoFormat) -> bytes:
    """ Generate a logo in PNG format."""
    return _bitmap_formatter(logodata, logoformat, device="png")


def png_print_formatter(logodata: LogoData, logoformat: LogoFormat) -> bytes:
    """ Generate a logo in PNG format with print quality (600 DPI) resolution."""
    logoformat.resolution = 600
    return _bitmap_formatter(logodata, logoformat, device="png")


def txt_formatter(logodata: LogoData, logoformat: LogoFormat) -> bytes:
    """ Create a text representation of the logo data.
    """
    return str(logodata).encode()


def eps_formatter(logodata: LogoData, logoformat: LogoFormat) -> bytes:
    """ Generate a logo in Encapsulated Postscript (EPS)"""
    substitutions = {}
    from_format = [
        "creation_date", "logo_width", "logo_height",
        "lines_per_logo", "line_width", "line_height",
        "line_margin_right", "line_margin_left", "line_margin_bottom",
        "line_margin_top", "title_height", "xaxis_label_height",
        "creator_text", "logo_title", "logo_margin",
        "stroke_width", "tic_length",
        "stacks_per_line", "stack_margin",
        "yaxis_label", "yaxis_tic_interval", "yaxis_minor_tic_interval",
        "xaxis_label", "xaxis_tic_interval", "number_interval",
        "fineprint", "shrink_fraction", "errorbar_fraction",
        "errorbar_width_fraction",
        "errorbar_gray", "small_fontsize", "fontsize",
        "title_fontsize", "number_fontsize", "text_font",
        "logo_font", "title_font",
        "logo_label", "yaxis_scale", "end_type",
        "debug", "show_title", "show_xaxis",
        "show_xaxis_label", "show_yaxis", "show_yaxis_label",
        "show_boxes", "show_errorbars", "show_fineprint",
        "rotate_numbers", "show_ends", "stack_height",
        "stack_width"
    ]

    for sf in from_format:
        substitutions[sf] = getattr(logoformat, sf)

    substitutions["shrink"] = str(logoformat.show_boxes).lower()

    def format_color(color: Color) -> str:    # (no fold)
        return " ".join(("[", str(color.red), str(color.green),
                         str(color.blue), "]"))

    substitutions["default_color"] = format_color(logoformat.default_color)

    data = []

    # Unit conversion. 'None' for probability units
    conv_factor = std_units[logoformat.unit_name]

    data.append("StartLine")

    seq_from = logoformat.logo_start - logoformat.first_index
    seq_to = logoformat.logo_end - logoformat.first_index + 1

    # seq_index : zero based index into sequence data
    # logo_index : User visible coordinate, first_index based
    # stack_index : zero based index of visible stacks
    for seq_index in range(seq_from, seq_to):
        # logo_index = seq_index + logoformat.first_index
        stack_index = seq_index - seq_from

        if stack_index != 0 and (stack_index % logoformat.stacks_per_line) == 0:
            data.append("")
            data.append("EndLine")
            data.append("StartLine")
            data.append("")

        data.append("(%s) StartStack" % logoformat.annotate[seq_index])

        if conv_factor:
            stack_height = logodata.entropy[seq_index] * std_units[logoformat.unit_name]
        else:
            stack_height = 1.0  # probability   # pragma: no cover

        # Sort by frequency. If equal frequency then reverse alphabetic
        # (So sort reverse alphabetic first, then frequencty)
        # TODO: doublecheck this actual works
        s = list(zip(logodata.counts[seq_index], logodata.alphabet))
        s.sort(key=lambda x: x[1])
        s.reverse()
        s.sort(key=lambda x: x[0])

        if not logoformat.reverse_stacks:
            s.reverse()         # pragma: no cover

        C = float(sum(logodata.counts[seq_index]))
        if C > 0.0:
            fraction_width = 1.0
            if logoformat.scale_width:
                fraction_width = logodata.weight[seq_index]
                # print(fraction_width, file=sys.stderr)
            for rank, c in enumerate(s):
                color = logoformat.color_scheme.symbol_color(seq_index, c[1], rank)
                data.append(" %f %f %s (%s) ShowSymbol" % (
                    fraction_width, c[0] * stack_height / C,
                    format_color(color), c[1]))

        # Draw error bar on top of logo. Replaced by DrawErrorbarFirst above.
        if logodata.entropy_interval is not None and conv_factor and C > 0.0:

            low, high = logodata.entropy_interval[seq_index]
            center = logodata.entropy[seq_index]
            low *= conv_factor
            high *= conv_factor
            center *= conv_factor
            if high > logoformat.yaxis_scale:
                high = logoformat.yaxis_scale       # pragma: no cover

            down = (center - low)
            up = (high - center)
            data.append(" %f %f DrawErrorbar" % (down, up))

        data.append("EndStack")
        data.append("")

    data.append("EndLine")
    substitutions["logo_data"] = "\n".join(data)

    # Create and output logo
    template = resource_string(__name__, 'template.eps', __file__).decode()
    logo = Template(template).substitute(substitutions)

    return logo.encode()


formatters = {
    'eps': eps_formatter,
    'pdf': pdf_formatter,
    'png': png_formatter,
    'png_print': png_print_formatter,
    'jpeg': jpeg_formatter,
    'svg': svg_formatter,
    'logodata': txt_formatter,
}
"""Map between output format names and corresponing logo formatter"""


default_formatter = eps_formatter
"""The default logo formatter."""


class GhostscriptAPI(object):
    """Interface to the command line program Ghostscript ('gs')"""

    formats = ('png', 'pdf', 'jpeg')

    def __init__(self, path: os.PathLike = None) -> None:
        """
        Raises:
            EnvironmentError: If cannot find Ghostscript executable on
                path
        """
        command = shutil.which('gs', path=path)
        if command is None:
            command = shutil.which('gswin64c.exe', path=path)   # pragma: no cover
        if command is None:
            command = shutil.which('gswin32c.exe', path=path)   # pragma: no cover
        if command is None:
            raise EnvironmentError("Could not find Ghostscript on path. "
                                   "There should be either a gs executable or a gswin32c.exe on "
                                   "your system's path")    # pragma: no cover

        self.command = command

    def version(self) -> str:
        """Returms: The ghostscript version string"""
        args = [self.command, '--version']
        try:
            p = Popen(args, stdout=PIPE)
            (out, err) = p.communicate()
        except OSError:     # pragma: no cover
            raise RuntimeError("Cannot communicate with ghostscript.")  # pragma: no cover
        return out.strip()

    def convert(self,
                format: str,
                postscript: str,
                width: int,
                height: int,
                resolution: int = 300) -> bytes:
        """Convert a string of postscript into a different graphical format

        Supported foramts are 'png', 'pdf', and 'jpeg'.

        Raises:
            ValueError: For an unregonized format.
        """
        device_map = {'png': 'png16m', 'pdf': 'pdfwrite', 'jpeg': 'jpeg'}

        try:
            device = device_map[format]
        except KeyError:                                # pragma: no cover
            raise ValueError("Unsupported format.")

        args = [self.command,
                "-sDEVICE=%s" % device,
                "-dPDFSETTINGS=/printer",
                # "-q",   # Quite: Do not dump messages to stdout.
                "-sstdout=%stderr",  # Redirect messages and errors to stderr
                # fix issue 36, problems with ghostscript 9.10
                "-dColorConversionStrategy=/LeaveColorUnchanged",
                "-sOutputFile=-",  # Stdout
                "-dDEVICEWIDTHPOINTS=%s" % str(width),
                "-dDEVICEHEIGHTPOINTS=%s" % str(height),
                "-dSAFER",  # For added security
                "-dNOPAUSE", ]

        if device != 'pdf':
            args.append("-r%s" % str(resolution))
            if resolution < 300:  # Antialias if resolution is Less than 300 DPI
                args.append("-dGraphicsAlphaBits=4")
                args.append("-dTextAlphaBits=4")
                args.append("-dAlignToPixels=0")

        args.append("-")  # Read from stdin. Must be last argument.

        error_msg = "Unrecoverable error : Ghostscript conversion failed " \
                    "(Invalid postscript?). %s" % " ".join(args)

        try:
            p = Popen(args, stdin=PIPE, stdout=PIPE, stderr=PIPE)
            (out, err) = p.communicate(postscript.encode())
        except OSError:                         # pragma: no cover
            raise RuntimeError(error_msg)       # pragma: no cover

        if p.returncode != 0:                   # pragma: no cover
            error_msg += '\nReturn code: %i\n' % p.returncode
            if err is not None:
                error_msg += err
            raise RuntimeError(error_msg)

        return out

# end class Ghostscript
