"""
Logo formatting. Each formatter is a function f(data, format) that draws
a representation of the logo. The main graphical formatter is pdf_formatter. A mapping
'formatters' containing all available formatters. Each formatter returns binary data.
"""

import os
import shutil
import tempfile
from subprocess import PIPE, Popen
from .logo import LogoData, LogoFormat

__all__ = [
    "pdf_formatter",
    "jpeg_formatter",
    "svg_formatter",
    "png_formatter",
    "txt_formatter",
    "formatters",
    "default_formatter",
]


def pdf_formatter(logodata: LogoData, logoformat: LogoFormat) -> bytes:
    """Generate a logo in PDF format."""
    from .pdf_formatter import native_pdf_formatter
    return native_pdf_formatter(logodata, logoformat)


def _bitmap_formatter(logodata: LogoData, logoformat: LogoFormat, device: str) -> bytes:
    """Convert native PDF to a bitmap format using Ghostscript."""
    pdf = pdf_formatter(logodata, logoformat)

    command = shutil.which("gs")
    if command is None:
        command = shutil.which("gswin64c.exe")  # pragma: no cover
    if command is None:
        command = shutil.which("gswin32c.exe")  # pragma: no cover
    if command is None:
        raise EnvironmentError(
            "Could not find Ghostscript on path. "
            "There should be either a gs executable or a gswin32c.exe on "
            "your system's path"
        )  # pragma: no cover

    device_map = {"png": "png16m", "jpeg": "jpeg"}
    gs_device = device_map[device]

    resolution = logoformat.resolution
    assert resolution is not None

    _, fname_pdf = tempfile.mkstemp(suffix=".pdf")
    try:
        with open(fname_pdf, "wb") as f:
            f.write(pdf)

        args = [
            command,
            "-sDEVICE=%s" % gs_device,
            "-sstdout=%stderr",
            "-sOutputFile=-",
            "-dDEVICEWIDTHPOINTS=%s" % str(logoformat.logo_width),
            "-dDEVICEHEIGHTPOINTS=%s" % str(logoformat.logo_height),
            "-dSAFER",
            "-dBATCH",
            "-dNOPAUSE",
            "-r%s" % str(resolution),
        ]

        if resolution < 300:
            args.append("-dGraphicsAlphaBits=4")
            args.append("-dTextAlphaBits=4")
            args.append("-dAlignToPixels=0")

        args.append(fname_pdf)

        error_msg = "Ghostscript conversion failed. %s" % " ".join(args)

        try:
            p = Popen(args, stdout=PIPE, stderr=PIPE)
            (out, err) = p.communicate()
        except OSError:  # pragma: no cover
            raise RuntimeError(error_msg)  # pragma: no cover

        if p.returncode != 0:  # pragma: no cover
            error_msg += "\nReturn code: %i\n" % p.returncode
            if err is not None:
                error_msg += str(err)
            raise RuntimeError(error_msg)

        return out
    finally:
        os.remove(fname_pdf)


def jpeg_formatter(logodata: LogoData, logoformat: LogoFormat) -> bytes:
    """Generate a logo in JPEG format."""
    return _bitmap_formatter(logodata, logoformat, device="jpeg")


def svg_formatter(logodata: LogoData, logoformat: LogoFormat) -> bytes:
    """Generate a logo in Scalable Vector Graphics (SVG) format.

    Converts native PDF output to SVG using the external tool 'pdf2svg'.
    """
    pdf = pdf_formatter(logodata, logoformat)

    command = shutil.which("pdf2svg")
    if command is None:
        raise EnvironmentError(
            "Scalable Vector Graphics (SVG) format requires the program 'pdf2svg'."
        )

    _, fname_pdf = tempfile.mkstemp(suffix=".pdf")
    _, fname_svg = tempfile.mkstemp(suffix=".svg")
    try:
        with open(fname_pdf, "wb") as fpdf:
            fpdf.write(pdf)

        args = [command, fname_pdf, fname_svg]
        p = Popen(args)
        p.communicate()

        with open(fname_svg) as fsvg:
            return fsvg.read().encode()
    finally:
        os.remove(fname_pdf)
        os.remove(fname_svg)


def png_formatter(logodata: LogoData, logoformat: LogoFormat) -> bytes:
    """Generate a logo in PNG format."""
    return _bitmap_formatter(logodata, logoformat, device="png")


def txt_formatter(logodata: LogoData, logoformat: LogoFormat) -> bytes:
    """Create a text representation of the logo data."""
    return str(logodata).encode()


def csv_formatter(logodata: LogoData, logoformat: LogoFormat) -> bytes:
    """Create a csv representation of the logo data."""
    return logodata.csv().encode()


formatters = {
    "pdf": pdf_formatter,
    "png": png_formatter,
    "jpeg": jpeg_formatter,
    "svg": svg_formatter,
    "logodata": txt_formatter,
    "csv": csv_formatter,
}
"""Map between output format names and corresponding logo formatter"""


default_formatter = pdf_formatter
"""The default logo formatter."""
