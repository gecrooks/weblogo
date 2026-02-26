import argparse
import shutil
import sys
from io import StringIO
from subprocess import PIPE, Popen
from typing import List, Optional, TextIO
from unittest.mock import MagicMock, patch

import pytest

from . import data_ref
from weblogo._cli import (
    _build_argument_parser,
    _build_logodata,
    _build_logoformat,
    _lookup,
    _parse_bool,
)


def _exec(
    args: List[str],
    outputtext: List[str],
    returncode: int = 0,
    stdin: Optional[TextIO] = None,
    binary: bool = False,
) -> None:
    if not stdin:
        stdin = data_ref("cap.fa").open()
    args = ["weblogo"] + args
    p = Popen(args, stdin=stdin, stdout=PIPE, stderr=PIPE)
    (out, err) = p.communicate()
    if returncode == 0 and p.returncode > 0:
        print(err)
    assert returncode == p.returncode
    if returncode == 0:
        assert len(err) == 0

    if outputtext:
        if binary:
            for item in outputtext:
                assert item.encode("latin-1") in out
        else:
            outs = out.decode()
            for item in outputtext:
                assert item in outs

    stdin.close()


def test_malformed_options() -> None:
    _exec(["--notarealoption"], [], 2)
    _exec(["extrajunk"], [], 2)
    _exec(["-I"], [], 2)


def test_help_option() -> None:
    _exec(["-h"], ["options"])
    _exec(["--help"], ["options"])


# def test_version_option():
#     _exec(['--version'], weblogo.__version__[0:5])


def test_default_build() -> None:
    _exec([], ["%PDF-1.4"], binary=True)


# Format options
def test_width() -> None:
    _exec(["-W", "1234"], [])
    _exec(["--stack-width", "1234"], [])


def test_height() -> None:
    _exec(["-W", "1000"], [])
    _exec(["-W", "1000", "--aspect-ratio", "2"], [])


def test_stacks_per_line() -> None:
    _exec(["-n", "7"], [])
    _exec(["--stacks-per-line", "7"], [])


def test_title() -> None:
    _exec(["-t", "3456"], [])
    _exec(["-t", ""], [])
    _exec(["--title", "3456"], [])


def test_annotate() -> None:
    _exec(["--annotate", "1,2,3,4"], [], 2)
    _exec(["--annotate", "1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,,"], [])


def test_color() -> None:
    _exec(["--color", "black", "AG", "Purine"], [])
    _exec(["--color", "not_a_color", "AG", "Purine"], [], 2)


def test_reverse_complement() -> None:
    _exec(["--complement"], [])
    _exec(["--reverse"], [])
    _exec(["--revcomp"], [])


def test_formats() -> None:
    _exec(["--format", "pdf"], [])
    _exec(["--format", "png"], [])
    _exec(["--format", "jpeg"], [])

    _exec(["--format", "logodata"], [])
    _exec(["--format", "csv"], [])


@pytest.mark.skipif(shutil.which("pdf2svg") is None, reason="requires pdf2svg")
def test_formats_svg() -> None:
    _exec(["--format", "svg"], [])


# ====================== In-process tests for _cli.py ======================


class TestParseBool:
    def test_true_values(self) -> None:
        assert _parse_bool("yes") is True
        assert _parse_bool("true") is True
        assert _parse_bool("1") is True
        assert _parse_bool("YES") is True
        assert _parse_bool("True") is True

    def test_false_values(self) -> None:
        assert _parse_bool("no") is False
        assert _parse_bool("false") is False
        assert _parse_bool("0") is False
        assert _parse_bool("NO") is False
        assert _parse_bool("False") is False

    def test_invalid_value(self) -> None:
        with pytest.raises(argparse.ArgumentTypeError, match="invalid choice"):
            _parse_bool("maybe")


class TestLookup:
    def test_valid_lookup(self) -> None:
        choices = {"foo": 1, "bar": 2}
        parse = _lookup(choices, "test")
        assert parse("foo") == 1
        assert parse("BAR") == 2

    def test_invalid_lookup(self) -> None:
        choices = {"foo": 1}
        parse = _lookup(choices, "test")
        with pytest.raises(argparse.ArgumentTypeError, match="invalid choice"):
            parse("missing")

    def test_name_set(self) -> None:
        parse = _lookup({}, "mylabel")
        assert parse.__name__ == "mylabel"


class TestBuildArgumentParser:
    def test_defaults(self) -> None:
        parser = _build_argument_parser()
        opts = parser.parse_args([])
        assert opts.serve is False
        assert opts.port == 8080
        assert opts.logo_title == ""
        assert opts.show_xaxis is True
        assert opts.show_yaxis is True
        assert opts.small_sample_correction is True

    def test_format_option(self) -> None:
        parser = _build_argument_parser()
        opts = parser.parse_args(["-F", "csv"])
        assert callable(opts.formatter)

    def test_sequence_type(self) -> None:
        parser = _build_argument_parser()
        opts = parser.parse_args(["-A", "dna"])
        assert opts.alphabet is not None

    def test_units(self) -> None:
        parser = _build_argument_parser()
        opts = parser.parse_args(["-U", "bits"])
        assert opts.unit_name == "bits"

    def test_title(self) -> None:
        parser = _build_argument_parser()
        opts = parser.parse_args(["-t", "My Title"])
        assert opts.logo_title == "My Title"

    def test_bool_options(self) -> None:
        parser = _build_argument_parser()
        opts = parser.parse_args(["--show-xaxis", "no", "--show-yaxis", "no"])
        assert opts.show_xaxis is False
        assert opts.show_yaxis is False
        opts = parser.parse_args(["--small-sample-correction", "no"])
        assert opts.small_sample_correction is False

    def test_color_option(self) -> None:
        parser = _build_argument_parser()
        opts = parser.parse_args(
            ["--color", "red", "AG", "Purine", "--color", "blue", "TC", "Pyrimidine"]
        )
        assert len(opts.colors) == 2

    def test_server_options(self) -> None:
        parser = _build_argument_parser()
        opts = parser.parse_args(["--serve", "--port", "9090"])
        assert opts.serve is True
        assert opts.port == 9090

    def test_transform_options(self) -> None:
        parser = _build_argument_parser()
        opts = parser.parse_args(["--reverse", "--complement"])
        assert opts.reverse is True
        assert opts.complement is True

    def test_advanced_options(self) -> None:
        parser = _build_argument_parser()
        opts = parser.parse_args(
            ["--aspect-ratio", "3.0", "--box", "yes", "--resolution", "300"]
        )
        assert opts.stack_aspect_ratio == 3.0
        assert opts.show_boxes is True
        assert opts.resolution == 300


class TestBuildLogodata:
    def test_from_fasta(self) -> None:
        parser = _build_argument_parser()
        opts = parser.parse_args([])
        opts.fin = data_ref("cap.fa").open()
        data = _build_logodata(opts)
        assert data.length > 0
        opts.fin.close()

    def test_reverse(self) -> None:
        from weblogo import LogoData

        parser = _build_argument_parser()
        opts = parser.parse_args(["--reverse"])
        opts.fin = data_ref("cap.fa").open()
        data = _build_logodata(opts)
        assert isinstance(data, LogoData)
        opts.fin.close()

    def test_fin_and_upload_incompatible(self) -> None:
        parser = _build_argument_parser()
        opts = parser.parse_args([])
        opts.fin = StringIO(">seq1\nACGT\n")
        opts.upload = "http://example.com/data.fa"
        with pytest.raises(ValueError, match="incompatible"):
            _build_logodata(opts)

    def test_transfac_input(self) -> None:
        parser = _build_argument_parser()
        opts = parser.parse_args(["-D", "transfac"])
        opts.fin = data_ref("transfac_matrix.txt").open()
        data = _build_logodata(opts)
        assert data.length > 0
        opts.fin.close()

    def test_stdin_input(self) -> None:
        """Read from stdin when no --fin provided"""
        parser = _build_argument_parser()
        opts = parser.parse_args([])
        opts.fin = None
        fasta = ">s1\nACGTACGT\n>s2\nACGTACGT\n>s3\nACGTACGT\n>s4\nACGTACGT\n"
        with patch("weblogo._cli.sys") as mock_sys:
            mock_sys.stdin.read.return_value = fasta
            data = _build_logodata(opts)
        assert data.length > 0

    def test_upload_url(self) -> None:
        """Upload from URL when no --fin provided."""
        parser = _build_argument_parser()
        opts = parser.parse_args([])
        opts.fin = None
        opts.upload = "http://example.com/seqs.fa"
        fasta = ">s1\nACGTACGT\n>s2\nACGTACGT\n>s3\nACGTACGT\n>s4\nACGTACGT\n"
        with patch("weblogo.logo._from_URL_fileopen", return_value=StringIO(fasta)):
            data = _build_logodata(opts)
        assert data.length > 0

    def test_transfac_parse_error_reraise(self) -> None:
        """Explicit transfac format with non-transfac data re-raises."""
        parser = _build_argument_parser()
        opts = parser.parse_args(["-D", "transfac"])
        opts.fin = StringIO(">s1\nACGT\n>s2\nACGT\n")
        with pytest.raises(ValueError):
            _build_logodata(opts)

    def test_transfac_ignore_lower_case_error(self) -> None:
        """Transfac with --ignore-lower-case raises."""
        parser = _build_argument_parser()
        opts = parser.parse_args(["--ignore-lower-case"])
        opts.fin = data_ref("transfac_matrix.txt").open()
        with pytest.raises(ValueError, match="ignore-lower-case"):
            _build_logodata(opts)
        opts.fin.close()

    def test_transfac_reverse(self) -> None:
        """Transfac with --reverse."""
        parser = _build_argument_parser()
        opts = parser.parse_args(["--reverse"])
        opts.fin = data_ref("transfac_matrix.txt").open()
        data = _build_logodata(opts)
        assert data.length > 0
        opts.fin.close()

    def test_transfac_complement(self) -> None:
        """Transfac with --complement."""
        parser = _build_argument_parser()
        opts = parser.parse_args(["--complement"])
        opts.fin = data_ref("transfac_matrix.txt").open()
        data = _build_logodata(opts)
        assert data.length > 0
        opts.fin.close()

    def test_transfac_no_small_sample_correction(self) -> None:
        """Transfac with small sample correction disabled."""
        parser = _build_argument_parser()
        opts = parser.parse_args(["--small-sample-correction", "no"])
        opts.fin = data_ref("transfac_matrix.txt").open()
        data = _build_logodata(opts)
        assert data.length > 0
        opts.fin.close()

    def test_complement_sequences(self) -> None:
        """Complement on DNA sequence data."""
        parser = _build_argument_parser()
        opts = parser.parse_args(["--complement", "-A", "dna"])
        opts.fin = data_ref("cap.fa").open()
        data = _build_logodata(opts)
        assert data.length > 0
        opts.fin.close()

    def test_sequences_no_small_sample_correction(self) -> None:
        """Sequence data with small sample correction disabled."""
        parser = _build_argument_parser()
        opts = parser.parse_args(["--small-sample-correction", "no"])
        opts.fin = data_ref("cap.fa").open()
        data = _build_logodata(opts)
        assert data.length > 0
        opts.fin.close()

    def test_complement_protein_error(self) -> None:
        """Complement on protein data raises ValueError."""
        parser = _build_argument_parser()
        opts = parser.parse_args(["--complement", "-A", "protein"])
        opts.fin = data_ref("Rv3829c.fasta").open()
        with pytest.raises(ValueError, match="non-nucleic"):
            _build_logodata(opts)
        opts.fin.close()


class TestHttpdServe:
    def test_httpd_serve_forever(self) -> None:
        """Test server setup and KeyboardInterrupt exit."""
        from weblogo._cli import httpd_serve_forever

        captured_handler_class = {}

        class FakeHTTPServer:
            def __init__(self, addr, handler_cls):  # type: ignore[no-untyped-def]  # noqa: ARG002
                captured_handler_class["cls"] = handler_cls

            def serve_forever(self) -> None:
                raise KeyboardInterrupt

        with (
            patch("weblogo._cli.importlib_resources") as mock_res,
            patch("weblogo._cli.os.chdir"),
            patch("http.server.HTTPServer", FakeHTTPServer),
            pytest.raises(SystemExit) as exc_info,
        ):
            import tempfile

            tmpdir = tempfile.mkdtemp()
            mock_ref = MagicMock()
            mock_res.files.return_value.__truediv__ = MagicMock(return_value=mock_ref)
            mock_res.as_file.return_value.__enter__ = MagicMock(return_value=tmpdir)
            mock_res.as_file.return_value.__exit__ = MagicMock(return_value=False)
            httpd_serve_forever(port=19876)

        assert exc_info.value.code == 0

        # Now test the handler class methods
        handler_cls = captured_handler_class["cls"]
        handler = handler_cls.__new__(handler_cls)

        # is_cgi with /create.cgi path
        handler.path = "/create.cgi"
        assert handler.is_cgi() is True
        assert handler.cgi_info == ("", "create.cgi")
        assert handler.have_fork is False

        # is_cgi with other path
        handler.path = "/index.html"
        assert handler.is_cgi() is False

        # is_python always returns True
        assert handler.is_python("/anything") is True


class TestBuildLogoformat:
    def test_basic(self) -> None:
        from weblogo import LogoFormat

        parser = _build_argument_parser()
        opts = parser.parse_args(["-t", "Test Title"])
        opts.fin = data_ref("cap.fa").open()
        data = _build_logodata(opts)
        fmt = _build_logoformat(data, opts)
        assert isinstance(fmt, LogoFormat)
        assert fmt.logo_title == "Test Title"
        opts.fin.close()

    def test_with_custom_colors(self) -> None:
        from weblogo import LogoFormat

        parser = _build_argument_parser()
        opts = parser.parse_args(["--color", "red", "AG", "Purine"])
        opts.fin = data_ref("cap.fa").open()
        data = _build_logodata(opts)
        fmt = _build_logoformat(data, opts)
        assert isinstance(fmt, LogoFormat)
        opts.fin.close()

    def test_with_annotate(self) -> None:
        from weblogo import LogoFormat

        parser = _build_argument_parser()
        annot = ",".join(str(i) for i in range(1, 23))
        opts = parser.parse_args(["--annotate", annot])
        opts.fin = data_ref("cap.fa").open()
        data = _build_logodata(opts)
        fmt = _build_logoformat(data, opts)
        assert isinstance(fmt, LogoFormat)
        opts.fin.close()


class TestMain:
    def test_main_pdf(self, tmp_path) -> None:  # type: ignore[no-untyped-def]
        from weblogo._cli import main

        outfile = tmp_path / "out.pdf"
        args = [
            "weblogo",
            "-f",
            str(data_ref("cap.fa")),
            "-o",
            str(outfile),
            "-F",
            "pdf",
        ]
        with patch.object(sys, "argv", args):
            main()
        content = outfile.read_bytes()
        assert content[:5] == b"%PDF-"

    def test_main_csv(self, tmp_path) -> None:  # type: ignore[no-untyped-def]
        from weblogo._cli import main

        outfile = tmp_path / "out.csv"
        args = [
            "weblogo",
            "-f",
            str(data_ref("cap.fa")),
            "-o",
            str(outfile),
            "-F",
            "csv",
        ]
        with patch.object(sys, "argv", args):
            main()
        content = outfile.read_text()
        assert len(content) > 0

    def test_main_logodata(self, tmp_path) -> None:  # type: ignore[no-untyped-def]
        from weblogo._cli import main

        outfile = tmp_path / "out.txt"
        args = [
            "weblogo",
            "-f",
            str(data_ref("cap.fa")),
            "-o",
            str(outfile),
            "-F",
            "logodata",
        ]
        with patch.object(sys, "argv", args):
            main()
        content = outfile.read_text()
        assert len(content) > 0

    def test_main_value_error(self, tmp_path) -> None:  # type: ignore[no-untyped-def]
        """Test that ValueError is caught and exits with code 2."""
        from weblogo._cli import main

        outfile = tmp_path / "out.pdf"
        # Provide an invalid color to trigger ValueError
        args = [
            "weblogo",
            "-f",
            str(data_ref("cap.fa")),
            "-o",
            str(outfile),
            "--color",
            "not_a_color",
            "AG",
            "Purine",
        ]
        with patch.object(sys, "argv", args):
            with pytest.raises(SystemExit) as exc_info:
                main()
            assert exc_info.value.code == 2
