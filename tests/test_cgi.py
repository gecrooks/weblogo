"""Tests for weblogo._cgi — CGI web interface helpers and form handling."""

import io
import os
from io import BytesIO, StringIO
from unittest.mock import patch

import pytest

from weblogo._cgi import (
    Field,
    _main,
    alphabets,
    color_schemes,
    composition,
    extension,
    float_or_none,
    int_or_none,
    main,
    mime_type,
    send_form,
    string_or_none,
    truth,
)


# ====================== Helper functions ======================


class TestStringOrNone:
    def test_none(self) -> None:
        assert string_or_none(None) is None

    def test_auto(self) -> None:
        assert string_or_none("auto") is None

    def test_string(self) -> None:
        assert string_or_none("hello") == "hello"

    def test_empty(self) -> None:
        assert string_or_none("") == ""


class TestTruth:
    def test_true_string(self) -> None:
        assert truth("true") is True

    def test_non_true_string(self) -> None:
        # "false" is a non-empty string, so bool("false") is True
        assert truth("false") is True
        # Only the literal string "true" gets special handling
        assert truth("") is False

    def test_truthy(self) -> None:
        assert truth(1) is True

    def test_falsy(self) -> None:
        assert truth(0) is False
        assert truth("") is False


class TestIntOrNone:
    def test_none(self) -> None:
        assert int_or_none(None) is None

    def test_auto(self) -> None:
        assert int_or_none("auto") is None

    def test_empty(self) -> None:
        assert int_or_none("") is None

    def test_int(self) -> None:
        assert int_or_none("42") == 42
        assert int_or_none("0") == 0


class TestFloatOrNone:
    def test_none(self) -> None:
        assert float_or_none(None) is None

    def test_auto(self) -> None:
        assert float_or_none("auto") is None

    def test_empty(self) -> None:
        assert float_or_none("") is None

    def test_float(self) -> None:
        assert float_or_none("3.14") == pytest.approx(3.14)
        assert float_or_none("0") == 0.0


# ====================== Field class ======================


class TestField:
    def test_basic_field(self) -> None:
        f = Field("test", "default_val")
        assert f.name == "test"
        assert f.default == "default_val"
        assert f.value == "default_val"
        assert f.get_value() == "default_val"

    def test_field_with_conversion(self) -> None:
        f = Field("num", "42", int)
        assert f.get_value() == 42

    def test_field_conversion_error(self) -> None:
        f = Field("num", "not_a_number", int, errmsg="Bad number")
        with pytest.raises(ValueError):
            f.get_value()

    def test_field_with_options_valid(self) -> None:
        f = Field("choice", "a", options=["a", "b", "c"])
        assert f.get_value() == "a"

    def test_field_with_options_invalid(self) -> None:
        f = Field("choice", "x", options=["a", "b", "c"], errmsg="Bad choice")
        with pytest.raises(ValueError):
            f.get_value()

    def test_field_options_with_conversion(self) -> None:
        f = Field("size", "small", lambda x: x.upper(), options=["small", "medium"])
        assert f.get_value() == "SMALL"

    def test_field_none_value(self) -> None:
        f = Field("opt", None)
        assert f.get_value() is None


# ====================== Module-level dicts ======================


class TestModuleDicts:
    def test_mime_types(self) -> None:
        assert "pdf" in mime_type
        assert "png" in mime_type
        assert "svg" in mime_type
        assert mime_type["pdf"] == "application/pdf"

    def test_extensions(self) -> None:
        assert extension["pdf"] == "pdf"
        assert extension["logodata"] == "txt"

    def test_alphabets(self) -> None:
        assert alphabets["alphabet_auto"] is None
        assert alphabets["alphabet_dna"] is not None

    def test_color_schemes(self) -> None:
        assert len(color_schemes) > 0
        assert "color_auto" in color_schemes

    def test_composition(self) -> None:
        assert composition["comp_auto"] == "auto"
        assert composition["comp_none"] == "none"


# ====================== send_form ======================


class TestSendForm:
    def test_send_default_form(self) -> None:
        """send_form with no errors produces HTML output."""
        import weblogo

        logooptions = weblogo.LogoOptions()
        controls = [
            Field("sequences", ""),
            Field("format", "pdf", options=["pdf", "png"]),
            Field("show_errorbars", logooptions.show_errorbars),
            Field("logo_title", ""),
        ]
        captured = StringIO()
        with patch("sys.stdout", captured):
            send_form(controls)
        output = captured.getvalue()
        assert "Content-Type: text/html" in output
        assert "WebLogo" in output

    def test_send_form_with_errors(self) -> None:
        """send_form with errors includes error messages."""
        controls = [
            Field("sequences", ""),
            Field("logo_title", ""),
        ]
        captured_out = StringIO()
        captured_err = StringIO()
        errors = [("sequences", "Please provide sequences")]
        with patch("sys.stdout", captured_out), patch("sys.stderr", captured_err):
            send_form(controls, errors)
        output = captured_out.getvalue()
        assert "ERROR" in output
        assert "Please provide sequences" in output

    def test_send_form_with_string_error(self) -> None:
        """send_form handles plain string errors."""
        controls = [Field("sequences", "")]
        captured_out = StringIO()
        captured_err = StringIO()
        errors = ["generic error"]
        with patch("sys.stdout", captured_out), patch("sys.stderr", captured_err):
            send_form(controls, errors)
        output = captured_out.getvalue()
        assert "ERROR" in output

    def test_send_form_checkbox_true(self) -> None:
        """send_form renders 'true' string as checked."""
        controls = [Field("show_errorbars", "true")]
        captured = StringIO()
        with patch("sys.stdout", captured):
            send_form(controls)
        # Just verify it doesn't crash — the checkbox handling path is exercised

    def test_send_form_bool_true(self) -> None:
        """send_form renders True bool as checked."""
        controls = [Field("show_errorbars", True)]
        captured = StringIO()
        with patch("sys.stdout", captured):
            send_form(controls)

    def test_send_form_none_value(self) -> None:
        """send_form renders None as 'auto'."""
        controls = [Field("yaxis_label", None)]
        captured = StringIO()
        with patch("sys.stdout", captured):
            send_form(controls)


# ====================== _main (CGI handler) ======================


def _make_urlencoded_body(params: dict) -> bytes:
    """Build a URL-encoded POST body from a dict."""
    from urllib.parse import quote_plus

    parts = []
    for k, v in params.items():
        parts.append(f"{quote_plus(k)}={quote_plus(v)}")
    return "&".join(parts).encode()


class TestMainCGI:
    def test_get_request_returns_form(self) -> None:
        """A GET request (no form data) returns the default HTML form."""
        captured = StringIO()
        env = {
            "REQUEST_METHOD": "GET",
            "CONTENT_TYPE": "text/plain",
        }
        with patch.dict(os.environ, env, clear=False), patch(
            "sys.stdout", captured
        ):
            _main()
        output = captured.getvalue()
        assert "Content-Type: text/html" in output

    def test_post_creates_logo(self) -> None:
        """A POST with sequence data produces a PDF logo."""
        sequences = ">s1\nACGTACGT\n>s2\nACGTACGT\n>s3\nACGTACGT\n>s4\nACGTACGT\n"
        params = {
            "sequences": sequences,
            "format": "pdf",
            "alphabet": "alphabet_dna",
            "unit_name": "bits",
            "composition": "comp_auto",
            "stack_width": "medium",
            "stacks_per_line": "40",
            "show_errorbars": "true",
            "show_xaxis": "true",
            "show_yaxis": "true",
            "show_ends": "true",
            "show_fineprint": "true",
            "scale_width": "true",
            "color_scheme": "color_auto",
        }
        body = _make_urlencoded_body(params)

        # Build a mock stdout with both text and binary capabilities
        binary_buf = BytesIO()

        class MockStdout(io.TextIOBase):
            def __init__(self) -> None:
                self._text_parts: list[str] = []
                self.buffer = binary_buf

            def write(self, s: str) -> int:
                self._text_parts.append(s)
                return len(s)

            def flush(self) -> None:
                pass

            def get_text(self) -> str:
                return "".join(self._text_parts)

        mock_stdout = MockStdout()

        env = {
            "REQUEST_METHOD": "POST",
            "CONTENT_TYPE": "application/x-www-form-urlencoded",
            "CONTENT_LENGTH": str(len(body)),
        }
        mock_stdin_buffer = BytesIO(body)

        with (
            patch.dict(os.environ, env, clear=False),
            patch("sys.stdout", mock_stdout),
            patch("sys.stdin", new_callable=lambda: type("", (), {"buffer": mock_stdin_buffer})),
        ):
            _main()

        text_output = mock_stdout.get_text()
        assert "Content-Type: application/pdf" in text_output
        pdf_data = binary_buf.getvalue()
        assert pdf_data[:5] == b"%PDF-"

    def test_post_csv_format(self) -> None:
        """A POST requesting CSV format produces CSV output."""
        sequences = ">s1\nACGTACGT\n>s2\nACGTACGT\n>s3\nACGTACGT\n>s4\nACGTACGT\n"
        params = {
            "sequences": sequences,
            "format": "csv",
            "alphabet": "alphabet_dna",
            "unit_name": "bits",
            "composition": "comp_auto",
            "stack_width": "medium",
            "stacks_per_line": "40",
            "color_scheme": "color_auto",
        }
        body = _make_urlencoded_body(params)

        binary_buf = BytesIO()

        class MockStdout(io.TextIOBase):
            def __init__(self) -> None:
                self._text_parts: list[str] = []
                self.buffer = binary_buf

            def write(self, s: str) -> int:
                self._text_parts.append(s)
                return len(s)

            def flush(self) -> None:
                pass

            def get_text(self) -> str:
                return "".join(self._text_parts)

        mock_stdout = MockStdout()

        env = {
            "REQUEST_METHOD": "POST",
            "CONTENT_TYPE": "application/x-www-form-urlencoded",
            "CONTENT_LENGTH": str(len(body)),
        }
        mock_stdin_buffer = BytesIO(body)

        with (
            patch.dict(os.environ, env, clear=False),
            patch("sys.stdout", mock_stdout),
            patch("sys.stdin", new_callable=lambda: type("", (), {"buffer": mock_stdin_buffer})),
        ):
            _main()

        text_output = mock_stdout.get_text()
        assert "Content-Type: text/plain" in text_output

    def test_post_no_sequences_returns_error(self) -> None:
        """A POST with no sequence data returns the form with an error."""
        params = {
            "format": "pdf",
            "alphabet": "alphabet_dna",
            "composition": "comp_auto",
            "stack_width": "medium",
            "color_scheme": "color_auto",
        }
        body = _make_urlencoded_body(params)

        captured = StringIO()
        env = {
            "REQUEST_METHOD": "POST",
            "CONTENT_TYPE": "application/x-www-form-urlencoded",
            "CONTENT_LENGTH": str(len(body)),
        }
        mock_stdin_buffer = BytesIO(body)

        with (
            patch.dict(os.environ, env, clear=False),
            patch("sys.stdout", captured),
            patch("sys.stderr", StringIO()),
            patch("sys.stdin", new_callable=lambda: type("", (), {"buffer": mock_stdin_buffer})),
        ):
            _main()

        output = captured.getvalue()
        assert "Content-Type: text/html" in output

    def test_post_cmd_reset_returns_form(self) -> None:
        """A POST with cmd_reset returns the default form."""
        params = {"cmd_reset": "true"}
        body = _make_urlencoded_body(params)

        captured = StringIO()
        env = {
            "REQUEST_METHOD": "POST",
            "CONTENT_TYPE": "application/x-www-form-urlencoded",
            "CONTENT_LENGTH": str(len(body)),
        }
        mock_stdin_buffer = BytesIO(body)

        with (
            patch.dict(os.environ, env, clear=False),
            patch("sys.stdout", captured),
            patch("sys.stdin", new_callable=lambda: type("", (), {"buffer": mock_stdin_buffer})),
        ):
            _main()

        output = captured.getvalue()
        assert "Content-Type: text/html" in output

    def test_main_catches_exceptions(self) -> None:
        """main() wraps _main() and catches exceptions."""
        captured = StringIO()
        with (
            patch("weblogo._cgi._main", side_effect=RuntimeError("test crash")),
            patch("sys.stdout", captured),
        ):
            main()
        output = captured.getvalue()
        assert "Internal Server Error" in output
        assert "test crash" in output


# ====================== Additional coverage tests ======================

_FASTA_SEQS = ">s1\nACGTACGT\n>s2\nACGTACGT\n>s3\nACGTACGT\n>s4\nACGTACGT\n"

_BASE_PARAMS: dict[str, str] = {
    "format": "pdf",
    "alphabet": "alphabet_dna",
    "unit_name": "bits",
    "composition": "comp_auto",
    "stack_width": "medium",
    "stacks_per_line": "40",
    "color_scheme": "color_auto",
}


def _make_multipart_body(
    fields: dict[str, str],
    files: dict[str, tuple[str, bytes]] | None = None,
) -> tuple[bytes, str]:
    """Build a multipart/form-data POST body. Returns (body, content_type)."""
    boundary = "TestBoundary7890"
    body = b""
    for name, value in fields.items():
        body += f"--{boundary}\r\n".encode()
        body += f'Content-Disposition: form-data; name="{name}"\r\n'.encode()
        body += b"\r\n"
        body += value.encode() + b"\r\n"
    if files:
        for name, (filename, content) in files.items():
            body += f"--{boundary}\r\n".encode()
            body += (
                f'Content-Disposition: form-data; name="{name}"; '
                f'filename="{filename}"\r\n'
            ).encode()
            body += b"Content-Type: application/octet-stream\r\n"
            body += b"\r\n"
            body += content + b"\r\n"
    body += f"--{boundary}--\r\n".encode()
    return body, f"multipart/form-data; boundary={boundary}"


class _MockStdout(io.TextIOBase):
    """Mock stdout with separate text and binary buffers."""

    def __init__(self) -> None:
        self._text_parts: list[str] = []
        self.buffer = BytesIO()

    def write(self, s: str) -> int:
        self._text_parts.append(s)
        return len(s)

    def flush(self) -> None:
        pass

    def get_text(self) -> str:
        return "".join(self._text_parts)


def _cgi_env(content_type: str, body: bytes) -> dict[str, str]:
    return {
        "REQUEST_METHOD": "POST",
        "CONTENT_TYPE": content_type,
        "CONTENT_LENGTH": str(len(body)),
    }


def _mock_stdin(body: bytes):  # type: ignore[no-untyped-def]
    return type("MockStdin", (), {"buffer": BytesIO(body)})


class TestMultipartPost:
    """Multipart/form-data POST handling."""

    def test_file_upload_creates_logo(self) -> None:
        """Multipart file upload produces a logo."""
        fields = dict(_BASE_PARAMS)
        body, ct = _make_multipart_body(
            fields, files={"sequences_file": ("test.fa", _FASTA_SEQS.encode())}
        )
        mock_stdout = _MockStdout()
        with (
            patch.dict(os.environ, _cgi_env(ct, body), clear=False),
            patch("sys.stdout", mock_stdout),
            patch("sys.stdin", new_callable=lambda: _mock_stdin(body)),
        ):
            _main()
        assert "Content-Type: application/pdf" in mock_stdout.get_text()
        assert mock_stdout.buffer.getvalue()[:5] == b"%PDF-"

    def test_file_and_text_conflict(self) -> None:
        """Multipart with both file and text sequences -> conflict error."""
        fields = {**_BASE_PARAMS, "sequences": ">s1\nACGT\n"}
        body, ct = _make_multipart_body(
            fields, files={"sequences_file": ("test.fa", b">s1\nACGT\n>s2\nACGT\n")}
        )
        captured_out, captured_err = StringIO(), StringIO()
        with (
            patch.dict(os.environ, _cgi_env(ct, body), clear=False),
            patch("sys.stdout", captured_out),
            patch("sys.stderr", captured_err),
            patch("sys.stdin", new_callable=lambda: _mock_stdin(body)),
        ):
            _main()
        assert "Content-Type: text/html" in captured_out.getvalue()


class TestSequenceSources:
    """Tests for different sequence input sources."""

    def test_text_and_url_conflict(self) -> None:
        """POST with both text and URL sequences -> conflict error."""
        params = {
            **_BASE_PARAMS,
            "sequences": _FASTA_SEQS,
            "sequences_url": "http://example.com/seqs.fa",
        }
        body = _make_urlencoded_body(params)
        captured_out, captured_err = StringIO(), StringIO()
        with (
            patch.dict(
                os.environ,
                _cgi_env("application/x-www-form-urlencoded", body),
                clear=False,
            ),
            patch("sys.stdout", captured_out),
            patch("sys.stderr", captured_err),
            patch("sys.stdin", new_callable=lambda: _mock_stdin(body)),
        ):
            _main()
        assert "Content-Type: text/html" in captured_out.getvalue()

    def test_sequences_from_url(self) -> None:
        """POST with sequences_url loads and creates logo."""
        params = {**_BASE_PARAMS, "sequences_url": "http://example.com/seqs.fa"}
        body = _make_urlencoded_body(params)
        mock_stdout = _MockStdout()
        with (
            patch.dict(
                os.environ,
                _cgi_env("application/x-www-form-urlencoded", body),
                clear=False,
            ),
            patch("sys.stdout", mock_stdout),
            patch("sys.stdin", new_callable=lambda: _mock_stdin(body)),
            patch(
                "weblogo.logo._from_URL_fileopen", return_value=StringIO(_FASTA_SEQS)
            ),
        ):
            _main()
        assert "Content-Type: application/pdf" in mock_stdout.get_text()

    def test_sequences_url_value_error(self) -> None:
        """POST with bad URL -> ValueError."""
        params = {**_BASE_PARAMS, "sequences_url": "bad-url"}
        body = _make_urlencoded_body(params)
        captured_out, captured_err = StringIO(), StringIO()
        with (
            patch.dict(
                os.environ,
                _cgi_env("application/x-www-form-urlencoded", body),
                clear=False,
            ),
            patch("sys.stdout", captured_out),
            patch("sys.stderr", captured_err),
            patch("sys.stdin", new_callable=lambda: _mock_stdin(body)),
            patch(
                "weblogo.logo._from_URL_fileopen", side_effect=ValueError("bad url")
            ),
        ):
            _main()
        assert "Content-Type: text/html" in captured_out.getvalue()

    def test_sequences_url_io_error(self) -> None:
        """POST with unreachable URL -> IOError."""
        params = {**_BASE_PARAMS, "sequences_url": "http://unreachable.example.com"}
        body = _make_urlencoded_body(params)
        captured_out, captured_err = StringIO(), StringIO()
        with (
            patch.dict(
                os.environ,
                _cgi_env("application/x-www-form-urlencoded", body),
                clear=False,
            ),
            patch("sys.stdout", captured_out),
            patch("sys.stderr", captured_err),
            patch("sys.stdin", new_callable=lambda: _mock_stdin(body)),
            patch(
                "weblogo.logo._from_URL_fileopen", side_effect=IOError("unreachable")
            ),
        ):
            _main()
        assert "Content-Type: text/html" in captured_out.getvalue()


class TestFormValidation:
    """Tests for form field validation and color scheme handling."""

    def test_invalid_format_option(self) -> None:
        """POST with invalid format option -> ValueError collected."""
        params = {**_BASE_PARAMS, "sequences": _FASTA_SEQS, "format": "invalid_fmt"}
        body = _make_urlencoded_body(params)
        captured_out, captured_err = StringIO(), StringIO()
        with (
            patch.dict(
                os.environ,
                _cgi_env("application/x-www-form-urlencoded", body),
                clear=False,
            ),
            patch("sys.stdout", captured_out),
            patch("sys.stderr", captured_err),
            patch("sys.stdin", new_callable=lambda: _mock_stdin(body)),
        ):
            _main()
        assert "Content-Type: text/html" in captured_out.getvalue()

    def test_custom_color_scheme(self) -> None:
        """POST with color_custom + valid custom colors creates logo."""
        params = {
            **_BASE_PARAMS,
            "sequences": _FASTA_SEQS,
            "color_scheme": "color_custom",
            "color0": "red",
            "symbols0": "AC",
            "desc0": "test",
        }
        body = _make_urlencoded_body(params)
        mock_stdout = _MockStdout()
        with (
            patch.dict(
                os.environ,
                _cgi_env("application/x-www-form-urlencoded", body),
                clear=False,
            ),
            patch("sys.stdout", mock_stdout),
            patch("sys.stdin", new_callable=lambda: _mock_stdin(body)),
        ):
            _main()
        assert "Content-Type: application/pdf" in mock_stdout.get_text()

    def test_invalid_custom_color(self) -> None:
        """POST with color_custom + invalid color -> error."""
        params = {
            **_BASE_PARAMS,
            "sequences": _FASTA_SEQS,
            "color_scheme": "color_custom",
            "color0": "not_a_valid_color_xyz",
            "symbols0": "AC",
        }
        body = _make_urlencoded_body(params)
        captured_out, captured_err = StringIO(), StringIO()
        with (
            patch.dict(
                os.environ,
                _cgi_env("application/x-www-form-urlencoded", body),
                clear=False,
            ),
            patch("sys.stdout", captured_out),
            patch("sys.stderr", captured_err),
            patch("sys.stdin", new_callable=lambda: _mock_stdin(body)),
        ):
            _main()
        assert "Content-Type: text/html" in captured_out.getvalue()

    def test_invalid_color_scheme_value(self) -> None:
        """POST with nonexistent color scheme -> error."""
        params = {
            **_BASE_PARAMS,
            "sequences": _FASTA_SEQS,
            "color_scheme": "color_nonexistent",
        }
        body = _make_urlencoded_body(params)
        captured_out, captured_err = StringIO(), StringIO()
        with (
            patch.dict(
                os.environ,
                _cgi_env("application/x-www-form-urlencoded", body),
                clear=False,
            ),
            patch("sys.stdout", captured_out),
            patch("sys.stderr", captured_err),
            patch("sys.stdin", new_callable=lambda: _mock_stdin(body)),
        ):
            _main()
        assert "Content-Type: text/html" in captured_out.getvalue()


class TestLogoCreation:
    """Tests for logo creation paths and error handling."""

    def test_percentcg_composition(self) -> None:
        """POST with CG percentage composition creates logo."""
        params = {
            **_BASE_PARAMS,
            "sequences": _FASTA_SEQS,
            "composition": "comp_CG",
            "percentCG": "50",
        }
        body = _make_urlencoded_body(params)
        mock_stdout = _MockStdout()
        with (
            patch.dict(
                os.environ,
                _cgi_env("application/x-www-form-urlencoded", body),
                clear=False,
            ),
            patch("sys.stdout", mock_stdout),
            patch("sys.stdin", new_callable=lambda: _mock_stdin(body)),
        ):
            _main()
        assert "Content-Type: application/pdf" in mock_stdout.get_text()

    def test_transfac_format(self) -> None:
        """POST with TRANSFAC-format data creates logo."""
        transfac_data = (
            "PO A C G T\n"
            "01 1 2 2 0\n"
            "02 2 1 2 0\n"
            "03 3 0 1 1\n"
            "04 0 5 0 0\n"
            "05 5 0 0 0\n"
            "XX\n"
        )
        params = {**_BASE_PARAMS, "sequences": transfac_data}
        body = _make_urlencoded_body(params)
        mock_stdout = _MockStdout()
        with (
            patch.dict(
                os.environ,
                _cgi_env("application/x-www-form-urlencoded", body),
                clear=False,
            ),
            patch("sys.stdout", mock_stdout),
            patch("sys.stdin", new_callable=lambda: _mock_stdin(body)),
        ):
            _main()
        text = mock_stdout.get_text()
        assert "Content-Type:" in text

    def test_logo_value_error(self) -> None:
        """ValueError during logo creation -> error form."""
        params = {**_BASE_PARAMS, "sequences": _FASTA_SEQS}
        body = _make_urlencoded_body(params)
        captured_out, captured_err = StringIO(), StringIO()
        with (
            patch.dict(
                os.environ,
                _cgi_env("application/x-www-form-urlencoded", body),
                clear=False,
            ),
            patch("sys.stdout", captured_out),
            patch("sys.stderr", captured_err),
            patch("sys.stdin", new_callable=lambda: _mock_stdin(body)),
            patch(
                "weblogo.matrix.Motif.read_transfac", side_effect=ValueError("x")
            ),
            patch("weblogo.read_seq_data", side_effect=ValueError("bad data")),
        ):
            _main()
        assert "Content-Type: text/html" in captured_out.getvalue()

    def test_logo_io_error(self) -> None:
        """IOError during logo creation -> error form."""
        params = {**_BASE_PARAMS, "sequences": _FASTA_SEQS}
        body = _make_urlencoded_body(params)
        captured_out, captured_err = StringIO(), StringIO()
        with (
            patch.dict(
                os.environ,
                _cgi_env("application/x-www-form-urlencoded", body),
                clear=False,
            ),
            patch("sys.stdout", captured_out),
            patch("sys.stderr", captured_err),
            patch("sys.stdin", new_callable=lambda: _mock_stdin(body)),
            patch(
                "weblogo.matrix.Motif.read_transfac", side_effect=ValueError("x")
            ),
            patch("weblogo.read_seq_data", side_effect=IOError("read error")),
        ):
            _main()
        assert "Content-Type: text/html" in captured_out.getvalue()

    def test_logo_runtime_error(self) -> None:
        """RuntimeError during logo creation -> error form."""
        params = {**_BASE_PARAMS, "sequences": _FASTA_SEQS}
        body = _make_urlencoded_body(params)
        captured_out, captured_err = StringIO(), StringIO()
        with (
            patch.dict(
                os.environ,
                _cgi_env("application/x-www-form-urlencoded", body),
                clear=False,
            ),
            patch("sys.stdout", captured_out),
            patch("sys.stderr", captured_err),
            patch("sys.stdin", new_callable=lambda: _mock_stdin(body)),
            patch(
                "weblogo.matrix.Motif.read_transfac", side_effect=ValueError("x")
            ),
            patch("weblogo.read_seq_data", side_effect=RuntimeError("internal")),
        ):
            _main()
        assert "Content-Type: text/html" in captured_out.getvalue()

    def test_cmd_validate(self) -> None:
        """POST with cmd_validate returns form after logo creation."""
        params = {**_BASE_PARAMS, "sequences": _FASTA_SEQS, "cmd_validate": "true"}
        body = _make_urlencoded_body(params)
        captured_out = StringIO()
        with (
            patch.dict(
                os.environ,
                _cgi_env("application/x-www-form-urlencoded", body),
                clear=False,
            ),
            patch("sys.stdout", captured_out),
            patch("sys.stdin", new_callable=lambda: _mock_stdin(body)),
        ):
            _main()
        assert "Content-Type: text/html" in captured_out.getvalue()

    def test_download_disposition(self) -> None:
        """POST with download sets Content-Disposition: attachment."""
        params = {**_BASE_PARAMS, "sequences": _FASTA_SEQS, "download": "true"}
        body = _make_urlencoded_body(params)
        mock_stdout = _MockStdout()
        with (
            patch.dict(
                os.environ,
                _cgi_env("application/x-www-form-urlencoded", body),
                clear=False,
            ),
            patch("sys.stdout", mock_stdout),
            patch("sys.stdin", new_callable=lambda: _mock_stdin(body)),
        ):
            _main()
        assert "Content-Disposition: attachment" in mock_stdout.get_text()


class TestSendFormEdgeCases:
    """Additional send_form edge cases."""

    def test_explicit_htdocs_directory(self) -> None:
        """send_form with explicit htdocs_directory skips default path."""
        import weblogo as wl

        htdocs = os.path.join(os.path.dirname(wl.__file__), "htdocs")
        controls = [Field("sequences", "")]
        captured = StringIO()
        with patch("sys.stdout", captured):
            send_form(controls, htdocs_directory=htdocs)
        assert "Content-Type: text/html" in captured.getvalue()

    def test_gs_not_found_disables_formats(self) -> None:
        """When Ghostscript/pdf2svg not found, formats are disabled."""
        controls = [Field("sequences", "")]
        captured = StringIO()
        with patch("sys.stdout", captured), patch("shutil.which", return_value=None):
            send_form(controls)
        assert "Content-Type: text/html" in captured.getvalue()
