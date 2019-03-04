
import shutil
from subprocess import Popen, PIPE

import pytest

# import weblogo

from . import data_stream


def _exec(args, outputtext, returncode=0, stdin=None):
    if not stdin:
        stdin = data_stream("cap.fa")
    args = ["weblogo"] + args
    p = Popen(args, stdin=stdin, stdout=PIPE, stderr=PIPE)
    (out, err) = p.communicate()
    if returncode == 0 and p.returncode > 0:
        print(err)
    assert returncode == p.returncode
    if returncode == 0:
        assert len(err) == 0

    if outputtext:
        out = out.decode()
        for item in outputtext:
            assert item in out

    stdin.close()


def test_malformed_options():
    _exec(["--notarealoption"], [], 2)
    _exec(["extrajunk"], [], 2)
    _exec(["-I"], [], 2)


def test_help_option():
    _exec(["-h"], ["options"])
    _exec(["--help"], ["options"])


# def test_version_option():
#     _exec(['--version'], weblogo.__version__[0:5])


def test_default_build():
    _exec([], ["%%Title:        Sequence Logo:"])


# Format options
def test_width():
    _exec(['-W', '1234'], ["/stack_width         1234"])
    _exec(['--stack-width', '1234'], ["/stack_width         1234"])


def test_height():
    _exec(['-W', '1000'], ["/stack_height        5000"])
    _exec(['-W', '1000', '--aspect-ratio', '2'], ["/stack_height        2000"])


def test_stacks_per_line():
    _exec(['-n', '7'], ["/stacks_per_line     7 def"])
    _exec(['--stacks-per-line', '7'], ["/stacks_per_line     7 def"])


def test_title():
    _exec(['-t', '3456'], ['/logo_title         (3456) def', '/show_title         True def'])
    _exec(['-t', ''], ['/logo_title         () def', '/show_title         False def'])
    _exec(['--title', '3456'], ['/logo_title         (3456) def', '/show_title         True def'])


def test_annotate():
    _exec(["--annotate", "1,2,3,4"], [], 2)
    _exec(["--annotate", "1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,,"], [])


def test_color():
    _exec(['--color', 'black', 'AG', 'Purine'], [])
    _exec(['--color', 'not_a_color', 'AG', 'Purine'], [], 2)


def test_reverse_complement():
    _exec(['--complement'], [])
    _exec(['--reverse'], [])
    _exec(['--revcomp'], [])


def test_formats():
    _exec(['--format', 'eps'], [])
    _exec(['--format', 'png'], [])
    _exec(['--format', 'png_print'], [])
    _exec(['--format', 'pdf'], [])
    _exec(['--format', 'jpeg'], [])

    _exec(['--format', 'logodata'], [])


@pytest.mark.skipif(shutil.which('pdf2svg') is None,
                    reason="requires pdf2svg")
def test_formats_svg():
    _exec(['--format', 'svg'], [])
