
from subprocess import Popen, PIPE

import weblogo._transformseq

from . import data_stream


def _exec(args, outputtext, returncode=0, stdin=None):
    if not stdin:
        stdin = data_stream("cap.fa")
    args = ["transformseq"] + args
    p = Popen(args, stdin=stdin, stdout=PIPE, stderr=PIPE)
    (out, err) = p.communicate()
    if returncode == 0 and p.returncode > 0:
        print(err)
    assert returncode == p.returncode
    if returncode == 0:
        assert len(err) == 0

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


def test_version_option():
    _exec(['--version'], weblogo._transformseq.__version__)


def test_clustal():
    _exec(['-F', 'clustal'], ['TCTTGTGATGTGGTTAACCAAT'])


def test_reverse():
    _exec(['--reverse'], ['TAACCAATTGGTGTAGTGTTCT'])


def test_complement():
    _exec(['--complement'], ['AGAACACTACACCAATTGGTTA'])


def test_seg():
    _exec(['--seg'], ['XXXXXXXXXXXXXXXXXXXXXX'])


def test_subsample():
    _exec(['--subsample', '0.4'], [])
