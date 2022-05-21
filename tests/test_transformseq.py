from subprocess import PIPE, Popen

import weblogo._transformseq

from . import data_stream


def _exec(args, outputtext, returncode=0, stdin=None):  # type: ignore
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


def test_malformed_options() -> None:
    _exec(["--notarealoption"], [], 2)
    _exec(["extrajunk"], [], 2)
    _exec(["-I"], [], 2)


def test_help_option() -> None:
    _exec(["-h"], ["options"])
    _exec(["--help"], ["options"])


def test_version_option() -> None:
    _exec(["--version"], weblogo._transformseq.__version__)


def test_clustal() -> None:
    _exec(["-F", "clustal"], ["TCTTGTGATGTGGTTAACCAAT"])


def test_reverse() -> None:
    _exec(["--reverse"], ["TAACCAATTGGTGTAGTGTTCT"])


def test_complement() -> None:
    _exec(["--complement"], ["AGAACACTACACCAATTGGTTA"])


def test_seg() -> None:
    _exec(["--seg"], ["XXXXXXXXXXXXXXXXXXXXXX"])


def test_subsample() -> None:
    _exec(["--subsample", "0.4"], [])
