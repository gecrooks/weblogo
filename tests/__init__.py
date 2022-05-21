# TODO: replace with direct calls to pkg_resources

from typing import TextIO

from weblogo.utils import resource_filename, resource_stream, resource_string


def data_string(name: str) -> bytes:
    return resource_string(__name__, "data/" + name, __file__)


def data_stream(name: str) -> TextIO:
    return resource_stream(__name__, "data/" + name, __file__)


def data_filename(name: str) -> str:
    return resource_filename(__name__, "data/" + name, __file__)
