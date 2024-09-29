from io import StringIO
from typing import TextIO

import importlib_resources


def data_string(name: str) -> bytes:
    ref = data_ref(name)
    with ref.open() as f:
        data = f.read()
    return data


def data_stream(name: str) -> TextIO:
    return StringIO(data_string(name))


def data_ref(name: str):
    ref = importlib_resources.files(__name__).joinpath("data").joinpath(name)
    return ref
