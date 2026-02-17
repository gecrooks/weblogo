#!/usr/bin/env python

import pytest

from weblogo.data import amino_acid_composition


def test_data_amino_acid_composition() -> None:
    cl = [amino_acid_composition[k] for k in "ARNDCQEGHILKMFPSTWYV"]
    assert sum(cl) == pytest.approx(1)
