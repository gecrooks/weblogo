#!/usr/bin/env python

from weblogo._ext.nexus import Nexus

from . import data_stream


def test_create() -> None:
    n = Nexus()
    assert n is not None


def test_parse_f0() -> None:
    f = data_stream("nexus/test_Nexus_input.nex")
    n = Nexus(f)
    # self.output_basics(n)

    expected = [
        "t1",
        "t2 the name",
        "isn'that [a] strange name?",
        "one should be punished, for (that)!",
        "t5",
        "t6",
        "t7",
        "t8",
        "t9",
    ]
    taxa = n.taxlabels
    assert taxa == expected
    f.close()


def test_parse_protein() -> None:
    f = data_stream("nexus/protein.nex")
    Nexus(f)
    f.close()


def test_parse_dna() -> None:
    f = data_stream("nexus/dna.nex")
    n = Nexus(f)

    taxa = n.taxlabels
    taxa.sort()
    assert len(taxa) == 10
    assert taxa[0] == "Carp"
    assert taxa[-1] == "Whale"
    f.close()


def test_TreeTest1() -> None:
    """Test Tree module."""
    f = data_stream("nexus/test_Nexus_input.nex")
    n = Nexus(f)
    t3 = n.trees[2]
    n.trees[2]
    t3.root_with_outgroup(["t1", "t5"])

    # Return node_id of common ancestor if
    # taxon_list is monophyletic, -1 otherwise.
    assert t3.is_monophyletic(["t1", "t5"]) == 13

    t3.split(parent_id=t3.search_taxon("t9"))
    f.close()
