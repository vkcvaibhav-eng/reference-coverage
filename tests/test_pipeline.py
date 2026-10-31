"""Sanity tests for the coverage deficit calculation."""
import os
import sys

import numpy as np
import pytest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "scripts"))

from deficit import (WEIGHTS, coverage_and_deficit, density_on_grid, pca2,
                     weights_for)


@pytest.fixture
def grid():
    gx, gy = np.meshgrid(np.linspace(-3, 3, 60), np.linspace(-3, 3, 60))
    return gx, gy


@pytest.fixture
def clouds():
    rng = np.random.default_rng(7)
    occurrence = rng.normal(0, 1.0, size=(400, 2))
    # Sequences cover only the left half of the occurrence cloud.
    sequences = rng.normal([-1.0, 0.0], 0.5, size=(120, 2))
    return occurrence, sequences


def test_pca2_shapes_and_variance():
    rng = np.random.default_rng(1)
    x = rng.normal(size=(200, 19))
    scores, var, _ = pca2(x)
    assert scores.shape == (200, 2)
    assert len(var) == 2
    assert 0 < var[1] <= var[0] < 1


def test_pca2_tolerates_constant_column():
    rng = np.random.default_rng(2)
    x = rng.normal(size=(50, 19))
    x[:, 3] = 4.0
    scores, _, _ = pca2(x)
    assert np.isfinite(scores).all()


def test_density_returns_positive_kernel_area(grid, clouds):
    gx, gy = grid
    _, sequences = clouds
    density, area = density_on_grid(sequences, gx, gy)
    assert density.shape == gx.shape
    assert (density >= 0).all()
    assert area > 0


def test_density_degrades_gracefully_on_tiny_input(grid):
    gx, gy = grid
    density, area = density_on_grid([[0.0, 0.0]], gx, gy)
    assert not density.any()
    assert area > 0


def test_coverage_within_bounds(grid, clouds):
    gx, gy = grid
    occurrence, sequences = clouds
    occ, _ = density_on_grid(occurrence, gx, gy)
    seq, area = density_on_grid(sequences, gx, gy)
    coverage, deficit, index = coverage_and_deficit(occ, seq, len(sequences), area)
    assert ((coverage >= 0) & (coverage <= 1)).all()
    assert (deficit >= 0).all()
    assert 0 <= index <= 1


def test_more_sequences_reduce_the_deficit(grid, clouds):
    """The index must fall as the effective library grows."""
    gx, gy = grid
    occurrence, sequences = clouds
    occ, _ = density_on_grid(occurrence, gx, gy)
    seq, area = density_on_grid(sequences, gx, gy)
    indices = [coverage_and_deficit(occ, seq, n, area)[2]
               for n in (10, 100, 1000)]
    assert indices[0] > indices[1] > indices[2]


def test_downweighting_raises_the_deficit(grid, clouds):
    """The central claim: doubting the labels must widen the measured gap."""
    gx, gy = grid
    occurrence, sequences = clouds
    occ, _ = density_on_grid(occurrence, gx, gy)
    seq, area = density_on_grid(sequences, gx, gy)
    naive = coverage_and_deficit(occ, seq, len(sequences), area)[2]
    weighted = coverage_and_deficit(occ, seq, len(sequences) * 0.3, area)[2]
    assert weighted > naive


def test_higher_threshold_raises_the_deficit(grid, clouds):
    gx, gy = grid
    occurrence, sequences = clouds
    occ, _ = density_on_grid(occurrence, gx, gy)
    seq, area = density_on_grid(sequences, gx, gy)
    low = coverage_and_deficit(occ, seq, len(sequences), area, 1.0)[2]
    high = coverage_and_deficit(occ, seq, len(sequences), area, 20.0)[2]
    assert high > low


def test_weights_map_verdicts_and_count_unmatched():
    rows = [{"record_id": "a"}, {"record_id": "b"}, {"record_id": "ghost"}]
    verdicts = {"a": "AGREES", "b": "DISAGREES"}
    w, unmatched = weights_for(rows, verdicts, WEIGHTS, 0.25)
    assert w[0] == 1.0
    assert w[1] == 0.0
    assert w[2] == 0.25
    assert unmatched == 1


def test_unknown_verdict_falls_back_to_the_unmatched_weight():
    rows = [{"record_id": "a"}]
    w, unmatched = weights_for(rows, {"a": "SOMETHING_NEW"}, WEIGHTS, 0.25)
    assert w[0] == 0.25
    assert unmatched == 0


from climate_extract import marker_label
from deficit import assign_roles


@pytest.mark.parametrize("text,expected", [
    ("COI-5P", "COI-5P"),
    ("COI-3P", "COI-3P"),
    ("coi-5p", "COI-5P"),
    ("cytochrome c oxidase subunit I", "COI-5P"),
    ("cytochrome oxidase subunit 1", "COI-5P"),
    ("COXI", "COI-5P"),
    ("COX1", "COI-5P"),
    ("cytochrome b", "CYTB"),
    ("CYTB", "CYTB"),
    ("18S-5P", "18S"),
    ("ITS2", "ITS2"),
    ("", ""),
    (None, ""),
])
def test_marker_label(text, expected):
    assert marker_label(text) == expected


def test_cytochrome_b_is_not_read_as_coi():
    """CYTB must be matched before the cytochrome-oxidase patterns."""
    assert marker_label("cytochrome b") != "COI-5P"


def test_non_barcode_sequences_become_occurrences():
    rows = [
        {"source": "GenBank", "marker": "COI-5P"},
        {"source": "BOLD", "marker": "COI-3P"},
        {"source": "BOLD", "marker": "CYTB"},
        {"source": "BOLD", "marker": ""},
        {"source": "GBIF", "marker": ""},
    ]
    occ, seq, moved = assign_roles(rows)
    assert len(seq) == 2                      # the COI-5P record and the blank-marker one
    assert len(occ) == 3                      # GBIF, plus COI-3P and CYTB moved across
    assert moved == {"COI-3P": 1, "CYTB": 1}


def test_tables_without_a_marker_column_keep_every_sequence_record():
    rows = [{"source": "GenBank"}, {"source": "BOLD"}, {"source": "GBIF"}]
    occ, seq, moved = assign_roles(rows)
    assert len(seq) == 2
    assert len(occ) == 1
    assert moved == {}


from deficit import FALLBACK_SHARED, SHARED_CLUSTER, label_accuracy


def test_shared_cluster_weight_is_one_over_k():
    """A label that fits k anchored species equally is worth 1/k, not a guess."""
    rows = [{"record_id": "two"}, {"record_id": "three"}, {"record_id": "four"}]
    verdicts = {r["record_id"]: SHARED_CLUSTER for r in rows}
    counts = {"two": 2, "three": 3, "four": 4}
    w, _ = weights_for(rows, verdicts, WEIGHTS, 0.25, counts)
    assert w[0] == pytest.approx(0.5)
    assert w[1] == pytest.approx(1 / 3)
    assert w[2] == pytest.approx(0.25)


def test_shared_cluster_falls_back_without_a_count():
    rows = [{"record_id": "unknown"}]
    w, _ = weights_for(rows, {"unknown": SHARED_CLUSTER}, WEIGHTS, 0.25, {})
    assert w[0] == FALLBACK_SHARED


def test_label_accuracy_uses_only_checkable_records():
    verdicts = {
        "a": "AGREES", "b": "AGREES", "c": SHARED_CLUSTER,
        "d": "DISAGREES", "e": "NO_ANCHOR", "f": "UNNAMED_SUGGESTED",
    }
    accuracy, checkable, wrong = label_accuracy(verdicts)
    assert checkable == 4          # NO_ANCHOR and UNNAMED are not checkable
    assert wrong == 1
    assert accuracy == pytest.approx(0.75)


def test_label_accuracy_is_none_without_checkable_records():
    accuracy, checkable, wrong = label_accuracy({"a": "NO_ANCHOR"})
    assert accuracy is None
    assert checkable == 0
