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
