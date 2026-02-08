import os
import numpy as np
import h5py
import pytest
from dynabench_loader import load_dynabench_data, _parse_file_indices


@pytest.fixture
def mock_dynabench_structure(tmp_path):
    """
    Creates a temporary directory structure mimicking Dynabench.
    Returns the base_path as a string.
    """
    base_path = tmp_path / "data"
    # Structure: data/advection/grid/low/
    target_dir = base_path / "advection" / "grid" / "low"
    target_dir.mkdir(parents=True)
    return str(base_path)


def create_h5_file(base_path, filename, num_sims, include_points=True):
    """Helper to create a dummy HDF5 file in the correct subdirectory."""
    target_dir = os.path.join(base_path, "advection", "grid", "low")
    fpath = os.path.join(target_dir, filename)

    with h5py.File(fpath, "w") as f:
        # Data shape: (N, Time, C, X, Y)
        # Using small dimensions for speed: (N, 2, 1, 5, 5)
        data = np.zeros((num_sims, 2, 1, 5, 5), dtype=np.float32)

        # Fill with index to identify simulation later (sim 0 -> 0.0, sim 1 -> 1.0)
        for i in range(num_sims):
            data[i] = i

        f.create_dataset("data", data=data)

        if include_points:
            # Points shape: (N, X, Y, 2)
            points = np.zeros((num_sims, 5, 5, 2), dtype=np.float32)
            f.create_dataset("points", data=points)

    return fpath


def test_parse_file_indices():
    """Test the filename parsing helper."""
    assert _parse_file_indices("advection_train_grid_low_0_499.h5") == (0, 499)
    assert _parse_file_indices("data_100_200.h5") == (100, 200)
    assert _parse_file_indices("invalid_name.h5") is None
    assert _parse_file_indices("no_indices.h5") is None


def test_load_data_success(mock_dynabench_structure):
    """Test loading a valid simulation index from a file."""
    base_path = mock_dynabench_structure

    # Create file for indices 0-9
    create_h5_file(base_path, "advection_train_grid_low_0_9.h5", num_sims=10)

    # Load index 5
    u, points = load_dynabench_data(base_path=base_path, sim_index=5)

    # Check shape
    assert u.shape == (2, 1, 5, 5)
    # Check data content (we filled with index, so value should be 5.0)
    assert np.all(u == 5)
    assert points.shape == (5, 5, 2)


def test_load_data_multiple_files(mock_dynabench_structure):
    """Test loading from the second file in a sequence."""
    base_path = mock_dynabench_structure

    # File 1: 0-9
    create_h5_file(base_path, "advection_train_grid_low_0_9.h5", num_sims=10)
    # File 2: 10-19
    create_h5_file(base_path, "advection_train_grid_low_10_19.h5", num_sims=10)

    # Load index 12 (should be index 2 in the second file)
    u, points = load_dynabench_data(base_path=base_path, sim_index=12)

    # In create_h5_file, we fill data with local index 0..9
    # So index 12 corresponds to local index 2 in the second file.
    # The data value should be 2.0.
    assert np.all(u == 2)


def test_implicit_points_generation(mock_dynabench_structure):
    """Test that points are generated if missing in HDF5."""
    base_path = mock_dynabench_structure

    # Create file without points
    create_h5_file(base_path, "advection_train_grid_low_0_9.h5", num_sims=10, include_points=False)

    u, points = load_dynabench_data(base_path=base_path, sim_index=0)

    # Check points shape and values
    assert points.shape == (5, 5, 2)
    # Check corners of generated grid [0,1]x[0,1]
    assert np.allclose(points[0, 0], [0.0, 0.0])
    assert np.allclose(points[-1, -1], [1.0, 1.0])


def test_fallback_behavior(mock_dynabench_structure, capsys):
    """Test fallback to first file if index is not found."""
    base_path = mock_dynabench_structure

    create_h5_file(base_path, "advection_train_grid_low_0_9.h5", num_sims=10)

    # Request index 100 (does not exist)
    u, points = load_dynabench_data(base_path=base_path, sim_index=100)

    # Should fallback to file 0, index 0
    assert np.all(u == 0)

    # Check stdout for warning
    captured = capsys.readouterr()
    assert "Warning: Could not locate sim_index=100" in captured.out
    assert "Defaulting to first file" in captured.out
