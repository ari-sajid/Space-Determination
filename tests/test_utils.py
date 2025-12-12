"""
Tests for utils module.
"""

import pytest
import math
from datetime import datetime
import numpy as np
import utils
import constants

class TestUtils:
    """Test cases for utility functions."""

    def test_altitude_conversion_round_trip(self):
        """Test altitude to semi-major axis and back."""
        altitude = 500.0  # km

        a = utils.altitude_to_semi_major_axis(altitude)
        altitude_back = utils.semi_major_axis_to_altitude(a)

        assert abs(altitude - altitude_back) < 1e-6

    def test_orbital_positions_circular(self):
        """Test orbital position generation for circular orbit."""
        a = 7000.0
        e = 0.0
        i = 0.0
        raan = 0.0
        arg_p = 0.0
        num_points = 100

        positions = utils.orbital_positions(a, e, i, raan, arg_p, num_points)

        assert len(positions) == num_points
        assert positions.shape == (num_points, 3)

        # All positions should have same magnitude for circular orbit
        magnitudes = np.linalg.norm(positions, axis=1)
        assert np.std(magnitudes) < 1.0  # Low standard deviation

    def test_orbital_positions_elliptical(self):
        """Test orbital position generation for elliptical orbit."""
        a = 8000.0
        e = 0.3
        i = 0.5
        raan = 0.0
        arg_p = 0.0
        num_points = 100

        positions = utils.orbital_positions(a, e, i, raan, arg_p, num_points)

        assert len(positions) == num_points

        # For elliptical orbit, max and min radius should differ
        magnitudes = np.linalg.norm(positions, axis=1)
        r_max = np.max(magnitudes)
        r_min = np.min(magnitudes)

        # Difference should be significant for e=0.3
        assert r_max - r_min > 1000  # At least 1000 km difference