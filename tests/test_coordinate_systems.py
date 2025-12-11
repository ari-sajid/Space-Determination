"""
Tests for coordinate_systems module.
"""

import pytest
import math
import numpy as np
from datetime import datetime
from coordinate_systems import CoordinateTransform
import constants


class TestCoordinateTransforms:
    """Test cases for CoordinateTransforms class."""

    def test_rotation_matrix_x(self):
        """Test X-axis rotation matrix."""
        coord = CoordinateTransform()
        angle = math.pi / 4  # 45 degrees

        R = coord.R_x(angle)

        # Check it's a proper rotation matrix
        assert R.shape == (3, 3)
        # Determinant should be 1
        assert abs(np.linalg.det(R) - 1.0) < 1e-10
        # Should be orthogonal: R * R^T = I
        I = np.dot(R, R.T)
        np.testing.assert_array_almost_equal(I, np.eye(3), decimal=10)

    def test_rotation_matrix_y(self):
        """Test Y-axis rotation matrix."""
        coord = CoordinateTransform()
        angle = math.pi / 3  # 60 degrees

        R = coord.R_y(angle)

        assert R.shape == (3, 3)
        assert abs(np.linalg.det(R) - 1.0) < 1e-10
        I = np.dot(R, R.T)
        np.testing.assert_array_almost_equal(I, np.eye(3), decimal=10)

    def test_rotation_matrix_z(self):
        """Test Z-axis rotation matrix."""
        coord = CoordinateTransform()
        angle = math.pi / 6  # 30 degrees

        R = coord.R_z(angle)

        assert R.shape == (3, 3)
        assert abs(np.linalg.det(R) - 1.0) < 1e-10
        I = np.dot(R, R.T)
        np.testing.assert_array_almost_equal(I, np.eye(3), decimal=10)

    def test_zero_rotation(self):
        """Test that zero rotation gives identity matrix."""
        coord = CoordinateTransform()

        Rx = coord.R_x(0.0)
        Ry = coord.R_y(0.0)
        Rz = coord.R_z(0.0)

        np.testing.assert_array_almost_equal(Rx, np.eye(3), decimal=10)
        np.testing.assert_array_almost_equal(Ry, np.eye(3), decimal=10)
        np.testing.assert_array_almost_equal(Rz, np.eye(3), decimal=10)

    def test_full_rotation(self):
        """Test that 2*pi rotation returns to original."""
        coord = CoordinateTransform()

        R = coord.R_z(2 * math.pi)
        np.testing.assert_array_almost_equal(R, np.eye(3), decimal=10)

    def test_perifocal_to_eci(self):
        """Test perifocal to ECI transformation."""
        coord = CoordinateTransform()

        # Simple test case: zero angles
        pos_pf = np.array([7000.0, 0.0, 0.0])
        vel_pf = np.array([0.0, 7.5, 0.0])

        pos_eci, vel_eci = coord.perifocal_to_eci(
            pos_pf, vel_pf,
            inclination=0.0,
            raan=0.0,
            arg_perigee=0.0
        )

        # With zero angles, should be identity transformation
        np.testing.assert_array_almost_equal(pos_eci, pos_pf, decimal=10)
        np.testing.assert_array_almost_equal(vel_eci, vel_pf, decimal=10)

    def test_eci_to_ecef_at_epoch(self):
        """Test ECI to ECEF transformation."""
        coord = CoordinateTransform()

        # Position on X-axis in ECI
        pos_eci = np.array([7000.0, 0.0, 0.0])
        epoch = datetime(2000, 1, 1, 12, 0, 0)

        pos_ecef = coord.eci_to_ecef(pos_eci, epoch)

        # ECEF position should have same magnitude
        mag_eci = np.linalg.norm(pos_eci)
        mag_ecef = np.linalg.norm(pos_ecef)
        assert abs(mag_eci - mag_ecef) < 1e-6

    def test_ecef_to_geodetic_equator(self):
        """Test ECEF to geodetic at equator."""
        coord = CoordinateTransform()

        # Point on equator at Earth's surface
        pos_ecef = np.array([constants.EARTH_RADIUS, 0.0, 0.0])

        lat, lon, alt = coord.ecef_to_geodetic(pos_ecef)

        # Should be at equator (lat=0), prime meridian (lon=0), zero altitude
        assert abs(lat) < 1e-6  # Near zero latitude
        assert abs(lon) < 1e-6  # Near zero longitude
        assert abs(alt) < 10.0  # Near zero altitude (within 10 km due to approximations)

    def test_ecef_to_geodetic_north_pole(self):
        """Test ECEF to geodetic at north pole."""
        coord = CoordinateTransform()

        # Point at north pole
        pos_ecef = np.array([0.0, 0.0, constants.EARTH_RADIUS])

        lat, lon, alt = coord.ecef_to_geodetic(pos_ecef)

        # Should be at north pole (lat=90 degrees)
        expected_lat = 90.0 * constants.DEG_TO_RAD
        assert abs(lat - expected_lat) < 0.01  # Within ~0.6 degrees

    def test_ecef_to_geodetic_altitude(self):
        """Test that altitude is calculated correctly."""
        coord = CoordinateTransform()

        # Point 1000 km above equator
        altitude = 1000.0
        pos_ecef = np.array([constants.EARTH_RADIUS + altitude, 0.0, 0.0])

        lat, lon, alt = coord.ecef_to_geodetic(pos_ecef)

        # Altitude should be approximately 1000 km
        assert abs(alt - altitude) < 50.0  # Within 50 km

    def test_round_trip_transformation(self):
        """Test that transforming back and forth preserves vectors."""
        coord = CoordinateTransform()

        # Original perifocal vectors
        pos_pf_original = np.array([8000.0, 1000.0, 0.0])
        vel_pf_original = np.array([0.5, 7.0, 0.0])

        # Some orbital angles
        i = 0.5
        raan = 1.0
        arg_p = 0.3

        # Transform to ECI
        pos_eci, vel_eci = coord.perifocal_to_eci(
            pos_pf_original, vel_pf_original, i, raan, arg_p
        )

        # Magnitudes should be preserved
        mag_pf = np.linalg.norm(pos_pf_original)
        mag_eci = np.linalg.norm(pos_eci)
        assert abs(mag_pf - mag_eci) < 1e-6

    def test_position_magnitude_invariant(self):
        """Test that coordinate transformations preserve vector magnitude."""
        coord = CoordinateTransform()

        pos_pf = np.array([7000.0, 500.0, 100.0])
        vel_pf = np.array([1.0, 7.0, 0.5])

        pos_eci, vel_eci = coord.perifocal_to_eci(
            pos_pf, vel_pf,
            inclination=0.7,
            raan=1.2,
            arg_perigee=0.5
        )

        # Check magnitudes are preserved
        assert abs(np.linalg.norm(pos_pf) - np.linalg.norm(pos_eci)) < 1e-9
        assert abs(np.linalg.norm(vel_pf) - np.linalg.norm(vel_eci)) < 1e-9