"""
Tests for orbital_elements module.
"""

import pytest
import math
from datetime import datetime
from orbital_elements import OrbitalElements
import constants


class TestOrbitalElements:
    """Test cases for OrbitalElements class."""

    def test_circular_orbit_creation(self):
        """Test creating a circular orbit."""
        elements = OrbitalElements(
            a=6878.137,
            e=0.0,
            i=0.0,
            raan=0.0,
            arg_perigee=0.0,
            mean_anomaly=0.0,
            epoch=datetime.now()
        )

        assert elements.a == 6878.137
        assert elements.e == 0.0
        assert elements.orbit_type == "circular"

    def test_elliptical_orbit_creation(self):
        """Test creating an elliptical orbit."""
        elements = OrbitalElements(
            a=10000.0,
            e=0.3,
            i=math.pi / 4,
            raan=0.0,
            arg_perigee=0.0,
            mean_anomaly=0.0,
            epoch=datetime.now()
        )

        assert elements.a == 10000.0
        assert elements.e == 0.3
        assert elements.orbit_type == "elliptical"

    def test_invalid_semi_major_axis(self):
        """Test that negative semi-major axis raises ValueError."""
        with pytest.raises(ValueError):
            OrbitalElements(
                a=-1000.0,
                e=0.0,
                i=0.0,
                raan=0.0,
                arg_perigee=0.0,
                mean_anomaly=0.0,
                epoch=datetime.now()
            )

    def test_invalid_eccentricity(self):
        """Test that negative eccentricity raises ValueError."""
        with pytest.raises(ValueError):
            OrbitalElements(
                a=7000.0,
                e=-0.1,
                i=0.0,
                raan=0.0,
                arg_perigee=0.0,
                mean_anomaly=0.0,
                epoch=datetime.now()
            )

    def test_invalid_inclination(self):
        """Test that invalid inclination raises ValueError."""
        with pytest.raises(ValueError):
            OrbitalElements(
                a=7000.0,
                e=0.0,
                i=4.0,  # > pi
                raan=0.0,
                arg_perigee=0.0,
                mean_anomaly=0.0,
                epoch=datetime.now()
            )

    def test_period_calculation(self):
        """Test orbital period calculation."""
        elements = OrbitalElements(
            a=6878.137,
            e=0.0,
            i=0.0,
            raan=0.0,
            arg_perigee=0.0,
            mean_anomaly=0.0,
            epoch=datetime.now()
        )

        # Expected period for LEO ~500km altitude
        expected_period = 2 * math.pi * math.sqrt(elements.a**3 / constants.EARTH_MU)
        assert abs(elements.period - expected_period) < 1.0  # Within 1 second

    def test_mean_motion_calculation(self):
        """Test mean motion calculation."""
        elements = OrbitalElements(
            a=6878.137,
            e=0.0,
            i=0.0,
            raan=0.0,
            arg_perigee=0.0,
            mean_anomaly=0.0,
            epoch=datetime.now()
        )

        expected_n = math.sqrt(constants.EARTH_MU / elements.a**3)
        assert abs(elements.mean_motion - expected_n) < 1e-10

    def test_apogee_perigee_calculation(self):
        """Test apogee and perigee calculation for elliptical orbit."""
        a = 8000.0
        e = 0.2

        elements = OrbitalElements(
            a=a,
            e=e,
            i=0.0,
            raan=0.0,
            arg_perigee=0.0,
            mean_anomaly=0.0,
            epoch=datetime.now()
        )

        expected_apogee = a * (1 + e)
        expected_perigee = a * (1 - e)

        assert abs(elements.apogee - expected_apogee) < 1e-6
        assert abs(elements.perigee - expected_perigee) < 1e-6

    def test_circular_orbit_apogee_perigee(self):
        """Test that circular orbit has equal apogee and perigee."""
        elements = OrbitalElements(
            a=7000.0,
            e=0.0,
            i=0.0,
            raan=0.0,
            arg_perigee=0.0,
            mean_anomaly=0.0,
            epoch=datetime.now()
        )

        assert abs(elements.apogee - elements.perigee) < 1e-6

    def test_to_dict_conversion(self):
        """Test conversion to dictionary."""
        elements = OrbitalElements(
            a=7000.0,
            e=0.1,
            i=0.5,
            raan=1.0,
            arg_perigee=0.5,
            mean_anomaly=0.0,
            epoch=datetime.now()
        )

        data_dict = elements.to_dict()

        assert 'a_km' in data_dict
        assert 'e' in data_dict
        assert 'i_deg' in data_dict
        assert data_dict['a_km'] == 7000.0

    def test_angle_normalization(self):
        """Test that angles are normalized to [-pi, pi]."""
        elements = OrbitalElements(
            a=7000.0,
            e=0.0,
            i=0.0,
            raan=7.0,  # > 2*pi
            arg_perigee=0.0,
            mean_anomaly=0.0,
            epoch=datetime.now()
        )

        # RAAN should be normalized
        assert -math.pi <= elements.raan <= math.pi

    def test_iss_like_orbit(self):
        """Test creating ISS-like orbital elements."""
        # ISS typical parameters
        elements = OrbitalElements(
            a=6796.0,
            e=0.0001,
            i=51.6 * constants.DEG_TO_RAD,
            raan=0.0,
            arg_perigee=0.0,
            mean_anomaly=0.0,
            epoch=datetime.now()
        )

        assert elements.orbit_type in ["circular", "elliptical"]
        period_minutes = elements.period / 60
        assert 90 < period_minutes < 95  # ISS period is ~93 minutes