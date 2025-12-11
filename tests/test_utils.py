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

    def test_julian_date_calculation(self):
        """Test Julian date calculation."""
        # J2000 epoch: January 1, 2000, 12:00:00 TT
        dt = datetime(2000, 1, 1, 12, 0, 0)
        jd = utils.julian_date(dt)

        # J2000 epoch JD is 2451545.0
        assert abs(jd - 2451545.0) < 1.0  # Within 1 day

    def test_gmst_calculation(self):
        """Test Greenwich Mean Sidereal Time calculation."""
        dt = datetime(2000, 1, 1, 0, 0, 0)
        gmst_rad = utils.gmst(dt)

        # GMST should be between 0 and 2*pi
        assert 0 <= gmst_rad < 2 * math.pi

    def test_orbital_period_iss(self):
        """Test orbital period calculation for ISS-like orbit."""
        # ISS semi-major axis
        a = 6796.0  # km

        period = utils.orbital_period(a)

        # ISS period is about 90-95 minutes
        period_minutes = period / 60
        assert 88 < period_minutes < 96

    def test_orbital_period_geo(self):
        """Test orbital period for geostationary orbit."""
        # GEO altitude is about 35,786 km
        altitude = 35786.0
        a = utils.altitude_to_semi_major_axis(altitude)

        period = utils.orbital_period(a)

        # GEO period should be 24 hours
        period_hours = period / 3600
        assert abs(period_hours - 24.0) < 0.5  # Within 30 minutes

    def test_altitude_conversion_round_trip(self):
        """Test altitude to semi-major axis and back."""
        altitude = 500.0  # km

        a = utils.altitude_to_semi_major_axis(altitude)
        altitude_back = utils.semi_major_axis_to_altitude(a)

        assert abs(altitude - altitude_back) < 1e-6

    def test_escape_velocity_earth_surface(self):
        """Test escape velocity at Earth's surface."""
        v_esc = utils.escape_velocity(constants.EARTH_RADIUS)

        # Earth's escape velocity is about 11.2 km/s
        assert abs(v_esc - 11.2) < 0.5

    def test_circular_velocity_leo(self):
        """Test circular velocity for LEO."""
        # 400 km altitude
        r = constants.EARTH_RADIUS + 400

        v_circ = utils.circular_velocity(r)

        # LEO circular velocity is about 7.7 km/s
        assert abs(v_circ - 7.7) < 0.3

    def test_hohmann_transfer(self):
        """Test Hohmann transfer calculation."""
        # LEO to GEO transfer
        r1 = constants.EARTH_RADIUS + 400   # LEO
        r2 = constants.EARTH_RADIUS + 35786  # GEO

        transfer = utils.hohmann_transfer(r1, r2)

        # Check return dictionary has required keys
        assert 'delta_v1' in transfer
        assert 'delta_v2' in transfer
        assert 'transfer_time' in transfer
        assert 'total_delta_v' in transfer

        # Total delta-v for LEO to GEO is about 3.9-4.0 km/s
        assert 3.5 < transfer['total_delta_v'] < 4.5

    def test_sun_synchronous_inclination(self):
        """Test sun-synchronous inclination calculation."""
        altitude = 800.0  # km
        i = utils.sun_synchronous_inclination(altitude)

        # Sun-synchronous orbits are typically 97-99 degrees
        i_deg = i * constants.RAD_TO_DEG
        assert 95 < i_deg < 100

    def test_atmospheric_density_sea_level(self):
        """Test atmospheric density at sea level."""
        rho = utils.atmospheric_density(0.0)

        # Sea level density is about 1.225 kg/m³
        assert 1.0 < rho < 1.5

    def test_atmospheric_density_decreases_with_altitude(self):
        """Test that density decreases with altitude."""
        rho_0 = utils.atmospheric_density(0.0)
        rho_100 = utils.atmospheric_density(100.0)
        rho_500 = utils.atmospheric_density(500.0)

        assert rho_0 > rho_100 > rho_500

    def test_visibility_angles(self):
        """Test satellite visibility angle calculation."""
        result = utils.visibility_angles(observer_altitude=0.0, satellite_altitude=500.0)

        assert 'max_range_km' in result
        assert 'horizon_angle_deg' in result
        assert result['max_range_km'] > 0
        assert 0 < result['horizon_angle_deg'] < 90

    def test_groundtrack_shift(self):
        """Test ground track shift calculation."""
        # ISS period
        period = 5560.0  # seconds

        shift = utils.groundtrack_shift_per_orbit(period)

        # Earth rotates about 22-23 degrees during one ISS orbit
        assert 20 < shift < 25

    def test_repeating_ground_track(self):
        """Test repeating ground track calculation."""
        revs = 15
        days = 1  # 15 revolutions per day

        semi_major_axis, period = utils.repeating_ground_track(revs, days)

        # Check that it gives reasonable values for LEO
        assert 6500 < semi_major_axis < 7500  # Reasonable LEO range
        assert period > 0

        # Check that it actually repeats
        orbits_per_day = 86400 / period
        assert abs(orbits_per_day - revs) < 0.1

    def test_eclipse_duration(self):
        """Test eclipse duration calculation."""
        a = 6878.0  # LEO
        e = 0.0

        duration = utils.eclipse_duration(a, e)

        # LEO eclipse is typically 30-40 minutes
        duration_minutes = duration / 60
        assert 20 < duration_minutes < 50

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

    def test_hohmann_same_orbit_zero_delta_v(self):
        """Test that Hohmann transfer to same orbit gives zero delta-v."""
        r = 7000.0

        transfer = utils.hohmann_transfer(r, r)

        # Delta-v should be effectively zero
        assert transfer['total_delta_v'] < 1e-6
        assert transfer['delta_v1'] < 1e-6
        assert transfer['delta_v2'] < 1e-6