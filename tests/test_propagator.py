"""
Tests for propagator module.
"""

import pytest
import math
import numpy as np
from datetime import datetime, timedelta
from orbital_elements import OrbitalElements
from propagator import TwoBodyPropagator, PropagationState
import constants

class TestTwoBodyPropagator:
    """Test cases for TwoBodyPropagator class."""

    def test_propagator_initialization(self):
        """Test creating a propagator."""
        elements = OrbitalElements(
            a=7000.0,
            e=0.0,
            i=0.0,
            raan=0.0,
            arg_perigee=0.0,
            mean_anomaly=0.0,
            epoch=datetime.now()
        )

        propagator = TwoBodyPropagator(elements)
        assert propagator.initial_elements == elements

    def test_zero_time_propagation(self):
        """Test that propagating for zero time returns initial state."""
        elements = OrbitalElements(
            a=7000.0,
            e=0.0,
            i=0.0,
            raan=0.0,
            arg_perigee=0.0,
            mean_anomaly=0.0,
            epoch=datetime.now()
        )

        propagator = TwoBodyPropagator(elements)
        state = propagator.propagate(0.0)

        # Position magnitude should equal semi-major axis for circular orbit at MA=0
        pos_mag = np.linalg.norm(state.position)
        assert abs(pos_mag - elements.a) < 1.0  # Within 1 km

    def test_circular_orbit_period(self):
        """Test that circular orbit returns to initial position after one period."""
        elements = OrbitalElements(
            a=7000.0,
            e=0.0,
            i=0.0,
            raan=0.0,
            arg_perigee=0.0,
            mean_anomaly=0.0,
            epoch=datetime.now()
        )

        propagator = TwoBodyPropagator(elements)

        # Propagate for one complete period
        state = propagator.propagate(elements.period)

        # Mean anomaly should be back to ~0 (or 2*pi)
        final_ma = state.orbital_elements.mean_anomaly
        # Normalize to [0, 2*pi]
        final_ma = final_ma % (2 * math.pi)

        assert abs(final_ma) < 0.01 or abs(final_ma - 2*math.pi) < 0.01

    def test_energy_conservation(self):
        """Test that orbital energy is conserved during propagation."""
        elements = OrbitalElements(
            a=7000.0,
            e=0.1,
            i=0.5,
            raan=0.0,
            arg_perigee=0.0,
            mean_anomaly=0.0,
            epoch=datetime.now()
        )

        propagator = TwoBodyPropagator(elements)

        # Calculate initial energy
        state0 = propagator.propagate(0.0)
        r0 = np.linalg.norm(state0.position)
        v0 = np.linalg.norm(state0.velocity)
        energy0 = 0.5 * v0**2 - constants.EARTH_MU / r0

        # Propagate forward
        state1 = propagator.propagate(3600.0)
        r1 = np.linalg.norm(state1.position)
        v1 = np.linalg.norm(state1.velocity)
        energy1 = 0.5 * v1**2 - constants.EARTH_MU / r1

        # Energy should be conserved
        assert abs(energy1 - energy0) < 0.01  # Within 0.01 km²/s²

    def test_semi_major_axis_conservation(self):
        """Test that semi-major axis is conserved for two-body problem."""
        elements = OrbitalElements(
            a=8000.0,
            e=0.2,
            i=0.3,
            raan=0.5,
            arg_perigee=1.0,
            mean_anomaly=0.0,
            epoch=datetime.now()
        )

        propagator = TwoBodyPropagator(elements)

        # Propagate for half an orbit
        state = propagator.propagate(elements.period / 2)

        # Semi-major axis should remain constant
        assert abs(state.orbital_elements.a - elements.a) < 1.0  # Within 1 km

    def test_inclination_conservation(self):
        """Test that inclination is conserved."""
        elements = OrbitalElements(
            a=7000.0,
            e=0.1,
            i=1.0,  # ~57 degrees
            raan=0.0,
            arg_perigee=0.0,
            mean_anomaly=0.0,
            epoch=datetime.now()
        )

        propagator = TwoBodyPropagator(elements)
        state = propagator.propagate(5000.0)

        # Inclination should remain constant
        assert abs(state.orbital_elements.i - elements.i) < 0.01  # Within 0.01 rad

    def test_multiple_propagations(self):
        """Test propagating to same time gives same result."""
        elements = OrbitalElements(
            a=7000.0,
            e=0.0,
            i=0.0,
            raan=0.0,
            arg_perigee=0.0,
            mean_anomaly=0.0,
            epoch=datetime.now()
        )

        propagator = TwoBodyPropagator(elements)

        state1 = propagator.propagate(3600.0)
        state2 = propagator.propagate(3600.0)

        # Results should be identical
        np.testing.assert_array_almost_equal(state1.position, state2.position, decimal=6)
        np.testing.assert_array_almost_equal(state1.velocity, state2.velocity, decimal=6)

    def test_backward_propagation(self):
        """Test propagating backward in time."""
        elements = OrbitalElements(
            a=7000.0,
            e=0.0,
            i=0.0,
            raan=0.0,
            arg_perigee=0.0,
            mean_anomaly=0.0,
            epoch=datetime.now()
        )

        propagator = TwoBodyPropagator(elements)

        # Propagate forward then backward
        state_forward = propagator.propagate(1800.0)
        state_back = propagator.propagate(-1800.0)

        # Forward and backward should be different
        assert not np.allclose(state_forward.position, state_back.position)

    def test_propagation_state_attributes(self):
        """Test that PropagationState has required attributes."""
        elements = OrbitalElements(
            a=7000.0,
            e=0.0,
            i=0.0,
            raan=0.0,
            arg_perigee=0.0,
            mean_anomaly=0.0,
            epoch=datetime.now()
        )

        propagator = TwoBodyPropagator(elements)
        state = propagator.propagate(1000.0)

        assert hasattr(state, 'time')
        assert hasattr(state, 'position')
        assert hasattr(state, 'velocity')
        assert hasattr(state, 'orbital_elements')
        assert isinstance(state.position, np.ndarray)
        assert isinstance(state.velocity, np.ndarray)
        assert len(state.position) == 3
        assert len(state.velocity) == 3

    def test_eccentric_orbit_propagation(self):
        """Test propagation of highly eccentric orbit."""
        elements = OrbitalElements(
            a=20000.0,
            e=0.7,  # Highly eccentric
            i=0.5,
            raan=0.0,
            arg_perigee=0.0,
            mean_anomaly=0.0,
            epoch=datetime.now()
        )

        propagator = TwoBodyPropagator(elements)
        state = propagator.propagate(elements.period / 4)

        # Should successfully propagate without errors
        assert state is not None
        pos_mag = np.linalg.norm(state.position)
        assert pos_mag > 0