"""
Example: Basic orbit propagation from custom orbital elements.

This example demonstrates:
- Creating orbital elements manually
- Initializing the propagator
- Propagating an orbit forward in time
- Analyzing the orbital state

Run this after installing the package:
    python examples/example_propagation.py
"""

import sys
import os
from datetime import datetime

# Add parent directory to path for imports
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from orbital_elements import OrbitalElements
from propagator import TwoBodyPropagator
import constants
import numpy as np


def main():
    """Demonstrate basic orbit propagation."""

    print("=" * 70)
    print("Basic Orbit Propagation Example")
    print("=" * 70)
    print()

    # Step 1: Define orbital elements for a sample Low Earth Orbit satellite
    print("Creating orbital elements for a sample LEO satellite...")
    print("(Similar to a typical Earth observation satellite)")
    print()

    # Define orbital parameters
    a = 6878.137  # Semi-major axis (km) - about 500 km altitude
    e = 0.001     # Eccentricity (nearly circular)
    i = 98.0 * constants.DEG_TO_RAD  # Inclination (radians) - sun-synchronous
    raan = 10.0 * constants.DEG_TO_RAD  # RAAN (radians)
    arg_perigee = 0.0  # Argument of perigee (radians)
    mean_anomaly = 0.0  # Mean anomaly (radians)
    epoch = datetime.now()

    orbital_elements = OrbitalElements(
        a=a,
        e=e,
        i=i,
        raan=raan,
        arg_perigee=arg_perigee,
        mean_anomaly=mean_anomaly,
        epoch=epoch
    )

    print("-" * 70)
    print("Initial Orbital Elements")
    print("-" * 70)
    print(f"Semi-major axis:    {orbital_elements.a:.3f} km")
    print(f"Eccentricity:       {orbital_elements.e:.6f}")
    print(f"Inclination:        {orbital_elements.i * constants.RAD_TO_DEG:.3f} degrees")
    print(f"RAAN:               {orbital_elements.raan * constants.RAD_TO_DEG:.3f} degrees")
    print(f"Arg of Perigee:     {orbital_elements.arg_perigee * constants.RAD_TO_DEG:.3f} degrees")
    print(f"Mean Anomaly:       {orbital_elements.mean_anomaly * constants.RAD_TO_DEG:.3f} degrees")
    print(f"Orbital Period:     {orbital_elements.period / 60:.2f} minutes")
    altitude = orbital_elements.a - constants.EARTH_RADIUS
    print(f"Average Altitude:   {altitude:.1f} km")
    print(f"Orbit Type:         {orbital_elements.orbit_type}")
    print("-" * 70)
    print()

    # Step 2: Initialize propagator
    print("Initializing two-body propagator...")
    propagator = TwoBodyPropagator(orbital_elements)
    print()

    # Step 3: Propagate for multiple orbital periods
    print("Propagating orbit for 3 complete orbits...")
    num_orbits = 3
    propagation_time = orbital_elements.period * num_orbits

    print()
    print("-" * 70)
    print("Propagation Progress")
    print("-" * 70)
    print(f"{'Orbit':<8} {'Time (min)':<12} {'Position Mag (km)':<18} {'Velocity Mag (km/s)':<20}")
    print("-" * 70)

    # Propagate and display state at each orbit
    for orbit_num in range(num_orbits + 1):
        t = orbital_elements.period * orbit_num
        state = propagator.propagate(t)

        pos_mag = np.linalg.norm(state.position)
        vel_mag = np.linalg.norm(state.velocity)

        print(f"{orbit_num:<8} {t/60:>11.2f}  {pos_mag:>17.3f}  {vel_mag:>19.6f}")

    print("-" * 70)
    print()

    # Step 4: Detailed final state
    final_state = propagator.propagate(propagation_time)

    print("-" * 70)
    print("Final State after 3 Orbits")
    print("-" * 70)
    print(f"Time:     {final_state.time}")
    print(f"Position: [{final_state.position[0]:.3f}, {final_state.position[1]:.3f}, {final_state.position[2]:.3f}] km")
    print(f"Velocity: [{final_state.velocity[0]:.3f}, {final_state.velocity[1]:.3f}, {final_state.velocity[2]:.3f}] km/s")
    print()
    print("Final Orbital Elements:")
    print(f"  Semi-major axis: {final_state.orbital_elements.a:.3f} km")
    print(f"  Eccentricity:    {final_state.orbital_elements.e:.6f}")
    print(f"  Inclination:     {final_state.orbital_elements.i * constants.RAD_TO_DEG:.3f} degrees")
    print(f"  Mean Anomaly:    {final_state.orbital_elements.mean_anomaly * constants.RAD_TO_DEG:.3f} degrees")
    print("-" * 70)
    print()

    print("Propagation complete!")
    print()
    print("Note: This uses two-body dynamics (no perturbations).")
    print("For real satellites, include J2, drag, and other perturbations.")
    print()

    return 0


if __name__ == "__main__":
    sys.exit(main())