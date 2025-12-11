"""
Example: Track the International Space Station using live TLE data.

This example demonstrates:
- Fetching live TLE data from orbit.ing-now.com
- Propagating the ISS orbit forward in time
- Calculating ground tracks
- Displaying orbital parameters

Run this after installing the package:
    python examples/example_iss_tracking.py
"""

import sys
import os

# Add parent directory to path for imports
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from data_fetcher import LiveDataFetcher
from propagator import TwoBodyPropagator
from coordinate_systems import CoordinateTransforms
import numpy as np


def main():
    """Track the ISS for 24 hours and display ground track."""

    print("=" * 70)
    print("ISS Tracking Example - Live TLE Data")
    print("=" * 70)
    print()

    # Step 1: Fetch live TLE data for ISS
    print("Fetching live TLE data for the International Space Station...")
    fetcher = LiveDataFetcher()
    orbital_elements = fetcher.get_orbital_elements('iss')

    if not orbital_elements:
        print("Error: Could not fetch ISS data. Check your internet connection.")
        return 1

    print()
    print("-" * 70)
    print("Current ISS Orbital Parameters")
    print("-" * 70)
    print(f"Semi-major axis:    {orbital_elements.a:.3f} km")
    print(f"Eccentricity:       {orbital_elements.e:.6f}")
    print(f"Inclination:        {orbital_elements.i * 57.2958:.3f} degrees")
    print(f"RAAN:               {orbital_elements.raan * 57.2958:.3f} degrees")
    print(f"Arg of Perigee:     {orbital_elements.arg_perigee * 57.2958:.3f} degrees")
    print(f"Mean Anomaly:       {orbital_elements.mean_anomaly * 57.2958:.3f} degrees")
    print(f"Orbital Period:     {orbital_elements.period / 60:.2f} minutes")
    altitude = orbital_elements.a - 6378.137
    print(f"Average Altitude:   {altitude:.1f} km")
    print("-" * 70)
    print()

    # Step 2: Create propagator
    print("Initializing orbit propagator...")
    propagator = TwoBodyPropagator(orbital_elements)

    # Step 3: Propagate for 24 hours
    propagation_time = 86400  # 24 hours in seconds
    print(f"Propagating orbit for 24 hours...")
    final_state = propagator.propagate(propagation_time)

    print()
    print("-" * 70)
    print("Final State after 24 Hours")
    print("-" * 70)
    print(f"Position: [{final_state.position[0]:.3f}, {final_state.position[1]:.3f}, {final_state.position[2]:.3f}] km")
    print(f"Velocity: [{final_state.velocity[0]:.3f}, {final_state.velocity[1]:.3f}, {final_state.velocity[2]:.3f}] km/s")
    pos_mag = np.linalg.norm(final_state.position)
    vel_mag = np.linalg.norm(final_state.velocity)
    print(f"Position magnitude: {pos_mag:.3f} km")
    print(f"Velocity magnitude: {vel_mag:.3f} km/s")
    print("-" * 70)
    print()

    # Step 4: Calculate ground track
    print("Calculating ground track...")
    num_points = 50
    time_step = propagation_time / num_points

    ground_track = []
    for i in range(num_points):
        t = i * time_step
        state = propagator.propagate(t)

        # Convert ECI to geodetic coordinates
        coord_transform = CoordinateTransforms()
        lat, lon, alt = coord_transform.eci_to_geodetic(
            state.position,
            state.time
        )
        ground_track.append((lat, lon, alt))

    # Display sample ground track points
    print()
    print("-" * 70)
    print("Ground Track Sample (First 10 Points)")
    print("-" * 70)
    print(f"{'Time (min)':<12} {'Latitude (deg)':<18} {'Longitude (deg)':<18} {'Altitude (km)':<15}")
    print("-" * 70)

    for i in range(min(10, len(ground_track))):
        t_minutes = (i * time_step) / 60
        lat, lon, alt = ground_track[i]
        print(f"{t_minutes:>11.2f}  {lat:>17.3f}  {lon:>17.3f}  {alt:>14.3f}")

    print("-" * 70)
    print()
    print("Tracking complete!")
    print()
    print("To save this data to a file, use the main.py script:")
    print("    python main.py --live iss --time 86400 --output iss_24hr.csv --format csv")
    print()

    return 0


if __name__ == "__main__":
    sys.exit(main())