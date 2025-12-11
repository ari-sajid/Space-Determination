"""
Example: Compare multiple satellites using live TLE data.

This example demonstrates:
- Fetching TLE data for multiple satellites
- Comparing orbital parameters
- Analyzing different orbit types

Run this after installing the package:
    python examples/example_satellite_comparison.py
"""

import sys
import os

# Add parent directory to path for imports
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from data_fetcher import LiveDataFetcher
import constants


def main():
    """Compare orbital parameters of multiple satellites."""

    print("=" * 80)
    print("Multi-Satellite Comparison Example")
    print("=" * 80)
    print()

    # Satellites to compare
    satellites_to_compare = ['iss', 'hubble', 'goes-16', 'terra']

    print("Fetching live TLE data for multiple satellites...")
    print("This may take a moment...")
    print()

    fetcher = LiveDataFetcher()
    satellite_data = {}

    # Fetch data for each satellite
    for sat_name in satellites_to_compare:
        print(f"Fetching {sat_name.upper()}...", end=" ")
        elements = fetcher.get_orbital_elements(sat_name)
        if elements:
            satellite_data[sat_name] = elements
            print("OK")
        else:
            print("FAILED")

    print()

    if not satellite_data:
        print("Error: Could not fetch data for any satellites.")
        return 1

    # Display comparison table
    print("=" * 80)
    print("Orbital Parameters Comparison")
    print("=" * 80)
    print()

    # Header
    print(f"{'Satellite':<15} {'Alt (km)':<12} {'Period (min)':<14} {'Inc (deg)':<12} {'Ecc':<10}")
    print("-" * 80)

    # Data rows
    for sat_name, elements in satellite_data.items():
        altitude = elements.a - constants.EARTH_RADIUS
        period_min = elements.period / 60 if elements.period else float('inf')
        inclination_deg = elements.i * constants.RAD_TO_DEG
        eccentricity = elements.e

        print(f"{sat_name.upper():<15} {altitude:>11.1f} {period_min:>13.2f} {inclination_deg:>11.3f} {eccentricity:>9.6f}")

    print("-" * 80)
    print()

    # Detailed analysis
    print("=" * 80)
    print("Detailed Analysis")
    print("=" * 80)
    print()

    for sat_name, elements in satellite_data.items():
        print(f"--- {sat_name.upper()} ---")
        print(f"  Orbit Type:        {elements.orbit_type}")
        print(f"  Semi-major axis:   {elements.a:.3f} km")
        print(f"  Altitude:          {elements.a - constants.EARTH_RADIUS:.1f} km")
        print(f"  Eccentricity:      {elements.e:.6f}")
        print(f"  Inclination:       {elements.i * constants.RAD_TO_DEG:.3f} degrees")
        print(f"  RAAN:              {elements.raan * constants.RAD_TO_DEG:.3f} degrees")

        if elements.period:
            print(f"  Orbital Period:    {elements.period / 60:.2f} minutes")
            print(f"  Mean Motion:       {elements.mean_motion * constants.RAD_TO_DEG * 60:.6f} deg/min")

            # Orbits per day
            orbits_per_day = 86400 / elements.period
            print(f"  Orbits per day:    {orbits_per_day:.2f}")

        if elements.apogee and elements.perigee:
            print(f"  Apogee altitude:   {elements.apogee - constants.EARTH_RADIUS:.1f} km")
            print(f"  Perigee altitude:  {elements.perigee - constants.EARTH_RADIUS:.1f} km")

        print()

    print("=" * 80)
    print("Orbit Classification")
    print("=" * 80)
    print()

    # Classify orbits by altitude
    for sat_name, elements in satellite_data.items():
        altitude = elements.a - constants.EARTH_RADIUS

        if altitude < 2000:
            orbit_class = "Low Earth Orbit (LEO)"
        elif altitude < 35786:
            orbit_class = "Medium Earth Orbit (MEO)"
        else:
            orbit_class = "Geostationary/High Earth Orbit"

        print(f"{sat_name.upper():<15} : {orbit_class}")

    print()
    print("=" * 80)
    print()

    print("Comparison complete!")
    print()
    print("To propagate any of these satellites, use:")
    print("    python main.py --live <satellite_name> --time <seconds>")
    print()

    return 0


if __name__ == "__main__":
    sys.exit(main())