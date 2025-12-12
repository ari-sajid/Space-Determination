"""
Main entry point for the Orbit Propagation & Prediction Tool.
APC 524 Project - Now with live TLE data support!
"""

import argparse
import sys
from datetime import datetime, timedelta
from orbital_elements import OrbitalElements
from input_parser import InputParser
from propagator import TwoBodyPropagator
from output_generator import OutputGenerator
from data_fetcher import LiveDataFetcher


def list_satellites():
    """Display available satellites."""
    fetcher = LiveDataFetcher()
    satellites = fetcher.list_available_satellites()

    # Format it to be beautiful
    print("\n=== Available Satellites ===\n")
    print(f"{'Key':<15} {'Full Name':<40}")
    print("-" * 55)
    for key, name in sorted(satellites.items()):
        print(f"{key:<15} {name:<40}")
    print("\nYou can also use any NORAD ID directly.")
    print("Example: python main.py --live iss --time 7200")
    sys.exit(0)


def main():
    """Main function to run the orbit propagation tool."""

    # Set up command line arguments and make it as easy to use as possible
    parser = argparse.ArgumentParser(
        description='Orbit Propagation & Prediction Tool - Now with Live TLE Data!',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Use live data for ISS, propagate for 2 hours
  python main.py --live iss --time 7200

  # Use live data for Hubble, propagate for 1 day
  python main.py --live hubble --time 86400 --output hubble_orbit.csv --format csv

  # Use a local TLE file
  python main.py tle_data.txt --time 3600

  # List all available satellites
  python main.py --list-satellites
        """
    )

    # Define arguments and help messages to guide users
    parser.add_argument(
        'input_file',
        nargs='?',
        help='Path to input file (TLE format or observation data). Not required with --live flag.'
    )
    parser.add_argument(
        '--live',
        type=str,
        metavar='SATELLITE',
        help='Fetch live TLE data for satellite (e.g., iss, hubble, tiangong, or NORAD ID)'
    )
    parser.add_argument(
        '--list-satellites',
        action='store_true',
        help='List all available satellites and exit'
    )
    parser.add_argument(
        '--time',
        type=float,
        default=3600,
        metavar='SECONDS',
        help='Propagation time in seconds (default: 3600). Use 86400 for 1 day, 604800 for 1 week.'
    )
    parser.add_argument(
        '--output',
        default='orbit_output.txt',
        metavar='FILE',
        help='Output file name (default: orbit_output.txt)'
    )
    parser.add_argument(
        '--format',
        choices=['txt', 'csv', 'json'],
        default='txt',
        help='Output format (default: txt)'
    )
    parser.add_argument(
        '--source',
        choices=['auto', 'orbiting-now', 'celestrak'],
        default='auto',
        help='Data source for live TLE (default: auto - tries both)'
    )

    args = parser.parse_args()

    # Handle list satellites request
    if args.list_satellites:
        list_satellites()

    # Validate arguments
    if not args.live and not args.input_file:
        parser.error("Either 'input_file' or '--live SATELLITE' must be provided")

    try:
        # Step 1: Get orbital elements (either from live data or file)
        if args.live:
            print(f"\n{'='*60}")
            print(f"Fetching live TLE data for: {args.live}")
            print(f"{'='*60}\n")

            fetcher = LiveDataFetcher()
            orbital_elements = fetcher.get_orbital_elements(
                satellite_name=args.live,
                source=args.source
            )

            if not orbital_elements:
                print(f"\nError: Could not fetch live data for '{args.live}'")
                print("Try using --list-satellites to see available options.")
                sys.exit(1)

            satellite_name = args.live.upper()
        else:
            print(f"\n{'='*60}")
            print(f"Reading input from: {args.input_file}")
            print(f"{'='*60}\n")

            input_parser = InputParser()
            orbital_elements = input_parser.parse_file(args.input_file)
            satellite_name = "SATELLITE"

        print("✓ Orbital elements loaded successfully!\n")

        # Display orbital information
        print(f"{'='*60}")
        print(f"Initial Orbital Elements")
        print(f"{'='*60}")
        print(f"Semi-major axis:    {orbital_elements.a:.3f} km")
        print(f"Eccentricity:       {orbital_elements.e:.6f}")
        print(f"Inclination:        {orbital_elements.i * 57.2958:.3f}°")
        print(f"RAAN:               {orbital_elements.raan * 57.2958:.3f}°")
        print(f"Arg of Perigee:     {orbital_elements.arg_perigee * 57.2958:.3f}°")
        print(f"Mean Anomaly:       {orbital_elements.mean_anomaly * 57.2958:.3f}°")
        print(f"Orbital Period:     {orbital_elements.period/60:.2f} minutes")
        altitude = orbital_elements.a - 6378.137
        print(f"Average Altitude:   {altitude:.1f} km")
        print(f"{'='*60}\n")

        # Step 2: Initialize propagator
        print("Initializing orbital propagator...")
        propagator = TwoBodyPropagator(orbital_elements)

        # Step 3: Propagate orbit
        hours = args.time / 3600
        print(f"Propagating orbit for {args.time:.0f} seconds ({hours:.2f} hours)...\n")
        final_state = propagator.propagate(args.time)

        # Step 4: Generate output
        print(f"Generating output report...")
        output_gen = OutputGenerator()
        output_gen.generate_report(
            initial_elements=orbital_elements,
            final_state=final_state,
            propagation_time=args.time,
            output_file=args.output,
            format=args.format
        )

        print(f"\n{'='*60}")
        print(f"SUCCESS! Results saved to: {args.output}")
        print(f"{'='*60}\n")

        # Display summary of final state
        print("Final State Summary:")
        print(f"  Position magnitude: {sum(x**2 for x in final_state.position)**0.5:.3f} km")
        print(f"  Velocity magnitude: {sum(x**2 for x in final_state.velocity)**0.5:.3f} km/s")
        final_altitude = final_state.orbital_elements.a - 6378.137
        print(f"  Final altitude:     {final_altitude:.1f} km")
        print()

    except FileNotFoundError:
        print(f"\nError: Input file '{args.input_file}' not found.")
        sys.exit(1)
    except ValueError as e:
        print(f"\nError parsing input: {e}")
        sys.exit(1)
    except KeyboardInterrupt:
        print("\n\nOperation cancelled by user.")
        sys.exit(0)
    except Exception as e:
        print(f"\nUnexpected error: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)


if __name__ == "__main__":
    main()