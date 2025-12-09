"""
Main entry point for the Orbit Propagation & Prediction Tool.
APC 524 Project
"""

import argparse
import sys
from datetime import datetime, timedelta
from orbital_elements import OrbitalElements
from input_parser import InputParser
from propagator import TwoBodyPropagator
from output_generator import OutputGenerator


def main():
    """Main function to run the orbit propagation tool."""
    
    # Set up command line arguments
    parser = argparse.ArgumentParser(
        description='Orbit Propagation & Prediction Tool'
    )
    parser.add_argument(
        'input_file',
        help='Path to input file (TLE format or observation data)'
    )
    parser.add_argument(
        '--time',
        type=float,
        default=3600,
        help='Propagation time in seconds (default: 3600)'
    )
    parser.add_argument(
        '--output',
        default='orbit_output.txt',
        help='Output file name (default: orbit_output.txt)'
    )
    parser.add_argument(
        '--format',
        choices=['txt', 'csv', 'json'],
        default='txt',
        help='Output format (default: txt)'
    )
    parser.add_argument(
        '--live',
        type=str,
        help='Fetch live data for satellite (e.g., iss, hubble, or NORAD ID)'
    )
    
    args = parser.parse_args()
    
    try:
        # Step 1: Parse input data
        print(f"Reading input from {args.input_file}...")
        parser = InputParser()
        orbital_elements = parser.parse_file(args.input_file)
        print("Input parsed successfully!")
        
        # Step 2: Initialize propagator
        print("Initializing orbital propagator...")
        propagator = TwoBodyPropagator(orbital_elements)
        
        # Step 3: Propagate orbit
        print(f"Propagating orbit for {args.time} seconds...")
        final_state = propagator.propagate(args.time)
        
        # Step 4: Generate output
        print(f"Generating output to {args.output}...")
        output_gen = OutputGenerator()
        output_gen.generate_report(
            initial_elements=orbital_elements,
            final_state=final_state,
            propagation_time=args.time,
            output_file=args.output,
            format=args.format
        )
        
        print(f"Success! Results saved to {args.output}")
        
    except FileNotFoundError:
        print(f"Error: Input file '{args.input_file}' not found.")
        sys.exit(1)
    except ValueError as e:
        print(f"Error parsing input: {e}")
        sys.exit(1)
    except Exception as e:
        print(f"Unexpected error: {e}")
        sys.exit(1)


if __name__ == "__main__":
    main()
