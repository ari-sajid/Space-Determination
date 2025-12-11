# Examples

This directory contains example scripts demonstrating the functionality of our project.

## Running Examples

After installing the package, you can run these examples directly:

```bash
python examples/example_iss_tracking.py
python examples/example_propagation.py
python examples/example_satellite_comparison.py
```

## Available Examples

### 1. ISS Tracking (example_iss_tracking.py)

Demonstrates real-time tracking of the International Space Station:
- Fetches live TLE data from orbit.ing-now.com
- Propagates the orbit for 24 hours
- Calculates ground track coordinates
- Displays orbital parameters

**What you'll learn:**
- How to fetch live satellite data
- Basic orbit propagation
- Ground track calculation
- Coordinate system transformations

### 2. Basic Propagation (example_propagation.py)

Shows how to create and propagate custom orbital elements:
- Manually defines orbital parameters
- Creates a sun-synchronous LEO orbit
- Propagates for multiple orbital periods
- Analyzes state evolution

**What you'll learn:**
- Creating OrbitalElements objects
- Understanding orbital parameters
- Multi-orbit propagation
- State vector analysis

### 3. Satellite Comparison (example_satellite_comparison.py)

Compares multiple satellites with different orbit types:
- Fetches data for ISS, Hubble, GOES-16, and Terra
- Creates comparison tables
- Classifies orbit types (LEO, MEO, GEO)
- Analyzes orbital characteristics

**What you'll learn:**
- Comparing different satellites
- Understanding orbit classifications
- Analyzing orbital characteristics
- Working with multiple data sources

## Installation

Before running examples, ensure the package is installed:

```bash
# Install in development mode
pip install -e .

# Or install with dev dependencies
pip install -e ".[dev]"
```

## Expected Output

Each example produces formatted text output showing:
- Satellite orbital parameters
- Propagation results
- Analysis and comparisons
- Helpful usage tips

## Creating Your Own Examples

Use these examples as templates for your own analysis:

```python
from data_fetcher import LiveDataFetcher
from propagator import TwoBodyPropagator

# Fetch satellite data
fetcher = LiveDataFetcher()
elements = fetcher.get_orbital_elements('iss')

# Create propagator and propagate
propagator = TwoBodyPagator(elements)
final_state = propagator.propagate(3600)  # 1 hour

# Access results
print(f"Position: {final_state.position}")
print(f"Velocity: {final_state.velocity}")
```

## Troubleshooting

### "ModuleNotFoundError"
Install the package first: `pip install -e .`

### "Could not fetch data"
Check internet connection and try again. The live data fetcher has automatic fallback to Celestrak.

### Import Errors
Make sure you're running from the project root directory.

## Next Steps

After running these examples:
1. Try modifying the propagation times
2. Experiment with different satellites
3. Create your own analysis scripts
4. Export data to files using the main CLI tool

For more information, see the main README.md in the project root.