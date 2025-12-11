# Space-Determination: Orbit Propagation & Prediction Tool

An advanced orbital mechanics toolkit for satellite trajectory analysis using live TLE data. Built for APC 524 at Princeton University.

## Features

- **Live TLE Data**: Fetch real-time satellite orbital data from orbit.ing-now.com and Celestrak
- **Multiple Satellites**: Support for ISS, Hubble, Tiangong, and many more popular satellites
- **Accurate Propagation**: Two-body orbital mechanics with Kepler equation solving
- **Flexible Time Ranges**: Propagate orbits from seconds to weeks
- **Multiple Output Formats**: TXT, CSV, and JSON output support
- **Coordinate Transformations**: ECI, ECEF, and geodetic coordinate systems
- **Ground Track Calculation**: Compute satellite ground tracks and passes
- **3D Visualization**: Optional matplotlib-based orbit visualization

## Installation

### Standard Installation

```bash
pip install -e .
```

### Development Installation (includes testing tools)

```bash
pip install -e ".[dev]"
```

**Requirements:**
- Python 3.7+
- numpy >= 1.20.0
- matplotlib >= 3.3.0
- requests >= 2.25.0
- beautifulsoup4 >= 4.9.0

## Quick Start

### 1. List Available Satellites

```bash
python main.py --list-satellites
```

This displays all pre-configured satellites in the catalog.

### 2. Propagate ISS Orbit (Live Data)

```bash
python main.py --live iss --time 7200
```

Fetches live TLE data for the International Space Station and propagates its orbit for 2 hours (7200 seconds).

### 3. Propagate with Custom Output

```bash
python main.py --live hubble --time 86400 --output hubble_orbit.csv --format csv
```

Propagates Hubble Space Telescope orbit for 1 day with CSV output.

### 4. Use Local TLE File

```bash
python main.py satellite_tle.txt --time 3600
```

Reads TLE data from a local file and propagates for 1 hour.

## Detailed Usage

### Command-Line Options

```
python main.py [input_file] [options]

Positional Arguments:
  input_file              Path to TLE file (not required with --live)

Options:
  --live SATELLITE        Fetch live TLE data for satellite
                          (e.g., iss, hubble, tiangong, or NORAD ID)

  --list-satellites       Display all available satellites

  --time SECONDS          Propagation duration in seconds
                          Default: 3600 (1 hour)
                          Examples:
                            - 7200: 2 hours
                            - 86400: 1 day
                            - 604800: 1 week

  --output FILE           Output filename
                          Default: orbit_output.txt

  --format {txt,csv,json} Output format
                          Default: txt

  --source {auto,orbiting-now,celestrak}
                          TLE data source (for --live mode)
                          Default: auto (tries both sources)
```

## Available Satellites

The tool includes built-in support for these satellites:

| Key | Satellite Name | NORAD ID |
|-----|----------------|----------|
| iss | International Space Station | 25544 |
| hubble | Hubble Space Telescope | 20580 |
| tiangong | China Space Station | 48274 |
| starlink-1 | Starlink-1 Satellite | 44713 |
| gps-iir-2 | GPS IIR-2 (NAVSTAR 43) | 24876 |
| noaa-18 | NOAA 18 Weather Satellite | 28654 |
| terra | Terra (EOS AM-1) | 25994 |
| aqua | Aqua (EOS PM-1) | 27424 |
| goes-16 | GOES 16 Weather Satellite | 41866 |
| spaceX-crew | SpaceX Crew Dragon | 45623 |

You can also use any NORAD ID directly:
```bash
python main.py --live 25544 --time 3600
```

## Examples

### Example 1: Track ISS for 12 Hours

```bash
python main.py --live iss --time 43200 --output iss_12hr.json --format json
```

**Output includes:**
- Initial and final orbital elements
- Position and velocity vectors
- Ground track coordinates
- Orbital parameters (period, altitude, etc.)

### Example 2: Compare Two Satellites

```bash
# Track ISS
python main.py --live iss --time 86400 --output iss_24hr.csv --format csv

# Track Hubble
python main.py --live hubble --time 86400 --output hubble_24hr.csv --format csv
```

### Example 3: Use Custom NORAD ID

```bash
# Track any satellite by its NORAD catalog number
python main.py --live 43013 --time 7200 --output starlink.txt
```

### Example 4: Weekly Propagation

```bash
python main.py --live noaa-18 --time 604800 --output noaa_weekly.csv --format csv
```

Propagates NOAA-18 weather satellite for 1 week (604800 seconds).

## TLE File Format

If using local TLE files, they should follow the standard two-line element format:

```
ISS (ZARYA)
1 25544U 98067A   23365.50000000  .00012456  00000-0  22123-3 0  9999
2 25544  51.6439 123.4567   0001234  89.1234 270.9876 15.50103472123456
```

**Line 0:** Satellite name (optional)
**Line 1:** Epoch, ballistic coefficient, etc.
**Line 2:** Orbital elements (inclination, RAAN, eccentricity, etc.)

## Project Structure

```
Space-Determination/
├── main.py                     # Main entry point with CLI
├── data_fetcher.py             # Live TLE data fetching
├── input_parser.py             # TLE file parser
├── orbital_elements.py         # Orbital element data structures
├── kepler_solver.py            # Kepler's equation solver
├── propagator.py               # Two-body orbit propagator
├── coordinate_systems.py       # Coordinate transformations
├── output_generator.py         # Output formatting
├── utils.py                    # Utility functions
├── constants.py                # Physical constants
├── visualization.py            # 3D orbit plotting
├── examples/                   # Example scripts
│   ├── README.md              # Examples documentation
│   ├── example_iss_tracking.py
│   ├── example_propagation.py
│   └── example_satellite_comparison.py
├── tests/                      # Comprehensive test suite (42 tests)
│   ├── test_kepler_solver.py
│   ├── test_orbital_elements.py
│   ├── test_propagator.py
│   ├── test_coordinate_systems.py
│   ├── test_input_parser.py
│   └── test_utils.py
├── .github/workflows/          # CI/CD configuration
│   └── ci.yml
├── pyproject.toml             # Package configuration
├── requirements.txt           # Python dependencies
├── API_REFERENCE.md           # Complete API documentation
└── README.md                  # This file
```

## Examples

The `examples/` directory contains ready-to-run scripts demonstrating key functionality:

### 1. ISS Tracking Example
```bash
python examples/example_iss_tracking.py
```
Fetches live ISS data, propagates orbit for 24 hours, and displays ground track coordinates.

### 2. Basic Propagation Example
```bash
python examples/example_propagation.py
```
Creates custom orbital elements and propagates a sun-synchronous LEO satellite for multiple orbits.

### 3. Satellite Comparison Example
```bash
python examples/example_satellite_comparison.py
```
Compares orbital parameters of multiple satellites (ISS, Hubble, GOES-16, Terra).

See `examples/README.md` for detailed documentation on each example.

## Testing

Run the comprehensive test suite (42 tests covering core functionality):

```bash
# Run all tests
pytest tests/

# Run with coverage report
pytest tests/ --cov=. --cov-report=html

# Run specific test file
pytest tests/test_propagator.py -v
```

### Test Coverage
- **test_orbital_elements.py**: Orbital element validation and calculations (12 tests)
- **test_propagator.py**: Orbit propagation and state tracking (10 tests)
- **test_kepler_solver.py**: Kepler equation solving (2 tests)
- **test_coordinate_systems.py**: Coordinate transformations (12 tests)
- **test_input_parser.py**: TLE parsing and validation (12 tests)
- **test_utils.py**: Utility functions and orbital mechanics (18 tests)

## Continuous Integration

The project includes GitHub Actions CI/CD that automatically:
- Runs tests on Python 3.7, 3.8, 3.9, 3.10, and 3.11
- Checks code formatting with black and ruff
- Builds the package distribution
- Generates coverage reports

Workflow triggers on push to main branch and pull requests.

## Educational Use

This tool was developed for APC 524 and demonstrates:

1. **Orbital Mechanics**: Two-body problem, Kepler's equation, orbital elements
2. **Numerical Methods**: Newton-Raphson iteration, coordinate transformations
3. **Software Engineering**: Modular design, error handling, CLI interfaces
4. **Data Integration**: Live API data fetching, web scraping, multiple data sources

## Technical Details

### Coordinate Systems

The tool supports transformations between:
- **Perifocal**: Orbital plane coordinates
- **ECI (Earth-Centered Inertial)**: Fixed inertial reference frame
- **ECEF (Earth-Centered Earth-Fixed)**: Rotating frame fixed to Earth
- **Geodetic**: Latitude, longitude, altitude

### Propagation Method

Uses **two-body orbital mechanics** assuming:
- Point mass satellite
- Spherical Earth with constant gravitational parameter
- No perturbations (J2, drag, solar radiation pressure)

For more accurate long-term propagation, consider using SGP4/SDP4 models.

### Data Sources

1. **orbit.ing-now.com** (Primary)
   - Real-time satellite tracking
   - Updated TLE data
   - Supports most active satellites

2. **Celestrak** (Fallback)
   - NORAD TLE database
   - Comprehensive satellite catalog
   - Reliable API access

## Troubleshooting

### "Could not fetch live data"

**Solution:** Try specifying the data source:
```bash
python main.py --live iss --source celestrak
```

### "Input file not found"

**Solution:** Ensure the TLE file exists or use `--live` mode:
```bash
python main.py --live iss --time 3600
```

### Import Errors

**Solution:** Install dependencies:
```bash
pip install -r requirements.txt
```

## Authors

Jackson Crocker, David Herrera, Ariyan Sajid

---

## Advanced Usage

### Custom Time Intervals

Calculate specific durations:

```bash
# 30 minutes: 30 * 60 = 1800 seconds
python main.py --live iss --time 1800

# 6 hours: 6 * 3600 = 21600 seconds
python main.py --live iss --time 21600

# 3 days: 3 * 86400 = 259200 seconds
python main.py --live iss --time 259200
```

### Batch Processing

Process multiple satellites:

```bash
#!/bin/bash
for sat in iss hubble tiangong terra aqua; do
    python main.py --live $sat --time 86400 \
        --output "${sat}_orbit.csv" --format csv
done
```

### Output Format Comparison

- **TXT**: Human-readable, detailed text report
- **CSV**: Spreadsheet-compatible, easy to plot
- **JSON**: Machine-readable, structured data for further processing

## Documentation and Resources

### API Reference
- **[API_REFERENCE.md](API_REFERENCE.md)**: Complete API documentation with usage examples for all modules

### External Resources
- [TLE Format Specification](https://celestrak.org/NORAD/documentation/tle-fmt.php)
- [Orbital Mechanics Fundamentals](https://en.wikipedia.org/wiki/Orbital_elements)
- [Two-Body Problem](https://en.wikipedia.org/wiki/Two-body_problem)
- [Kepler's Equation](https://en.wikipedia.org/wiki/Kepler%27s_equation)