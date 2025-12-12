"""
Physical and mathematical constants for calculations.
All units are in SI unless otherwise specified.
"""

import math

# Earth parameters
EARTH_MU = 398600.4418  # km³/s² - Earth's gravitational parameter
EARTH_RADIUS = 6378.137  # km - Earth's equatorial radius
EARTH_J2 = 1.08262668e-3  # J2 perturbation coefficient - not currently used :(
EARTH_ROTATION_RATE = 7.2921159e-5  # rad/s - Earth's rotation rate

# Mathematical constants
TWO_PI = 2 * math.pi
DEG_TO_RAD = math.pi / 180.0  # Convert degrees to radians
RAD_TO_DEG = 180.0 / math.pi  # Convert radians to degrees

# Numerical tolerances
KEPLER_TOLERANCE = 1e-10  # Tolerance for Kepler equation solver
MAX_KEPLER_ITERATIONS = 100  # Maximum iterations for Kepler solver

# Time constants
SECONDS_PER_DAY = 86400.0
MINUTES_PER_DAY = 1440.0

# Orbit type thresholds
CIRCULAR_ECC_THRESHOLD = 1e-4  # Eccentricity below this is considered circular
PARABOLIC_ECC_THRESHOLD = 0.9999  # Eccentricity thresholds for orbit types
HYPERBOLIC_ECC_THRESHOLD = 1.0001

# TLE constants (Two-Line Element format)
TLE_LINE_LENGTH = 69  # Standard TLE line length
TLE_EPOCH_YEAR_CUTOFF = 57  # Year cutoff for century determination (1957-2056)

# Default values
DEFAULT_INCLINATION = 0.0  # degrees
DEFAULT_RAAN = 0.0  # degrees - Right Ascension of Ascending Node
DEFAULT_ARG_PERIGEE = 0.0  # degrees - Argument of Perigee


def get_mu(body='earth'):
    """
    Get gravitational parameter for different celestial bodies.
    
    Args:
        body (str): Name of celestial body
        
    Returns:
        float: Gravitational parameter in km³/s²
    """
    bodies = {
        'earth': EARTH_MU,
        'moon': 4902.8,
        'sun': 132712440018.0,
        'mars': 42828.37,
        'jupiter': 126686534.0
    }
    
    return bodies.get(body.lower(), EARTH_MU)


def normalize_angle(angle, center=0.0):
    """
    Normalize an angle to be within a range centered at 'center'.
    Default range is [-π, π].
    
    Args:
        angle (float): Angle in radians
        center (float): Center of the range (default: 0)
        
    Returns:
        float: Normalized angle in radians
    """
    return angle - TWO_PI * math.floor((angle + math.pi - center) / TWO_PI)