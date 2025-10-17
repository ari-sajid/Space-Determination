"""
Utility functions for the orbit propagation tool.
Common helper functions used across multiple modules.
"""

import math
import numpy as np
from datetime import datetime, timedelta
from typing import Tuple, List, Optional
import constants


def julian_date(dt: datetime) -> float:
    """
    Convert datetime to Julian Date.
    
    Args:
        dt: DateTime object (UTC)
        
    Returns:
        Julian Date as float
    """
    # Reference: J2000.0 = January 1, 2000, 12:00 UTC = JD 2451545.0
    j2000 = datetime(2000, 1, 1, 12, 0, 0)
    delta = dt - j2000
    jd = 2451545.0 + delta.total_seconds() / 86400.0
    return jd


def gmst(dt: datetime) -> float:
    """
    Calculate Greenwich Mean Sidereal Time.
    
    Args:
        dt: DateTime object (UTC)
        
    Returns:
        GMST in radians
    """
    # Simplified GMST calculation
    jd = julian_date(dt)
    t_centuries = (jd - 2451545.0) / 36525.0
    
    # GMST at 0h UTC
    gmst_0h = 280.46061837 + 360.98564736629 * (jd - 2451545.0) + \
              0.000387933 * t_centuries**2 - t_centuries**3 / 38710000.0
    
    # Add the time since 0h UTC
    hour_angle = dt.hour + dt.minute/60.0 + dt.second/3600.0
    gmst_now = gmst_0h + hour_angle * 15.0  # 15 degrees per hour
    
    # Normalize to [0, 360] degrees and convert to radians
    gmst_rad = (gmst_now % 360.0) * constants.DEG_TO_RAD
    
    return gmst_rad


def orbital_period(semi_major_axis: float, mu: float = constants.EARTH_MU) -> float:
    """
    Calculate orbital period from semi-major axis.
    
    Args:
        semi_major_axis: Semi-major axis in km
        mu: Gravitational parameter in km³/s²
        
    Returns:
        Orbital period in seconds
    """
    return constants.TWO_PI * math.sqrt(semi_major_axis**3 / mu)


def altitude_to_semi_major_axis(altitude: float) -> float:
    """
    Convert altitude above Earth's surface to semi-major axis.
    
    Args:
        altitude: Altitude in km
        
    Returns:
        Semi-major axis in km
    """
    return constants.EARTH_RADIUS + altitude


def semi_major_axis_to_altitude(semi_major_axis: float) -> float:
    """
    Convert semi-major axis to mean altitude.
    
    Args:
        semi_major_axis: Semi-major axis in km
        
    Returns:
        Mean altitude in km
    """
    return semi_major_axis - constants.EARTH_RADIUS


def escape_velocity(radius: float, mu: float = constants.EARTH_MU) -> float:
    """
    Calculate escape velocity at a given radius.
    
    Args:
        radius: Distance from Earth's center in km
        mu: Gravitational parameter in km³/s²
        
    Returns:
        Escape velocity in km/s
    """
    return math.sqrt(2 * mu / radius)


def circular_velocity(radius: float, mu: float = constants.EARTH_MU) -> float:
    """
    Calculate circular orbital velocity at a given radius.
    
    Args:
        radius: Distance from Earth's center in km
        mu: Gravitational parameter in km³/s²
        
    Returns:
        Circular velocity in km/s
    """
    return math.sqrt(mu / radius)


def hohmann_transfer(r1: float, r2: float, mu: float = constants.EARTH_MU) -> dict:
    """
    Calculate Hohmann transfer orbit parameters.
    
    Args:
        r1: Initial orbit radius in km
        r2: Final orbit radius in km
        mu: Gravitational parameter in km³/s²
        
    Returns:
        Dictionary with transfer parameters
    """
    # Semi-major axis of transfer orbit
    a_transfer = (r1 + r2) / 2
    
    # Velocities
    v1_circular = circular_velocity(r1, mu)
    v2_circular = circular_velocity(r2, mu)
    
    # Transfer orbit velocities at periapsis and apoapsis
    v_transfer_peri = math.sqrt(mu * (2/r1 - 1/a_transfer))
    v_transfer_apo = math.sqrt(mu * (2/r2 - 1/a_transfer))
    
    # Delta-v requirements
    delta_v1 = abs(v_transfer_peri - v1_circular)
    delta_v2 = abs(v2_circular - v_transfer_apo)
    delta_v_total = delta_v1 + delta_v2
    
    # Transfer time (half period of transfer orbit)
    transfer_time = math.pi * math.sqrt(a_transfer**3 / mu)
    
    return {
        'delta_v1_km_s': delta_v1,
        'delta_v2_km_s': delta_v2,
        'delta_v_total_km_s': delta_v_total,
        'transfer_time_seconds': transfer_time,
        'transfer_time_hours': transfer_time / 3600,
        'transfer_semi_major_axis_km': a_transfer
    }


def sun_synchronous_inclination(altitude: float, eccentricity: float = 0.0) -> float:
    """
    Calculate required inclination for sun-synchronous orbit.
    
    Args:
        altitude: Altitude above Earth's surface in km
        eccentricity: Orbital eccentricity
        
    Returns:
        Required inclination in radians
    """
    a = altitude_to_semi_major_axis(altitude)
    
    # Approximate formula for sun-synchronous inclination
    # Requires precession rate = Earth's mean motion around Sun
    n = math.sqrt(constants.EARTH_MU / a**3)  # Mean motion
    
    # Required precession rate (approximately 360 degrees per year)
    precession_rate = constants.TWO_PI / (365.25 * 86400)  # rad/s
    
    # Calculate inclination using J2 perturbation
    cos_i = -precession_rate * 2 * a**2 * (1 - eccentricity**2)**2 / \
            (3 * n * constants.EARTH_J2 * constants.EARTH_RADIUS**2)
    
    # Ensure valid range
    if abs(cos_i) > 1.0:
        raise ValueError(f"Sun-synchronous orbit not possible at altitude {altitude} km")
    
    return math.acos(cos_i)


def atmospheric_density(altitude: float) -> float:
    """
    Simple exponential atmosphere model for drag calculations.
    
    Args:
        altitude: Altitude above Earth's surface in km
        
    Returns:
        Atmospheric density in kg/m³
    """
    # Reference densities and scale heights for different altitude ranges
    # This is a very simplified model
    
    if altitude < 0:
        return 1.225  # Sea level density
    elif altitude < 100:
        rho_0 = 1.225
        H = 8.5  # Scale height in km
        return rho_0 * math.exp(-altitude / H)
    elif altitude < 200:
        rho_0 = 2.5e-9
        H = 20.0
        return rho_0 * math.exp(-(altitude - 100) / H)
    elif altitude < 500:
        rho_0 = 5.0e-12
        H = 40.0
        return rho_0 * math.exp(-(altitude - 200) / H)
    else:
        # Very thin atmosphere above 500 km
        rho_0 = 1e-14
        H = 100.0
        return rho_0 * math.exp(-(altitude - 500) / H)


def visibility_angles(observer_altitude: float = 0.0, 
                     minimum_elevation: float = 0.0) -> float:
    """
    Calculate maximum angle from nadir for satellite visibility.
    
    Args:
        observer_altitude: Observer altitude above sea level in km
        minimum_elevation: Minimum elevation angle in radians
        
    Returns:
        Maximum angle from nadir in radians
    """
    r_earth = constants.EARTH_RADIUS
    r_observer = r_earth + observer_altitude
    
    # Maximum Earth central angle for visibility
    if minimum_elevation <= 0:
        # Geometric horizon
        return math.acos(r_earth / r_observer)
    else:
        # With minimum elevation constraint
        return math.acos(r_earth / r_observer) - minimum_elevation


def groundtrack_shift_per_orbit(period: float) -> float:
    """
    Calculate westward shift of ground track per orbit due to Earth rotation.
    
    Args:
        period: Orbital period in seconds
        
    Returns:
        Westward shift in degrees longitude per orbit
    """
    earth_rotation_deg_per_sec = 360.0 / constants.SECONDS_PER_DAY
    shift = earth_rotation_deg_per_sec * period
    return shift


def repeating_ground_track(revolutions: int, days: int) -> Tuple[float, float]:
    """
    Calculate orbit parameters for repeating ground track.
    
    Args:
        revolutions: Number of satellite revolutions
        days: Number of days for repeat cycle
        
    Returns:
        Tuple of (semi_major_axis, altitude) in km
    """
    # Period for repeating ground track
    period = days * constants.SECONDS_PER_DAY / revolutions
    
    # Semi-major axis from period
    a = (constants.EARTH_MU * period**2 / (4 * math.pi**2)) ** (1/3)
    
    # Altitude
    altitude = a - constants.EARTH_RADIUS
    
    return a, altitude


def eclipse_duration(semi_major_axis: float, eccentricity: float = 0.0) -> float:
    """
    Estimate maximum eclipse duration for a satellite.
    Simplified calculation assuming circular orbit and cylindrical Earth shadow.
    
    Args:
        semi_major_axis: Semi-major axis in km
        eccentricity: Orbital eccentricity
        
    Returns:
        Maximum eclipse duration in seconds
    """
    if eccentricity > 0.1:
        # For highly elliptical orbits, this is more complex
        print("Warning: Eclipse calculation simplified for low eccentricity orbits")
    
    # Angle subtended by Earth as seen from orbit
    beta = math.asin(constants.EARTH_RADIUS / semi_major_axis)
    
    # Eclipse occurs over approximately this angle (simplified)
    eclipse_angle = 2 * beta
    
    # Period
    period = orbital_period(semi_major_axis)
    
    # Eclipse duration
    eclipse_time = (eclipse_angle / constants.TWO_PI) * period
    
    return eclipse_time


def orbital_positions(semi_major_axis: float,
                      inclination: float = 0.0,
                      raan: float = 0.0,
                      arg_perigee: float = 0.0,
                      num_points: int = 360) -> np.ndarray:
    """
    Generate ECI position samples for a (near-)circular orbit given classical elements.
    Assumes semi_major_axis in km and angles in radians. Returns an (N,3) array of positions (km).

    Args:
        semi_major_axis: semi-major axis (km). For circular orbits this is orbital radius.
        inclination: inclination (radians)
        raan: right ascension of the ascending node (radians)
        arg_perigee: argument of perigee (radians)
        num_points: number of sample points around the orbit

    Returns:
        positions: numpy.ndarray shape (num_points, 3) of ECI coordinates (km)
    """
    # Circular radius
    r = float(semi_major_axis)

    # True anomaly samples
    nu = np.linspace(0.0, 2.0 * math.pi, num_points, endpoint=False)

    # Positions in perifocal frame (p-frame)
    x_pf = r * np.cos(nu)
    y_pf = r * np.sin(nu)
    z_pf = np.zeros_like(nu)
    r_pf = np.vstack((x_pf, y_pf, z_pf))  # shape (3, N)

    # Rotation matrices
    def R_z(theta: float) -> np.ndarray:
        c = math.cos(theta); s = math.sin(theta)
        return np.array([[c, -s, 0.0],
                         [s,  c, 0.0],
                         [0.0, 0.0, 1.0]])

    def R_x(theta: float) -> np.ndarray:
        c = math.cos(theta); s = math.sin(theta)
        return np.array([[1.0, 0.0, 0.0],
                         [0.0,  c, -s],
                         [0.0,  s,  c]])

    # Perifocal -> ECI: R = R_z(raan) * R_x(inclination) * R_z(arg_perigee)
    R = R_z(raan) @ R_x(inclination) @ R_z(arg_perigee)

    # Rotate all samples into ECI
    r_eci = R @ r_pf  # shape (3, N)
    return r_eci.T  # shape (N, 3)


def test_utils():
    """Test utility functions."""
    print("Testing Utility Functions...")
    print("-" * 40)
    
    # Test 1: Orbital period
    print("Test 1: Orbital Period")
    iss_altitude = 400  # km
    iss_sma = altitude_to_semi_major_axis(iss_altitude)
    period = orbital_period(iss_sma)
    print(f"  ISS altitude: {iss_altitude} km")
    print(f"  Semi-major axis: {iss_sma:.1f} km")
    print(f"  Period: {period/60:.1f} minutes")
    
    # Test 2: Hohmann transfer
    print("\nTest 2: Hohmann Transfer (LEO to GEO)")
    r_leo = altitude_to_semi_major_axis(400)  # 400 km altitude
    r_geo = altitude_to_semi_major_axis(35786)  # GEO altitude
    transfer = hohmann_transfer(r_leo, r_geo)
    print(f"  Delta-V total: {transfer['delta_v_total_km_s']:.3f} km/s")
    print(f"  Transfer time: {transfer['transfer_time_hours']:.1f} hours")
    
    # Test 3: Sun-synchronous orbit
    print("\nTest 3: Sun-Synchronous Orbit")
    sso_altitude = 800  # km
    try:
        sso_inclination = sun_synchronous_inclination(sso_altitude)
        print(f"  Altitude: {sso_altitude} km")
        print(f"  Required inclination: {math.degrees(sso_inclination):.2f}°")
    except ValueError as e:
        print(f"  Error: {e}")
    
    # Test 4: Atmospheric density
    print("\nTest 4: Atmospheric Density")
    altitudes = [0, 100, 200, 400]
    for alt in altitudes:
        density = atmospheric_density(alt)
        print(f"  Altitude {alt:3d} km: {density:.2e} kg/m³")
    
    # Test 5: Repeating ground track
    print("\nTest 5: Repeating Ground Track")
    revs = 15  # revolutions
    days = 1   # day
    sma, altitude = repeating_ground_track(revs, days)
    print(f"  {revs} revolutions in {days} day(s)")
    print(f"  Required altitude: {altitude:.1f} km")
    print(f"  Semi-major axis: {sma:.1f} km")
    
    # Test 6: Eclipse duration
    print("\nTest 6: Eclipse Duration")
    leo_sma = altitude_to_semi_major_axis(400)
    eclipse_time = eclipse_duration(leo_sma)
    print(f"  LEO (400 km altitude)")
    print(f"  Max eclipse duration: {eclipse_time/60:.1f} minutes")
    
    # Test 7: Orbital positions
    print("\nTest 7: Orbital Positions")
    sma_test = altitude_to_semi_major_axis(500)  # 500 km altitude
    incl_test = math.radians(98.7)  # Retrograde
    raan_test = 0.0
    argp_test = 0.0
    positions = orbital_positions(sma_test, incl_test, raan_test, argp_test, num_points=10)
    for i, pos in enumerate(positions):
        print(f"  Position {i+1}: X={pos[0]:.1f} km, Y={pos[1]:.1f} km, Z={pos[2]:.1f} km")
    
    print("-" * 40)
    print("Tests complete!")


if __name__ == "__main__":
    test_utils()
