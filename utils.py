"""
Utility functions for the orbit propagation tool
"""

import math
import numpy as np
from typing import Tuple, List, Optional
import constants


def altitude_to_semi_major_axis(altitude: float) -> float:
    """
    Convert altitude above Earth's surface to semi-major axis.
    
    Args:
        altitude: Altitude in km
        
    Returns:
        Semi-major axis in km
    """
    return constants.EARTH_RADIUS + altitude

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

