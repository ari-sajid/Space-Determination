"""
Simple 3D orbit visualizer.

Usage:
    python visualization.py
This will open a matplotlib 3D plot showing Earth and an example orbit.
"""
import numpy as np
import math

from mpl_toolkits.mplot3d import Axes3D  # noqa: F401 (needed by mpl)
import matplotlib.pyplot as plt

from constants import EARTH_RADIUS
from utils import altitude_to_semi_major_axis, orbital_positions


def _draw_earth(ax, radius_km=EARTH_RADIUS, color=(0.2, 0.4, 0.8), alpha=0.6):
    """Draw a textured-like sphere for Earth (simple colored sphere)."""
    u = np.linspace(0, 2 * math.pi, 60)
    v = np.linspace(0, math.pi, 30)
    x = radius_km * np.outer(np.cos(u), np.sin(v))
    y = radius_km * np.outer(np.sin(u), np.sin(v))
    z = radius_km * np.outer(np.ones_like(u), np.cos(v))
    ax.plot_surface(x, y, z, rstride=1, cstride=1, color=color, alpha=alpha, linewidth=0)


def plot_orbit_from_altitude(altitude_km: float,
                             inclination_deg: float = 0.0,
                             raan_deg: float = 0.0,
                             arg_perigee_deg: float = 0.0,
                             num_points: int = 720,
                             show=True,
                             save_path: str = None):
    """
    Plot a circular orbit defined by altitude and orientation.

    Args:
        altitude_km: altitude above Earth's surface (km)
        inclination_deg: inclination in degrees
        raan_deg: RAAN in degrees
        arg_perigee_deg: argument of perigee in degrees
        num_points: number of points to render along orbit
        show: whether to call plt.show()
        save_path: if provided, save the figure to this path
    """
    sma = altitude_to_semi_major_axis(altitude_km)

    # Convert angles to radians
    incl = math.radians(inclination_deg)
    raan = math.radians(raan_deg)
    argp = math.radians(arg_perigee_deg)

    # Compute ECI positions
    pos = orbital_positions(sma, inclination=incl, raan=raan, arg_perigee=argp, num_points=num_points)

    # Plot
    fig = plt.figure(figsize=(8, 8))
    ax = fig.add_subplot(111, projection='3d')
    _draw_earth(ax, radius_km=EARTH_RADIUS)

    # Orbit line
    ax.plot(pos[:, 0], pos[:, 1], pos[:, 2], color='red', linewidth=1.5, label='Orbit')

    # Mark a few points (ascending node, etc.)
    ax.scatter([pos[0, 0]], [pos[0, 1]], [pos[0, 2]], color='black', s=30, label='Start')

    # Set equal aspect (approximation)
    max_range = np.max(np.abs(pos)) * 1.1
    for axis in (ax.set_xlim, ax.set_ylim, ax.set_zlim):
        axis([-max_range, max_range])

    ax.set_xlabel('ECI X (km)')
    ax.set_ylabel('ECI Y (km)')
    ax.set_zlabel('ECI Z (km)')
    ax.set_title(f'Orbit: alt {altitude_km} km, incl {inclination_deg}°')
    ax.legend()

    if save_path:
        plt.savefig(save_path, dpi=200, bbox_inches='tight')

    if show:
        plt.show()


if __name__ == "__main__":
    # Example: ISS-like orbit
    iss_altitude = 400.0
    iss_incl = 51.6
    plot_orbit_from_altitude(iss_altitude, inclination_deg=iss_incl, raan_deg=0.0, arg_perigee_deg=0.0)
