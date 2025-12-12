"""
Coordinate system transformations for orbital mechanics.
Handles conversions between different reference frames:
- Orbital plane coordinates
- Earth-Centered Inertial (ECI)
- Earth-Centered Earth-Fixed (ECEF)
"""

import math
import numpy as np
from typing import Tuple, List
from datetime import datetime
import constants


class CoordinateTransform:
    """Handle coordinate transformations for orbital mechanics."""
    
    @staticmethod
    def rotation_matrix_x(angle: float) -> np.ndarray:
        """
        Create rotation matrix about X-axis.
        
        Args:
            angle: Rotation angle in radians
            
        Returns:
            3x3 rotation matrix
        """
        c = math.cos(angle)
        s = math.sin(angle)
        return np.array([
            [1, 0, 0],
            [0, c, -s],
            [0, s, c]
        ])
    
    @staticmethod
    def rotation_matrix_y(angle: float) -> np.ndarray:
        """
        Create rotation matrix about Y-axis.
        
        Args:
            angle: Rotation angle in radians
            
        Returns:
            3x3 rotation matrix
        """
        c = math.cos(angle)
        s = math.sin(angle)
        return np.array([
            [c, 0, s],
            [0, 1, 0],
            [-s, 0, c]
        ])
    
    @staticmethod
    def rotation_matrix_z(angle: float) -> np.ndarray:
        """
        Create rotation matrix about Z-axis.
        
        Args:
            angle: Rotation angle in radians
            
        Returns:
            3x3 rotation matrix
        """
        c = math.cos(angle)
        s = math.sin(angle)
        return np.array([
            [c, -s, 0],
            [s, c, 0],
            [0, 0, 1]
        ])
    
    @classmethod
    def orbital_to_eci(cls, r_orbital: np.ndarray, v_orbital: np.ndarray,
                       i: float, raan: float, arg_perigee: float) -> Tuple[np.ndarray, np.ndarray]:
        """
        Transform from orbital plane coordinates to Earth-Centered Inertial (ECI).
        
        The orbital plane coordinate system has:
        - X-axis pointing to periapsis
        - Y-axis in the orbital plane, 90° ahead
        - Z-axis perpendicular to orbital plane (angular momentum direction)
        
        Args:
            r_orbital: Position vector in orbital plane [x, y, z] (km)
            v_orbital: Velocity vector in orbital plane [vx, vy, vz] (km/s)
            i: Inclination (radians)
            raan: Right Ascension of Ascending Node (radians)
            arg_perigee: Argument of perigee (radians)
            
        Returns:
            Tuple of (r_eci, v_eci): Position and velocity in ECI frame
        """
        # Build the transformation matrix from orbital to ECI
        # This is done through three rotations:
        # 1. Rotate by -arg_perigee about Z (moves periapsis to ascending node)
        # 2. Rotate by -i about X (aligns with equatorial plane)
        # 3. Rotate by -raan about Z (aligns with vernal equinox)
        
        R1 = cls.rotation_matrix_z(-arg_perigee)
        R2 = cls.rotation_matrix_x(-i)
        R3 = cls.rotation_matrix_z(-raan)
        
        # Combined transformation matrix
        R_total = R3 @ R2 @ R1
        
        # Transform position and velocity
        r_eci = R_total @ r_orbital
        v_eci = R_total @ v_orbital
        
        return r_eci, v_eci
    
    @classmethod
    def eci_to_orbital(cls, r_eci: np.ndarray, v_eci: np.ndarray,
                       i: float, raan: float, arg_perigee: float) -> Tuple[np.ndarray, np.ndarray]:
        """
        Transform from Earth-Centered Inertial (ECI) to orbital plane coordinates.
        
        Args:
            r_eci: Position vector in ECI frame [x, y, z] (km)
            v_eci: Velocity vector in ECI frame [vx, vy, vz] (km/s)
            i: Inclination (radians)
            raan: Right Ascension of Ascending Node (radians)
            arg_perigee: Argument of perigee (radians)
            
        Returns:
            Tuple of (r_orbital, v_orbital): Position and velocity in orbital frame
        """
        # Inverse transformation (transpose of rotation matrices)
        R1 = cls.rotation_matrix_z(arg_perigee)
        R2 = cls.rotation_matrix_x(i)
        R3 = cls.rotation_matrix_z(raan)
        
        # Combined transformation matrix (reverse order)
        R_total = R1 @ R2 @ R3
        
        # Transform position and velocity
        r_orbital = R_total @ r_eci
        v_orbital = R_total @ v_eci
        
        return r_orbital, v_orbital
    
    @staticmethod
    def perifocal_to_orbital(r: float, true_anomaly: float, 
                           eccentricity: float, 
                           semi_major_axis: float,
                           mu: float = constants.EARTH_MU) -> Tuple[np.ndarray, np.ndarray]:
        """
        Calculate position and velocity in perifocal (orbital plane) coordinates.
        
        Perifocal coordinates have origin at the focus (Earth center) with:
        - X-axis pointing to periapsis
        - Y-axis in orbital plane, 90° ahead
        - Z-axis perpendicular to orbital plane
        
        Args:
            r: Orbital radius (km)
            true_anomaly: True anomaly (radians)
            eccentricity: Orbital eccentricity
            semi_major_axis: Semi-major axis (km)
            mu: Gravitational parameter (km³/s²)
            
        Returns:
            Tuple of (r_vec, v_vec): Position and velocity vectors in orbital plane
        """
        nu = true_anomaly
        
        # Position in perifocal coordinates
        r_vec = np.array([
            r * math.cos(nu),
            r * math.sin(nu),
            0.0
        ])
        
        # Velocity in perifocal coordinates
        # Using vis-viva equation components
        p = semi_major_axis * (1 - eccentricity**2)  # Semi-latus rectum
        h = math.sqrt(mu * p)  # Specific angular momentum magnitude
        
        v_vec = np.array([
            -mu / h * math.sin(nu),
            mu / h * (eccentricity + math.cos(nu)),
            0.0
        ])
        
        return r_vec, v_vec
    
    @staticmethod
    def eci_to_ecef(r_eci: np.ndarray, time: datetime) -> np.ndarray:
        """
        Transform from Earth-Centered Inertial to Earth-Centered Earth-Fixed.
        
        ECEF rotates with the Earth, while ECI is fixed in space.
        
        Args:
            r_eci: Position vector in ECI frame [x, y, z] (km)
            time: UTC time for the transformation
            
        Returns:
            Position vector in ECEF frame [x, y, z] (km)
        """
        # Calculate Greenwich Mean Sidereal Time (GMST)
        # Rudimentary calculation - write a more precise version if perturbations are added
        
        # Days since J2000.0 epoch
        j2000 = datetime(2000, 1, 1, 12, 0, 0)
        delta_t = (time - j2000).total_seconds() / constants.SECONDS_PER_DAY
        
        # GMST in degrees (simplified formula)
        gmst = 280.46061837 + 360.98564736629 * delta_t
        gmst_rad = (gmst % 360) * constants.DEG_TO_RAD
        
        # Rotation matrix about Z-axis
        c = math.cos(gmst_rad)
        s = math.sin(gmst_rad)
        
        r_ecef = np.array([
            c * r_eci[0] + s * r_eci[1],
            -s * r_eci[0] + c * r_eci[1],
            r_eci[2]
        ])
        
        return r_ecef
    
    @staticmethod
    def ecef_to_geodetic(r_ecef: np.ndarray) -> Tuple[float, float, float]:
        """
        Convert ECEF coordinates to geodetic (latitude, longitude, altitude).
        
        Uses iterative method for latitude calculation.
        
        Args:
            r_ecef: Position vector in ECEF frame [x, y, z] (km)
            
        Returns:
            Tuple of (latitude, longitude, altitude):
                latitude in radians (-π/2 to π/2)
                longitude in radians (-π to π)
                altitude in km above ellipsoid
        """
        x, y, z = r_ecef
        
        # Ellipsoid parameters
        a = constants.EARTH_RADIUS  # Equatorial radius (km)
        f = 1/298.257223563  # Flattening
        e2 = 2*f - f**2  # First eccentricity squared
        
        # Longitude is straightforward
        lon = math.atan2(y, x)
        
        # Distance from Z-axis
        p = math.sqrt(x**2 + y**2)
        
        # Initial guess for latitude
        lat = math.atan2(z, p * (1 - e2))
        
        # Iterative refinement for latitude and altitude
        for _ in range(5):  # We found it converges in 1-3 iterations usually
            sin_lat = math.sin(lat)
            N = a / math.sqrt(1 - e2 * sin_lat**2)
            alt = p / math.cos(lat) - N
            lat = math.atan2(z, p * (1 - e2 * N / (N + alt)))
        
        # Final altitude calculation
        sin_lat = math.sin(lat)
        N = a / math.sqrt(1 - e2 * sin_lat**2)
        if abs(lat) < math.pi/4:  # Near equator
            alt = p / math.cos(lat) - N
        else:  # Near poles
            alt = z / sin_lat - N * (1 - e2)
        
        return lat, lon, alt
    
    @staticmethod
    def calculate_ground_track(r_eci_history: List[np.ndarray], 
                              times: List[datetime]) -> List[Tuple[float, float]]:
        """
        Calculate ground track (sub-satellite points) from ECI position history.
        
        Args:
            r_eci_history: List of ECI position vectors
            times: Corresponding UTC times
            
        Returns:
            List of (latitude, longitude) tuples in degrees
        """
        ground_track = []
        transform = CoordinateTransform()
        
        for r_eci, time in zip(r_eci_history, times):
            # Convert to ECEF
            r_ecef = transform.eci_to_ecef(r_eci, time)
            
            # Convert to geodetic
            lat, lon, _ = transform.ecef_to_geodetic(r_ecef)
            
            # Convert to degrees
            ground_track.append((
                lat * constants.RAD_TO_DEG,
                lon * constants.RAD_TO_DEG
            ))
        
        return ground_track
    
    @staticmethod
    def calculate_look_angles(observer_lat: float, observer_lon: float, observer_alt: float,
                             satellite_eci: np.ndarray, time: datetime) -> Tuple[float, float, float]:
        """
        Calculate look angles from ground observer to satellite.
        
        Args:
            observer_lat: Observer latitude (radians)
            observer_lon: Observer longitude (radians)
            observer_alt: Observer altitude (km)
            satellite_eci: Satellite position in ECI (km)
            time: UTC time for calculation
            
        Returns:
            Tuple of (azimuth, elevation, range):
                azimuth in radians (0 to 2π, from North)
                elevation in radians (-π/2 to π/2)
                range in km
        """
        transform = CoordinateTransform()
        
        # Convert satellite position to ECEF
        sat_ecef = transform.eci_to_ecef(satellite_eci, time)
        
        # Observer position in ECEF
        a = constants.EARTH_RADIUS
        N = a  # Assume spherical Earth - future work with perturbations
        obs_ecef = np.array([
            (N + observer_alt) * math.cos(observer_lat) * math.cos(observer_lon),
            (N + observer_alt) * math.cos(observer_lat) * math.sin(observer_lon),
            (N + observer_alt) * math.sin(observer_lat)
        ])
        
        # Range vector from observer to satellite
        range_vec = sat_ecef - obs_ecef
        range_km = np.linalg.norm(range_vec)
        
        # Transform to topocentric horizon coordinates (East-North-Up)
        sin_lat = math.sin(observer_lat)
        cos_lat = math.cos(observer_lat)
        sin_lon = math.sin(observer_lon)
        cos_lon = math.cos(observer_lon)
        
        # Rotation matrix from ECEF to ENU
        R = np.array([
            [-sin_lon, cos_lon, 0],
            [-sin_lat*cos_lon, -sin_lat*sin_lon, cos_lat],
            [cos_lat*cos_lon, cos_lat*sin_lon, sin_lat]
        ])
        
        enu = R @ range_vec
        e, n, u = enu
        
        # Calculate azimuth and elevation
        azimuth = math.atan2(e, n)  # From North, clockwise
        if azimuth < 0:
            azimuth += constants.TWO_PI
        
        horizontal_range = math.sqrt(e**2 + n**2)
        elevation = math.atan2(u, horizontal_range)
        
        return azimuth, elevation, range_km