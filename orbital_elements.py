"""
Orbital elements data structures and conversions.
Represents the six classical orbital elements that uniquely define an orbit.
"""

import math
from dataclasses import dataclass
from datetime import datetime
from typing import Optional, Tuple
import constants


@dataclass
class OrbitalElements:
    """
    Classical orbital elements representation.
    
    Attributes:
        a: Semi-major axis (km)
        e: Eccentricity (dimensionless)
        i: Inclination (radians)
        raan: Right Ascension of Ascending Node (radians)
        arg_perigee: Argument of perigee (radians)
        mean_anomaly: Mean anomaly at epoch (radians)
        epoch: Reference epoch for the elements (datetime)
        mu: Gravitational parameter (km³/s²)
    """
    
    a: float  # Semi-major axis (km)
    e: float  # Eccentricity
    i: float  # Inclination (radians)
    raan: float  # Right Ascension of Ascending Node (radians)
    arg_perigee: float  # Argument of perigee (radians)
    mean_anomaly: float  # Mean anomaly at epoch (radians)
    epoch: datetime  # Reference epoch
    mu: float = constants.EARTH_MU  # Gravitational parameter
    
    def __post_init__(self):
        """
        Validate orbital elements after initialization.
        
        Raises
        ------
        ValueError
            If orbital elements are invalid 
        """
        self.validate()
    
    def validate(self):
        """
        Validate that orbital elements are physically meaningful.
        
        Raises:
            ValueError: If any orbital element is invalid
        """
        if self.a <= 0:
            raise ValueError(f"Semi-major axis must be positive, got {self.a}")
        
        if self.e < 0:
            raise ValueError(f"Eccentricity must be non-negative, got {self.e}")
        
        if self.i < 0 or self.i > math.pi:
            raise ValueError(f"Inclination must be between 0 and π, got {self.i}")
        
        # Normalize angles to [0, 2π]
        self.raan = constants.normalize_angle(self.raan, math.pi)
        self.arg_perigee = constants.normalize_angle(self.arg_perigee, math.pi)
        self.mean_anomaly = constants.normalize_angle(self.mean_anomaly, math.pi)
    
    @property
    def period(self) -> float:
        """
        Calculate orbital period.
        
        Returns:
            float: Orbital period in seconds
            None for non-elliptical orbits
        """
        if self.e >= 1.0:  # Parabolic or hyperbolic orbit
            return None
        return constants.TWO_PI * math.sqrt(self.a**3 / self.mu)
    
    @property
    def mean_motion(self) -> float:
        """
        Calculate mean motion n.
        
        Returns:
            float: Mean motion in rad/s
        """
        return math.sqrt(self.mu / self.a**3)
    
    @property
    def apogee(self) -> Optional[float]:
        """
        Calculate apogee distance (furthest point from Earth).
        
        Returns:
            float: Apogee distance in km 
            None for non-elliptical orbits
        """
        if self.e >= 1.0:
            return None
        return self.a * (1 + self.e)
    
    @property
    def perigee(self) -> float:
        """
        Calculate perigee distance (closest point to Earth).
        
        Returns:
            float: Perigee distance in km
        """
        return self.a * (1 - self.e)
    
    @property
    def orbit_type(self) -> str:
        """
        Determine the type of orbit based on eccentricity.
        
        Returns:
            str: Type of orbit (circular, elliptical, parabolic, or hyperbolic)
        """
        if self.e < constants.CIRCULAR_ECC_THRESHOLD:
            return "circular"
        elif self.e < constants.PARABOLIC_ECC_THRESHOLD:
            return "elliptical"
        elif self.e < constants.HYPERBOLIC_ECC_THRESHOLD:
            return "parabolic"
        else:
            return "hyperbolic"
    
    @classmethod
    def from_state_vector(cls, r: Tuple[float, float, float], 
                         v: Tuple[float, float, float],
                         epoch: datetime,
                         mu: float = constants.EARTH_MU):
        """
        Create orbital elements from position and velocity vectors.
        
        Args:
            r: Position vector (x, y, z) in km
            v: Velocity vector (vx, vy, vz) in km/s
            epoch: Epoch for the state vector
            mu: Gravitational parameter (km³/s²)
            
        Returns:
            OrbitalElements: Computed orbital elements

        Raises
        ------
        ValueError
            If position or velocity vector is not of length 3

        Example
        --------
        >>> from datetime import datetime
        >>> r = (7000, 0, 0)
        >>> v = (0, 7.5, 0)
        >>> epoch = datetime.utcnow()
        >>> oe = OrbitalElements.from_state_vector(r, v, epoch)
        >>> print(oe.a)
        7000.0
        """
        # This is a placeholder - implement the full conversion
        # This involves computing h (angular momentum), n (node vector),
        # and using various vector operations
        
        # For now, return dummy values
        print("Warning: from_state_vector not fully implemented yet")
        return cls(
            a=7000.0,  # dummy value
            e=0.001,   # dummy value
            i=0.0,     # dummy value
            raan=0.0,  # dummy value
            arg_perigee=0.0,  # dummy value
            mean_anomaly=0.0,  # dummy value
            epoch=epoch,
            mu=mu
        )
    
    def to_dict(self) -> dict:
        """
        Convert orbital elements to dictionary.
        
        Returns:
            dict: Dictionary representation of orbital elements, shown below 

        dict
            Dictionary with keys:
            - 'semi_major_axis_km'
            - 'eccentricity'
            - 'inclination_deg'
            - 'raan_deg'
            - 'arg_perigee_deg'
            - 'mean_anomaly_deg'
            - 'epoch'
            - 'period_hours'
            - 'orbit_typ

        """
        return {
            'semi_major_axis_km': self.a,
            'eccentricity': self.e,
            'inclination_deg': self.i * constants.RAD_TO_DEG,
            'raan_deg': self.raan * constants.RAD_TO_DEG,
            'arg_perigee_deg': self.arg_perigee * constants.RAD_TO_DEG,
            'mean_anomaly_deg': self.mean_anomaly * constants.RAD_TO_DEG,
            'epoch': self.epoch.isoformat(),
            'period_hours': self.period / 3600 if self.period else None,
            'orbit_type': self.orbit_type
        }
    
    def __str__(self) -> str:
        """String representation of orbital elements.
        
        Returns
        -------
        str
            human-readable formatted string of orbital elements
            
        """
        return (
            f"Orbital Elements (Epoch: {self.epoch.isoformat()})\n"
            f"  Semi-major axis: {self.a:.3f} km\n"
            f"  Eccentricity: {self.e:.6f}\n"
            f"  Inclination: {self.i * constants.RAD_TO_DEG:.3f}°\n"
            f"  RAAN: {self.raan * constants.RAD_TO_DEG:.3f}°\n"
            f"  Arg of Perigee: {self.arg_perigee * constants.RAD_TO_DEG:.3f}°\n"
            f"  Mean Anomaly: {self.mean_anomaly * constants.RAD_TO_DEG:.3f}°\n"
            f"  Orbit Type: {self.orbit_type}\n"
            f"  Period: {self.period/3600:.2f} hours" if self.period else ""
        )
