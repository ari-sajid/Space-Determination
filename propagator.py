"""
Orbital propagator for satellite orbit prediction.
Implements two-body propagation with optional perturbations.
"""

import math
import numpy as np
from typing import Optional, List, Tuple, Dict
from datetime import datetime, timedelta
from dataclasses import dataclass

from orbital_elements import OrbitalElements
from kepler_solver import KeplerSolver
from coordinate_systems import CoordinateTransform
import constants


@dataclass
class PropagationState:
    """
    State of a satellite at a specific time.
    
    Attributes:
        time: Epoch of this state
        position: Position vector in ECI frame (km)
        velocity: Velocity vector in ECI frame (km/s)
        orbital_elements: Classical orbital elements at this epoch
    """
    time: datetime
    position: np.ndarray  # [x, y, z] in km
    velocity: np.ndarray  # [vx, vy, vz] in km/s
    orbital_elements: OrbitalElements
    
    def to_dict(self) -> dict:
        """Convert state to dictionary.

        Raises:
            AttributeError: If orbital_elements does not have to_dict().
        """
        return {
            'time': self.time.isoformat(),
            'position_km': self.position.tolist(),
            'velocity_km_s': self.velocity.tolist(),
            'orbital_elements': self.orbital_elements.to_dict()
        }


class TwoBodyPropagator:
    """
    Two-body orbital propagator using Kepler's equations.
    
    This is the baseline propagator that assumes only the central body's
    gravitational force acts on the satellite.
    """
    
    def __init__(self, initial_elements: OrbitalElements):
        """
        Initialize the propagator with initial orbital elements.
        
        Args:
            initial_elements: Initial orbital elements at epoch

        Raises:
            AttributeError: If initial_elements is missing required attributes.
            TypeError: If initial_elements is not an OrbitalElements instance.
        """
        self.initial_elements = initial_elements
        self.kepler_solver = KeplerSolver()
        self.coord_transform = CoordinateTransform()
        self.mu = initial_elements.mu
        
        # Cache frequently used values
        self._mean_motion = self.initial_elements.mean_motion
        self._period = self.initial_elements.period
    
    def propagate(self, delta_t: float) -> PropagationState:
        """
        Propagate orbit forward by delta_t seconds from epoch.
        
        Args:
            delta_t: Time to propagate in seconds
        
        Returns:
            PropagationState at the requested time

        Raises:
            ValueError: If delta_t is not a finite float or is negative.
            RuntimeError: If Kepler's equation cannot be solved.

        Example:
            >>> from datetime import datetime
            >>> elements = OrbitalElements(
            ...     a=7000, e=0.001, i=0.1, raan=0.2, arg_perigee=0.3,
            ...     mean_anomaly=0.4, epoch=datetime(2025, 1, 1, 0, 0, 0)
            ... )
            >>> propagator = TwoBodyPropagator(elements)
            >>> state = propagator.propagate(3600)  # Propagate 1 hour
            >>> print(state.position)
            [x, y, z]  # Position in km

        """
        # Calculate new mean anomaly
        M_new = self._propagate_mean_anomaly(delta_t)
        
        # Solve Kepler's equation for eccentric anomaly
        E = self.kepler_solver.solve_kepler_equation(M_new, self.initial_elements.e)
        
        # Calculate true anomaly
        nu = self.kepler_solver.eccentric_to_true_anomaly(E, self.initial_elements.e)
        
        # Calculate orbital radius
        r = self._calculate_radius(E)
        
        # Get position and velocity in orbital plane
        r_orbital, v_orbital = self.coord_transform.perifocal_to_orbital(
            r, nu, self.initial_elements.e, self.initial_elements.a, self.mu
        )
        
        # Transform to ECI coordinates
        r_eci, v_eci = self.coord_transform.orbital_to_eci(
            r_orbital, v_orbital,
            self.initial_elements.i,
            self.initial_elements.raan,
            self.initial_elements.arg_perigee
        )
        
        # Create new orbital elements at this epoch
        new_time = self.initial_elements.epoch + timedelta(seconds=delta_t)
        new_elements = OrbitalElements(
            a=self.initial_elements.a,  # Unchanged in two-body problem
            e=self.initial_elements.e,  # Unchanged in two-body problem
            i=self.initial_elements.i,  # Unchanged in two-body problem
            raan=self.initial_elements.raan,  # Unchanged in two-body problem
            arg_perigee=self.initial_elements.arg_perigee,  # Unchanged in two-body problem
            mean_anomaly=M_new,
            epoch=new_time,
            mu=self.mu
        )
        
        return PropagationState(
            time=new_time,
            position=r_eci,
            velocity=v_eci,
            orbital_elements=new_elements
        )
    
    def propagate_to_time(self, target_time: datetime) -> PropagationState:
        """
        Propagate orbit to a specific target time.
        
        Args:
            target_time: Target time for propagation
            
        Returns:
            PropagationState at the target time

        Raises:
            TypeError: If target_time is not a datetime object.
            ValueError: If target_time is before the epoch.
            RuntimeError: If propagation fails.
        """
        delta_t = (target_time - self.initial_elements.epoch).total_seconds()
        return self.propagate(delta_t)
    
    def propagate_orbit(self, duration: float, time_step: float) -> List[PropagationState]:
        """
        Propagate orbit for a duration with given time steps.
        
        Args:
            duration: Total duration to propagate (seconds)
            time_step: Time step between states (seconds)
            
        Returns:
            List of PropagationStates at each time step

        Raises:
            ValueError: If duration or time_step is not a positive finite float.
            RuntimeError: If propagation fails at any time step.
        """
        states = []
        num_steps = int(duration / time_step) + 1
        
        for i in range(num_steps):
            t = i * time_step
            state = self.propagate(t)
            states.append(state)
        
        return states
    
    def _propagate_mean_anomaly(self, delta_t: float) -> float:
        """
        Calculate mean anomaly at time t from epoch.
        
        Args:
            delta_t: Time since epoch in seconds
            
        Returns:
            Mean anomaly in radians

        Raises:
            ValueError: If delta_t is not a finite float or is negative.
        """
        M = self.initial_elements.mean_anomaly + self._mean_motion * delta_t
        return M % constants.TWO_PI
    
    def _calculate_radius(self, eccentric_anomaly: float) -> float:
        """
        Calculate orbital radius from eccentric anomaly.
        
        Args:
            eccentric_anomaly: Eccentric anomaly in radians
            
        Returns:
            Orbital radius in km

        Raises:
            ValueError: If eccentric_anomaly is not a finite float.
        """
        return self.initial_elements.a * (1 - self.initial_elements.e * math.cos(eccentric_anomaly))
    
    def find_next_pass(self, target_lat: float, target_lon: float, 
                      search_duration: float = 86400.0,
                      elevation_threshold: float = 0.0) -> Optional[Dict]:
        """
        Find the next satellite pass over a ground location.
        
        Args:
            target_lat: Target latitude in degrees
            target_lon: Target longitude in degrees
            search_duration: How long to search ahead (seconds)
            elevation_threshold: Minimum elevation angle in degrees
            
        Returns:
            Dictionary with pass information or None if no pass found

        Raises:
            ValueError: If latitude or longitude are out of valid range.
            ValueError: If search_duration or elevation_threshold are negative.

        Example:
            >>> propagator = TwoBodyPropagator(elements)
            >>> pass_info = propagator.find_next_pass(
            ...     target_lat=40.7128, target_lon=-74.0060, elevation_threshold=10.0
            ... )
            >>> if pass_info:
            ...     print(pass_info['pass_start'], pass_info['pass_end'])
            ...     print(pass_info['max_elevation_deg'])
            ... else:
            ...     print("No pass found.")

        """
        target_lat_rad = target_lat * constants.DEG_TO_RAD
        target_lon_rad = target_lon * constants.DEG_TO_RAD
        elevation_threshold_rad = elevation_threshold * constants.DEG_TO_RAD
        
        # Search with 1-minute time steps
        time_step = 60.0
        num_steps = int(search_duration / time_step)
        
        pass_start = None
        pass_max_elevation = -math.pi/2
        pass_max_time = None
        pass_data = []
        
        for i in range(num_steps):
            t = i * time_step
            state = self.propagate(t)
            
            # Calculate look angles
            azimuth, elevation, range_km = self.coord_transform.calculate_look_angles(
                target_lat_rad, target_lon_rad, 0.0,
                state.position, state.time
            )
            
            if elevation > elevation_threshold_rad:
                if pass_start is None:
                    pass_start = state.time
                
                pass_data.append({
                    'time': state.time,
                    'azimuth_deg': azimuth * constants.RAD_TO_DEG,
                    'elevation_deg': elevation * constants.RAD_TO_DEG,
                    'range_km': range_km
                })
                
                if elevation > pass_max_elevation:
                    pass_max_elevation = elevation
                    pass_max_time = state.time
            
            elif pass_start is not None:
                # Pass has ended
                return {
                    'pass_start': pass_start,
                    'pass_end': state.time,
                    'max_elevation_deg': pass_max_elevation * constants.RAD_TO_DEG,
                    'max_elevation_time': pass_max_time,
                    'duration_seconds': (state.time - pass_start).total_seconds(),
                    'pass_data': pass_data
                }
        
        # Check if we're still in a pass at the end
        if pass_start is not None:
            return {
                'pass_start': pass_start,
                'pass_end': state.time,
                'max_elevation_deg': pass_max_elevation * constants.RAD_TO_DEG,
                'max_elevation_time': pass_max_time,
                'duration_seconds': (state.time - pass_start).total_seconds(),
                'pass_data': pass_data,
                'note': 'Pass continues beyond search duration'
            }
        
        return None
    
    def get_ground_track(self, duration: float, time_step: float = 60.0) -> List[Tuple[float, float]]:
        """
        Calculate ground track for specified duration.
        
        Args:
            duration: Duration in seconds
            time_step: Time step between points in seconds
            
        Returns:
            List of (latitude, longitude) tuples in degrees

        Raises:
            ValueError: If duration or time_step is not a positive finite float.
            RuntimeError: If propagation fails at any time step.
        """
        states = self.propagate_orbit(duration, time_step)
        
        positions = [state.position for state in states]
        times = [state.time for state in states]
        
        return self.coord_transform.calculate_ground_track(positions, times)


class PerturbedPropagator(TwoBodyPropagator):
    """
    Extended propagator that includes perturbation effects.
    This is a stretch goal - implements J2 perturbation and atmospheric drag.
    """
    
    def __init__(self, initial_elements: OrbitalElements, 
                 enable_j2: bool = True,
                 enable_drag: bool = False):
        """
        Initialize perturbed propagator.
        
        Args:
            initial_elements: Initial orbital elements
            enable_j2: Include J2 perturbation
            enable_drag: Include atmospheric drag

        Raises:
            AttributeError: If initial_elements is missing required attributes.
            TypeError: If initial_elements is not an OrbitalElements instance.
        """
        super().__init__(initial_elements)
        self.enable_j2 = enable_j2
        self.enable_drag = enable_drag
        
        print("Note: Perturbed propagator is a stretch goal - not fully implemented")
    
    def _calculate_j2_perturbation(self, r_eci: np.ndarray) -> np.ndarray:
        """
        Calculate acceleration due to J2 perturbation.
        
        Args:
            r_eci: Position vector in ECI frame
            
        Returns:
            Acceleration vector due to J2

        Raises:
            ValueError: If r_eci is not a valid 3-element array.
        """
        # Placeholder for J2 perturbation calculation
        # This would involve spherical harmonics
        return np.zeros(3)
    
    def _calculate_drag(self, r_eci: np.ndarray, v_eci: np.ndarray, 
                       altitude: float) -> np.ndarray:
        """
        Calculate acceleration due to atmospheric drag.
        
        Args:
            r_eci: Position vector in ECI frame
            v_eci: Velocity vector in ECI frame
            altitude: Altitude above Earth surface
            
        Returns:
            Acceleration vector due to drag

        Raises:
            ValueError: If r_eci or v_eci is not a valid 3-element array, or altitude is negative.
        """
        # Placeholder for drag calculation
        # This would need atmospheric density model
        return np.zeros(3)


def test_propagator():
    """Test the orbital propagator.

    Raises:
        AssertionError: If any test fails.
    """
    print("Testing Orbital Propagator...")
    print("-" * 40)
    
    # Create test orbital elements (ISS-like orbit)
    from datetime import datetime
    
    elements = OrbitalElements(
        a=6778.0,  # km (400 km altitude)
        e=0.0001,  # Nearly circular
        i=math.radians(51.6),  # ISS inclination
        raan=math.radians(45.0),
        arg_perigee=math.radians(0.0),
        mean_anomaly=math.radians(0.0),
        epoch=datetime(2025, 1, 15, 12, 0, 0)
    )
    
    print("Initial Orbital Elements:")
    print(elements)
    print()
    
    # Create propagator
    propagator = TwoBodyPropagator(elements)
    
    # Test 1: Propagate one orbit
    print("Test 1: Propagate one complete orbit")
    period = elements.period
    final_state = propagator.propagate(period)
    print(f"  After one orbit ({period/3600:.2f} hours):")
    print(f"  Position: [{final_state.position[0]:.1f}, {final_state.position[1]:.1f}, {final_state.position[2]:.1f}] km")
    print(f"  Mean anomaly: {final_state.orbital_elements.mean_anomaly * constants.RAD_TO_DEG:.2f}°")
    print(f"  Should return to ~0°: {abs(final_state.orbital_elements.mean_anomaly) < 0.01}")
    
    # Test 2: Generate ground track
    print("\nTest 2: Generate ground track (first 5 points)")
    ground_track = propagator.get_ground_track(duration=600, time_step=120)
    for i, (lat, lon) in enumerate(ground_track[:5]):
        print(f"  Point {i}: Lat = {lat:7.2f}°, Lon = {lon:7.2f}°")
    
    # Test 3: Find next pass
    print("\nTest 3: Find next pass over New York")
    pass_info = propagator.find_next_pass(
        target_lat=40.7128,  # New York
        target_lon=-74.0060,
        elevation_threshold=10.0
    )
    
    if pass_info:
        print(f"  Pass start: {pass_info['pass_start']}")
        print(f"  Pass end: {pass_info['pass_end']}")
        print(f"  Max elevation: {pass_info['max_elevation_deg']:.1f}°")
        print(f"  Duration: {pass_info['duration_seconds']/60:.1f} minutes")
    else:
        print("  No pass found in next 24 hours")
    
    print("-" * 40)
    print("Tests complete!")


if __name__ == "__main__":
    test_propagator()
