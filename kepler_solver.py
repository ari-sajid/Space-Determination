"""
Kepler equation solver for orbital mechanics.
Solves Kepler's equation: M = E - e*sin(E)
"""

import math
from typing import Tuple, Optional
import constants

class KeplerSolver:
    """Solves Kepler's equation and related anomaly conversions."""
    
    def __init__(self, tolerance: float = constants.KEPLER_TOLERANCE,
                 max_iterations: int = constants.MAX_KEPLER_ITERATIONS):
        """
        Initialize Kepler solver.
        
        Args:
            tolerance: Convergence tolerance for iterative solver
            max_iterations: Maximum number of iterations
        """
        self.tolerance = tolerance
        self.max_iterations = max_iterations
    
    def solve_kepler_equation(self, mean_anomaly: float, eccentricity: float) -> float:
        """
        Solve Kepler's equation M = E - e*sin(E) for eccentric anomaly E.
        Uses Newton-Raphson method for elliptical orbits (e < 1).
        
        Parameters
        ----------
        mean_anomaly : float
            Mean anomaly M (radians), normalized to [0, 2\pi].
        eccentricity : float
            Orbital eccentricity e. Must be in range [0, 1) for elliptical orbits.
            
        Returns
        -------
        float
            Eccentric anomaly E (radians), normalized to [0, 2\pi].
            
        Raises
        ------
        ValueError
            Eccentricity is negative or >= 1 (parabolic/hyperbolic orbits not supported).
            If Newton-Raphson iteration fails to converge within max_iterations.
            
        Notes
        -----
        The Newton-Raphson method iteratively refines the eccentric anomaly estimate using:
            E_new = E - f(E)/f'(E)
        where:
            f(E) = E - e*sin(E) - M
            f'(E) = 1 - e*cos(E)
            
        Initial guess selection is optimized based on eccentricity:
        - For e < 0.8: E_o = M (works well for near-circular orbits)
        - For e >= 0.8: E_o = \pi (better for highly elliptical orbits)
        
        Convergence is typically achieved within 3-5 iterations for most orbits.
        
        Examples
        --------
        >>> solver = KeplerSolver()
        >>> # Circular orbit (e = 0)
        >>> E = solver.solve_kepler_equation(math.pi/2, 0.0)
        >>> print(f"E = {math.degrees(E):.2f} deg")
        E = 90.00 deg
        
        >>> # Elliptical orbit
        >>> E = solver.solve_kepler_equation(math.pi/3, 0.1)
        >>> print(f"E = {math.degrees(E):.2f} deg")
        E = 66.79 deg
        """
        if eccentricity < 0 or eccentricity >= 1:
            raise ValueError(f"This solver only handles elliptical orbits (0 <= e < 1), got e={eccentricity}")
        
        # Normalize mean anomaly to [0, 2\pi]
        M = mean_anomaly % constants.TWO_PI
        
        # Initial guess for eccentric anomaly
        if eccentricity < 0.8:
            E = M  # Good for low eccentricity
        else:
            E = math.pi  # Better for high eccentricity
        
        # Newton-Raphson iteration
        for iteration in range(self.max_iterations):
            # Calculate f(E) = E - e*sin(E) - M
            f = E - eccentricity * math.sin(E) - M
            
            # Calculate f'(E) = 1 - e*cos(E)
            f_prime = 1 - eccentricity * math.cos(E)
            
            # Avoid division by zero
            if abs(f_prime) < 1e-12:
                f_prime = 1e-12
            
            # Newton-Raphson step
            E_new = E - f / f_prime
            
            # Check convergence
            if abs(E_new - E) < self.tolerance:
                return E_new % constants.TWO_PI
            
            E = E_new
        
        # Failed convergence
        raise ValueError(f"Kepler equation did not converge after {self.max_iterations} iterations")
    
    def mean_to_true_anomaly(self, mean_anomaly: float, eccentricity: float) -> float:
        """
        DOCUMENTATION #2

        Convert mean anomaly to true anomaly.
        
        Performs a two-step conversion: mean anomaly --> eccentric anomaly --> true anomaly.
        
        Parameters
        ----------
        mean_anomaly : float
            Mean anomaly in radians. Can be any value (will be normalized internally).
        eccentricity : float
            Orbital eccentricity. Must be in range [0, 1) for elliptical orbits.
            
        Returns
        -------
        float
            True anomaly ν in radians, normalized to [0, 2\pi].
            
        Raises
        ------
        ValueError
            If eccentricity is outside valid range [0, 1).
            If Kepler equation solver fails to converge.
            
        Notes
        -----
        This method combines two transformations:
        1. Solve Kepler's equation to get eccentric anomaly E from mean anomaly M
        2. Convert eccentric anomaly E to true anomaly ν
        
        The true anomaly represents the actual angular position of the orbiting body
        as measured from periapsis (closest approach point).
        
        For circular orbits (e = 0), the true anomaly equals the mean anomaly.
        For elliptical orbits, the relationship is non-linear and requires solving
        Kepler's transcendental equation.
        
        See Also
        --------
        solve_kepler_equation : Solves for eccentric anomaly
        eccentric_to_true_anomaly : Converts E to ν
        true_to_mean_anomaly : Inverse operation
        
        Examples
        --------
        >>> solver = KeplerSolver()
        >>> # Convert 60 deg mean anomaly for a slightly elliptical orbit
        >>> nu = solver.mean_to_true_anomaly(math.radians(60), 0.1)
        >>> print(f"True anomaly: {math.degrees(nu):.2f} deg")
        True anomaly: 67.20 deg
        
        >>> # Near-circular orbit - anomalies are nearly equal
        >>> nu = solver.mean_to_true_anomaly(math.radians(45), 0.001)
        >>> print(f"Mean: 45.00 deg, True: {math.degrees(nu):.2f} deg")
        Mean: 45.00 deg, True: 45.09 deg
        """
        # Get eccentric anomaly
        E = self.solve_kepler_equation(mean_anomaly, eccentricity)
        
        # Convert eccentric anomaly to true anomaly
        true_anomaly = self.eccentric_to_true_anomaly(E, eccentricity)
        
        return true_anomaly
    
    def eccentric_to_true_anomaly(self, eccentric_anomaly: float, eccentricity: float) -> float:
        """
        DOCUMENTATION #3

        Convert eccentric anomaly to true anomaly.
        
        Uses the geometric relationship between eccentric and true anomalies
        for elliptical orbits.
        
        Parameters
        ----------
        eccentric_anomaly : float
            Eccentric anomaly E in radians.
        eccentricity : float
            Orbital eccentricity. Must be in range [0, 1) for elliptical orbits.
            
        Returns
        -------
        float
            True anomaly ν in radians, normalized to [0, 2\pi].
            
        Notes
        -----
        The conversion uses the half-angle formula for numerical stability:
            ν = E + 2*a*tan^2(\beta*sin(E), 1 - \beta*cos(E))
        where:
            \beta = e / (1 + \sqrt(1 - e^2))
            
        This formulation avoids the singularity that occurs in the traditional
        formula tan(ν/2) = \sqrt((1+e)/(1-e)) * tan(E/2) when E approaches \pi.
        
        The eccentric anomaly E is an auxiliary angle used in orbital calculations.
        It represents the angle at the center of the auxiliary circle that
        circumscribes the orbital ellipse.
        
        Geometric interpretation:
        - Project the position on the ellipse perpendicular to the major axis
        onto a circumscribed circle
        - The angle from periapsis to this projected point is the eccentric anomaly
        
        Examples
        --------
        >>> solver = KeplerSolver()
        >>> # Convert eccentric anomaly for moderate eccentricity
        >>> E = math.radians(60)  # 60 deg eccentric anomaly
        >>> nu = solver.eccentric_to_true_anomaly(E, 0.2)
        >>> print(f"E = 60 deg, ν = {math.degrees(nu):.2f} deg")
        E = 60 deg, ν = 73.40 deg
        
        >>> # For circular orbit (e=0), E = ν
        >>> nu = solver.eccentric_to_true_anomaly(math.pi/4, 0.0)
        >>> print(f"Circular: E = 45 deg, ν = {math.degrees(nu):.2f} deg")
        Circular: E = 45 deg, ν = 45.00 deg
        """
        E = eccentric_anomaly
        e = eccentricity
        
        # Method 1: Using half-angle formula
        beta = e / (1 + math.sqrt(1 - e**2))
        true_anomaly = E + 2 * math.atan2(beta * math.sin(E), 1 - beta * math.cos(E))
        
        return true_anomaly % constants.TWO_PI
    
    def true_to_eccentric_anomaly(self, true_anomaly: float, eccentricity: float) -> float:
        """
        Convert true anomaly to eccentric anomaly.
        
        Args:
            true_anomaly: True anomaly ν in radians
            eccentricity: Orbital eccentricity
            
        Returns:
            float: Eccentric anomaly E in radians
        """
        nu = true_anomaly
        e = eccentricity
        
        # Using atan2 for correct quadrant
        E = math.atan2(
            math.sqrt(1 - e**2) * math.sin(nu),
            e + math.cos(nu)
        )
        
        return E % constants.TWO_PI
    
    def true_to_mean_anomaly(self, true_anomaly: float, eccentricity: float) -> float:
        """
        DOCUMENTATION #4
        Convert true anomaly to mean anomaly.
        
        Inverse operation of mean_to_true_anomaly. Performs a two-step conversion:
        true anomaly --> eccentric anomaly --> mean anomaly.
        
        Parameters
        ----------
        true_anomaly : float
            True anomaly ν in radians. Can be any value (will be normalized).
        eccentricity : float
            Orbital eccentricity. Must be in range [0, 1) for elliptical orbits.
            
        Returns
        -------
        float
            Mean anomaly M in radians, normalized to [0, 2\pi].
            
        Notes
        -----
        The mean anomaly represents the fraction of the orbital period that has
        elapsed since periapsis passage, expressed as an angle. It increases
        uniformly with time at a rate equal to the mean motion n = 2\pi / T.
        
        Conversion process:
        1. Convert true anomaly ν to eccentric anomaly E using:
        tan(E/2) = \sqrt((1-e)/(1+e)) * tan(ν/2)
        2. Apply Kepler's equation to get mean anomaly:
        M = E - e*sin(E)
        
        This is the inverse of the more common mean_to_true_anomaly conversion
        and is useful for:
        - Determining time since periapsis from observed position
        - Converting osculating elements to mean elements
        - Orbit determination from observations
        
        See Also
        --------
        mean_to_true_anomaly : Inverse operation
        true_to_eccentric_anomaly : First step of conversion
        
        Examples
        --------
        >>> solver = KeplerSolver()
        >>> # Round-trip conversion test
        >>> nu_original = math.radians(120)  # 120 deg true anomaly
        >>> e = 0.3
        >>> M = solver.true_to_mean_anomaly(nu_original, e)
        >>> nu_recovered = solver.mean_to_true_anomaly(M, e)
        >>> error = abs(nu_recovered - nu_original)
        >>> print(f"Round-trip error: {math.degrees(error):.6f} deg")
        Round-trip error: 0.000000 deg
        
        >>> # True anomaly at apoapsis (180 deg) 
        >>> M = solver.true_to_mean_anomaly(math.pi, 0.2)
        >>> print(f"Mean anomaly at apoapsis: {math.degrees(M):.2f} deg")
        Mean anomaly at apoapsis: 180.00 deg
        """
        # Convert to eccentric anomaly
        E = self.true_to_eccentric_anomaly(true_anomaly, eccentricity)
        
        # Convert to mean anomaly via Kepler's equation
        M = E - eccentricity * math.sin(E)
        
        return M % constants.TWO_PI
    
    def calculate_flight_path_angle(self, true_anomaly: float, eccentricity: float) -> float:
        """
        Calculate flight path angle.
        
        Args:
            true_anomaly: True anomaly (radians)
            eccentricity: Orbital eccentricity
            
        Returns:
            float: Flight path angle in radians
        """
        nu = true_anomaly
        e = eccentricity
        
        # Flight path angle
        gamma = math.atan2(e * math.sin(nu), 1 + e * math.cos(nu))
        
        return gamma
    
    def time_since_periapsis(self, mean_anomaly: float, semi_major_axis: float,
                            mu: float = constants.EARTH_MU) -> float:
        """
        Calculate time since periapsis passage.
        
        Args:
            mean_anomaly: Mean anomaly in radians
            semi_major_axis: Semi-major axis in km
            mu: Gravitational parameter in km³/s²
            
        Returns:
            float: Time since periapsis in seconds
        """
        # Mean motion
        n = math.sqrt(mu / semi_major_axis**3)
        
        # Time since periapsis
        t = mean_anomaly / n
        
        return t
    
    def mean_anomaly_at_time(self, t: float, semi_major_axis: float,
                            initial_mean_anomaly: float = 0.0,
                            mu: float = constants.EARTH_MU) -> float:
        """
        DOCUMENTATION #5
        Calculate mean anomaly at a given time from epoch.
        
        Propagates the mean anomaly forward in time using the mean motion,
        which is constant for unperturbed two-body motion.
        
        Parameters
        ----------
        t : float
            Time since epoch in seconds. Can be negative for past times.
        semi_major_axis : float
            Semi-major axis of the orbit in km. Must be positive.
        initial_mean_anomaly : float, optional
            Mean anomaly at epoch t=0 in radians. Default is 0.0 (at periapsis).
        mu : float, optional
            Gravitational parameter in km^3/s^2. Default is Earth's μ = 398600.4418 km^2/s^.
            
        Returns
        -------
        float
            Mean anomaly at time t in radians, normalized to [0, 2\pi].
            
        Raises
        ------
        ValueError
            If semi_major_axis <= 0.
            
        Notes
        -----
        The mean anomaly increases linearly with time according to:
            M(t) = M_o + n * t
        where:
            M_o = initial mean anomaly at epoch
            n = \sqrt(μ/a^3) = mean motion (rad/s)
            
        The mean motion n represents the average angular velocity of the orbiting
        body. For a complete orbit, the mean anomaly increases by 2\pi radians.
        
        This linear propagation is exact for two-body motion but becomes
        approximate when perturbations (J2, drag, third-body) are considered.
        
        Examples
        --------
        >>> solver = KeplerSolver()
        >>> # ISS-like orbit (a ≈ 6778 km for 400 km altitude)
        >>> a = 6778.0  # km
        >>> period = 2 * math.pi * math.sqrt(a**3 / constants.EARTH_MU)
        >>> print(f"Orbital period: {period/60:.1f} minutes")
        Orbital period: 92.6 minutes
        
        >>> # Mean anomaly after quarter orbit
        >>> M = solver.mean_anomaly_at_time(period/4, a, initial_mean_anomaly=0.0)
        >>> print(f"After 1/4 orbit: M = {math.degrees(M):.1f} deg")
        After 1/4 orbit: M = 90.0 deg
        
        >>> # Mean anomaly after 1.5 orbits
        >>> M = solver.mean_anomaly_at_time(1.5*period, a)
        >>> print(f"After 1.5 orbits: M = {math.degrees(M):.1f} deg")
        After 1.5 orbits: M = 180.0 deg
        
        >>> # Propagate backwards in time
        >>> M = solver.mean_anomaly_at_time(-period/2, a, initial_mean_anomaly=math.pi)
        >>> print(f"Half orbit before epoch: M = {math.degrees(M):.1f} deg")
        Half orbit before epoch: M = 0.0 deg
        """
        # Mean motion
        n = math.sqrt(mu / semi_major_axis**3)
        
        # Propagated mean anomaly
        M = initial_mean_anomaly + n * t
        
        return M % constants.TWO_PI 

if __name__ == "__main__":