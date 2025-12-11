"""
Input parser for various orbital data formats.
Primarily handles Two-Line Element (TLE) format from NORAD/Celestrak.
"""

import re
from datetime import datetime, timedelta
from typing import List, Tuple
from orbital_elements import OrbitalElements
import constants


class InputParser:
    """Parser for various orbital data input formats."""
    
    def __init__(self):
        """Initialize the input parser."""
        self.tle_pattern = re.compile(r'^[12] ')  # TLE lines start with 1 or 2
        
    def parse_file(self, filename: str) -> OrbitalElements:
        """
        Parse an input file and return orbital elements.
        
        Args:
            filename: Path to input file
            
        Returns:
            OrbitalElements: Parsed orbital elements
            
        Raises:
            FileNotFoundError: If file doesn't exist
            ValueError: If file format is invalid
        """
        with open(filename, 'r') as f:
            lines = f.readlines()
        
        # Clean lines (strip whitespace)
        lines = list(filter(None, map(str.strip, lines)))

        
        # Try to detect format
        if self._is_tle_format(lines):
            return self.parse_tle(lines)
        else:
            # Try to parse as simple orbital elements format
            return self.parse_simple_format(lines)
    
    def _is_tle_format(self, lines: List[str]) -> bool:
        """
        Check if lines appear to be in TLE format.
        
        Args:
            lines: List of file lines
            
        Returns:
            bool: True if TLE format detected
        """
        # TLE has either 2 lines (no name) or 3 lines (with name)
        if len(lines) < 2:
            return False
        
        # Check if we have lines starting with "1 " and "2 "
        return {'1 ', '2 '} <= {line[:2] for line in lines}

    
    def parse_tle(self, lines: List[str]) -> OrbitalElements:
        """
        Parse Two-Line Element (TLE) format.
        
        TLE Format:
        Line 0 (optional): Satellite name
        Line 1: Contains epoch, ballistic coefficient, etc.
        Line 2: Contains orbital elements
        
        Args:
            lines: TLE lines (2 or 3 lines)
            
        Returns:
            OrbitalElements: Parsed orbital elements
            
        Raises:
            ValueError: If TLE format is invalid
        """
        # Find the actual TLE lines (starting with 1 and 2)
        line1 = None
        line2 = None
        satellite_name = None
        
        for i, line in enumerate(lines):
            if line.startswith('1 '):
                line1 = line
                # Check if previous line is the satellite name
                if i > 0 and not lines[i-1].startswith(('1 ', '2 ')):
                    satellite_name = lines[i-1]
            elif line.startswith('2 '):
                line2 = line
        
        if not line1 or not line2:
            raise ValueError("Invalid TLE format: missing line 1 or line 2")
        
        # Validate line lengths
        if len(line1) != constants.TLE_LINE_LENGTH or len(line2) != constants.TLE_LINE_LENGTH:
            raise ValueError(f"Invalid TLE line length (expected {constants.TLE_LINE_LENGTH})")
        
        try:
            # Parse Line 1
            # Columns 19-32: Epoch (YYDDD.DDDDDDDD)
            epoch_year = int(line1[18:20])
            epoch_day = float(line1[20:32])
            
            # Determine century
            if epoch_year < constants.TLE_EPOCH_YEAR_CUTOFF:
                epoch_year += 2000
            else:
                epoch_year += 1900
            
            # Convert to datetime
            epoch = datetime(epoch_year, 1, 1) + timedelta(days=epoch_day - 1)
            
            # Parse Line 2 - Orbital Elements
            # Columns 9-16: Inclination (degrees)
            inclination = float(line2[8:16])
            
            # Columns 18-25: RAAN (degrees)
            raan = float(line2[17:25])
            
            # Columns 27-33: Eccentricity (decimal point assumed)
            eccentricity = float('0.' + line2[26:33])
            
            # Columns 35-42: Argument of Perigee (degrees)
            arg_perigee = float(line2[34:42])
            
            # Columns 44-51: Mean Anomaly (degrees)
            mean_anomaly = float(line2[43:51])
            
            # Columns 53-63: Mean Motion (revolutions per day)
            mean_motion_rpd = float(line2[52:63])
            
            # Convert mean motion to rad/s and derive semi-major axis
            mean_motion = mean_motion_rpd * constants.TWO_PI / constants.SECONDS_PER_DAY
            semi_major_axis = (constants.EARTH_MU / (mean_motion ** 2)) ** (1/3)
            
            # Convert angles to radians
            inclination_rad = inclination * constants.DEG_TO_RAD
            raan_rad = raan * constants.DEG_TO_RAD
            arg_perigee_rad = arg_perigee * constants.DEG_TO_RAD
            mean_anomaly_rad = mean_anomaly * constants.DEG_TO_RAD
            
            # Create and return orbital elements
            elements = OrbitalElements(
                a=semi_major_axis,
                e=eccentricity,
                i=inclination_rad,
                raan=raan_rad,
                arg_perigee=arg_perigee_rad,
                mean_anomaly=mean_anomaly_rad,
                epoch=epoch
            )
            
            print(f"Successfully parsed TLE for: {satellite_name or 'Unknown satellite'}")
            return elements
            
        except (ValueError, IndexError) as e:
            raise ValueError(f"Error parsing TLE format: {e}")
    
    def parse_simple_format(self, lines: List[str]) -> OrbitalElements:
        """
        Parse a simple text format with orbital elements.
        
        Expected format (one value per line or key-value pairs):
        semi_major_axis: 7000.0  # km
        eccentricity: 0.001
        inclination: 45.0  # degrees
        raan: 30.0  # degrees
        arg_perigee: 60.0  # degrees
        mean_anomaly: 0.0  # degrees
        epoch: 2025-01-01T00:00:00
        
        Args:
            lines: Lines containing orbital elements
            
        Returns:
            OrbitalElements: Parsed orbital elements
        """
        elements_dict = {}
        
        for line in lines:
            # Skip comments and empty lines
            if line.startswith('#') or not line.strip():
                continue
            
            # Try to parse key-value pairs
            if ':' in line:
                key, value = line.split(':', 1)
                key = key.strip().lower()
                value = value.split('#')[0].strip()  # Remove inline comments
                
                # Map keys to values
                if 'semi' in key or 'axis' in key:
                    elements_dict['a'] = float(value)
                elif 'ecc' in key:
                    elements_dict['e'] = float(value)
                elif 'inc' in key:
                    elements_dict['i'] = float(value) * constants.DEG_TO_RAD
                elif 'raan' in key or 'ascending' in key:
                    elements_dict['raan'] = float(value) * constants.DEG_TO_RAD
                elif 'arg' in key or 'perigee' in key:
                    elements_dict['arg_perigee'] = float(value) * constants.DEG_TO_RAD
                elif 'mean' in key or 'anomaly' in key:
                    elements_dict['mean_anomaly'] = float(value) * constants.DEG_TO_RAD
                elif 'epoch' in key:
                    elements_dict['epoch'] = datetime.fromisoformat(value)
        
        # Set defaults for missing values
        if 'epoch' not in elements_dict:
            elements_dict['epoch'] = datetime.now()
        
        # Check for required elements
        required = ['a', 'e', 'i', 'raan', 'arg_perigee', 'mean_anomaly']
        missing = [k for k in required if k not in elements_dict]
        
        if missing:
            raise ValueError(f"Missing required orbital elements: {missing}")
        
        return OrbitalElements(**elements_dict)
    
    def parse_observation_data(self, observations: List[Tuple[datetime, float, float, float]]) -> OrbitalElements:
        """
        Parse ground-based observation data to determine orbital elements.
        
        Args:
            observations: List of (time, azimuth, elevation, range) tuples
            
        Returns:
            OrbitalElements: Estimated orbital elements
            
        Note:
            This is a placeholder for orbit determination from observations.
            Full implementation would require methods like Gauss or Laplace.
        """
        print("Warning: Orbit determination from observations not yet implemented")
        print(f"Received {len(observations)} observations")
        
        # Return dummy orbital elements for now
        return OrbitalElements(
            a=7000.0,
            e=0.01,
            i=45.0 * constants.DEG_TO_RAD,
            raan=0.0,
            arg_perigee=0.0,
            mean_anomaly=0.0,
            epoch=observations[0][0] if observations else datetime.now()
        )


def create_sample_tle_file(filename: str = "sample_tle.txt"):
    """
    Create a sample TLE file for testing.
    
    Args:
        filename: Output filename
    """
    tle_content = """ISS (ZARYA)
1 25544U 98067A   25015.50000000  .00003456  00000-0  12345-3 0  9999
2 25544  51.6439  123.4567  0001234  234.5678  125.4321  15.50123456789012"""
    
    with open(filename, 'w') as f:
        f.write(tle_content)
    
    print(f"Sample TLE file created: {filename}")


def create_sample_simple_file(filename: str = "sample_elements.txt"):
    """
    Create a sample simple format file for testing.
    
    Args:
        filename: Output filename
    """
    content = """# Sample orbital elements file
semi_major_axis: 7000.0  # km
eccentricity: 0.001
inclination: 45.0  # degrees
raan: 30.0  # degrees
arg_perigee: 60.0  # degrees  
mean_anomaly: 0.0  # degrees
epoch: 2025-01-15T12:00:00"""
    
    with open(filename, 'w') as f:
        f.write(content)
    
    print(f"Sample orbital elements file created: {filename}")
