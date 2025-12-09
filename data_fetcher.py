"""
Live orbital data fetcher for satellite TLE information.
Fetches current TLE data from orbit.ing-now.com and other sources.
"""

import requests
import re
from datetime import datetime
from typing import Optional, Dict, List
from orbital_elements import OrbitalElements
from input_parser import InputParser


class LiveDataFetcher:
    """Fetch live satellite orbital data from various sources."""
    
    def __init__(self):
        """Initialize the data fetcher with common satellites."""
        self.base_url = "https://orbit.ing-now.com"
        self.parser = InputParser()
        
        # Common satellite IDs for quick access
        self.satellite_catalog = {
            'iss': {'id': '25544', 'designator': '1998-067a', 'name': 'iss'},
            'hubble': {'id': '20580', 'designator': '1990-037b', 'name': 'hst'},
            'tiangong': {'id': '48274', 'designator': '2021-035a', 'name': 'css'},
        }
    
    def fetch_tle_from_orbiting_now(self, satellite_name: str = 'iss') -> Optional[str]:
        """
        Fetch TLE data from orbit.ing-now.com for a specific satellite.
        
        Parameters
        ----------
        satellite_name : str
            Name of the satellite ('iss', 'hubble', 'tiangong') or NORAD ID
            
        Returns
        -------
        str or None
            TLE data as a string, or None if fetch fails
            
        Examples
        --------
        >>> fetcher = LiveDataFetcher()
        >>> tle_data = fetcher.fetch_tle_from_orbiting_now('iss')
        >>> print(tle_data)
        ISS (ZARYA)
        1 25544U 98067A   ...
        2 25544  51.6439  ...
        """
        try:
            # Check if it's a known satellite name
            if satellite_name.lower() in self.satellite_catalog:
                sat_info = self.satellite_catalog[satellite_name.lower()]
                url = f"{self.base_url}/satellite/{sat_info['id']}/{sat_info['designator']}/{sat_info['name']}/"
            else:
                # Assume it's a NORAD ID
                url = f"{self.base_url}/satellite/{satellite_name}/"
            
            # Fetch the webpage
            response = requests.get(url, timeout=10)
            response.raise_for_status()
            
            # Extract TLE from the page
            # The website typically has TLE data in a specific format
            tle_data = self._extract_tle_from_html(response.text)
            
            return tle_data
            
        except requests.RequestException as e:
            print(f"Error fetching data from orbit.ing-now.com: {e}")
            return None
    
    def _extract_tle_from_html(self, html_content: str) -> Optional[str]:
        """
        Extract TLE data from the HTML content of orbit.ing-now.com.
        
        Parameters
        ----------
        html_content : str
            HTML content from the webpage
            
        Returns
        -------
        str or None
            Extracted TLE data
        """
        # Look for TLE pattern in the HTML
        # TLE lines have a very specific format
        tle_pattern = r'([A-Z0-9\s\(\)]+)\n(1 \d{5}[A-Z\s]+.{62})\n(2 \d{5}.{63})'
        
        # Search for pre-formatted text blocks that might contain TLE
        pre_pattern = r'<pre[^>]*>(.*?)</pre>'
        pre_matches = re.findall(pre_pattern, html_content, re.DOTALL)
        
        for pre_content in pre_matches:
            # Clean HTML entities
            pre_content = pre_content.replace('&nbsp;', ' ')
            pre_content = re.sub(r'<[^>]+>', '', pre_content)  # Remove any nested tags
            
            # Look for TLE pattern
            tle_match = re.search(tle_pattern, pre_content)
            if tle_match:
                return '\n'.join(tle_match.groups())
        
        # Alternative: Look for TLE lines directly in the page
        line1_pattern = r'1 \d{5}[A-Z\s]+[\d\s\.\-\+]+\s+\d+'
        line2_pattern = r'2 \d{5}[\d\s\.]+\s+\d+'
        
        line1_matches = re.findall(line1_pattern, html_content)
        line2_matches = re.findall(line2_pattern, html_content)
        
        if line1_matches and line2_matches:
            # Try to find satellite name
            name_pattern = r'([A-Z][A-Z0-9\s\(\)]+)(?=\s*1 \d{5})'
            name_match = re.search(name_pattern, html_content)
            
            tle_lines = []
            if name_match:
                tle_lines.append(name_match.group(1).strip())
            tle_lines.append(line1_matches[0].strip())
            tle_lines.append(line2_matches[0].strip())
            
            return '\n'.join(tle_lines)
        
        return None
    
    def fetch_tle_from_celestrak(self, satellite_name: str) -> Optional[str]:
        """
        Fetch TLE data from Celestrak.
        
        Parameters
        ----------
        satellite_name : str
            Name or NORAD ID of the satellite
            
        Returns
        -------
        str or None
            TLE data as string
        """
        celestrak_urls = {
            'iss': 'https://celestrak.org/NORAD/elements/gp.php?CATNR=25544&FORMAT=TLE',
            'stations': 'https://celestrak.org/NORAD/elements/gp.php?GROUP=stations&FORMAT=TLE',
        }
        
        try:
            if satellite_name.lower() in celestrak_urls:
                url = celestrak_urls[satellite_name.lower()]
            elif satellite_name.isdigit():
                # It's a NORAD ID
                url = f'https://celestrak.org/NORAD/elements/gp.php?CATNR={satellite_name}&FORMAT=TLE'
            else:
                # Try as part of active satellites
                url = 'https://celestrak.org/NORAD/elements/gp.php?GROUP=active&FORMAT=TLE'
            
            response = requests.get(url, timeout=10)
            response.raise_for_status()
            
            # For group queries, find the specific satellite
            if 'GROUP=' in url:
                lines = response.text.strip().split('\n')
                for i in range(0, len(lines), 3):
                    if satellite_name.upper() in lines[i].upper():
                        return '\n'.join(lines[i:i+3])
            
            return response.text.strip()
            
        except requests.RequestException as e:
            print(f"Error fetching from Celestrak: {e}")
            return None
    
    def get_orbital_elements(self, satellite_name: str = 'iss', 
                            source: str = 'auto') -> Optional[OrbitalElements]:
        """
        Get orbital elements for a satellite from live sources.
        
        Parameters
        ----------
        satellite_name : str
            Name or NORAD ID of the satellite
        source : str
            Data source ('orbiting-now', 'celestrak', or 'auto')
            
        Returns
        -------
        OrbitalElements or None
            Parsed orbital elements object
            
        Examples
        --------
        >>> fetcher = LiveDataFetcher()
        >>> elements = fetcher.get_orbital_elements('iss')
        >>> print(f"ISS altitude: {elements.a - 6378.137:.1f} km")
        """
        tle_data = None
        
        if source == 'auto' or source == 'orbiting-now':
            print(f"Fetching from orbit.ing-now.com...")
            tle_data = self.fetch_tle_from_orbiting_now(satellite_name)
        
        if not tle_data and (source == 'auto' or source == 'celestrak'):
            print(f"Fetching from Celestrak...")
            tle_data = self.fetch_tle_from_celestrak(satellite_name)
        
        if tle_data:
            # Save to temporary file for parsing
            temp_file = f'temp_{satellite_name}_tle.txt'
            with open(temp_file, 'w') as f:
                f.write(tle_data)
            
            # Parse using existing InputParser
            try:
                elements = self.parser.parse_file(temp_file)
                print(f"Successfully fetched live data for {satellite_name}")
                return elements
            except Exception as e:
                print(f"Error parsing TLE data: {e}")
                return None
        
        print(f"Could not fetch data for {satellite_name}")
        return None
    
    def track_satellite_realtime(self, satellite_name: str = 'iss', 
                                update_interval: int = 60):
        """
        Track a satellite with periodic updates.
        
        Parameters
        ----------
        satellite_name : str
            Name of satellite to track
        update_interval : int
            How often to fetch new data (seconds)
        """
        import time
        from propagator import TwoBodyPropagator
        
        while True:
            try:
                # Get latest orbital elements
                elements = self.get_orbital_elements(satellite_name)
                if elements:
                    # Create propagator
                    propagator = TwoBodyPropagator(elements)
                    
                    # Get current position
                    current_state = propagator.propagate(0)
                    
                    print(f"\n{datetime.now().isoformat()}")
                    print(f"Current position of {satellite_name}:")
                    print(f"  Position: {current_state.position} km")
                    print(f"  Altitude: {np.linalg.norm(current_state.position) - 6378.137:.1f} km")
                    
                    # Calculate ground track
                    ground_track = propagator.get_ground_track(300, 60)
                    if ground_track:
                        lat, lon = ground_track[0]
                        print(f"  Ground position: {lat:.2f}°, {lon:.2f}°")
                
                # Wait before next update
                time.sleep(update_interval)
                
            except KeyboardInterrupt:
                print("\nStopping real-time tracking.")
                break
            except Exception as e:
                print(f"Error in tracking: {e}")
                time.sleep(update_interval)


def test_live_fetcher():
    """Test the live data fetcher."""
    import numpy as np
    
    print("Testing Live Data Fetcher...")
    print("-" * 40)
    
    fetcher = LiveDataFetcher()
    
    # Test 1: Fetch ISS data
    print("Test 1: Fetching ISS orbital data")
    elements = fetcher.get_orbital_elements('iss')
    if elements:
        print(elements)
        altitude = elements.a - 6378.137
        print(f"  Current altitude: {altitude:.1f} km")
        print(f"  Period: {elements.period/60:.1f} minutes")
    
    print("-" * 40)
    print("Tests complete!")


if __name__ == "__main__":
    test_live_fetcher()