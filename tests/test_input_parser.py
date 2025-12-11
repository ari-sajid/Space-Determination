"""
Tests for input_parser module.
"""

import pytest
import os
import tempfile
from datetime import datetime
from input_parser import InputParser
from orbital_elements import OrbitalElements


class TestInputParser:
    """Test cases for InputParser class."""

    def test_parser_initialization(self):
        """Test creating an input parser."""
        parser = InputParser()
        assert parser is not None

    def test_parse_valid_tle(self):
        """Test parsing a valid TLE file."""
        # Create a temporary TLE file
        tle_content = """ISS (ZARYA)
            1 25544U 98067A   08264.51782528 -.00002182  00000-0 -11606-4 0  2927
            2 25544  51.6416 247.4627 0006703 130.5360 325.0288 15.72125391563537"""

        with tempfile.NamedTemporaryFile(mode='w', delete=False, suffix='.txt') as f:
            f.write(tle_content)
            temp_file = f.name

        try:
            parser = InputParser()
            elements = parser.parse_file(temp_file)

            assert isinstance(elements, OrbitalElements)
            assert elements.a > 6000  # Reasonable LEO altitude
            assert elements.e >= 0 and elements.e < 1  # Valid eccentricity
            assert elements.i > 0  # ISS has inclination ~51.6 degrees
        finally:
            os.unlink(temp_file)

    def test_parse_tle_without_name(self):
        """Test parsing TLE without satellite name line."""
        tle_content = """1 25544U 98067A   08264.51782528 -.00002182  00000-0 -11606-4 0  2927
            2 25544  51.6416 247.4627 0006703 130.5360 325.0288 15.72125391563537"""

        with tempfile.NamedTemporaryFile(mode='w', delete=False, suffix='.txt') as f:
            f.write(tle_content)
            temp_file = f.name

        try:
            parser = InputParser()
            elements = parser.parse_file(temp_file)
            assert isinstance(elements, OrbitalElements)
        finally:
            os.unlink(temp_file)

    def test_parse_invalid_tle_short_lines(self):
        """Test that invalid TLE (short lines) raises error."""
        tle_content = """ISS
            1 25544U
            2 25544"""

        with tempfile.NamedTemporaryFile(mode='w', delete=False, suffix='.txt') as f:
            f.write(tle_content)
            temp_file = f.name

        try:
            parser = InputParser()
            with pytest.raises(ValueError):
                parser.parse_file(temp_file)
        finally:
            os.unlink(temp_file)

    def test_parse_empty_file(self):
        """Test that empty file raises error."""
        with tempfile.NamedTemporaryFile(mode='w', delete=False, suffix='.txt') as f:
            f.write("")
            temp_file = f.name

        try:
            parser = InputParser()
            with pytest.raises(ValueError):
                parser.parse_file(temp_file)
        finally:
            os.unlink(temp_file)

    def test_parse_nonexistent_file(self):
        """Test that nonexistent file raises FileNotFoundError."""
        parser = InputParser()
        with pytest.raises(FileNotFoundError):
            parser.parse_file("nonexistent_file.txt")

    def test_tle_checksum_validation(self):
        """Test TLE checksum validation."""
        parser = InputParser()

        # Valid TLE line
        line1 = "1 25544U 98067A   08264.51782528 -.00002182  00000-0 -11606-4 0  2927"
        assert parser.validate_tle_checksum(line1) == True

        # Invalid checksum (last digit changed)
        line1_invalid = "1 25544U 98067A   08264.51782528 -.00002182  00000-0 -11606-4 0  2920"
        assert parser.validate_tle_checksum(line1_invalid) == False

    def test_epoch_conversion(self):
        """Test TLE epoch conversion to datetime."""
        parser = InputParser()

        # Year 2008, day 264.51782528
        epoch_str = "08264.51782528"
        epoch = parser.parse_tle_epoch(epoch_str)

        assert isinstance(epoch, datetime)
        assert epoch.year == 2008
        # Day 264 is in September
        assert epoch.month == 9

    def test_epoch_year_2000_boundary(self):
        """Test epoch conversion around year 2000."""
        parser = InputParser()

        # Year 1999 (99)
        epoch_99 = parser.parse_tle_epoch("99001.0")
        assert epoch_99.year == 1999

        # Year 2000 (00)
        epoch_00 = parser.parse_tle_epoch("00001.0")
        assert epoch_00.year == 2000

        # Year 2057 (cutoff year)
        epoch_57 = parser.parse_tle_epoch("57001.0")
        assert epoch_57.year == 2057

        # Year 1958 (should be 2058)
        epoch_58 = parser.parse_tle_epoch("58001.0")
        assert epoch_58.year == 1958

    def test_mean_motion_to_semi_major_axis(self):
        """Test conversion from mean motion to semi-major axis."""
        parser = InputParser()

        # ISS mean motion ~ 15.7 revs/day
        mean_motion_rev_per_day = 15.7

        a = parser.mean_motion_to_semi_major_axis(mean_motion_rev_per_day)

        # ISS semi-major axis should be around 6770-6800 km
        assert 6700 < a < 6900

    def test_simple_format_parsing(self):
        """Test parsing simple key-value format."""
        simple_content = """a: 7000.0
            e: 0.1
            i: 51.6
            raan: 10.0
            arg_perigee: 20.0
            mean_anomaly: 30.0"""

        with tempfile.NamedTemporaryFile(mode='w', delete=False, suffix='.txt') as f:
            f.write(simple_content)
            temp_file = f.name

        try:
            parser = InputParser()
            elements = parser.parse_file(temp_file)

            assert isinstance(elements, OrbitalElements)
            assert elements.a == 7000.0
            assert abs(elements.e - 0.1) < 1e-6
        finally:
            os.unlink(temp_file)

    def test_hubble_tle_parsing(self):
        """Test parsing Hubble Space Telescope TLE."""
        tle_content = """HUBBLE SPACE TELESCOPE
            1 20580U 90037B   23364.50000000  .00001234  00000-0  12345-4 0  9999
            2 20580  28.5000  45.0000 0001000  90.0000 270.0000 15.00000000123456"""

        with tempfile.NamedTemporaryFile(mode='w', delete=False, suffix='.txt') as f:
            f.write(tle_content)
            temp_file = f.name

        try:
            parser = InputParser()
            elements = parser.parse_file(temp_file)

            # Hubble is in higher orbit than ISS
            assert elements.a > 6900  # Higher altitude
            # Hubble has lower inclination
            assert elements.i < 1.0  # Less than ~57 degrees
        finally:
            os.unlink(temp_file)