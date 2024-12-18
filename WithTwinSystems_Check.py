# -*- coding: utf-8 -*-
"""
Created on Sun Dec 15 14:40:08 2024

@author: repha
"""

import math
import numpy as np
import pandas as pd
from tkinter import Tk, filedialog, simpledialog


class Constants:
    """Holds all constants for calculations."""
    def __init__(self, ca_ratio, a_lattice, c_lattice):
        self.ca_ratio = ca_ratio
        self.a = a_lattice
        self.c = c_lattice
        self.d = math.sqrt(1 + 0 + self.c ** 2)  # Example derived constant


class FileManager:
    """Handles file operations."""
    @staticmethod
    def select_file(prompt):
        """Prompt user to select a file."""
        root = Tk()
        root.withdraw()
        root.attributes("-topmost", True)
        filepath = filedialog.askopenfilename(title=prompt)
        root.destroy()
        return filepath

    @staticmethod
    def read_files():
        """Select and read required CSV files."""
        grain_file = FileManager.select_file("Select the Grain Data CSV")
        band_file = FileManager.select_file("Select the Band Data CSV")
        mapping_file = FileManager.select_file("Select the Band-Grain Mapping CSV")
        
        grain_data = pd.read_csv(grain_file)
        band_data = pd.read_csv(band_file)
        mapping_data = pd.read_csv(mapping_file)

        return grain_data, band_data, mapping_data


def get_user_input():
    """Prompt user for lattice parameters."""
    root = Tk()
    root.withdraw()
    root.attributes("-topmost", True)

    ca_ratio = simpledialog.askfloat("Input", "Enter c/a ratio for the HCP lattice:")
    a_lattice = simpledialog.askfloat("Input", "Enter lattice parameter 'a':")
    c_lattice = ca_ratio * a_lattice

    root.destroy()
    return ca_ratio, a_lattice, c_lattice


class OrientationAnalysis:
    """Handles orientation, trace matching, and shear/burgers calculations."""
    def __init__(self, constants):
        self.constants = constants
        self.slip_systems = self.generate_slip_systems()
        self.twin_systems = self.generate_twin_systems()

    def generate_slip_systems(self):
        """Generate theoretical slip systems for HCP."""
        # Format: (Plane Normal, Slip Direction)
        return [
            ((0, 0, 0, 1), (1, 0, -1, 0)),  # Basal
            ((1, 0, -1, 0), (1, -1, 0, 0)),  # Prismatic
            ((1, 0, -1, 1), (1, 0, -1, 0))   # Pyramidal
        ]

    def generate_twin_systems(self):
        """Generate theoretical twin systems for HCP."""
        # Format: (Twin Plane Normal, Shear Direction)
        return [
            ((1, 0, -1, 2), (0, 1, 0, 0)),  # Tensile Twin
            ((1, 0, -1, -2), (0, -1, 0, 0))  # Compression Twin
        ]

    def calculate_rotation_matrix(self, euler_angles):
        """Calculate rotation matrix from Euler angles."""
        phi1, phi, phi2 = [math.radians(angle) for angle in euler_angles]

        return [
            [math.cos(phi1) * math.cos(phi2) - math.sin(phi1) * math.cos(phi) * math.sin(phi2),
             -math.cos(phi1) * math.sin(phi2) - math.sin(phi1) * math.cos(phi) * math.cos(phi2),
             math.sin(phi1) * math.sin(phi)],
            [math.sin(phi1) * math.cos(phi2) + math.cos(phi1) * math.cos(phi) * math.sin(phi2),
             -math.sin(phi1) * math.sin(phi2) + math.cos(phi1) * math.cos(phi) * math.cos(phi2),
             -math.cos(phi1) * math.sin(phi)],
            [math.sin(phi2) * math.sin(phi),
             math.sin(phi) * math.cos(phi2),
             math.cos(phi)]
        ]

    def calculate_trace_angle(self, normal, rotation_matrix):
        """Calculate trace angle of a system."""
        rotated_normal = np.dot(rotation_matrix, normal[:3])
        trace_angle = math.degrees(math.atan2(rotated_normal[1], rotated_normal[0]))
        return trace_angle % 360

    def calculate_burgers_vector(self, slip_direction):
        """Calculate Burgers vector for slip systems."""
        b_x = slip_direction[0] * self.constants.a
        b_y = slip_direction[1] * self.constants.a
        b_z = slip_direction[2] * self.constants.c
        b_mag = math.sqrt(b_x ** 2 + b_y ** 2 + b_z ** 2)
        return b_x, b_y, b_z, b_mag

    def calculate_twinning_shear(self, twin_plane_normal, shear_direction):
        """Calculate twinning shear for twin systems."""
        # Example formula: gamma = (c/a ratio) * cos(angle between plane and shear direction)
        return self.constants.ca_ratio * abs(np.dot(twin_plane_normal, shear_direction))

    def match_trace_angle(self, measured_angle, rotation_matrix):
        """Match measured trace angle with theoretical systems."""
        results = []
        for i, (normal, direction) in enumerate(self.slip_systems):
            angle = self.calculate_trace_angle(normal, rotation_matrix)
            b_x, b_y, b_z, b_mag = self.calculate_burgers_vector(direction)
            angle_diff = abs(angle - measured_angle)
            results.append(("Slip", i + 1, angle, angle_diff, b_x, b_y, b_z, b_mag))
        
        for i, (normal, direction) in enumerate(self.twin_systems):
            angle = self.calculate_trace_angle(normal, rotation_matrix)
            shear = self.calculate_twinning_shear(normal, direction)
            angle_diff = abs(angle - measured_angle)
            results.append(("Twin", i + 1, angle, angle_diff, shear, None, None, None))

        return sorted(results, key=lambda x: x[3])  # Sort by angle difference


def main():
    # Get user inputs
    ca_ratio, a_lattice, c_lattice = get_user_input()
    constants = Constants(ca_ratio, a_lattice, c_lattice)

    # File operations
    grain_data, band_data, mapping_data = FileManager.read_files()

    analysis = OrientationAnalysis(constants)

    summary_results = []
    detailed_results = []

    for _, band in band_data.iterrows():
        band_id = band['Band ID']
        measured_angle = band['Measured Angle']

        # Find grain associated with the band
        grain_id = mapping_data.loc[mapping_data['Band ID'] == band_id, 'Grain ID'].values[0]
        euler_angles = grain_data.loc[grain_data['Grain ID'] == grain_id, ['Phi1', 'Phi', 'Phi2']].values[0]

        rotation_matrix = analysis.calculate_rotation_matrix(euler_angles)
        matches = analysis.match_trace_angle(measured_angle, rotation_matrix)

        summary_results.append((band_id, grain_id, measured_angle, *matches[0][:4]))
        for match in matches:
            detailed_results.append((band_id, grain_id, measured_angle, *match))

    # Save results
    pd.DataFrame(summary_results, columns=[
        "Band ID", "Grain ID", "Measured Angle", "Matched Type",
        "System ID", "System Angle", "Angle Difference"
    ]).to_csv("summary_results.csv", index=False)

    pd.DataFrame(detailed_results, columns=[
        "Band ID", "Grain ID", "Measured Angle", "Type",
        "System ID", "System Angle", "Angle Difference",
        "Burgers X", "Burgers Y", "Burgers Z", "Burgers Magnitude"
    ]).to_csv("detailed_results.csv", index=False)


if __name__ == "__main__":
    main()
