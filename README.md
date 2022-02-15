# DFT-Turbomole

## Description

This Wano plots spectra (IR, UV/Vis) based on the output of DFT calculations for one or several molecules.

## Server setup

The code is based on python and the necessary virtual environment on the server is provided by Simstack.

## Required input

The required input consists DFT spectrum calculations in the ```results.yml``` format of the [DFT-Turbomole WaNo](https://github.com/KIT-Workflows/DFT-Turbomole) (results of other codes in the same format work as well).

## WaNo Settings

Choose which spectra to plot and the corresponding yml files (the calculated transitions must be present in the yml file).

## Output

The output of this WaNo consists of (if chosen) one IR and one U/Vis spectrum (absorption and emission) per molecule.
