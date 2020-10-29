# LE4PD-XYZ
Codes and data for Frontiers in Molecular Biosciences article on LE4PD-XYZ model

# Dependencies
The following "non-standard" Python libraries MUST be installed for this code to work correctly:
* physt (https://github.com/janpipek/physt)

# Installation
Run the 'install.sh' file using the shell to install.

# Requirements
To run the analysis, all that is required is a trajectory file in .g96 format and a .pdb file with the structure of the protein of interest. The .g96 file should have been processed to remove rigid body rotation and translational motions. There is a 'process.sh' file included in this repository that will perform the necessary corrections to the trajectory when run. 

These codes are set up to analyze MD trajectories of proteins generated using GROMACS (http://www.gromacs.org/). If you have a trajectory generated from another MD simulation engine, you can use a piece of software such as Open Babel (http://openbabel.org) to convert to .g96 format or 2) write your own code to convert to .g96 format. Future implementations might have to ability to accept other types of MD trajectory files.

# Examples
Examples for using LE4PD are located in LE4PD/examples. Right now there is only a single notebook showing how to run the analysis for a 1-ns trajectory of a monomer of the HIV-1 protease (from PDB file 1EBW). There is an example in the 'ensemble' directory, but it is a dead end (the codes are not there to support it).

# Computational Performance

Running the full LE4PD analysis (including save the model to file) takes about 1 minute per 50 000 frames of trajectory data and scales approximately linearly up to 1 500 000 frames of trajectory data. These benchmarking data were generated running the code on the Comet supercomputer at the San Diego Supercomputing Center (https://www.sdsc.edu/support/user_guides/comet.html). 

# Caveats
This implementation of the LE4PD code is still a work in progress. Please report any issues that arise with usage of the codes, and hoepfully someone will respond in a reasonable amount of time.
