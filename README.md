# g4bnb_ana: A lightweight analysis framework for G4BNB

## OVERVIEW

g4bnb_ana is a small, modular C++ analysis tool designed for processing output from G4BNB, a GEANT4-based simulation of the Booster Neutrino Beamline (BNB) at Fermilab.
It provides a simple interface for reading simulation output files, running custom analysis routines, and generating basic diagnostics or physics-driven summaries. This repo relays on the alg4bnb docker container that returns the output files.

# REPOSITORY STRUCTURE

```
├── CMakeLists.txt (Build system configuration)
├── main.cpp (Program entry point)
├── input.txt (Example configuration / input file list)
├── setup_g41061_detached.sh (Environment script for detached G4/BNB builds)
├── include/ (Header files)
├── src/ (Implementation source files)
└── README.txt (This documentation)
```
# FEATURES

Minimal and easy-to-understand C++17 codebase

Clean separation between header and source files

CMake-based build system

Configurable input file lists (via input.txt)

Compatible with detached GEANT4/G4BNB setups

Easy to extend with ROOT, JSON configs, or custom analysis modules

## REQUIREMENTS

C++17 compiler (GCC, Clang, or MSVC)

CMake version 3.10 or later

GEANT4 installation (matching the version used to generate the simulation outputs)

Optional: ROOT for histogramming or file output

If GEANT4 is not automatically detected, set:
```
-DGeant4_DIR=/path/to/geant4/lib/cmake/Geant4
```
## BUILD INSTRUCTIONS

Clone the repository:

```
git clone https://github.com/marvlad/g4bnb_ana.git

cd g4bnb_ana
```

Build:
```
mkdir build && cd build
cmake ..
make -j4
```
(Optional) Install:

```sudo make install```

If using a detached G4BNB setup, run:

```source setup_g41061_detached.sh```

## RUNNING THE ANALYSIS

Option 1: Using an input file list:

```./g4bnb_ana```



Open a pull request

QUESTIONS

If you need assistance extending the analysis or integrating new components, feel free to ask.
