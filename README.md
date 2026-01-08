# MHD Nanofluid Couette Flow - Analytical Validation Suite

## Overview
MATLAB implementation of a Spectral Quasi-Linearization Method (SQLM) for solving Magnetohydrodynamic (MHD) nanofluid Couette flow with comprehensive analytical validation. This repository contains validation codes for various MHD nanofluid flow configurations, including parametric studies of key dimensionless parameters.

## Repository Contents
This repository includes MATLAB codes for:
- MHD nanofluid Couette flow validation (Current main code)
- Parametric studies of various dimensionless numbers (forthcoming)
- Analytical solutions for validation cases
- Visualization tools for publication-quality figures

## Current Files
- mhd_couette_validation.m - Main validation suite for four fundamental cases
- README.md - This documentation file
- LICENSE - MIT License for open academic use

## Quick Start
1. Clone or download this repository:
   git clone https://github.com/MrMosala/mhd-nanofluid-couette-validation.git
2. Open MATLAB and navigate to the project directory
3. Run the main validation script:
   mhd_couette_validation
4. Follow the interactive menu to select validation mode

## Validation Cases Included
Case 1: Simple Couette Flow
- Conditions: Ha = 0, Ec = 0, G = 0
- Analytical solution: Linear velocity profile
- Purpose: Basic validation of numerical scheme

Case 2: MHD Couette Flow
- Conditions: Ha ≠ 0, Ec = 0, G = 0
- Analytical solution: Hyperbolic sine profile
- Purpose: Validation of magnetic field effects

Case 3: Viscous Dissipation
- Conditions: Ha = 0, Ec ≠ 0, G = 0
- Analytical solution: Parabolic temperature profile
- Purpose: Validation of viscous heating effects

Case 4: Pressure Gradient
- Conditions: Ha = 0, Ec = 0, G ≠ 0
- Analytical solution: Quadratic velocity profile
- Purpose: Validation of pressure-driven flow

## Parametric Studies (Forthcoming)
Additional codes will be added for comprehensive parametric analysis:

Dimensionless Parameters Analysis:
1. Reynolds Number (Re) - Flow regime effects
2. Hartmann Number (Ha) - Magnetic field strength effects
3. Eckert Number (Ec) - Viscous dissipation effects
4. Prandtl Number (Pr) - Thermal diffusivity effects
5. Biot Number (Bi) - Convection boundary condition effects
6. Pressure Gradient (G) - Driving force effects
7. Slip Parameter (λ) - Velocity slip effects

Nanoparticle Studies:
1. Volume Fraction (φ) - 0% to 20% variations
2. Nanoparticle Types - Cu, Al₂O₃, comparisons
3. Thermophysical Properties - Temperature-dependent properties

## Fluid Properties
Supported Fluid Types:
- Base Fluid (Water: A₁ = A₂ = A₃ = 1.0)
- Copper (Cu) Nanofluid - φ = 0.01 to 0.10
- Alumina (Al₂O₃) Nanofluid - φ = 0.01 to 0.10

Property Models:
- Viscosity: Brinkman model
- Thermal Conductivity: Maxwell model
- Electrical Conductivity: Maxwell-Garnett model
- Density & Heat Capacity: Simple mixing rule

## Outputs Generated
Automatic Figure Generation:
1. Velocity Analysis (Profiles & Skin Friction)
2. Temperature Analysis (Profiles & Nusselt Number)
3. Analytical Validation (All 4 cases comparison)
4. Validation Summary (Error analysis)
5. Individual Case Figures (Thesis-quality plots)

File Formats:
- .png - For quick viewing
- .pdf - For publications
- .fig - For MATLAB editing

Data Tables:
- Skin Friction Coefficient (C_f) at multiple η locations
- Nusselt Number (Nu) at multiple η locations
- Maximum validation errors for each case
- Parameter sensitivity tables

## Mathematical Formulation
Governing Equations:
Momentum: A₁·W'' - A₂·Ha²·W + G = 0
Energy: A₃·θ'' + A₁·Pr·Ec·(W')² + A₂·Pr·Ec·Ha²·W² = 0

Dimensionless Parameters:
Symbol | Parameter | Range
Re | Reynolds Number | 0.1 - 100
Ha | Hartmann Number | 0 - 10
Ec | Eckert Number | 0 - 1
Pr | Prandtl Number | 0.7 - 7
Bi | Biot Number | 0.1 - 10
λ | Slip Parameter | 0 - 1
G | Pressure Gradient | 0 - 5
φ | Volume Fraction | 0 - 0.2

Boundary Conditions:
W(0) = 0, W(1) = Re - λ·W'(1)
θ(0) = 0, θ'(1) + Bi·θ(1) = 0

## Dependencies
- MATLAB R2018a or newer
- No additional toolboxes required
- Pure MATLAB implementation using native functions

## Usage Examples
Single Case Validation:
>> mhd_couette_validation
Enter choice (1-3): 1
Select case (1-4): 2
Enter Ha: 2.0

Comprehensive Analysis:
>> mhd_couette_validation
Enter choice (1-3): 3

## Sample Validation Output
╔════════════════════════════════════════════════════════════════╗
║                 VALIDATION SUMMARY - ALL CASES                ║
╚════════════════════════════════════════════════════════════════╝

Case 1: Simple Couette (Ha=0, Ec=0, G=0):
  Max velocity error:      3.12e-14 ✓ PASS
  Max temperature error:   4.56e-15 ✓ PASS

Case 2: MHD Couette (Ha=2.0, Ec=0, G=0):
  Max velocity error:      2.89e-13 ✓ PASS
  Max temperature error:   3.45e-14 ✓ PASS

✓ ALL VALIDATION CASES PASSED! 🎉

## Code Structure for Future Additions
mhd-nanofluid-couette-validation/
├── validation/          # Core validation codes
│   ├── mhd_couette_validation.m
│   ├── case1_simple.m
│   └── ...
├── parametric/          # Parametric studies (forthcoming)
│   ├── re_study.m
│   ├── ha_study.m
│   ├── ec_study.m
│   └── ...
├── utilities/           # Helper functions
│   ├── cheb_diff_matrices.m
│   ├── nanofluid_properties.m
│   └── plot_utilities.m
├── data/               # Generated data
├── figures/            # Saved figures
└── docs/               # Documentation

## Academic Context
- Author: Mosala S.I. (mrmosalas.i@gmail.com)
- Supervisor: Prof. O.D. Makinde
- Institution: Nelson Mandela University, South Africa
- Department: Mathematical Sciences
- Sponsor: NITheCS (National Institute for Theoretical and Computational Sciences)
- Research Area: Computational Fluid Dynamics, MHD, Nanofluids

## Methodology
- Numerical Method: Spectral Quasi-Linearization Method (SQLM)
- Spatial Discretization: Chebyshev spectral collocation
- Validation: Exact analytical solutions
- Convergence: Machine precision (10⁻¹² tolerance)
- Grid Points: N = 100 (adjustable)

## References
1. Makinde, O.D. & Anwar Bég, O. (2010). Application of SQLM to MHD flows
2. Sheikholeslami, M. (2017). Nanofluid heat transfer analysis
3. Canuto, C., et al. (2006). Spectral Methods: Fundamentals
4. Buongiorno, J. (2006). Convective transport in nanofluids

## Contact Information
- Primary Email: mrmosalas.i@gmail.com
- Mobile: +27 72 952 0988
- GitHub: @MrMosala

## Collaboration
This research is open for collaboration. Please contact the author if:
- You want to extend the code for different geometries
- You need validation data for your MHD simulations
- You're interested in joint publications
- You have experimental data for comparison

## License
MIT License - See LICENSE file for details.

Permissions:
- Commercial use
- Modification
- Distribution
- Private use

Conditions:
- Copyright notice must be retained
- License text must be included

Limitations:
- No liability
- No warranty

## Related Repositories
- mind-nanofluid-simulation
- mind-fluid-app
- mind-unsteady-app

## Future Development Roadmap
1. Phase 1: Core validation suite (Current)
2. Phase 2: Parametric studies (Re, Ha, Ec, etc.) (Forthcoming)
3. Phase 3: Additional nanofluid types (Ag, TiO₂)
4. Phase 4: Temperature-dependent properties
5. Phase 5: Unsteady flow extension
6. Phase 6: GUI interface for easier use

## Issue Reporting
Found a bug or have a feature request? Please:
1. Check existing issues
2. Create a new issue with:
   - MATLAB version
   - Error message
   - Steps to reproduce
   - Expected vs actual behavior

Last Updated: 2026-01-08
MATLAB Code Version: 1.0.0
Validation Status: All cases pass with machine precision
Repository Status: Active development
