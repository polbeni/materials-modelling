# Materials modelling processing code

Compendium of scripts for preprocessing and postprocessing data from materials modelling simulations.
Here, you can find the different scripts that I have generated for the materials simulations in my PhD research. 

The scripts and available content: 
- **Bash scripts**: Basic bash scripts to accelerate some routine processes.
- **VASP INCAR files**: Files with the necessary tags to perform different calculations with [VASP](https://www.vasp.at/) software. The INCAR files available can be used to:
   - Compute energies and forces.
   - Perform ionic relaxations.
   - Determine phonon frequencies at the $\Gamma$-point.
   - Perform ab-initio molecular dynamics (AIMD).
   - Compute the dielectric constant
   - Compute the elastic tensor.
   - Compute the dielectric tensor for optical properties.
   - Determine the electronic structure.
   - Compute the energy bands (with and without hybrid functionals).
   - Compute the charge density.
   - Perform band alignment calculations.
   - VCA calculations.
   - DFT calculations within the GW approximation.
- **Quantum ESPRESSO input files**: Input files to perform DFT calculations with [Quantum ESPRESSO](https://www.quantum-espresso.org/) software. The input files available can be used to:
   - Single point energy calculations.
   - Ionic relaxations.
   - Phonon calculations.
   - Molecular dynamics calculations.
   - Determine the electron-phonon matrix elements.
- **Python scripts**: An assortment of Python codes for various purposes. These scripts can:
   - Generate KPOINTS and INCAR files for a study of energy convergence.
   - Generate POTCAR files.
   - Get the optical and electronic (band gap) properties.
   - Calculate effective electron and hole masses using [sumo](https://github.com/SMTG-Bham/sumo).
   - Plot phonon frequencies and phonon density of states from [phonopy](https://phonopy.github.io/phonopy/) results.
   - Plot phonon frequencies from QuantumESPRESSO calculations.
   - Obtain thermal properties from [phonopy](https://phonopy.github.io/phonopy/) results and study diffusion from molecular dynamics simulations.
   - Execute [DynaPhoPy](https://github.com/abelcarreras/DynaPhoPy) for anharmonic phonon renormalization.
   - Plot electronic density of states (eDOS) with and without thermal correction.
   - Plot energy bands.
   - Compute bond lengths.
   - Send VASP crystal structures found with [PyMCSP](https://github.com/polbeni/PyMCSP).
   - Obtain the resulting VASP energies from the determined structures with [PyMCSP](https://github.com/polbeni/PyMCSP).
   - Use [Materials Project](https://next-gen.materialsproject.org/) API to download data from their dataset.
   - Use [M3GNet](https://github.com/materialsvirtuallab/m3gnet) for single-point energy calculation, ionic relaxation, and molecular dynamics simulation.
   - Retrain M3GNet from DFT results.
   - Use various functionalities of [hiPhive](https://hiphive.materialsmodeling.org/).
   - Generate POSCAR files from a XDATCAR file.
   - Study the movement of ions from molecular dynamics simulations.
   - Change the atomic species in POSCAR files.
   - Find new possible structures reached with molecular dynamics simulations.
   - Search in Togo's phonon database [Phonondb](https://github.com/atztogo/phonondb)
   - Band alignments and slab calculations.
   - Generate solid solution structures.
   - Generate VCA structures.
   - Machine Learning methods applied to solid solutions and VCA calculations.
   - Generate catalysis structures.
   - Compute the Fröhlich coupling correction to the band gap for polar materials.
   - Distort an structure given a phonon mode eigenvectors and an amplitude.
   - Compute the dielectric constant for an anisotropic dielectric tensor and compute the thermal average for a given T.
   - Generate distorted structures with uniform or gaussian noise for each atom.
   - Get POSCAR from a phonopy.yaml file.
   - Plot a given band in the full Brillouin zone of a 2D material.
   - Prepare calculations to compute e-p matrix elements with [Quantum ESPRESSO](https://www.quantum-espresso.org/).
   - Cut XDATCAR file.
   - Impose rotational sum rules to phonon force constants.
   - Study phonon mode modulus difference above and below a phonon gap.
   - Generate distorted structures with Monte Carlo, from phonon eigenvectors at a given T.
   - Use [MACE](https://mace-docs.readthedocs.io/en/latest/index.html) for single shot calculations, relaxations and molecular dynamics simulations.
   - Compute phonon frequencies (harmonic and anharmonic) using phonopy (and dynaphopy for anharmonicity) from the forces computed using MACE.
   - Generates a pT-diagram from phonon calculations.
   - Determine if a phonon mode is polar or not.
   - Use phonopy with python scripts.


## Disclaimer

These codes are not intended for general use or as tools to assist other researchers, but rather for my own quick access. However, feel free to use them if you think they could be useful for your research :)

## Authors

This code and repository are being developed by:
- Pol Benítez Colominas (pol.benitez@upc.edu)
