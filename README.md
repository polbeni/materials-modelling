# Materials modelling processing code

Compendium of scripts for preprocessing and postprocessing data from materials modelling simulations.
Here can be found the different scripts that I generate for materials simulations of my PhD research. 

The scripts and available content: 
- **Bash scripts**: basic bash scripts to accelerate some routine processes.
- **VASP INCAR files**: files with the necessary tags to perform different calculations with [VASP](https://www.vasp.at/) software. The INCAR files that can be found can be used to:
   - Compute energies and forces.
   - Perform ionic relaxations.
   - Determine phonon frequencies in $\Gamma$-point.
   - Perform ab-initio molecular dynamics (aimd).
   - Compute the dielectric constant
   - Compute the elastic tensor.
   - Compute the dielectric tensor for optical properties.
   - Determine the electronic structure.
   - Compute the energy bands (with and without hybrid functionals).
- **Python scripts**: assortment of python codes for many different things. With these scripts is possible to:
   - Generate KPOINTS and INCAR file for a study of energy convergence.
   - Get the optical and electronic (band gap) properties.
   - Plot phonon frequencies and phonon density of states from [phonopy](https://phonopy.github.io/phonopy/) results.
   - Get thermal properties from [phonopy](https://phonopy.github.io/phonopy/) results and study diffusion from molecular dynamics simulations.
   - Plot electronic density of states (eDOS) with and without thermal correction.
   - Plot energy bands.
   - Compute bond lengths.
   - Functionalities to send to VASP crystal structures found with [PyMCSP](https://github.com/polbeni/PyMCSP).
   - Get the resulting VASP energies from the determined structures with [PyMCSP](https://github.com/polbeni/PyMCSP).
   - Use [Materials Project](https://next-gen.materialsproject.org/) API to download data from their dataset.
   - Use [M3GNet](https://github.com/materialsvirtuallab/m3gnet).
   - Use different functionalities of [hiPhive](https://hiphive.materialsmodeling.org/).
   - Generate POSCAR files from a XDATCAR file.
   - Study the movement of ions from molecular dynamics simulations.
   - Change the atomic species from POSCAR files.
   - Find new possible structure reached with molecular dynamics simulations.


## Disclaimer

These codes are not intended for a general use or tool to help others researchers, but for me to have fast access to the codes. However, feel free to use them if you think that they could be useful for your research:)

## Authors

This code and repository are being developed by:
- Pol Benítez Colominas (pol.benitez@upc.edu)
