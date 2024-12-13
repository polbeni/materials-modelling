# Pol Benítez Colominas, December 2024
# Universitat Politècnica de Catalunya

# Impose rotational sum rules to a FORCE_CONSTANTS file using hiphive software
# We need to provide the following files: POSCAR, SPOSCAR, FORCE_CONSTANTS
# We need to provide the following parameters:
#       supercell_matrix -> dimension of the supercell
#       cutoffs -> cutoff value for the cluster space, just to second order force constants

from phonopy import Phonopy
from phonopy.file_IO import write_force_constants_to_hdf5, parse_FORCE_CONSTANTS
from phonopy.interface.vasp import read_vasp

from hiphive import ClusterSpace, ForceConstants, ForceConstantPotential, enforce_rotational_sum_rules
from hiphive.utilities import extract_parameters

from ase.io import read


### Convert FORCE_CONSTANTS to hdf5 format
unitcell = read_vasp('POSCAR')
supercell_matrix = [[5, 0, 0], [0, 5, 0], [0, 0, 1]] 

phonon = Phonopy(unitcell, supercell_matrix)
force_constants = parse_FORCE_CONSTANTS('FORCE_CONSTANTS') 
phonon.force_constants = force_constants

write_force_constants_to_hdf5(phonon.force_constants, filename='fcs_not_corrected.hdf5')

### Create hiphive cluster space
prim = read('POSCAR')

cutoffs = [8.0]
cs = ClusterSpace(prim, cutoffs)

### Read the force constants
supercell = read('SPOSCAR')

fc = ForceConstants.read_phonopy(supercell, 'fcs_not_corrected.hdf5')

### Extract the parameters
parameters = extract_parameters(fc, cs)

### Determine the parameters when rotational sum rules are imposed
parameters_huang = enforce_rotational_sum_rules(cs, parameters, ['Huang']) # Huang invariance
parameters_bornhuang = enforce_rotational_sum_rules(cs, parameters, ['Born-Huang']) # Born+Huang invariance
parameters_rot = enforce_rotational_sum_rules(cs, parameters, ['Huang', 'Born-Huang']) # Huang and Born+Huang invariance
#! Here we are using the last correction

# Generate the corrected force constants and save it in the hdf5 format
fcp_rot = ForceConstantPotential(cs, parameters_rot)
fcs = fcp_rot.get_force_constants(supercell)

fcs.write_to_phonopy('force_constants.hdf5')
