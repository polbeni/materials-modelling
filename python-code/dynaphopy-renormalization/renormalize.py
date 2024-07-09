# Pol Benítez Colominas, July 2024
# Universitat Politècnica de Catalunya

# This code uses DynaPhoPy to renormalize harmonic phonon frequencies from the results of a aimd simulation 
# IMPORTANT: you should have in the same directory the results of the aimd with the name XDATCAR-file
# and also a input file with the other important data (take a look at input text file)

import subprocess

dynaphopy_command = 'dynaphopy input XDATCAR-file -ts 0.0015 -sfc FORCE_CONSTANTS'

result = subprocess.run(dynaphopy_command, shell=True, capture_output=True, text=True)