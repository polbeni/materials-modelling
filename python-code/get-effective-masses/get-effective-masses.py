# Pol Benítez Colominas, June 2024
# Universitat Politècnica de Catalunya

# This code uses sumo to get effective masses and optoelectronic properties and save it in a text file
# IMPORTANT: you should have in the same directory the following files: KPOINTS and vasprun.xml

import subprocess

sumo_command = 'sumo-bandstats'

result = subprocess.run(sumo_command, shell=True, capture_output=True, text=True)