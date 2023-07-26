# Pol Benítez Colominas, June 2023
# This code changes the atoms in POSCAR file


for x in range(20):

    file_path = 'new-phases-' + str(x+1) + '/POSCAR'

    with open(file_path,'r') as file:
        poscar = file.readlines()

    poscar[5] = 'atom1 atom2 atom3\n'

    with open(file_path, 'w') as file:
        file.writelines(poscar)
