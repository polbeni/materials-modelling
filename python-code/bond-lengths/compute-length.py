from pymatgen.core.structure import Structure

def length_bond(struc, atom1, atom2):
    """
    Computes the lenght of the bond between two atoms

    Inputs:
        struc -> unit cell structure in pymatgen format
        atom1 -> first atom index
        atom2 -> second atom index
    """

    length = struc.get_distance(atom1, atom2)

    return length

lengths_file = open('lengths_file.txt', 'w')
lengths_file.write('Bond length     u=0     u=0.2     u=0.4     u=0.6\n')

print('Length bonds Ag-S for Ag3SBr are: ')

path = 'Ag3SBr/'

poscar0 = path + 'POSCAR-0'
structure = Structure.from_file(poscar0)
dist01 = length_bond(structure, 0, 3)
dist02 = length_bond(structure, 1, 3)
dist03 = length_bond(structure, 2, 3)
print(f'For the non distorted structures distances (in angstrom) are: {dist01:.2f}, {dist02:.2f}, {dist03:.2f}')

poscar2 = path + 'POSCAR-2'
structure = Structure.from_file(poscar2)
dist21 = length_bond(structure, 0, 3)
dist22 = length_bond(structure, 1, 3)
dist23 = length_bond(structure, 2, 3)
print(f'For distorted structures with u=0.2 distances (in angstrom) are: {dist21:.2f}, {dist22:.2f}, {dist23:.2f}')

poscar4 = path + 'POSCAR-4'
structure = Structure.from_file(poscar4)
dist41 = length_bond(structure, 0, 3)
dist42 = length_bond(structure, 1, 3)
dist43 = length_bond(structure, 2, 3)
print(f'For distorted structures with u=0.4 distances (in angstrom) are: {dist41:.2f}, {dist42:.2f}, {dist43:.2f}')

poscar6 = path + 'POSCAR-6'
structure = Structure.from_file(poscar6)
dist61 = length_bond(structure, 0, 3)
dist62 = length_bond(structure, 1, 3)
dist63 = length_bond(structure, 2, 3)
print(f'For distorted structures with u=0.6 distances (in angstrom) are: {dist61:.2f}, {dist62:.2f}, {dist63:.2f}')

lengths_file.write(f'Ag-S            {dist01:.2f}    {dist21:.2f}      {dist41:.2f}      {dist61:.2f}\n')
lengths_file.write(f'Ag-S            {dist02:.2f}    {dist22:.2f}      {dist42:.2f}      {dist62:.2f}\n')
lengths_file.write(f'Ag-S            {dist03:.2f}    {dist23:.2f}      {dist43:.2f}      {dist63:.2f}\n')

print('')

print('Length bonds Pb-O for HfPbO3 are: ')

path = 'HfPbO3/'

poscar0 = path + 'POSCAR-0'
structure = Structure.from_file(poscar0)
dist01 = length_bond(structure, 0, 4)
dist02 = length_bond(structure, 1, 4)
dist03 = length_bond(structure, 2, 4)
print(f'For the non distorted structures distances (in angstrom) are: {dist01:.2f}, {dist02:.2f}, {dist03:.2f}')

poscar2 = path + 'POSCAR-2'
structure = Structure.from_file(poscar2)
dist21 = length_bond(structure, 0, 4)
dist22 = length_bond(structure, 1, 4)
dist23 = length_bond(structure, 2, 4)
print(f'For distorted structures with u=0.2 distances (in angstrom) are: {dist21:.2f}, {dist22:.2f}, {dist23:.2f}')

poscar4 = path + 'POSCAR-4'
structure = Structure.from_file(poscar4)
dist41 = length_bond(structure, 0, 4)
dist42 = length_bond(structure, 1, 4)
dist43 = length_bond(structure, 2, 4)
print(f'For distorted structures with u=0.4 distances (in angstrom) are: {dist41:.2f}, {dist42:.2f}, {dist43:.2f}')

poscar6 = path + 'POSCAR-6'
structure = Structure.from_file(poscar6)
dist61 = length_bond(structure, 0, 4)
dist62 = length_bond(structure, 1, 4)
dist63 = length_bond(structure, 2, 4)
print(f'For distorted structures with u=0.6 distances (in angstrom) are: {dist61:.2f}, {dist62:.2f}, {dist63:.2f}')

lengths_file.write(f'Pb-O            {dist01:.2f}    {dist21:.2f}      {dist41:.2f}      {dist61:.2f}\n')
lengths_file.write(f'Pb-O            {dist02:.2f}    {dist22:.2f}      {dist42:.2f}      {dist62:.2f}\n')
lengths_file.write(f'Pb-O            {dist03:.2f}    {dist23:.2f}      {dist43:.2f}      {dist63:.2f}\n')

print('')

print('Length bonds Ti-O for BaTiO3 are: ')

path = 'BaTiO3/'

poscar0 = path + 'POSCAR-0'
structure = Structure.from_file(poscar0)
dist01 = length_bond(structure, 0, 3)
dist02 = length_bond(structure, 1, 3)
dist03 = length_bond(structure, 2, 3)
print(f'For the non distorted structures distances (in angstrom) are: {dist01:.2f}, {dist02:.2f}, {dist03:.2f}')

poscar2 = path + 'POSCAR-2'
structure = Structure.from_file(poscar2)
dist21 = length_bond(structure, 0, 3)
dist22 = length_bond(structure, 1, 3)
dist23 = length_bond(structure, 2, 3)
print(f'For distorted structures with u=0.2 distances (in angstrom) are: {dist21:.2f}, {dist22:.2f}, {dist23:.2f}')

poscar4 = path + 'POSCAR-4'
structure = Structure.from_file(poscar4)
dist41 = length_bond(structure, 0, 3)
dist42 = length_bond(structure, 1, 3)
dist43 = length_bond(structure, 2, 3)
print(f'For distorted structures with u=0.4 distances (in angstrom) are: {dist41:.2f}, {dist42:.2f}, {dist43:.2f}')

poscar6 = path + 'POSCAR-6'
structure = Structure.from_file(poscar6)
dist61 = length_bond(structure, 0, 3)
dist62 = length_bond(structure, 1, 3)
dist63 = length_bond(structure, 2, 3)
print(f'For distorted structures with u=0.6 distances (in angstrom) are: {dist61:.2f}, {dist62:.2f}, {dist63:.2f}')

lengths_file.write(f'Ti-O            {dist01:.2f}    {dist21:.2f}      {dist41:.2f}      {dist61:.2f}\n')
lengths_file.write(f'Ti-O            {dist02:.2f}    {dist22:.2f}      {dist42:.2f}      {dist62:.2f}\n')
lengths_file.write(f'Ti-O            {dist03:.2f}    {dist23:.2f}      {dist43:.2f}      {dist63:.2f}\n')

lengths_file.close()