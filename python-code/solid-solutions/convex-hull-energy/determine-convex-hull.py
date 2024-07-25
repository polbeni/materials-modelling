import sympy as sp

from pymatgen.io.vasp import Poscar
from pymatgen.io.vasp import Vasprun

compounds = {
    # primary compounds of CAP materials, energy per atom (eV/1atom) times the number of atoms
    'Ag': -3.312115,
    'S': -4.37164,
    'Se': -3.789887187,
    'Br': -1.75290445,
    'I': -1.6907348,
    'Ag2S': -3.699058833*3,
    'Ag2Se': -3.542973417*3,
    'AgBr': -2.849598375*2,
    'AgI': -2.733631625*2,
    'Ag3SBr': -3.2734616*5,
    'Ag3SI': -3.2037046*5,
    'Ag3SeBr': -3.076793*5,
    'Ag3SeI': -3.0271424*5
}

# A:Ag, B:S, C:Se, D:Br, E:I, F:Ag2S, G:Ag2Se, H:'AgBr', I:'AgI', J:'Ag3SBr', K:'Ag3SI', L:'Ag3SeBr', M:'Ag3SeI'
chem_A = {'Ag': 1}
chem_B = {'S': 1}
chem_C = {'Se': 1}
chem_D = {'Br': 1}
chem_E = {'I': 1}
chem_F = {'Ag': 2, 'S': 1}
chem_G = {'Ag': 2, 'Se': 1}
chem_H = {'Ag': 1, 'Br': 1}
chem_I = {'Ag': 1, 'I': 1}
chem_J = {'Ag': 3, 'S': 1, 'Br': 1}
chem_K = {'Ag': 3, 'S': 1, 'I': 1}
chem_L = {'Ag': 3, 'Se': 1, 'Br': 1}
chem_M = {'Ag': 3, 'Se': 1, 'I': 1}

def compute_convex_hull_energy(path_POSCAR, path_INCAR, path_vasprun):
    """
    Computes the convex hull energy for a given solid solution compound

    Inputs:
        path_POSCAR: path to the POSCAR file
        path_INCAR: path to the INCAR file
        path_vasprun: path to the vasprun.xml
    """
    poscar = Poscar.from_file(path_POSCAR)
    structure = poscar.structure
    list_atoms = [site.species_string for site in structure.sites]
    unique_atom_list = list(dict.fromkeys(list_atoms))
    total_num_atoms = 5 # we have the same number of atoms

    coefs = []
    with open(path_INCAR, 'r') as file:
        for line in file:
            values = line.split()
    
    for x in range(len(values)):
        if (x != 0) and (x != 1):
            coefs.append(float(values[x]))

    vasprun = Vasprun(path_vasprun)
    energy = float(vasprun.final_energy)
    print(energy)
    energy_per_atom = energy/total_num_atoms

    atoms_coef = []
    for x in range(len(coefs)):
        if x == 0:
            atoms_coef.append(3.0*coefs[0])
        else:
            atoms_coef.append(coefs[x])

    chem_solid_solution = {}
    for x in range(len(atoms_coef)):
        chem_solid_solution[unique_atom_list[x]] = atoms_coef[x]

    list_primary_compounds = [chem_A]
    list_primary_compounds_energy = [compounds['Ag']]
    if 'S' in unique_atom_list:
        list_primary_compounds.append(chem_B)
        list_primary_compounds_energy.append(compounds['S'])
        list_primary_compounds.append(chem_F)
        list_primary_compounds_energy.append(compounds['Ag2S'])
    if 'Se' in unique_atom_list:
        list_primary_compounds.append(chem_C)
        list_primary_compounds_energy.append(compounds['Se'])
        list_primary_compounds.append(chem_G)
        list_primary_compounds_energy.append(compounds['Ag2Se'])
    if 'Br' in unique_atom_list:
        list_primary_compounds.append(chem_D)
        list_primary_compounds_energy.append(compounds['Br'])
        list_primary_compounds.append(chem_H)
        list_primary_compounds_energy.append(compounds['AgBr'])
    if 'I' in unique_atom_list:
        list_primary_compounds.append(chem_E)
        list_primary_compounds_energy.append(compounds['I'])
        list_primary_compounds.append(chem_I)
        list_primary_compounds_energy.append(compounds['AgI'])
    if ('S' in unique_atom_list) and ('Br' in unique_atom_list) and (len(unique_atom_list) >= 4):
        list_primary_compounds.append(chem_J)
        list_primary_compounds_energy.append(compounds['Ag3SBr'])
    if ('S' in unique_atom_list) and ('I' in unique_atom_list) and (len(unique_atom_list) >= 4):
        list_primary_compounds.append(chem_K)
        list_primary_compounds_energy.append(compounds['Ag3SI'])
    if ('Se' in unique_atom_list) and ('Br' in unique_atom_list) and (len(unique_atom_list) >= 4):
        list_primary_compounds.append(chem_L)
        list_primary_compounds_energy.append(compounds['Ag3SeBr'])
    if ('Se' in unique_atom_list) and ('I' in unique_atom_list) and (len(unique_atom_list) >= 4):
        list_primary_compounds.append(chem_M)
        list_primary_compounds_energy.append(compounds['Ag3SeI'])

    reactants = [chem_solid_solution]
    products = list_primary_compounds
    #print(reactants, products)

    coefficients = sp.symbols(f'a:{len(reactants) + len(products)}')

    elements = set()
    for compound in reactants + products:
        elements.update(compound.keys())

    equations = []
    for element in elements:
        equation = sp.Eq(
            sum(coefficients[i] * reactants[i].get(element, 0) for i in range(len(reactants))),
            sum(coefficients[i + len(reactants)] * products[i].get(element, 0) for i in range(len(products)))
        )
        equations.append(equation)
    
    solution = sp.solve(equations, coefficients, dict=True)
    general_solution = list(solution[0])
    
    independet_var = set(general_solution)
    all_var = set(coefficients)
    dependent_var = all_var - independet_var

    fixed_value = 1
    fixed_values = {var: fixed_value for var in dependent_var}

    results = {key: expr.subs(fixed_values) for key, expr in solution[0].items()}

    combined_dict = {**fixed_values, **results}

    def extract_number(symbol):
    # Convert symbolic variables to string and extract the number part
        return int(str(symbol)[1:])

    sorted_keys = sorted(combined_dict.keys(), key=extract_number)

    ordered_values = [combined_dict[key] for key in sorted_keys]

    convex_energy = energy_per_atom*5*ordered_values[0]
    #print(convex_energy)

    for x in range(len(list_primary_compounds_energy)):
        convex_energy = convex_energy - list_primary_compounds_energy[x]*ordered_values[x + 1]
        #print(convex_energy)

    #print(list_primary_compounds, list_primary_compounds_energy, ordered_values)

    return convex_energy

path = '../generate-VCA-grid/VCA_structures/vca-116/relaxation/'
convex_energy = compute_convex_hull_energy(path + 'POSCAR', path + 'INCAR', path + 'vasprun.xml')
print(convex_energy)

path = '../generate-VCA-grid/VCA_structures/vca-'
for x in range(121):
    path_ss = path + str(x + 1).zfill(3) + '/relaxation/'
    convex_energy = compute_convex_hull_energy(path_ss + 'POSCAR', path_ss + 'INCAR', path_ss + 'vasprun.xml')
    print(convex_energy)
