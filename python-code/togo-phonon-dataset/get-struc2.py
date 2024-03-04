# get materials with the desired bandgap

from mp_api.client import MPRester

file_names = open('valid_structures1.txt', 'r')
file_names.readline()
id_names = []
for x in range(42):
    line = file_names.readline()
    material_id = line.split()[0]
    if material_id != 'mp-864801-20180417': # this is a problematic phase
        id_names.append(material_id)

file_names.close()

structures_names = []
for path_mat in id_names:
    parts = path_mat.split('-')
    name_mp = '-'.join(parts[:2])
    structures_names.append(name_mp)

api_key = 'mp api key'

bandgaps = []
for struc in structures_names:
    with MPRester(api_key) as mpr:
        docs = mpr.materials.summary.search(material_ids=[struc], fields=["band_gap"])
        band_gap = docs[0].band_gap
        print(struc, band_gap)

    bandgaps.append(band_gap)

print(structures_names)
print(bandgaps)

valid_structures = open('valid_structures2.txt', 'w')
valid_structures.write(f'mp-id     bg (eV)\n')
for x in range(len(structures_names)):
    if (bandgaps[x] >= 1.5) and (bandgaps[x] <= 4):
        valid_structures.write(f'{structures_names[x]}     {bandgaps[x]}\n')
        
valid_structures.close()
