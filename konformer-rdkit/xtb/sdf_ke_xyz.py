from rdkit import Chem

# Baca semua konformer kuersetin dari file SDF hasil skrip RDKit.
# removeHs=False: atom H harus ikut, karena xTB menghitung semua atom.
pembaca = Chem.SDMolSupplier("../hasil/kuersetin.sdf", removeHs=False)

for nomor, mol in enumerate(pembaca, start=1):
    # Simpan tiap konformer sebagai file XYZ terpisah.
    Chem.MolToXYZFile(mol, f"konformer{nomor}.xyz")
    print(f"konformer{nomor}.xyz: {mol.GetNumAtoms()} atom")