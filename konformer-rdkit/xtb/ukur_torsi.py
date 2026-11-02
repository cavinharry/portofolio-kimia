import math
from rdkit import Chem
from rdkit.Chem import rdMolTransforms

HARTREE_KE_KKAL = 627.5095          # 1 Eh = 627,5095 kkal/mol
RT = 0.0019872 * 298.15             # kkal/mol pada 298 K
TORSI = (2, 4, 5, 6)                # atom C3, C2, C1', C6' (nomor urut atom di file)

# Energi total dari grep "TOTAL ENERGY" (satuan Eh)
energi = {"k1": -65.720508540427, "k2": -65.718963367141}

for nama in ("k1", "k2"):
    berkas = {"MMFF94": f"konformer{nama[1]}.xyz", "xTB": f"{nama}/xtbopt.xyz"}
    for metode, file in berkas.items():
        mol = Chem.MolFromXYZFile(file)
        t = rdMolTransforms.GetDihedralDeg(mol.GetConformer(), *TORSI)
        dari_bidang = min(abs(t), 180 - abs(t))
        print(f"{nama} {metode}: torsi {t:.1f}°, {dari_bidang:.1f}° dari sebidang")

dE = (energi["k2"] - energi["k1"]) * HARTREE_KE_KKAL
w = math.exp(-dE / RT)
print(f"Selisih energi xTB (k2 - k1): {dE:.2f} kkal/mol")
print(f"Populasi 298 K: k1 {100 / (1 + w):.1f}%, k2 {100 * w / (1 + w):.1f}%")