"""
konformer_rdkit.py
Pencarian konformer 3D untuk beberapa molekul memakai RDKit.

Alur kerja:
    SMILES -> tambah atom H -> bangkitkan banyak konformer (ETKDG)
    -> optimasi geometri tiap konformer (MMFF94) -> buang duplikat (RMSD)
    -> urutkan berdasarkan energi -> simpan SDF dan ringkasan JSON.

Catatan penting:
    Energi MMFF94 adalah perkiraan fase gas dari medan gaya (force field),
    bukan energi kuantum dan bukan energi dalam larutan. Pakai untuk
    membandingkan konformer dari SATU molekul yang sama, jangan untuk
    membandingkan molekul yang berbeda.
"""
import json
import math
from pathlib import Path

from rdkit import Chem
from rdkit.Chem import (AllChem, Crippen, Descriptors, Lipinski, rdMolAlign,
                        rdMolDescriptors, rdMolTransforms)

HASIL = Path(__file__).parent / "hasil"  # folder hasil di sebelah skrip ini
SUHU_K = 298.15
RT_KKAL = 0.0019872 * SUHU_K  # RT dalam kkal/mol (R = 1,9872 kal/mol/K)

# id: (nama tampilan, SMILES)
MOLEKUL = {
    "kuersetin": ("Kuersetin (quercetin)", "O=c1c(O)c(-c2ccc(O)c(O)c2)oc2cc(O)cc(O)c12"),
    # (1S,3R,4R,5R)-asam klorogenat, C16H18O9. SMILES ini diturunkan dari konstitusinya lalu
    # dipilih stereoisomer yang label CIP-nya cocok (bukan disalin dari ingatan).
    "asam_klorogenat": (
        "Asam klorogenat",
        "O=C(/C=C/c1ccc(O)c(O)c1)O[C@@H]1C[C@](O)(C(=O)O)C[C@@H](O)[C@H]1O",
    ),
    "asam_tereftalat": ("Asam tereftalat", "O=C(O)c1ccc(C(=O)O)cc1"),
}

# Torsi yang ingin dilaporkan per molekul: empat indeks atom (urutan atom pada SMILES).
# Kuersetin: C3-C2-C1'-C2' = seberapa miring cincin B terhadap inti kromon.
TORSI = {"kuersetin": (2, 4, 5, 6)}


def hitung_deskriptor(smiles):
    """Deskriptor sederhana yang lazim dipakai di kemoinformatika."""
    m = Chem.MolFromSmiles(smiles)
    return {
        "rumus": rdMolDescriptors.CalcMolFormula(m),
        "mr": round(Descriptors.MolWt(m), 2),
        "logp": round(Crippen.MolLogP(m), 2),
        "tpsa": round(rdMolDescriptors.CalcTPSA(m), 1),
        "donor_h": Lipinski.NumHDonors(m),
        "akseptor_h": Lipinski.NumHAcceptors(m),
        "ikatan_putar": Lipinski.NumRotatableBonds(m),
        "atom_berat": m.GetNumHeavyAtoms(),
    }


def bayangan_cermin(mol_satu_konformer):
    """Salinan molekul dengan koordinat x dibalik (cermin lewat bidang yz)."""
    cermin = Chem.Mol(mol_satu_konformer)
    konf = cermin.GetConformer()
    for i in range(cermin.GetNumAtoms()):
        p = konf.GetAtomPosition(i)
        konf.SetAtomPosition(i, (-p.x, p.y, p.z))
    return cermin


def cari_konformer(smiles, jumlah_percobaan=500, ambang_rmsd=0.5, seed=42):
    """Kembalikan (molekul_H, [(energi, id_konformer)] unik dan terurut, jumlah konvergen)."""
    molekul = Chem.AddHs(Chem.MolFromSmiles(smiles))
    if not AllChem.MMFFHasAllMoleculeParams(molekul):
        raise ValueError("Parameter MMFF94 tidak lengkap untuk molekul ini")

    # 1. Bangkitkan konformer awal. randomSeed membuat hasil bisa diulang persis.
    params = AllChem.ETKDGv3()
    params.randomSeed = seed
    AllChem.EmbedMultipleConfs(molekul, numConfs=jumlah_percobaan, params=params)

    # 2. Optimasi geometri. Tanpa langkah ini, energi hanya mencerminkan tebakan
    #    awal ETKDG, jadi urutan "paling stabil" belum bisa dipercaya.
    hasil = AllChem.MMFFOptimizeMoleculeConfs(molekul, maxIters=5000, numThreads=0)
    kandidat = []
    for konf, (kode, energi) in zip(molekul.GetConformers(), hasil):
        if kode == 0:  # 0 = optimasi konvergen
            kandidat.append((energi, konf.GetId()))
    kandidat.sort()

    # 3. Buang konformer duplikat. Banyak tebakan awal berakhir di lembah energi
    #    yang sama. Dibandingkan lewat RMSD atom berat (memperhitungkan simetri),
    #    jadi konformer yang hanya beda arah H pada gugus OH dianggap sama.
    #    Untuk molekul achiral, bayangan cermin sebuah konformer secara fisik
    #    adalah konformer yang sama, jadi ikut dianggap duplikat. Itu sebabnya
    #    energinya persis sama. (Untuk molekul kiral, cermin adalah molekul lain.)
    achiral = not Chem.FindMolChiralCenters(Chem.MolFromSmiles(smiles), includeUnassigned=True)
    tanpa_h = Chem.RemoveHs(molekul)
    unik = []
    for energi, cid in kandidat:
        calon = Chem.Mol(tanpa_h, False, cid)
        pembanding = [calon] + ([bayangan_cermin(calon)] if achiral else [])
        sama = False
        for _, cid_ref in unik:
            acuan = Chem.Mol(tanpa_h, False, cid_ref)
            if any(rdMolAlign.GetBestRMS(Chem.Mol(c), acuan) < ambang_rmsd for c in pembanding):
                sama = True
                break
        if not sama:
            unik.append((energi, cid))
    return molekul, unik, len(kandidat)


def susun_hasil(id_mol, nama, smiles, simpan_maks=12):
    molekul, unik, jumlah_konvergen = cari_konformer(smiles)

    # Konformer unik dikumpulkan ke satu molekul baru, dan disejajarkan ke yang
    # berenergi terendah supaya mudah dibandingkan secara visual.
    baru = Chem.Mol(molekul)
    baru.RemoveAllConformers()
    for urutan, (_, cid) in enumerate(unik):
        konf = Chem.Conformer(molekul.GetConformer(cid))
        konf.SetId(urutan)
        baru.AddConformer(konf, assignId=False)
    atom_berat = [a.GetIdx() for a in baru.GetAtoms() if a.GetAtomicNum() > 1]
    if len(unik) > 1:
        rdMolAlign.AlignMolConformers(baru, atomIds=atom_berat)

    energi = [e for e, _ in unik]
    delta = [e - energi[0] for e in energi]
    bobot = [math.exp(-d / RT_KKAL) for d in delta]
    populasi = [b / sum(bobot) * 100 for b in bobot]  # perkiraan Boltzmann, fase gas

    # Simpan SDF (semua konformer unik, lengkap dengan energinya)
    penulis = Chem.SDWriter(str(HASIL / f"{id_mol}.sdf"))
    for urutan in range(len(unik)):
        baru.SetProp("_Name", f"{nama} konformer {urutan}")
        baru.SetDoubleProp("MMFF94_energi_kkal_mol", round(energi[urutan], 3))
        baru.SetDoubleProp("delta_E_kkal_mol", round(delta[urutan], 3))
        penulis.write(baru, confId=urutan)
    penulis.close()

    konformer = []
    for i in range(min(simpan_maks, len(unik))):
        data = {
            "urutan": i,
            "energi": round(energi[i], 3),
            "delta": round(delta[i], 3),
            "populasi": round(populasi[i], 1),
            "molblock": Chem.MolToMolBlock(baru, confId=i),
        }
        if id_mol in TORSI:
            a, b, c, d = TORSI[id_mol]
            data["torsi"] = round(rdMolTransforms.GetDihedralDeg(baru.GetConformer(i), a, b, c, d), 1)
        konformer.append(data)

    return {
        "id": id_mol,
        "nama": nama,
        "smiles": smiles,
        "deskriptor": hitung_deskriptor(smiles),
        "jumlah_konvergen": jumlah_konvergen,
        "jumlah_unik": len(unik),
        "konformer": konformer,
    }


if __name__ == "__main__":
    HASIL.mkdir(exist_ok=True)
    semua = []
    for id_mol, (nama, smiles) in MOLEKUL.items():
        ringkasan = susun_hasil(id_mol, nama, smiles)
        semua.append(ringkasan)
        k = ringkasan["konformer"]
        print(
            f"{nama}: {ringkasan['jumlah_konvergen']} konvergen -> "
            f"{ringkasan['jumlah_unik']} unik | "
            f"E terendah {k[0]['energi']:.2f} kkal/mol | "
            f"rentang dE {k[-1]['delta']:.2f} kkal/mol"
        )
    with open(HASIL / "konformer.json", "w", encoding="utf-8") as f:
        json.dump({"molekul": semua}, f, ensure_ascii=False)
    print(f"Selesai: {HASIL / 'konformer.json'} dan file SDF per molekul.")
