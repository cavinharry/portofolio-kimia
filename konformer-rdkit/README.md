# Konformer 3D dengan RDKit

Latihan kimia komputasi: mencari bentuk 3D berenergi rendah (konformer) untuk tiga molekul, lalu membandingkan energi relatifnya. Ini latihan belajar, bukan penelitian.

## Alur kerja

1. SMILES diubah menjadi molekul, lalu atom H ditambahkan.
2. 500 konformer awal dibangkitkan dengan ETKDGv3 (benih acak 42, jadi hasilnya bisa diulang persis).
3. Setiap konformer dioptimasi dengan medan gaya MMFF94.
4. Duplikat dibuang: dua konformer dianggap sama bila RMSD atom beratnya di bawah 0,5 Å. Untuk molekul achiral, bayangan cermin juga dianggap sama.
5. Konformer diurutkan berdasarkan energi, lalu populasi Boltzmann pada 298 K dihitung.

## Hasil

| Molekul | Rumus | Konformer unik dari 500 | Selisih energi | Populasi konformer terendah |
|---|---|---|---|---|
| Kuersetin | C₁₅H₁₀O₇ | 2 | 0,15 kkal/mol | 56% |
| Asam klorogenat | C₁₆H₁₈O₉ | 77 | 2,90 kkal/mol di antara 12 terendah | 67% |
| Asam tereftalat | C₈H₆O₄ | 1 | – | 100% |

Asam tereftalat praktis kaku karena cincin benzenanya. Asam klorogenat jauh lebih fleksibel karena cincin sikloheksana dan ikatan esternya.

File SDF di `hasil/` bisa dibuka di Avogadro, PyMOL, atau py3Dmol.

## Batasan

- Energi MMFF94 adalah perkiraan fase gas dari medan gaya, bukan perhitungan kuantum dan bukan energi dalam larutan. Energi hanya dibandingkan antar konformer dari molekul yang sama.
- Torsi cincin B kuersetin di model ini sekitar 45° dari sebidang. Perhitungan DFT B3LYP fase gas memberi sekitar 0° dan struktur kristal sekitar 5° ([sumber](https://pmc.ncbi.nlm.nih.gov/articles/PMC2885091)). Jadi angka torsi ini adalah batasan medan gaya, bukan nilai kuersetin yang sebenarnya.
- Medan gaya fase gas cenderung melebih-lebihkan ikatan hidrogen di dalam molekul, jadi urutan konformer asam klorogenat dalam air bisa berbeda.

## Perubahan dari skrip latihan awal

Skrip awal saya membandingkan energi 50 konformer tanpa optimasi geometri, sehingga "pose paling stabil" belum bisa dipercaya. Versi ini menambahkan optimasi MMFF94, pembuangan duplikat, populasi Boltzmann, dan deskriptor molekul (log P, TPSA, donor/akseptor H).

## Menjalankan

```bash
pip install -r ../requirements.txt
python konformer_rdkit.py
```

Butuh sekitar 25 detik. Hasilnya `hasil/konformer.json` dan satu file SDF per molekul.

Skrip dirapikan dengan bantuan Claude (Anthropic) dari skrip latihan saya.
