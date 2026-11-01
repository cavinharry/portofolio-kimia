# Konformer 3D dengan RDKit

Latihan kimia komputasi: mencari bentuk 3D berenergi rendah (konformer) untuk tiga molekul, lalu membandingkan energi relatifnya. Ini latihan belajar, bukan penelitian.

## Alur kerja

1. SMILES diubah menjadi molekul, lalu atom H ditambahkan.
2. 500 konformer awal dibangkitkan dengan ETKDGv3 (benih acak 42, jadi hasilnya bisa diulang persis).
3. Setiap konformer dioptimasi dengan medan gaya MMFF94.
4. Duplikat dibuang: dua konformer dianggap sama bila RMSD atom beratnya di bawah 0,5 Å. Dua atom O pada gugus COOH tidak dianggap setara (lihat bagian perbaikan). Untuk molekul achiral, bayangan cermin juga dianggap sama.
5. Konformer diurutkan berdasarkan energi, lalu populasi Boltzmann pada 298 K dihitung.

## Hasil

| Molekul | Rumus | Konformer unik dari 500 | Selisih energi | Populasi konformer terendah |
|---|---|---|---|---|
| Kuersetin | C₁₅H₁₀O₇ | 2 | 0,15 kkal/mol | 56% |
| Asam klorogenat | C₁₆H₁₈O₉ | 121 | 2,80 kkal/mol di antara 12 terendah | 64% |
| Asam tereftalat | C₈H₆O₄ | 2 | 0,04 kkal/mol | 52% |

Asam tereftalat punya dua ikatan berputar (cincin–COOH). Kedua gugus karboksil tetap sebidang dengan cincin, jadi yang tersisa hanya dua bentuk: C=O kedua gugus berlawanan arah (trans) atau searah (cis). Selisih energinya 0,04 kkal/mol, jauh di bawah RT, jadi keduanya praktis sama banyak. Di 500 percobaan tidak ada yang berakhir dengan H karboksil menghadap menjauhi C=O (bentuk anti). Asam klorogenat jauh lebih fleksibel karena cincin sikloheksana dan ikatan esternya.

File SDF di `hasil/` bisa dibuka di Avogadro, PyMOL, atau py3Dmol.

## Batasan

- Energi MMFF94 adalah perkiraan fase gas dari medan gaya, bukan perhitungan kuantum dan bukan energi dalam larutan. Energi hanya dibandingkan antar konformer dari molekul yang sama.
- Torsi cincin B kuersetin (C3–C2–C1′–C2′) di model ini sekitar 45–48° dari sebidang. Di literatur, perhitungan DFT B3LYP fase gas memberi minimum sebidang, 0° ([Aparicio, 2010](https://pmc.ncbi.nlm.nih.gov/articles/PMC2885091)); studi DFT/MP2 lain menemukan 48 konformer stabil, separuh sebidang dan separuh terpilin 40,9–44,3° ([Brovarets' dan Hovorun, 2020](https://doi.org/10.3390/app10031147)); struktur padatan hampir sebidang, sekitar 5–7° (Aparicio, 2010; data Rossi dkk., 1986 dalam [Marković dkk., 2013](https://www.rsc.org/suppdata/cp/c3/c3cp44605k/c3cp44605k.pdf)). Jadi sudut MMFF94 dibaca sebagai hasil medan gaya yang perlu dicek dengan metode kuantum, bukan nilai pasti kuersetin.
- Medan gaya fase gas cenderung melebih-lebihkan ikatan hidrogen di dalam molekul, jadi urutan konformer asam klorogenat dalam air bisa berbeda.

## Perbaikan: asam tereftalat sempat terhitung 1 konformer

Versi sebelumnya melaporkan asam tereftalat hanya punya 1 konformer unik. Penyebabnya ada di langkah pembuangan duplikat: `rdMolAlign.GetBestRMS` secara bawaan (`symmetrizeConjugatedTerminalGroups=True`) menganggap dua atom O pada gugus karboksil setara. Itu benar untuk ion karboksilat (COO⁻), tetapi tidak untuk asam (C=O berbeda dengan C–OH). Karena RMSD dihitung tanpa atom H, rotamer cis dan trans punya posisi atom berat yang sama persis (hanya O=C dan O–H yang bertukar tempat), sehingga RMSD-nya 0,07 Å dan keduanya digabung. Setelah opsi itu dimatikan, asam tereftalat punya 2 konformer dan asam klorogenat (yang juga punya COOH) naik dari 77 menjadi 121. Kuersetin tidak berubah.

## Perubahan dari skrip latihan awal

Skrip awal saya membandingkan energi 50 konformer tanpa optimasi geometri, sehingga "pose paling stabil" belum bisa dipercaya. Versi ini menambahkan optimasi MMFF94, pembuangan duplikat, populasi Boltzmann, dan deskriptor molekul (log P, TPSA, donor/akseptor H).

## Menjalankan

```bash
pip install -r ../requirements.txt
python konformer_rdkit.py
```

Butuh sekitar 25 detik. Hasilnya `hasil/konformer.json` dan satu file SDF per molekul.

Skrip dirapikan dengan bantuan Claude (Anthropic) dari skrip latihan saya.
