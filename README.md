# Portofolio Nabilah Nur Fitri

Sarjana Pendidikan Kimia (S.Pd.), Universitas Halu Oleo. Berpengalaman dalam analisis fitokimia di laboratorium (spektrofotometri UV-Vis) dan sedang belajar kimia komputasi dengan Python.

*Chemistry Education graduate (Universitas Halu Oleo) with hands-on UV-Vis phytochemical analysis, currently learning computational chemistry in Python.*

**Kontak:** [nbilanf@gmail.com](mailto:nbilanf@gmail.com) · [github.com/cavinharry](https://github.com/cavinharry)

**Website:** [cavinharry.github.io/portofolio-kimia](https://cavinharry.github.io/portofolio-kimia/), versi interaktif dari isi repo ini: molekul kuersetin 3D dari koordinat RDKit, grafik kalibrasi dan DPPH beserta deret kuvet yang warnanya dihitung dari absorbansi, animasi 3D penampang cat Purivon, dan penampil konformer.

| Proyek | Isi | Status |
|---|---|---|
| [Air sadapan batang pisang merah](skripsi-red-dacca/) | Analisis ulang data skripsi dari absorbansi mentah: kurva kalibrasi, kadar total fenolik dan flavonoid, presisi, LOD/LOQ, aktivitas antioksidan DPPH, dan IC50 | Penelitian tugas akhir, 2026 |
| [Konformer 3D dengan RDKit](konformer-rdkit/) | Pencarian konformer (ETKDGv3 dan MMFF94) untuk kuersetin, asam klorogenat, dan asam tereftalat, beserta batasan medan gaya | Latihan belajar |
| [Purivon Technologies](https://cavinharry.github.io/portofolio-kimia/#purivon) | Studi kasus proposal startup cat fotokatalitik, finalis 10 besar PENDIKAR Startup Innovation 2025 | Konsep, belum ada data uji |

## Menjalankan

```bash
pip install -r requirements.txt
python skripsi-red-dacca/analisis.py
python konformer-rdkit/konformer_rdkit.py
python web/bangun_web.py    # membangun ulang docs/index.html dari data di atas
```

Website ada di `docs/index.html` dan dibangun dari `web/index.template.html`. Angka di website dihitung ulang di peramban dari data mentah yang sama. Model 3D dan animasi dibuat dengan Three.js langsung dari koordinat di `konformer-rdkit/hasil/konformer.json`; kalau 3D tidak bisa dimuat, grafik, tabel, dan teks tetap berfungsi.

Data, prosedur laboratorium, dan interpretasi berasal dari saya. Kode dan README dirapikan dengan bantuan Claude (Anthropic).
