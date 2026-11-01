"""Membangun docs/index.html (website GitHub Pages) dari templat dan data repo.

Data yang disuntikkan:
    - hasil konformer: konformer-rdkit/hasil/konformer.json
    - data mentah skripsi: skripsi-red-dacca/data/*.csv dan parameter di analisis.py
    - kontak di bawah (yang kosong tidak ditampilkan di halaman)

Jalankan dari folder repo:  python web/bangun_web.py
"""
import importlib.util
import json
from pathlib import Path

import pandas as pd
import rdkit

AKAR = Path(__file__).resolve().parent.parent
REPO = "https://github.com/cavinharry/portofolio-kimia"

# Kontak yang tampil di halaman. Isi linkedin dengan URL lengkap profil, misalnya
# "https://www.linkedin.com/in/nama-anda". Nilai kosong tidak ditampilkan.
KONTAK = {
    "email": "nbilanf@gmail.com",
    "linkedin": "",
    "github": "https://github.com/cavinharry",
}


def muat_analisis():
    spec = importlib.util.spec_from_file_location("analisis", AKAR / "skripsi-red-dacca" / "analisis.py")
    modul = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(modul)
    return modul


def data_skripsi(analisis):
    folder = AKAR / "skripsi-red-dacca" / "data"
    hasil = {}
    for nama in ("fenolik", "flavonoid"):
        uji = analisis.UJI[nama]
        kal = pd.read_csv(folder / f"kalibrasi_{nama}.csv")
        sampel = pd.read_csv(folder / f"sampel_{nama}.csv")
        hasil[nama] = {
            "standar": kal[["konsentrasi_ppm", "absorbansi"]].values.tolist(),
            "sampel": sampel.to_dict(orient="records"),
            "V_L": uji["V_L"], "FP": uji["FP"], "m_g": uji["m_g"],
            "excel": list(uji["koef_excel"]),
        }
    dpph = pd.read_csv(folder / "dpph_sampel.csv")
    kontrol = pd.read_csv(folder / "dpph_kontrol.csv").set_index("larutan")
    hasil["dpph"] = {
        "konsentrasi": dpph["konsentrasi_persen_vv"].tolist(),
        "absorbansi": dpph["absorbansi"].tolist(),
        "blanko": float(kontrol.loc["blanko DPPH", "absorbansi"]),
        "vitamin_c": {"konsentrasi_ppm": float(kontrol.loc["vitamin C", "konsentrasi"]),
                      "absorbansi": float(kontrol.loc["vitamin C", "absorbansi"])},
    }
    hasil["nilai_skripsi"] = analisis.NILAI_SKRIPSI
    return hasil


def aman_untuk_script(teks):
    # Mencegah "</script>" di dalam data menutup tag lebih awal.
    return teks.replace("<", "\\u003c")


if __name__ == "__main__":
    analisis = muat_analisis()
    data = json.loads((AKAR / "konformer-rdkit" / "hasil" / "konformer.json").read_text(encoding="utf-8"))
    data["skripsi"] = data_skripsi(analisis)
    data["kontak"] = KONTAK
    html = (Path(__file__).parent / "index.template.html").read_text(encoding="utf-8")
    html = html.replace("__DATA__", aman_untuk_script(json.dumps(data, ensure_ascii=False)))
    html = html.replace("__RDKIT__", rdkit.__version__)
    html = html.replace("__REPO__", REPO)
    for sisa in ("__DATA__", "__RDKIT__", "__REPO__"):
        assert sisa not in html, f"placeholder {sisa} belum terganti"
    keluaran = AKAR / "docs" / "index.html"
    keluaran.parent.mkdir(exist_ok=True)
    keluaran.write_text(html, encoding="utf-8")
    print(f"{keluaran.relative_to(AKAR)}: {len(html) / 1024:.0f} KB")
