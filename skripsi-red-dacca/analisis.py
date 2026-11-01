"""
analisis.py
Analisis ulang data skripsi: air sadapan batang semu pisang merah
(Musa acuminata Colla 'Red Dacca').

Yang dihitung dari data mentah (folder data/):
    1. Kurva kalibrasi asam galat dan kuersetin: regresi linear, R², residual,
       serta perkiraan LOD dan LOQ.
    2. Kadar total fenolik (KTF) dan kadar total flavonoid (KTFla):
       absorbansi -> konsentrasi dari kurva -> kadar per gram sampel,
       lalu rerata, SD, dan RSD dari tiga ulangan.
    3. Aktivitas antioksidan DPPH: % inhibisi dan IC50.

Setiap kadar dihitung dua kali:
    - "koefisien Excel": slope dan intersep yang dibulatkan seperti di Excel
      skripsi. Ini untuk membuktikan angka skripsi bisa direproduksi.
    - "koefisien penuh": slope dan intersep hasil regresi tanpa pembulatan.

Jalankan:  python analisis.py
Hasil:     grafik dan ringkasan di folder hasil/
"""
from pathlib import Path

import matplotlib
matplotlib.use("Agg")  # simpan ke file, tanpa membuka jendela
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.ticker import FuncFormatter

FOLDER = Path(__file__).parent
DATA = FOLDER / "data"
HASIL = FOLDER / "hasil"

# Warna grafik (sudah dicek aman untuk buta warna)
BIRU = "#1C4DB8"    # standar dan garis regresi
OKER = "#A85A0B"    # sampel
TINTA = "#1F2A28"   # teks
ABU = "#5B6B67"     # teks sekunder dan garis bantu

# Parameter rumus kadar = C x V x FP / m, diambil dari rumus di Excel skripsi:
#   fenolik   : =(0.001*C*20)/1      -> V = 0,001 L, FP = 20,   m = 1 g
#   flavonoid : =(0.003*C*3.33)/3    -> V = 0,003 L, FP = 3,33, m = 3 g
UJI = {
    "fenolik": {
        "judul": "Kadar total fenolik (Folin–Ciocalteu)",
        "standar": "asam galat",
        "satuan": "mg GAE/g",
        "V_L": 0.001, "FP": 20, "m_g": 1,
        "koef_excel": (0.0153, -0.1994),
    },
    "flavonoid": {
        "judul": "Kadar total flavonoid (AlCl₃)",
        "standar": "kuersetin",
        "satuan": "mg QE/g",
        "V_L": 0.003, "FP": 3.33, "m_g": 3,
        "koef_excel": (0.001, -0.0228),
    },
}

# Angka yang tertulis di Bab IV skripsi, untuk dibandingkan.
NILAI_SKRIPSI = {
    "fenolik": {"R2": 0.9914, "kadar": [1.187, 1.195, 1.190], "rerata": 1.191, "SD": 0.004, "RSD": 0.34},
    "flavonoid": {"R2": 0.9902, "kadar": [0.426, 0.432, 0.442], "rerata": 0.433, "SD": 0.008, "RSD": 1.87},
    "dpph": {"inhibisi": [42.82, 55.37, 62.34, 72.11, 74.90], "vitamin_c": 52.58},
}


def koma(x, desimal, tanda=False):
    """Format angka gaya Indonesia: 1.234 -> '1,234'. tanda=True menulis + atau - di depan."""
    teks = f"{x:+.{desimal}f}" if tanda else f"{x:.{desimal}f}"
    return teks.replace(".", ",")


def regresi_linear(x, y):
    """Regresi y = slope*x + intersep. Mengembalikan slope, intersep, R², residual."""
    slope, intersep = np.polyfit(x, y, 1)
    y_prediksi = slope * x + intersep
    residual = y - y_prediksi
    r2 = 1 - (residual ** 2).sum() / ((y - y.mean()) ** 2).sum()
    return slope, intersep, r2, residual


def lod_loq(slope, residual):
    """Perkiraan LOD dan LOQ dari kurva kalibrasi (pendekatan ICH Q2):
    LOD = 3,3 s/slope dan LOQ = 10 s/slope, dengan s = simpangan baku residual."""
    s = np.sqrt((residual ** 2).sum() / (len(residual) - 2))
    return 3.3 * s / slope, 10 * s / slope


def kadar_per_gram(absorbansi, slope, intersep, uji):
    """Absorbansi -> konsentrasi larutan ukur (mg/L) -> kadar (mg/g sampel)."""
    konsentrasi = (absorbansi - intersep) / slope
    kadar = konsentrasi * uji["V_L"] * uji["FP"] / uji["m_g"]
    return konsentrasi, kadar


def statistik(nilai):
    """Rerata, SD sampel (n-1, sama dengan STDEV di Excel), dan RSD (%)."""
    rerata = nilai.mean()
    sd = nilai.std(ddof=1)
    return rerata, sd, sd / rerata * 100


def gaya_sumbu(ax):
    """Grid dan sumbu dibuat samar supaya data yang menonjol."""
    ax.grid(True, color="#E3E8E6", linewidth=1)
    ax.set_axisbelow(True)
    for sisi in ("top", "right"):
        ax.spines[sisi].set_visible(False)
    for sisi in ("left", "bottom"):
        ax.spines[sisi].set_color(ABU)
    ax.tick_params(colors=ABU, labelcolor=TINTA)
    angka_koma = FuncFormatter(lambda v, _: f"{0 if abs(v) < 1e-12 else round(v, 6):g}".replace(".", ","))
    ax.xaxis.set_major_formatter(angka_koma)
    ax.yaxis.set_major_formatter(angka_koma)


def analisis_kadar(nama):
    uji = UJI[nama]
    kal = pd.read_csv(DATA / f"kalibrasi_{nama}.csv")
    sampel = pd.read_csv(DATA / f"sampel_{nama}.csv")
    x = kal["konsentrasi_ppm"].to_numpy(float)
    y = kal["absorbansi"].to_numpy(float)
    a_sampel = sampel["absorbansi"].to_numpy(float)

    slope, intersep, r2, residual = regresi_linear(x, y)
    lod, loq = lod_loq(slope, residual)

    hasil = {"nama": nama, "uji": uji, "slope": slope, "intersep": intersep, "R2": r2,
             "residual": residual, "LOD": lod, "LOQ": loq, "rentang": (x.min(), x.max())}
    for label, (s, b) in (("excel", uji["koef_excel"]), ("penuh", (slope, intersep))):
        konsentrasi, kadar = kadar_per_gram(a_sampel, s, b, uji)
        rerata, sd, rsd = statistik(kadar)
        hasil[label] = {"konsentrasi": konsentrasi, "kadar": kadar, "rerata": rerata, "SD": sd, "RSD": rsd}

    # Grafik: kurva kalibrasi di atas, residual di bawah
    fig, (ax, ax_res) = plt.subplots(2, 1, figsize=(7, 6.2), sharex=True,
                                     gridspec_kw={"height_ratios": [3, 1.2], "hspace": 0.08})
    x_garis = np.linspace(0, x.max() * 1.05, 200)
    ax.plot(x_garis, slope * x_garis + intersep, color=BIRU, linewidth=2, label="regresi linear")
    ax.plot(x, y, "o", color=BIRU, markersize=8, markeredgecolor="white", markeredgewidth=1.5,
            label=f"standar {uji['standar']}")
    ax.plot(hasil["penuh"]["konsentrasi"], a_sampel, "D", color=OKER, markersize=8,
            markeredgecolor="white", markeredgewidth=1.5, label="sampel (3 ulangan)")
    # Daerah abu-abu = di bawah LOQ, tempat konsentrasi belum bisa dikuantifikasi dengan andal
    ax.axvspan(0, loq, color=ABU, alpha=0.08, linewidth=0)
    ax.axvline(loq, color=ABU, linewidth=1, linestyle="--")
    ax.text(loq, 0.03, f"  LOQ ≈ {koma(loq, 0)} ppm", transform=ax.get_xaxis_transform(),
            ha="left", va="bottom", fontsize=9, color=ABU)
    ax.text(0.98, 0.28, f"A = {koma(slope, 6)}·C {koma(intersep, 4, tanda=True)}\nR² = {koma(r2, 4)}",
            transform=ax.transAxes, ha="right", va="bottom", fontsize=10, color=TINTA, family="monospace")
    ax.set_ylabel("Absorbansi", color=TINTA)
    ax.set_title(f"Kurva kalibrasi {uji['standar']}", loc="left", fontsize=12, color=TINTA, fontweight="bold")
    ax.legend(frameon=False, loc="upper left", labelcolor=TINTA)
    gaya_sumbu(ax)

    ax_res.axhline(0, color=ABU, linewidth=1)
    ax_res.plot(x, residual, "o", color=BIRU, markersize=7, markeredgecolor="white", markeredgewidth=1.2)
    ax_res.set_ylabel("Residual", color=TINTA)
    ax_res.set_xlabel(f"Konsentrasi {uji['standar']} (ppm)", color=TINTA)
    batas = np.abs(residual).max() * 1.6
    ax_res.set_ylim(-batas, batas)
    gaya_sumbu(ax_res)
    ax_res.set_xlim(0, x.max() * 1.05)

    fig.savefig(HASIL / f"kalibrasi_{nama}.png", dpi=150, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    return hasil


def analisis_dpph():
    sampel = pd.read_csv(DATA / "dpph_sampel.csv")
    kontrol = pd.read_csv(DATA / "dpph_kontrol.csv").set_index("larutan")
    a_blanko = kontrol.loc["blanko DPPH", "absorbansi"]
    a_vitc = kontrol.loc["vitamin C", "absorbansi"]

    konsentrasi = sampel["konsentrasi_persen_vv"].to_numpy(float)
    inhibisi = (a_blanko - sampel["absorbansi"].to_numpy(float)) / a_blanko * 100
    inhibisi_vitc = (a_blanko - a_vitc) / a_blanko * 100

    # IC50 cara umum: regresi linear % inhibisi terhadap konsentrasi, cari x saat y = 50
    slope, intersep, r2, _ = regresi_linear(konsentrasi, inhibisi)
    ic50 = (50 - intersep) / slope
    # Pembanding: regresi terhadap log10 konsentrasi (kurva dosis-respons biasanya melengkung)
    slope_log, intersep_log, r2_log, _ = regresi_linear(np.log10(konsentrasi), inhibisi)
    ic50_log = 10 ** ((50 - intersep_log) / slope_log)

    fig, ax = plt.subplots(figsize=(7, 4.6))
    x_garis = np.linspace(0, konsentrasi.max() * 1.08, 200)
    ax.plot(x_garis, slope * x_garis + intersep, color=BIRU, linewidth=2, label="regresi linear")
    ax.plot(konsentrasi, inhibisi, "o", color=OKER, markersize=8, markeredgecolor="white",
            markeredgewidth=1.5, label="air sadapan")
    ax.axhline(50, color=ABU, linewidth=1, linestyle="--")
    ax.plot([ic50, ic50], [0, 50], color=ABU, linewidth=1, linestyle="--")
    ax.text(ic50 + 0.002, 8, f"IC50 ≈ {koma(ic50, 3)}% v/v", color=TINTA, fontsize=10)
    ax.text(0.98, 0.06, f"y = {koma(slope, 1)}·x {koma(intersep, 2, tanda=True)}\nR² = {koma(r2, 4)}",
            transform=ax.transAxes, ha="right", va="bottom", fontsize=10, color=TINTA, family="monospace")
    ax.set_xlim(0, konsentrasi.max() * 1.08)
    ax.set_ylim(0, 100)
    ax.set_xlabel("Konsentrasi air sadapan (% v/v)", color=TINTA)
    ax.set_ylabel("% inhibisi DPPH", color=TINTA)
    ax.set_title("Aktivitas antioksidan DPPH", loc="left", fontsize=12, color=TINTA, fontweight="bold")
    ax.legend(frameon=False, loc="upper left", labelcolor=TINTA)
    gaya_sumbu(ax)
    fig.savefig(HASIL / "dpph.png", dpi=150, bbox_inches="tight", facecolor="white")
    plt.close(fig)

    return {"konsentrasi": konsentrasi, "inhibisi": inhibisi, "inhibisi_vitc": inhibisi_vitc,
            "slope": slope, "intersep": intersep, "R2": r2, "IC50": ic50,
            "R2_log": r2_log, "IC50_log": ic50_log}


def tulis_ringkasan(fenolik, flavonoid, dpph):
    baris = ["# Ringkasan hasil analisis ulang", "",
             "Dihasilkan otomatis oleh `analisis.py` dari data di folder `data/`.", ""]

    baris += ["## Kurva kalibrasi", "",
              "| Standar | Persamaan (koefisien penuh) | R² | Koefisien di Excel | LOD (ppm) | LOQ (ppm) | Sampel terbaca (ppm) |",
              "|---|---|---|---|---|---|---|"]
    for h in (fenolik, flavonoid):
        s_ex, b_ex = h["uji"]["koef_excel"]
        c = h["penuh"]["konsentrasi"]
        baris.append(f"| {h['uji']['standar']} | A = {koma(h['slope'], 6)}·C {koma(h['intersep'], 4, tanda=True)} "
                     f"| {koma(h['R2'], 4)} | A = {str(s_ex).replace('.', ',')}·C {koma(b_ex, 4, tanda=True)} "
                     f"| {koma(h['LOD'], 1)} | {koma(h['LOQ'], 1)} "
                     f"| {koma(c.min(), 1)}–{koma(c.max(), 1)} |")

    baris += ["", "## Kadar total", "",
              "| Uji | Skripsi (Bab IV) | Hitung ulang, koefisien Excel | Hitung ulang, koefisien penuh |",
              "|---|---|---|---|"]
    for h in (fenolik, flavonoid):
        sk = NILAI_SKRIPSI[h["nama"]]
        ex, pn = h["excel"], h["penuh"]
        sat = h["uji"]["satuan"]
        baris.append(f"| {h['uji']['judul']} | {koma(sk['rerata'], 3)} ± {koma(sk['SD'], 3)} {sat} (RSD {koma(sk['RSD'], 2)}%) "
                     f"| {koma(ex['rerata'], 3)} ± {koma(ex['SD'], 3)} (RSD {koma(ex['RSD'], 2)}%) "
                     f"| {koma(pn['rerata'], 3)} ± {koma(pn['SD'], 3)} (RSD {koma(pn['RSD'], 2)}%) |")
    rasio = NILAI_SKRIPSI["flavonoid"]["rerata"] / NILAI_SKRIPSI["fenolik"]["rerata"] * 100
    baris += ["", f"Rasio KTFla/KTF = {koma(rasio, 2)}%. Ini perbandingan ekuivalen kuersetin terhadap ekuivalen "
                  "asam galat dari dua reaksi berbeda, jadi dibaca sebagai rasio, bukan persen massa."]

    baris += ["", "## Aktivitas antioksidan DPPH", "",
              "| Konsentrasi (% v/v) | % inhibisi (hitung ulang) | % inhibisi (skripsi) |", "|---|---|---|"]
    for c, i, s in zip(dpph["konsentrasi"], dpph["inhibisi"], NILAI_SKRIPSI["dpph"]["inhibisi"]):
        baris.append(f"| {koma(c, 3)} | {koma(i, 2)} | {koma(s, 2)} |")
    baris.append(f"| vitamin C 3 ppm | {koma(dpph['inhibisi_vitc'], 2)} | {koma(NILAI_SKRIPSI['dpph']['vitamin_c'], 2)} |")
    baris += ["",
              f"IC50 (regresi linear, R² = {koma(dpph['R2'], 4)}): **{koma(dpph['IC50'], 4)}% v/v**. "
              f"Pembanding regresi terhadap log konsentrasi (R² = {koma(dpph['R2_log'], 4)}): {koma(dpph['IC50_log'], 4)}% v/v. "
              "IC50 tidak dilaporkan di skripsi; nilai ini dihitung ulang dari data yang sama.", ""]
    (HASIL / "ringkasan.md").write_text("\n".join(baris), encoding="utf-8")


if __name__ == "__main__":
    HASIL.mkdir(exist_ok=True)
    fenolik = analisis_kadar("fenolik")
    flavonoid = analisis_kadar("flavonoid")
    dpph = analisis_dpph()
    tulis_ringkasan(fenolik, flavonoid, dpph)

    for h in (fenolik, flavonoid):
        sat = h["uji"]["satuan"]
        print(f"{h['uji']['judul']}")
        print(f"  kurva: A = {h['slope']:.6f} C {h['intersep']:+.4f}, R² = {h['R2']:.4f}, "
              f"LOD {h['LOD']:.1f} ppm, LOQ {h['LOQ']:.1f} ppm")
        for label in ("excel", "penuh"):
            r = h[label]
            print(f"  koefisien {label:5s}: {r['rerata']:.3f} ± {r['SD']:.3f} {sat} (RSD {r['RSD']:.2f}%)")
    print(f"DPPH: IC50 {dpph['IC50']:.4f}% v/v (linear, R² {dpph['R2']:.4f}); "
          f"{dpph['IC50_log']:.4f}% v/v (log, R² {dpph['R2_log']:.4f}); vitamin C 3 ppm {dpph['inhibisi_vitc']:.2f}%")
    print(f"Grafik dan ringkasan tersimpan di {HASIL}")
