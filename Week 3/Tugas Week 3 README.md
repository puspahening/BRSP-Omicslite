# Analisis Ekspresi Gen Diferensial pada Kanker Payudara — Dataset GSE42568

**Oleh:** Puspa Hening  
**Mata Kuliah / Kegiatan:** Bioinformatika — Laporan Minggu 3  
**Platform Data:** Gene Expression Omnibus (GEO), NCBI

---

## 📌 Deskripsi

Laporan ini menyajikan analisis ekspresi gen diferensial (*Differentially Expressed Genes*/DEGs) menggunakan dataset publik **GSE42568** dari GEO. Tujuannya adalah mengidentifikasi gen-gen yang mengalami perubahan aktivitas pada jaringan kanker payudara dibandingkan jaringan normal, serta memahami makna biologisnya melalui analisis enrichment.

---

## 📊 Dataset

| Atribut | Detail |
|---|---|
| **GEO Accession** | GSE42568 |
| **Total Sampel** | 121 sampel |
| **Sampel Kanker** | 104 jaringan kanker payudara invasif |
| **Sampel Normal** | 17 jaringan payudara normal |
| **Platform Microarray** | Affymetrix Human Genome U133 Plus 2.0 (GPL570) |

---

## 🔬 Metode

| Tahapan | Tools / Paket R |
|---|---|
| Download & manajemen data | `GEOquery`, `Biobase` |
| Normalisasi (skala log₂) | Bioconductor |
| Quality Control (boxplot, density plot) | `ggplot2` |
| Reduksi dimensi | UMAP |
| Identifikasi DEGs | `limma` (empirical Bayes moderated t-test) |
| Visualisasi (volcano plot, heatmap) | `ggplot2`, `pheatmap` |
| Analisis enrichment (GO & KEGG) | `clusterProfiler`, `enrichplot`, `org.Hs.eg.db` |

**Kriteria DEGs:**
- Adjusted p-value < 0.05 (koreksi Benjamini–Hochberg / FDR)
- |log₂ Fold Change| ≥ 1

---

## 📁 Struktur Laporan

```
├── 1. Pendahuluan
├── 2. Metode
│   ├── 2.1 Dataset
│   ├── 2.2 Pra-pemrosesan Data
│   ├── 2.3 Analisis DEGs (limma)
│   ├── 2.4 Visualisasi (Volcano Plot, Heatmap)
│   └── 2.5 Analisis Enrichment (GO-BP, GO-MF, KEGG)
├── 3. Hasil dan Interpretasi
│   ├── 3.1 Eksplorasi Kualitas Data (Boxplot, Density Plot, UMAP)
│   ├── 3.2 Volcano Plot
│   ├── 3.3 Heatmap 50 DEGs Teratas
│   └── 3.4 Analisis Enrichment
└── 4. Kesimpulan
```

---

## 🔑 Temuan Utama

### DEGs
- **Dominasi downregulation**: lebih banyak gen yang mengalami penekanan ekspresi pada kanker payudara dibanding peningkatan (logFC hingga −7 s/d −8 vs. maksimal ~+5/+6).
- Pemisahan transkriptomik kanker vs. normal dikonfirmasi secara visual melalui **UMAP**.

### 50 DEGs Teratas (Heatmap)
Gen-gen yang menonjol beserta relevansi biologisnya:

| Gen | Fungsi | Relevansi Kanker |
|---|---|---|
| FOXP2 | Faktor transkripsi perkembangan | Dediferensiasi sel |
| CA4 | Carbonic Anhydrase 4 | Adaptasi metabolik kanker |
| ACADL | Acyl-CoA Dehydrogenase | Penurunan oksidasi asam lemak (efek Warburg) |
| CIDEA | Cell death-inducing effector A | Resistensi apoptosis |
| GSN | Gelsolin (pengikatan aktin) | Invasivitas tumor |
| SLC16A7 | Monocarboxylate transporter 2 | Transpor laktat |
| ACSM5 | Acyl-CoA Synthetase | Perubahan metabolisme lipid |
| ADH1C | Alcohol/retinol dehydrogenase | Sering diturunkan pada kanker payudara |

### Analisis Enrichment

**GO – Biological Process (GO-BP):**
- Actin filament organization (EMT)
- Cell-substrate & cell-matrix adhesion
- Wound healing
- Jalur PI3K/AKT signaling

**GO – Molecular Function (GO-MF):**
- Actin binding
- GTPase regulator/activator activity (Rho/Ras)
- Extracellular matrix structural constituent
- Integrin binding & cadherin binding (penanda EMT)

**KEGG Pathway:**
- Focal adhesion & ECM-receptor interaction
- Regulation of actin cytoskeleton
- Proteoglycans in cancer
- FoxO signaling (tumor suppressor)
- Reprogramming metabolisme lipid (Fatty acid metabolism, TCA cycle)
- Insulin resistance & PPAR signaling

---

## 💡 Kesimpulan

Temuan analisis konsisten dengan *hallmarks* kanker payudara yang dikenal dalam literatur:
1. **Disregulasi adhesi sel** — perubahan luas pada interaksi sel-ECM
2. **Peningkatan motilitas & invasivitas** — melalui EMT dan reorganisasi sitoskeleton
3. **Aktivasi sinyal pro-survival** — jalur PI3K/AKT
4. **Reprogramming metabolik** — penurunan oksidasi asam lemak, pergeseran ke metabolisme glikolitik (efek Warburg)

Gen dan jalur yang teridentifikasi merupakan kandidat biomarker diagnostik maupun target terapeutik potensial untuk penelitian lanjutan.

---

## 💻 Script R

Script lengkap tersedia di GitHub:  
[Script R — Breast Cancer Analysis & Enrichment](https://github.com/puspahening/BRSP-Omicslite/blob/main/Script%20R%20-%20Breast%20Cancer%203%20and%20Enrichment%20R.R)

---

## 🛠️ Dependensi R

```r
# Instalasi paket yang diperlukan
BiocManager::install(c("GEOquery", "Biobase", "limma", 
                       "clusterProfiler", "enrichplot", 
                       "org.Hs.eg.db"))

install.packages(c("ggplot2", "pheatmap", "umap"))
```
