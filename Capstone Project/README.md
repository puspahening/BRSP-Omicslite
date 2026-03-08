# 🧬 Analisis Transkriptomik Tumor Prostat Pasca Terapi Poly-ICLC

> **Capstone Project — Program Bioteknologi**  
> Analisis ekspresi gen jaringan tumor kanker prostat sebelum dan sesudah terapi imunomodulator Poly-ICLC menggunakan dataset publik GSE288248.

---

## 📋 Deskripsi

Proyek ini mereplikasi dan memvalidasi temuan dari:

> Nair, S. et al. (2025). *Intratumoral poly-ICLC activates innate immunity and immune cell infiltration in localized prostate cancer.* **Med** (CellPress). https://doi.org/10.1016/j.medj.2025.100669

Poly-ICLC (Hiltonol®) adalah analog RNA rantai ganda sintetis yang bekerja sebagai agonis TLR3 dan MDA5 — meniru sinyal infeksi virus untuk mengaktifkan sistem imun di dalam tumor (*viral mimic*). Pada trial klinis ini, 12 pasien kanker prostat lokalisata mendapat injeksi Poly-ICLC intratumoral sebelum radical prostatectomy, menghasilkan pasangan sampel biopsi (pre) dan spesimen operasi (post) dari pasien yang sama.

**Pertanyaan utama:** Apakah Poly-ICLC mengubah profil ekspresi gen tumor prostat secara signifikan — khususnya pada jalur imun innate dan adaptif?

---

## 📊 Hasil Utama

| Analisis | Temuan |
|----------|--------|
| **DEG (Limma Paired)** | 17 gen signifikan (15 ↑, 2 ↓) · Top gen: *RGS1* (logFC=+1.457, FDR=6.2×10⁻⁸) |
| **GSEA Hallmark** | 26 dari 39 pathway signifikan (FDR<0.05) · semua NES positif |
| **Pathway dominan** | TNFα/NFκB (NES=2.752) · MYC Targets (NES=2.256) · Apoptosis (NES=1.974) |
| **Interferon** | IFN-γ signifikan (FDR=0.017) · IFN-α hanya tren (FDR=0.177) |
| **g:Profiler ORA** | 170+ term signifikan · top: *defense response to virus* (p=2.9×10⁻¹⁵) |
| **KEGG hsa05168** | IRF7, PKR/EIF2AK2, OAS1 aktif di jalur dsRNA sensing |

---

## 🗂️ Struktur Repositori

```
capstone-prostate-polyICLC/
│
├── README.md
│
├── script/
│   └── Script_Prostate_Paired_and_GSEA.R   # Script R lengkap (9 bagian)
│
├── data/
│   └── (data diunduh otomatis dari GEO via GEOquery)
│
├── output/
│   ├── Hasil_GSE288248_DEG_limma_paired.csv        # 17 DEG signifikan
│   ├── Hasil_GSE288248_all_genes_limma_paired.csv  # Semua gen (ranked)
│   ├── GSEA_Hallmark_GSE288248_post_vs_pre.csv     # Hasil GSEA lengkap
│   └── LeadingEdge_IFNa.csv                        # 29 leading-edge ISG
│
├── figures/
│   ├── Boxxplot_Paired_GSEA.png
│   ├── Density_Plot_Paired_GSEA.png
│   ├── UMAP_Plot_Paired_GSEA.png
│   ├── Volcano_Plot_Paired-GSEA.png
│   ├── Heatmap_Paired_GSEA.png
│   ├── GSEA_Hallmark.png
│   ├── Hallmark_Pathway_Enrichment.png
│   ├── GSEA_Hallmark_Inflammatory_Response.png
│   ├── Ridgeplot.png
│   └── Heatmap_Leading_Genes.png
│
└── laporan/
    └── Laporan_Capstone_FINAL_v4.docx
```

---

## ⚙️ Pipeline Analisis

```
GEO (GSE288248)
      │
      ▼
Preprocessing (log2 check, metadata paired)
      │
      ▼
Quality Control (Boxplot · Density Plot · UMAP)
      │
      ▼
Limma Paired DEG ──────────────────────────────► DEG CSV
(duplicateCorrelation + BH correction)
      │
      ▼
GSEA Hallmark MSigDB ──────────────────────────► GSEA CSV
(clusterProfiler, t-stat ranked list, set.seed=123)
      │
      ├──► Dotplot · Pathway Volcano · Enrichment Curve · Ridgeplot
      └──► Heatmap Leading-Edge IFN genes
      │
      ▼
Validasi
  ├── g:Profiler ORA (15 upregulated DEG → GO · KEGG · CORUM · TF)
  └── KEGG Pathway Map hsa05168 (HSV-1 Infection)
```

---

## 🛠️ Cara Menjalankan

### 1. Prasyarat

Pastikan sudah terinstal **R** (≥ 4.2) dan **RStudio**.

### 2. Install package

Package akan diinstal otomatis saat script pertama kali dijalankan. Untuk instalasi manual:

```r
# Bioconductor
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(c("GEOquery", "limma"))

# CRAN
install.packages(c("pheatmap", "ggplot2", "dplyr", "umap",
                   "clusterProfiler", "msigdbr", "enrichplot"))
```

### 3. Jalankan script

```r
source("script/Script_Prostate_Paired_and_GSEA.R")
```

Data GSE288248 akan **diunduh otomatis** dari NCBI GEO — tidak perlu unduh manual. Pastikan koneksi internet aktif saat pertama kali menjalankan.

### 4. Output yang dihasilkan

| File | Deskripsi |
|------|-----------|
| `Hasil_GSE288248_DEG_limma_paired.csv` | DEG signifikan (FDR<0.05, \|logFC\|≥0.5) |
| `Hasil_GSE288248_all_genes_limma_paired.csv` | Seluruh gen beserta t-statistik |
| `GSEA_Hallmark_GSE288248_post_vs_pre.csv` | Hasil GSEA (NES, FDR, leading-edge) |
| `LeadingEdge_IFNa.csv` | 29 gen leading-edge IFN-α |
| `*.png` | Semua visualisasi (10 gambar) |

---

## 📦 Package yang Digunakan

| Package | Fungsi | Sumber |
|---------|--------|--------|
| `GEOquery` | Unduh dataset dari NCBI GEO | Bioconductor |
| `limma` | Analisis ekspresi diferensial paired | Bioconductor |
| `clusterProfiler` | GSEA | Bioconductor |
| `msigdbr` | Koleksi Hallmark MSigDB | CRAN |
| `enrichplot` | Visualisasi GSEA | Bioconductor |
| `pheatmap` | Heatmap | CRAN |
| `ggplot2` | Visualisasi | CRAN |
| `umap` | UMAP dimensionality reduction | CRAN |
| `dplyr` | Manipulasi data | CRAN |

---

## 🔬 Detail Dataset

| | |
|---|---|
| **GEO Accession** | [GSE288248](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE288248) |
| **Platform** | GPL5175 · Affymetrix Human Exon 1.0 ST Array |
| **Jumlah sampel** | 24 (12 pre-treatment + 12 post-treatment) |
| **Tipe jaringan** | FFPE tumor prostat (biopsi vs spesimen RP) |
| **Desain** | Paired — 12 pasien kanker prostat lokalisata |
| **Intervensi** | Injeksi Poly-ICLC (Hiltonol®) intratumoral |

---

## 📈 Ringkasan Temuan Biologis

Poly-ICLC mengaktifkan program **antiviral innate immunity** di jaringan tumor prostat, yang ditandai oleh:

- **TNFα/NFκB** sebagai jalur paling dominan → konsisten dengan aktivasi TLR3/MDA5
- **IRF7, PKR, OAS1** aktif di jalur dsRNA sensing (KEGG hsa05168)
- **IFN-γ** signifikan; **IFN-α** menunjukkan heterogenitas respons antar pasien
- **Apoptosis** dan **Androgen Response** teraktivasi → selaras dengan Fig. 1F Nair et al.
- g:Profiler: top term *"defense response to virus"* (p=2.9×10⁻¹⁵) mengkonfirmasi mekanisme *viral mimic*

Seluruh temuan **konsisten dengan Nair et al. (2025)** yang dianalisis menggunakan dataset klinis internal.

---

## ⚠️ Keterbatasan

- Jumlah sampel kecil (n=12) — membatasi *statistical power*
- Perbedaan jenis jaringan: biopsi (pre) vs spesimen RP (post) → potensi bias ischemia/FFPE
- Belum terintegrasi dengan data NanoString darah dari dataset GSE294391

---

## 📄 Lisensi

Kode dalam repositori ini menggunakan lisensi **MIT**.  
Dataset GSE288248 tersedia secara publik di [NCBI GEO](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE288248).

---

## 📚 Referensi

1. Nair, S. et al. (2025). Intratumoral poly-ICLC activates innate immunity and immune cell infiltration in localized prostate cancer. *Med* (CellPress). https://doi.org/10.1016/j.medj.2025.100669

2. Ritchie, M.E. et al. (2015). limma powers differential expression analyses for RNA-sequencing and microarray studies. *Nucleic Acids Research*, 43(7), e47.

3. Wu, T. et al. (2021). clusterProfiler 4.0: A universal enrichment tool for interpreting omics data. *Innovation*, 2(3), 100141.

4. Liberzon, A. et al. (2015). The Molecular Signatures Database Hallmark Gene Set Collection. *Cell Systems*, 1(6), 417-425.

---
