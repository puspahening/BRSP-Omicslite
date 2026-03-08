############################################################
# ANALISIS TRANSCRIPTOMICS GSE288248 (Prostate Cancer)
# pre vs post poly-ICLC (PAIRED DESIGN) + GSEA (Hallmark)
# Platform: GPL5175 Affymetrix Human Exon 1.0 ST Array
############################################################

#############################
# PART A. SETUP & PACKAGE
#############################

# 1) Bioconductor manager
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

# 2) Install Bioconductor packages
BiocManager::install(c("GEOquery", "limma"), ask = FALSE, update = FALSE)

# 3) Install CRAN packages
pkgs_cran <- c("pheatmap", "ggplot2", "dplyr", "umap", "clusterProfiler", "msigdbr", "enrichplot")
to_install <- pkgs_cran[!pkgs_cran %in% rownames(installed.packages())]
if (length(to_install) > 0) install.packages(to_install)

# 4) Load libraries
library(GEOquery)
library(limma)
library(pheatmap)
library(ggplot2)
library(dplyr)
library(umap)
library(clusterProfiler)
library(msigdbr)
library(enrichplot)

#############################
# PART B. DOWNLOAD DATA GEO
#############################

gset <- getGEO("GSE288248", GSEMatrix = TRUE, AnnotGPL = TRUE)[[1]]

ex <- exprs(gset)
info <- pData(gset)
feature_data <- fData(gset)

cat("Dimensi matriks ekspresi (probe x sample):\n")
print(dim(ex))

#############################
# PART C. LOG2 TRANSFORM (IF NEEDED)
#############################

qx <- as.numeric(quantile(ex, c(0, 0.25, 0.5, 0.75, 0.99, 1), na.rm = TRUE))
LogTransform <- (qx[5] > 100) || (qx[6] - qx[1] > 50 && qx[2] > 0)

if (LogTransform) {
  ex[ex <= 0] <- NA
  ex <- log2(ex)
  cat("Log2 transform diterapkan.\n")
} else {
  cat("Log2 transform tidak diperlukan (data sudah log-like).\n")
}

#############################
# PART D. DEFINE GROUP (PRE/POST)
#############################

stopifnot("treatment time-point:ch1" %in% colnames(info))

group_info <- info$`treatment time-point:ch1`
cat("\nNilai unik treatment time-point:\n")
print(unique(group_info))

groups <- make.names(group_info)
group <- factor(groups, levels = c("pre_treat", "post_treat"))
stopifnot(all(!is.na(group)))

cat("\nDistribusi group:\n")
print(table(group))

#############################
# PART E. DEFINE PATIENT ID (PAIRED)
#############################

# Cari kolom ID pasien yang umum di GEO
candidate_cols <- c(
  "subject:ch1", "Subject:ch1",
  "patient:ch1", "Patient:ch1",
  "individual:ch1", "Individual:ch1",
  "donor:ch1", "Donor:ch1",
  "participant:ch1", "Participant:ch1"
)

patient_col <- candidate_cols[candidate_cols %in% colnames(info)][1]

if (!is.na(patient_col) && length(patient_col) == 1) {
  patient <- factor(info[[patient_col]])
  cat("\nPatient ID diambil dari kolom metadata: ", patient_col, "\n", sep = "")
} else {
  # fallback: ekstrak dari title (contoh: Subject-1 / Subject 1)
  ttl <- as.character(info$title)
  patient_extracted <- gsub(".*(Subject|SUBJECT|Patient|PATIENT)[ _-]*([0-9]+).*", "S\\2", ttl)
  
  if (any(grepl("^S[0-9]+$", patient_extracted))) {
    patient <- factor(patient_extracted)
    cat("\nPatient ID diekstrak dari kolom 'title'.\n")
  } else {
    stop(
      "Gagal menentukan patient ID.\n",
      "Silakan cek pData(gset) untuk kolom ID pasien, lalu set 'patient' secara manual."
    )
  }
}

cat("\nJumlah pasien unik:\n")
print(length(unique(patient)))

# cek pairing (harus kira-kira 2 sampel per pasien)
cat("\nCek pairing (jumlah sampel per pasien):\n")
print(table(patient))



#############################
# PART G. LIMMA PAIRED (duplicateCorrelation)
#############################

design <- model.matrix(~0 + group)
colnames(design) <- levels(group)

corfit <- duplicateCorrelation(ex, design, block = patient)

fit <- lmFit(ex, design, block = patient, correlation = corfit$consensus)

# Kontras: post - pre (positif = naik setelah terapi)
cm <- makeContrasts(post_treat - pre_treat, levels = design)

fit2 <- eBayes(contrasts.fit(fit, cm))

tt <- topTable(fit2, number = Inf, adjust.method = "BH", sort.by = "P")
tt$PROBEID <- rownames(tt)

cat("\nTop genes (unannotated):\n")
print(head(tt))

#############################
# PART H. ANNOTATE GENE SYMBOL
#############################

extract_symbol <- function(x) {
  x <- as.character(x)
  if (is.na(x) || x == "---" || x == "") return(NA_character_)
  parts <- strsplit(x, " // ", fixed = TRUE)[[1]]
  if (length(parts) >= 2) return(trimws(parts[2]))
  return(NA_character_)
}

if ("gene_assignment" %in% colnames(feature_data)) {
  feature_data$SYMBOL <- vapply(feature_data$gene_assignment, extract_symbol, FUN.VALUE = character(1))
} else if ("Gene Symbol" %in% colnames(feature_data)) {
  feature_data$SYMBOL <- as.character(feature_data$`Gene Symbol`)
} else if ("SYMBOL" %in% colnames(feature_data)) {
  feature_data$SYMBOL <- as.character(feature_data$SYMBOL)
} else {
  feature_data$SYMBOL <- NA_character_
}

tt$SYMBOL <- feature_data[tt$PROBEID, "SYMBOL"]

# Optional: label fallback jika SYMBOL kosong
tt$SYMBOL2 <- ifelse(is.na(tt$SYMBOL) | tt$SYMBOL=="", tt$PROBEID, tt$SYMBOL)

cat("\nTop genes (annotated):\n")
print(head(tt[, c("SYMBOL2","logFC","P.Value","adj.P.Val","t")]))

#############################
# PART I. DEFINE DEG (FLEXIBLE THRESHOLD)
#############################

padj_cut <- 0.05
lfc_cut  <- 0.5

deg <- subset(tt, adj.P.Val < padj_cut & abs(logFC) >= lfc_cut)

cat("\nJumlah DEG (FDR<", padj_cut, " & |logFC|>=", lfc_cut, "): ", nrow(deg), "\n", sep = "")
if (nrow(deg) > 0) print(head(deg[, c("SYMBOL2","logFC","adj.P.Val")]))

#############################
# PART F. QC PLOTS (BOXPLOT & DENSITY)
#############################

# Boxplot per sample (urut group)
# Tambah margin kanan
par(mar = c(5, 4, 4, 8))   # bottom, left, top, right

ord <- order(group)

boxplot(ex[, ord],
        col = ifelse(group[ord] == "pre_treat", "black", "firebrick"),
        outline = FALSE,
        las = 2,
        ylab = "Expression (log2)",
        main = "Boxplot Distribusi Ekspresi per Sampel")

# Izinkan plotting di luar area
par(xpd = TRUE)

legend("topright",
       inset = c(-0.25, 0),   # geser keluar ke kanan
       legend = levels(group),
       fill   = c("black","firebrick"),
       bty    = "n",
       cex    = 0.9)

# Density plot
expr_long <- data.frame(
  Expression = as.vector(ex),
  Group = rep(group, each = nrow(ex))
)

ggplot(expr_long, aes(x = Expression, color = Group)) +
  geom_density(linewidth = 1) +
  theme_minimal() +
  labs(title="Distribusi Nilai Ekspresi Gen", x="Expression Value (log2)", y="Density")

#############################
# PART J. VOLCANO PLOT
#############################

logFC <- tt$logFC
padj  <- tt$adj.P.Val
padj[padj == 0] <- .Machine$double.xmin
y <- -log10(padj)

is_sig  <- padj < padj_cut
is_up   <- is_sig & (logFC >=  lfc_cut)
is_down <- is_sig & (logFC <= -lfc_cut)
is_ns   <- !(is_up | is_down)

xmax <- max(abs(logFC), na.rm = TRUE)
plot(logFC[is_ns], y[is_ns],
     pch=16, cex=0.55, col="black",
     xlim=c(-xmax, xmax),
     xlab="log2FC (post - pre)",
     ylab="-log10(FDR)",
     main="Differential Expression: Post vs Pre")

points(logFC[is_down], y[is_down], pch=16, cex=0.75, col="deepskyblue3")
points(logFC[is_up],   y[is_up],   pch=16, cex=0.75, col="red")

abline(v=c(-lfc_cut, lfc_cut), lty=2)
abline(h=-log10(padj_cut), lty=2)

par(xpd = TRUE)

legend("topright",
       inset = c(-0.28, 0),
       legend = c(
         paste0("Up (", sum(is_up), ")"),
         paste0("Down (", sum(is_down), ")"),
         paste0("FDR < ", padj_cut),
         paste0("|logFC| ≥ ", lfc_cut)
       ),
       col = c("#d7301f", "#2b8cbe", "black", "black"),
       pch = c(16, 16, NA, NA),
       lty = c(NA, NA, 2, 2),
       pt.cex = 0.8,
       bty = "n")

par(xpd = FALSE)


#############################
# PART K. HEATMAP TOP 50 (post - pre)
#############################

topN <- 50
tt_sorted <- tt[order(tt$adj.P.Val), ]
topN_tbl <- head(tt_sorted, topN)

mat_heatmap <- ex[topN_tbl$PROBEID, , drop = FALSE]

gene_label <- topN_tbl$SYMBOL2
rownames(mat_heatmap) <- gene_label

# remove NA rows and zero-variance rows
mat_heatmap <- mat_heatmap[rowSums(is.na(mat_heatmap)) == 0, , drop = FALSE]
gene_variance <- apply(mat_heatmap, 1, var)
mat_heatmap <- mat_heatmap[gene_variance > 0, , drop = FALSE]

annotation_col <- data.frame(Group = group)
rownames(annotation_col) <- colnames(mat_heatmap)

pheatmap(
  mat_heatmap,
  scale = "row",
  annotation_col = annotation_col,
  show_colnames = FALSE,
  show_rownames = TRUE,
  fontsize_row = 7,
  clustering_method = "complete",
  main = paste0("Top ", topN, " DE Genes (post - pre)")
)

#############################
# PART L. UMAP (OPTIONAL)
#############################

umap_input <- t(ex)
umap_result <- umap(umap_input)

umap_df <- data.frame(
  UMAP1 = umap_result$layout[, 1],
  UMAP2 = umap_result$layout[, 2],
  Group = group
)

ggplot(umap_df, aes(x = UMAP1, y = UMAP2, color = Group)) +
  geom_point(size = 4, alpha = 0.9) +
  theme_minimal(base_size = 13) +
  stat_ellipse(level = 0.8, linetype = "dashed") +
  labs(title="UMAP Plot — Prostate Cancer", subtitle="Pre vs Post Treatment", x="UMAP 1", y="UMAP 2")

#############################
# PART M. GSEA (HALLMARK) — SESUAI PAPER
#############################

# --- Buat ranked gene list (HARUS named numeric vector) ---
gene_list <- tt$t
names(gene_list) <- tt$SYMBOL

# buang NA/kosong
keep <- !is.na(names(gene_list)) & names(gene_list) != "" & !is.na(gene_list)
gene_list <- gene_list[keep]

# pastikan numeric vector (bukan array/matrix)
gene_list <- as.numeric(gene_list)
names(gene_list) <- tt$SYMBOL[keep]

# kalau ada duplikat SYMBOL, ambil |t| terbesar (wajib supaya names unik)
gene_list <- tapply(gene_list, names(gene_list), function(x) x[which.max(abs(x))])

# tapply menghasilkan "array" -> ubah jadi numeric vector
gene_list <- as.numeric(gene_list)
names(gene_list) <- names(tapply(tt$t[keep], tt$SYMBOL[keep], function(x) x[which.max(abs(x))]))

# rapikan dan sort decreasing
gene_list <- sort(gene_list, decreasing = TRUE)

# cek cepat: harus TRUE
is.numeric(gene_list)
is.null(dim(gene_list))
head(gene_list)

m_df <- msigdbr(species = "Homo sapiens", collection = "H")
hallmark <- m_df[, c("gs_name", "gene_symbol")]

set.seed(123)
gsea_h <- GSEA(
  geneList      = gene_list,
  TERM2GENE     = hallmark,
  pvalueCutoff  = 0.2,
  verbose       = FALSE
)

res_h <- gsea_h@result
res_h <- res_h[order(res_h$p.adjust), ]
head(res_h[, c("ID","NES","p.adjust")], 20)

############################################################
# PART 4 — PLOT 1: DOTPLOT (ringkasan top pathways)
############################################################

p_dot <- dotplot(gsea_h, showCategory = 15) +
  ggtitle("GSEA Hallmark — Top Enriched Pathways (post vs pre)") +
  theme_minimal()

print(p_dot)

############################################################
# PART 5 — PLOT 2: ENRICHMENT CURVES (paper main figures)
############################################################
# Ganti ID sesuai yang muncul pada res_h$ID bila perlu

targets <- c(
  "HALLMARK_INTERFERON_ALPHA_RESPONSE",
  "HALLMARK_INTERFERON_GAMMA_RESPONSE",
  "HALLMARK_TNFA_SIGNALING_VIA_NFKB",
  "HALLMARK_INFLAMMATORY_RESPONSE"
)

for (pid in targets) {
  if (pid %in% res_h$ID) {
    p <- gseaplot2(gsea_h, geneSetID = pid,
                   title = paste0("GSEA — ", pid, " (post vs pre)"))
    print(p)
  } else {
    message("Pathway tidak ditemukan: ", pid)
  }
}

############################################################
# PART 6 — PLOT 3: "VOLCANO" PATHWAY-LEVEL (NES vs -log10(FDR))
############################################################

# 1) Ambil hasil GSEA
res <- as.data.frame(gsea_h@result)

# 2) Parameter
fdr_cut <- 0.05
top_n   <- 12   # jumlah pathway yang dilabel (ubah 10–15 sesuai kebutuhan)

# 3) Buat kategori significance (biar legend tidak TRUE/FALSE)
res$Significance <- ifelse(res$p.adjust < fdr_cut,
                           "Significant",
                           "Not significant")

# 4) Siapkan kolom -log10(FDR)
res$neglogFDR <- -log10(res$p.adjust)

# 5) Data label: hanya TOP N yang signifikan (biar rapi)
lab_df <- subset(res, p.adjust < fdr_cut)
lab_df <- lab_df[order(lab_df$p.adjust), ]
lab_df <- head(lab_df, top_n)

# 6) Plot
p <- ggplot(res, aes(x = NES, y = neglogFDR)) +
  geom_point(aes(color = Significance), size = 3, alpha = 0.9) +
  scale_color_manual(
    values = c("Not significant" = "grey70",
               "Significant"     = "#1b9e77")
  ) +
  geom_hline(yintercept = -log10(fdr_cut), linetype = "dashed") +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_text_repel(
    data = lab_df,
    aes(label = Description),
    size = 3.2,             # kecilkan kalau masih padat (mis. 2.8)
    box.padding = 0.5,
    point.padding = 0.35,
    max.overlaps = Inf,
    min.segment.length = 0,
    segment.alpha = 0.6
  ) +
  scale_y_continuous(expand = expansion(mult = c(0.02, 0.15))) +
  coord_cartesian(clip = "off") +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(size = 18, face = "bold"),
    legend.title = element_text(size = 13),
    legend.text  = element_text(size = 12)
  ) +
  labs(
    title = "Hallmark Pathway Enrichment After Poly-ICLC",
    x = "NES (Post vs Pre)",
    y = "-log10(FDR)",
    color = "Pathway Status"
  )

print(p)


############################################################
# PART 7 — PLOT 4: RIDGEPLOT (opsional, ringkasan distribusi)
############################################################

p_ridge <- ridgeplot(gsea_h, showCategory = 15) +
  ggtitle("GSEA Hallmark — Ridgeplot (Top 15)") +
  theme_minimal()

print(p_ridge)

############################################################
# PART 8 — HEATMAP LEADING-EDGE GENES (GSEA-driven) [ADVANCED]
############################################################
# Ini butuh matriks ekspresi level gene (ex_gene), bukan probe.
# Kalau kamu belum punya ex_gene, jalankan blok "BUILD ex_gene" di bawah.
#
# Input: ex (probe x sample), feature_data$SYMBOL, group (pre/post)
############################################################

# ===== BUILD ex_gene (collapse probe -> gene by highest variance) =====
# SKIP jika kamu sudah punya ex_gene

if (!exists("ex_gene")) {
  stopifnot(exists("ex"), exists("feature_data"))
  probe2gene <- feature_data$SYMBOL
  names(probe2gene) <- rownames(feature_data)
  
  gene_vec <- probe2gene[rownames(ex)]
  keep2 <- !is.na(gene_vec) & gene_vec != ""
  
  ex2 <- ex[keep2, , drop = FALSE]
  gene_vec2 <- gene_vec[keep2]
  
  vars <- apply(ex2, 1, var, na.rm = TRUE)
  df <- data.frame(gene = gene_vec2, var = vars, row = seq_len(nrow(ex2)))
  
  idx <- df %>%
    group_by(gene) %>%
    slice_max(order_by = var, n = 1, with_ties = FALSE) %>%
    pull(row)
  
  ex_gene <- ex2[idx, , drop = FALSE]
  rownames(ex_gene) <- gene_vec2[idx]
  
  cat("ex_gene dibuat. Dim:", dim(ex_gene), "\n")
}

# ===== helper ambil leading-edge genes =====
get_leading_edge <- function(gsea_res, pathway_id) {
  row <- gsea_res[gsea_res$ID == pathway_id, , drop = FALSE]
  if (nrow(row) == 0) return(character(0))
  le <- row$core_enrichment[1]
  unlist(strsplit(le, "/"))
}

# Ambil leading-edge IFN alpha + IFN gamma (kalau ada)
le_ifna <- if ("HALLMARK_INTERFERON_ALPHA_RESPONSE" %in% res_h$ID)
  get_leading_edge(res_h, "HALLMARK_INTERFERON_ALPHA_RESPONSE") else character(0)

le_ifng <- if ("HALLMARK_INTERFERON_GAMMA_RESPONSE" %in% res_h$ID)
  get_leading_edge(res_h, "HALLMARK_INTERFERON_GAMMA_RESPONSE") else character(0)

le_genes <- unique(c(le_ifna, le_ifng))
le_genes <- le_genes[le_genes %in% rownames(ex_gene)]

cat("Leading-edge genes tersedia:", length(le_genes), "\n")

if (length(le_genes) >= 5) {
  # batasi agar heatmap enak (mis. 50 gen)
  if (length(le_genes) > 50) le_genes <- le_genes[1:50]
  
  mat_le <- ex_gene[le_genes, , drop = FALSE]
  
  annotation_col <- data.frame(Group = group)
  rownames(annotation_col) <- colnames(mat_le)
  
  pheatmap(mat_le,
           scale = "row",
           annotation_col = annotation_col,
           show_colnames = FALSE,
           fontsize_row = 7,
           main = "Perubahan Ekspresi Gen Jalur Interferon (Pre vs Post Poly-ICLC)")
} else {
  message("Leading-edge IFN genes terlalu sedikit / pathway IFN tidak ditemukan. Heatmap dilewati.")
}

############################################################
# PART 9 — EXPORT: leading-edge genes to CSV
############################################################

if (length(le_ifna) > 0) write.csv(data.frame(LeadingEdge_IFNa = le_ifna),
                                   "LeadingEdge_IFNa.csv", row.names = FALSE)
if (length(le_ifng) > 0) write.csv(data.frame(LeadingEdge_IFNg = le_ifng),
                                   "LeadingEdge_IFNg.csv", row.names = FALSE)

message("\nSelesai ✅ Plot sudah dibuat & output tersimpan.")



#############################
# PART N. SAVE OUTPUTS
#############################

write.csv(tt,  "Hasil_GSE288248_all_genes_limma_paired.csv", row.names = FALSE)
write.csv(deg, "Hasil_GSE288248_DEG_limma_paired.csv",       row.names = FALSE)
write.csv(res_h, "GSEA_Hallmark_GSE288248_post_vs_pre.csv",  row.names = FALSE)

message("\nSELESAI ✅ Output tersimpan:")
message("- Hasil_GSE288248_all_genes_limma_paired.csv")
message("- Hasil_GSE288248_DEG_limma_paired.csv")
message("- GSEA_Hallmark_GSE288248_post_vs_pre.csv")