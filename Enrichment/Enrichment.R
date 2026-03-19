# ==============================================================================
# File: Enrichment.R
# Author: Yuanwei
# Date: 2026-03-19
# Description:
#   AI Skill 核心执行脚本 —— 智能基因富集分析引擎。
#   架构分为两层：
#     Layer 1 (CLI Adapter): 解析命令行参数，自动识别输入格式，
#                            将各种来源的数据洗为标准内部结构。
#     Layer 2 (Core Engine): 接受标准基因向量，执行 ORA / GSEA 富集分析
#                            并生成高质量可视化图表。
#
# Usage:
#   Rscript Enrichment.R --help
#   Rscript Enrichment.R --input <file> --gene_col <col> [--metric_col <col>] [options]
# ==============================================================================


# ==============================================================================
# §0  包加载 + CLI 参数定义
# ==============================================================================

# --- CRAN 包自动安装 Bootstrap（首次运行自动处理缺失包）---------------------
.cran_pkgs <- c("optparse", "ggplot2", "ggridges", "dplyr",
                "forcats", "stringr", "data.table")
.missing   <- .cran_pkgs[!sapply(.cran_pkgs, requireNamespace, quietly = TRUE)]
if (length(.missing) > 0) {
  message(">>> [Bootstrap] Installing missing CRAN packages: ",
          paste(.missing, collapse = ", "))
  install.packages(.missing, repos = "https://cloud.r-project.org", quiet = TRUE)
}

# Bioconductor 包需要单独处理
.bioc_pkgs <- c("clusterProfiler", "enrichplot")
.missing_bioc <- .bioc_pkgs[!sapply(.bioc_pkgs, requireNamespace, quietly = TRUE)]
if (length(.missing_bioc) > 0) {
  message(">>> [Bootstrap] Installing missing Bioconductor packages: ",
          paste(.missing_bioc, collapse = ", "))
  if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager", repos = "https://cloud.r-project.org", quiet = TRUE)
  }
  BiocManager::install(.missing_bioc, ask = FALSE, update = FALSE)
}

suppressPackageStartupMessages({
  library(optparse)
  library(clusterProfiler)
  library(enrichplot)
  library(ggplot2)
  library(ggridges)
  library(dplyr)
  library(forcats)
  library(stringr)
  library(data.table)
})

options(timeout = 300)
options(clusterProfiler.kegg_mirror = "https://rest.kegg.cn/")

# --- CLI 参数定义 -----------------------------------------------------------
option_list <- list(
  # 输入
  make_option(c("-i", "--input"),
    type = "character", default = NULL,
    help = "[必填] 输入文件路径（支持 csv / tsv / txt，自动识别分隔符）"),
  make_option(c("-g", "--gene_col"),
    type = "character", default = NULL,
    help = "[可选] 基因名所在列名。未指定时自动探测含 'gene'/'symbol' 的列，或退回第1列。"),
  make_option(c("-m", "--metric_col"), type = "character", default = NULL,
              help = "[可选] 排序指标列名（如 log2FoldChange、logFC 等）。
           提供后且 --mode 为 auto/GSEA 时触发 GSEA。"),
  make_option(c("--mode"), type = "character", default = "auto",
              help = "[可选] 分析模式: auto, ORA, GSEA。默认 auto (有排序指标则 GSEA，否则 ORA)"),
  make_option(c("-p", "--pval_col"), type = "character", default = NULL,
    help = "[可选] ORA 模式下用于过滤的 p 值列（如 padj、adj.P.Val）。
           未提供则假设输入文件已经是显著基因集，不做额外过滤。"),
  make_option(c("--pval_cutoff"),
    type = "double", default = 0.05,
    help = "[可选] ORA p 值过滤阈值，默认 0.05 [default: %default]"),
  make_option(c("-b", "--background"),
    type = "character", default = NULL,
    help = "[可选] ORA 背景基因文件路径（单列 txt / csv 均可）"),
  make_option(c("--organism"),
    type = "character", default = "ssc",
    help = "[可选] 物种 KEGG 代码，默认 ssc (Sus scrofa / Pig) [default: %default]
           常用值: hsa (Human), mmu (Mouse), rno (Rat), bta (Cattle)"),
  make_option(c("--orgdb"),
    type = "character", default = "org.Ss.eg.db",
    help = "[可选] 对应 OrgDb R 包名，默认 org.Ss.eg.db [default: %default]
           常用值: org.Hs.eg.db, org.Mm.eg.db, org.Rn.eg.db, org.Bt.eg.db"),
  make_option(c("-o", "--outdir"),
    type = "character", default = "./enrichment_out",
    help = "[可选] 结果输出目录，默认 ./enrichment_out [default: %default]"),
  make_option(c("--task_id"),
    type = "character", default = "Enrichment",
    help = "[可选] 任务标识符，用于图表标题，默认 Enrichment [default: %default]"),
  # 分析参数
  make_option(c("--kegg_p_cutoff"),
    type = "double", default = 0.05,
    help = "[可选] KEGG 富集 p 值阈值 [default: %default]"),
  make_option(c("--go_p_cutoff"),
    type = "double", default = 0.05,
    help = "[可选] GO 富集 p 值阈值 [default: %default]"),
  make_option(c("--go_q_cutoff"),
    type = "double", default = 0.2,
    help = "[可选] GO 富集 q 值阈值 [default: %default]"),
  make_option(c("--show_n"),
    type = "integer", default = 15,
    help = "[可选] 各图表展示条目数 [default: %default]")
)

opt <- parse_args(OptionParser(
  option_list  = option_list,
  prog         = "Rscript Enrichment.R",
  description  = "智能基因富集分析 Skill — 支持 ORA 和 GSEA 双模式，多物种，自动 ID 探测转换。"
))

# 检查必填参数
if (is.null(opt$input)) {
  stop("ERROR: --input 参数为必填项！请使用 --help 查看用法。")
}
if (!file.exists(opt$input)) {
  stop(sprintf("ERROR: 文件不存在: %s", opt$input))
}


# ==============================================================================
# §1  全局绘图配置 + 辅助函数
# ==============================================================================

style_config <- list(
  color_sig       = "#A91601FF",
  color_ns        = "#4666DDFF",
  base_font_size  = 14,
  show_n          = opt$show_n
)

# --- 同步保存 PDF + PNG -------------------------------------------------------
save_plot_pair <- function(plot_obj, file_pdf, file_png, width = 8, height = 8) {
  ggplot2::ggsave(filename = file_pdf, plot = plot_obj, width = width, height = height)
  ggplot2::ggsave(filename = file_png, plot = plot_obj, width = width, height = height, dpi = 300)
}

# --- ID 类型自动探针 ----------------------------------------------------------
#' 通过正则规则检测 gene_list 中的 ID 类型
#' @param gene_vec 字符向量，基因 ID 样本
#' @return "ENTREZID" | "ENSEMBL" | "SYMBOL"
detect_id_type <- function(gene_vec) {
  s <- head(na.omit(as.character(gene_vec)), 30)
  if (all(grepl("^[0-9]+$", s)))                   return("ENTREZID")
  if (all(grepl("^ENS[A-Z]*G[0-9]+(\\.\\d+)?$", s)))  return("ENSEMBL")
  return("SYMBOL")
}

# --- 列名自动猜测 -------------------------------------------------------------
#' 当用户未指定 gene_col 时，从数据框中尝试推断
guess_gene_col <- function(df) {
  cols <- colnames(df)
  # 优先匹配包含 gene / symbol / id 的列名（忽略大小写）
  hit <- grep("gene|symbol", cols, ignore.case = TRUE, value = TRUE)
  if (length(hit) > 0) return(hit[1])
  # fallback 到第一列
  message("  [Adapter] gene_col not specified, falling back to column 1: ", cols[1])
  return(cols[1])
}


# ==============================================================================
# §2  ID 转换模块
# ==============================================================================

#' 将任意 ID 类型转换为 ENTREZID，同时构建 gene_map（原 ID ↔ ENTREZID 双向映射）
#' @param gene_vec  字符向量，原始基因 ID
#' @param id_type   "SYMBOL" | "ENSEMBL" | "ENTREZID"
#' @param orgdb_obj OrgDb 对象
#' @return list(entrez_ids = 字符向量, gene_map = data.frame)；失败则返回 NULL
convert_to_entrez <- function(gene_vec, id_type, orgdb_obj) {
  if (id_type == "ENTREZID") {
    # 已经是 Entrez，直接构建 map
    gene_map <- data.frame(INPUT = gene_vec, ENTREZID = gene_vec, stringsAsFactors = FALSE)
    return(list(entrez_ids = unique(gene_vec), gene_map = gene_map))
  }

  gene_map <- tryCatch({
    clusterProfiler::bitr(
      geneID   = gene_vec,
      fromType = id_type,
      toType   = "ENTREZID",
      OrgDb    = orgdb_obj
    )
  }, error = function(e) {
    warning("  [ID Convert] bitr 转换失败: ", e$message)
    NULL
  })

  if (is.null(gene_map) || nrow(gene_map) == 0) {
    return(NULL)
  }

  # 统一列名方便下游使用
  colnames(gene_map)[1] <- "INPUT"
  list(entrez_ids = unique(gene_map$ENTREZID), gene_map = gene_map)
}


# ==============================================================================
# §3  核心 ORA 分析层
# ==============================================================================

#' 执行 ORA (Over-Representation Analysis) 富集分析并输出图表
#' @param entrez_ids  字符向量，目标基因 Entrez ID
#' @param gene_map    data.frame，INPUT ↔ ENTREZID 映射（来自 convert_to_entrez）
#' @param fc_map      命名数值向量（names 为原始 gene symbol），有则画山脊图；NULL 则跳过
#' @param bg_entrez   字符向量，背景基因 Entrez ID；NULL 使用全库默认背景
#' @param orgdb_obj   OrgDb 对象
#' @param organism    KEGG 物种代码
#' @param out_dir     输出目录
#' @param task_id     任务标识符（图标题）
#' @param cfg         style_config 列表
#' @return list，包含各项富集结果的统计指标
run_ora <- function(entrez_ids, gene_map, fc_map = NULL,
                    bg_entrez = NULL, orgdb_obj, organism,
                    out_dir, task_id, cfg,
                    kegg_p_cutoff = 0.05,
                    go_p_cutoff = 0.05, go_q_cutoff = 0.2) {

  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  has_fc <- !is.null(fc_map) && length(fc_map) > 0

  kegg_stats <- list(count = 0)
  go_stats   <- list(count = 0, bp = 0, cc = 0, mf = 0)
  kegg_plot_df <- NULL
  go_plot_df   <- NULL

  # --- 3.1  KEGG ORA --------------------------------------------------------
  try({
    kegg_res <- clusterProfiler::enrichKEGG(
      gene          = entrez_ids,
      organism      = organism,
      keyType       = "kegg",
      pAdjustMethod = "none",
      pvalueCutoff  = kegg_p_cutoff,
      qvalueCutoff  = 1,
      universe      = bg_entrez
    )

    if (!is.null(kegg_res) && nrow(kegg_res) > 0) {
      kegg_res  <- clusterProfiler::setReadable(kegg_res, OrgDb = orgdb_obj, keyType = "ENTREZID")
      kegg_df   <- as.data.frame(kegg_res)
      kegg_df$log_p <- -log10(kegg_df$p.adjust + 1e-300)

      saveRDS(kegg_res, file.path(out_dir, "01.kegg_results.rds"))
      write.csv(kegg_df, file.path(out_dir, "01.kegg_results.csv"), row.names = FALSE)

      # -- 散点图 --
      kegg_plot_df <- kegg_df %>%
        dplyr::arrange(p.adjust) %>%
        dplyr::slice_head(n = cfg$show_n)

      p_kegg <- ggplot2::ggplot(kegg_plot_df,
                   ggplot2::aes(x = FoldEnrichment,
                                y = forcats::fct_reorder(Description, FoldEnrichment))) +
        ggplot2::geom_point(ggplot2::aes(size = Count, color = log_p)) +
        ggplot2::scale_color_gradient(low = cfg$color_ns, high = cfg$color_sig) +
        ggplot2::labs(color = "-log10(p.adjust)", size = "Count",
                      x = "Fold Enrichment", y = NULL,
                      title = paste0("KEGG ORA: ", task_id)) +
        ggplot2::theme_classic(base_size = cfg$base_font_size)

      save_plot_pair(p_kegg,
                     file.path(out_dir, "03.kegg_scatter.pdf"),
                     file.path(out_dir, "03.kegg_scatter.png"),
                     width = 10, height = 8)

      # -- 山脊散点联合图（需要 fc_map）--
      if (has_fc) {
        try({
          df_pathway_kegg <- kegg_plot_df %>%
            dplyr::mutate(pathway = factor(Description, levels = rev(Description)))

          df_genes_kegg <- do.call(rbind, lapply(seq_len(nrow(df_pathway_kegg)), function(i) {
            genes   <- unlist(strsplit(df_pathway_kegg$geneID[i], "/"))
            fc_vals <- fc_map[genes]
            fc_vals <- fc_vals[!is.na(fc_vals)]
            if (length(fc_vals) == 0) fc_vals <- 0
            data.frame(
              pathway  = df_pathway_kegg$pathway[i],
              log_p    = df_pathway_kegg$log_p[i],
              gene_ratio = fc_vals
            )
          }))

          if (nrow(df_genes_kegg) > 0) {
            min_x   <- min(df_genes_kegg$gene_ratio, na.rm = TRUE)
            max_x   <- max(df_genes_kegg$gene_ratio, na.rm = TRUE)
            x_range <- max(max_x - min_x, 1)
            scatter_x <- min_x - x_range * 0.15
            left_lim  <- min_x - x_range * 0.35

            p_kegg_ridge <- ggplot2::ggplot() +
              ggridges::geom_density_ridges(
                data = df_genes_kegg,
                ggplot2::aes(x = gene_ratio, y = pathway, fill = log_p),
                scale = 1.6, color = "white", linewidth = 0.5) +
              ggplot2::geom_point(
                data = df_pathway_kegg,
                ggplot2::aes(x = scatter_x, y = pathway,
                             color = FoldEnrichment, size = Count)) +
              ggplot2::scale_size_continuous(range = c(2, 6), guide = "none") +
              ggplot2::scale_fill_distiller(palette = "RdYlBu",
                                            name = "-Log10(p.adjust)", direction = -1) +
              ggplot2::scale_color_gradientn(
                colors = c("#7f489f", "#ca5ba1", "#ed506a", "#e82335"),
                name = "FoldEnrichment") +
              ggplot2::scale_x_continuous(
                limits = c(left_lim, max_x + x_range * 0.1)) +
              ggplot2::labs(x = "log2(FoldChange)", y = NULL,
                            title = paste0("KEGG ORA: ", task_id)) +
              ggplot2::theme_classic(base_size = cfg$base_font_size) +
              ggplot2::theme(
                plot.title   = ggplot2::element_text(hjust = 0.5, face = "bold"),
                panel.border = ggplot2::element_rect(color = "black", fill = NA, linewidth = 1),
                axis.line    = ggplot2::element_blank(),
                panel.grid.major.x = ggplot2::element_line(color = "grey90", linewidth = 0.3)
              )

            save_plot_pair(p_kegg_ridge,
                           file.path(out_dir, "04.kegg_ridge.pdf"),
                           file.path(out_dir, "04.kegg_ridge.png"),
                           width = 10, height = 8)
          }
        }, silent = TRUE)
      } else {
        message("  [Plot] KEGG Ridge skipped (no metric / fc_map provided).")
      }

      kegg_stats$count <- nrow(kegg_df)
    }
  })

  # --- 3.2  GO ORA ----------------------------------------------------------
  try({
    go_res <- clusterProfiler::enrichGO(
      gene          = entrez_ids,
      OrgDb         = orgdb_obj,
      keyType       = "ENTREZID",
      ont           = "ALL",
      pAdjustMethod = "BH",
      pvalueCutoff  = go_p_cutoff,
      qvalueCutoff  = go_q_cutoff,
      readable      = TRUE,
      universe      = bg_entrez
    )

    if (!is.null(go_res) && nrow(go_res) > 0) {
      go_df       <- as.data.frame(go_res)
      go_df$log_p <- -log10(go_df$p.adjust + 1e-300)

      saveRDS(go_res, file.path(out_dir, "02.go_results.rds"))
      write.csv(go_df, file.path(out_dir, "02.go_results.csv"), row.names = FALSE)

      # -- 散点图 --
      go_plot_df <- go_df %>%
        dplyr::filter(ONTOLOGY %in% c("BP", "CC", "MF")) %>%
        dplyr::group_by(ONTOLOGY) %>%
        dplyr::arrange(p.adjust) %>%
        dplyr::slice_head(n = cfg$show_n) %>%
        dplyr::ungroup()

      p_go <- ggplot2::ggplot(
          go_plot_df,
          ggplot2::aes(x = FoldEnrichment,
                       y = forcats::fct_reorder(Description, FoldEnrichment))) +
        ggplot2::geom_point(ggplot2::aes(size = Count, color = log_p, shape = ONTOLOGY)) +
        ggplot2::scale_color_gradient(low = cfg$color_ns, high = cfg$color_sig) +
        ggplot2::labs(color = "-log10(p.adjust)", size = "Count",
                      x = "Fold Enrichment", y = NULL,
                      title = paste0("GO ORA: ", task_id)) +
        ggplot2::theme_classic(base_size = cfg$base_font_size) +
        ggplot2::facet_grid(ONTOLOGY ~ ., scales = "free_y", space = "free")

      save_plot_pair(p_go,
                     file.path(out_dir, "05.go_scatter.pdf"),
                     file.path(out_dir, "05.go_scatter.png"),
                     width = 12, height = 10)

      # -- GO 山脊图（需要 fc_map）--
      if (has_fc) {
        try({
          go_plot_df_comb <- go_plot_df %>%
            dplyr::mutate(ONTOLOGY = factor(ONTOLOGY, levels = c("MF", "CC", "BP"))) %>%
            dplyr::arrange(ONTOLOGY, FoldEnrichment)
          go_plot_df_comb$y_num <- seq_len(nrow(go_plot_df_comb))

          df_genes_go <- do.call(rbind, lapply(seq_len(nrow(go_plot_df_comb)), function(i) {
            genes   <- unlist(strsplit(go_plot_df_comb$geneID[i], "/"))
            fc_vals <- fc_map[genes]
            fc_vals <- fc_vals[!is.na(fc_vals)]
            if (length(fc_vals) == 0) fc_vals <- 0
            data.frame(
              ontology   = go_plot_df_comb$ONTOLOGY[i],
              log_p      = go_plot_df_comb$log_p[i],
              gene_ratio = fc_vals,
              y_num      = go_plot_df_comb$y_num[i]
            )
          }))

          if (nrow(df_genes_go) > 0) {
            min_x_go  <- min(df_genes_go$gene_ratio, na.rm = TRUE)
            max_x_go  <- max(df_genes_go$gene_ratio, na.rm = TRUE)
            x_range_go <- max(max_x_go - min_x_go, 1)
            scatter_x_go <- min_x_go - x_range_go * 0.15
            text_x_go    <- min_x_go - x_range_go * 0.35
            line_x_go    <- min_x_go - x_range_go * 0.25
            left_lim_go  <- min_x_go - x_range_go * 0.45

            group_pos <- go_plot_df_comb %>%
              dplyr::group_by(ONTOLOGY) %>%
              dplyr::summarize(y_min  = min(y_num),
                               y_max  = max(y_num),
                               y_mean = mean(y_num),
                               .groups = "drop")
            y_cutoffs <- group_pos$y_max[-length(group_pos$y_max)] + 0.5

            p_go_ridge <- ggplot2::ggplot() +
              ggplot2::geom_hline(yintercept = y_cutoffs,
                                  linetype = "dashed", color = "grey80", linewidth = 0.5) +
              ggplot2::geom_segment(
                data = group_pos,
                ggplot2::aes(x = line_x_go, xend = line_x_go,
                             y = y_min - 0.4, yend = y_max + 0.4),
                linewidth = 1.5, color = "grey40", lineend = "round") +
              ggplot2::geom_text(
                data = group_pos,
                ggplot2::aes(x = text_x_go, y = y_mean, label = ONTOLOGY),
                angle = 90, fontface = "bold", size = 4.5, color = "grey20") +
              ggridges::geom_density_ridges(
                data = df_genes_go,
                ggplot2::aes(x = gene_ratio, y = y_num,
                             group = y_num, fill = log_p),
                scale = 1.5, color = "white", linewidth = 0.5) +
              ggplot2::geom_point(
                data = go_plot_df_comb,
                ggplot2::aes(x = scatter_x_go, y = y_num,
                             color = FoldEnrichment, size = Count)) +
              ggplot2::scale_size_continuous(range = c(2, 6), guide = "none") +
              ggplot2::scale_fill_distiller(palette = "RdYlBu",
                                            name = "-Log10(p.adjust)", direction = -1) +
              ggplot2::scale_color_gradientn(
                colors = c("#7f489f", "#ca5ba1", "#ed506a", "#e82335"),
                name = "FoldEnrichment") +
              ggplot2::scale_x_continuous(
                limits = c(left_lim_go, max_x_go + x_range_go * 0.1)) +
              ggplot2::scale_y_continuous(
                breaks = go_plot_df_comb$y_num,
                labels = go_plot_df_comb$Description) +
              ggplot2::labs(x = "log2(FoldChange)", y = NULL,
                            title = paste0("GO ORA: ", task_id)) +
              ggplot2::theme_classic(base_size = cfg$base_font_size) +
              ggplot2::theme(
                plot.title   = ggplot2::element_text(hjust = 0.5, face = "bold"),
                panel.border = ggplot2::element_rect(color = "black", fill = NA, linewidth = 1),
                axis.line    = ggplot2::element_blank(),
                panel.grid.major.x = ggplot2::element_line(color = "grey90", linewidth = 0.3)
              )

            save_plot_pair(p_go_ridge,
                           file.path(out_dir, "06.go_ridge.pdf"),
                           file.path(out_dir, "06.go_ridge.png"),
                           width = 11, height = 9)
          }
        }, silent = TRUE)
      } else {
        message("  [Plot] GO Ridge skipped (no metric / fc_map provided).")
      }

      go_stats$count <- nrow(go_df)
      go_stats$bp    <- sum(go_df$ONTOLOGY == "BP")
      go_stats$cc    <- sum(go_df$ONTOLOGY == "CC")
      go_stats$mf    <- sum(go_df$ONTOLOGY == "MF")
    }
  })

  # --- 3.3  多模态 Barcode 综合图 -------------------------------------------
  try({
    if (!is.null(kegg_plot_df) && !is.null(go_plot_df) &&
        nrow(kegg_plot_df) > 0 && nrow(go_plot_df) > 0) {

      kegg_sub <- kegg_plot_df %>%
        dplyr::mutate(category = "KEGG") %>%
        dplyr::select(category, pathway = Description,
                      log_p, count_num = Count, gene_id = geneID)

      go_sub <- go_plot_df %>%
        dplyr::rename(category = ONTOLOGY, pathway = Description,
                      count_num = Count, gene_id = geneID) %>%
        dplyr::select(category, pathway, log_p, count_num, gene_id)

      cat_levels <- c("BP", "CC", "MF", "KEGG")
      barcode_df <- dplyr::bind_rows(kegg_sub, go_sub) %>%
        dplyr::filter(category %in% cat_levels) %>%
        dplyr::group_by(category) %>%
        dplyr::arrange(dplyr::desc(log_p)) %>%
        dplyr::slice_head(n = 5) %>%
        dplyr::ungroup() %>%
        dplyr::mutate(category = factor(category, levels = cat_levels)) %>%
        dplyr::arrange(dplyr::desc(category), log_p)

      if (nrow(barcode_df) > 4) {
        barcode_df$genes_clean <- sapply(barcode_df$gene_id, function(g) {
          gl <- unlist(strsplit(g, "/"))
          if (length(gl) > 10) gl <- gl[1:10]
          paste(gl, collapse = "/")
        })

        n_total <- nrow(barcode_df)
        barcode_df$y_num <- rev(seq_len(n_total))
        barcode_df$y_cat <- factor(barcode_df$y_num,
                                   levels = sort(unique(barcode_df$y_num)))

        color_mapping <- c("BP" = "#3e528b", "CC" = "#089d7b",
                           "MF" = "#60bed5", "KEGG" = "#df4428")

        group_pos_bar <- barcode_df %>%
          dplyr::group_by(category) %>%
          dplyr::summarize(y_min  = min(y_num) - 0.45,
                           y_max  = max(y_num) + 0.45,
                           y_mean = mean(y_num),
                           .groups = "drop")

        max_logp_bar <- max(barcode_df$log_p, na.rm = TRUE)
        x_breaks     <- seq(0, ceiling(max_logp_bar / 5) * 5, by = 5)

        p_bar <- ggplot2::ggplot() +
          ggplot2::geom_segment(
            data = barcode_df,
            ggplot2::aes(x = 0, xend = log_p,
                         y  = as.numeric(y_cat), yend = as.numeric(y_cat),
                         color = category),
            linewidth = 6.5, lineend = "round", alpha = 0.85) +
          ggplot2::geom_rect(
            data = group_pos_bar,
            ggplot2::aes(ymin = y_min, ymax = y_max,
                         xmin = -6.5, xmax = -4.5, fill = category)) +
          ggplot2::geom_text(
            data = group_pos_bar,
            ggplot2::aes(y = y_mean, x = -5.5, label = category),
            angle = 90, fontface = "bold", size = 5.0) +
          ggplot2::geom_point(
            data = barcode_df,
            ggplot2::aes(y = as.numeric(y_cat), x = -2.5,
                         fill = category, size = count_num),
            color = "black", shape = 21, stroke = 1.0) +
          ggplot2::geom_text(
            data = barcode_df,
            ggplot2::aes(y = as.numeric(y_cat), x = -2.5, label = count_num),
            color = "black", size = 3.8) +
          ggplot2::geom_text(
            data = barcode_df,
            ggplot2::aes(y = as.numeric(y_cat), x = 0, label = pathway),
            hjust = 0, color = "black", size = 4.2) +
          ggplot2::geom_text(
            data = barcode_df,
            ggplot2::aes(y = as.numeric(y_cat) - 0.40, x = 0,
                         label = genes_clean, color = category),
            hjust = 0, fontface = "italic", size = 3.0) +
          ggplot2::geom_segment(
            ggplot2::aes(x = 0, xend = max(x_breaks),
                         y = -0.3, yend = -0.3),
            color = "black", linewidth = 1.5) +
          ggplot2::geom_segment(
            data = data.frame(x = x_breaks),
            ggplot2::aes(x = x, xend = x, y = -0.3, yend = -0.55),
            color = "black", linewidth = 1.5) +
          ggplot2::geom_text(
            data = data.frame(x = x_breaks),
            ggplot2::aes(x = x, y = -0.9, label = x),
            color = "black", size = 5.5, fontface = "bold") +
          ggplot2::geom_text(
            ggplot2::aes(x = max(x_breaks) / 2, y = -1.8,
                         label = "-log10(p.adjust)"),
            color = "black", size = 6, fontface = "bold") +
          ggplot2::scale_fill_manual(values = color_mapping) +
          ggplot2::scale_color_manual(values = color_mapping) +
          ggplot2::scale_size_continuous(range = c(6, 11), name = NULL) +
          ggplot2::scale_x_continuous(
            limits = c(-7, max(max(x_breaks) * 1.5, 20))) +
          ggplot2::scale_y_continuous(
            limits = c(-2.2, n_total + 0.5), expand = c(0, 0)) +
          ggplot2::theme_classic(base_size = 14) +
          ggplot2::theme(
            axis.title   = ggplot2::element_blank(),
            axis.text    = ggplot2::element_blank(),
            axis.ticks   = ggplot2::element_blank(),
            axis.line    = ggplot2::element_blank(),
            legend.position = "right",
            legend.title = ggplot2::element_blank()
          )

        save_plot_pair(p_bar,
                       file.path(out_dir, "07.enrichment_barcode.pdf"),
                       file.path(out_dir, "07.enrichment_barcode.png"),
                       width = 11, height = max(6, n_total * 0.45))
      }
    }
  }, silent = TRUE)

  list(
    mode         = "ORA",
    kegg_terms   = kegg_stats$count,
    go_terms     = go_stats$count,
    count_bp     = go_stats$bp,
    count_cc     = go_stats$cc,
    count_mf     = go_stats$mf
  )
}


# ==============================================================================
# §4  核心 GSEA 分析层
# ==============================================================================

#' 执行 GSEA (Gene Set Enrichment Analysis) 富集分析并输出图表
#' @param ranked_vec  按 metric（如 log2FC）降序排列的命名数值向量（names 为 ENTREZID）
#' @param orgdb_obj   OrgDb 对象
#' @param organism    KEGG 物种代码
#' @param out_dir     输出目录
#' @param task_id     任务标识符
#' @param cfg         style_config 列表
#' @param kegg_p_cutoff KEGG p 值阈值
#' @param go_p_cutoff   GO p 值阈值
#' @return list，包含各项富集结果统计
run_gsea <- function(ranked_vec, orgdb_obj, organism,
                     out_dir, task_id, cfg,
                     kegg_p_cutoff = 0.05, go_p_cutoff = 0.05) {

  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  kegg_stats <- list(count = 0)
  go_stats   <- list(count = 0)

  # --- 4.1  KEGG GSEA --------------------------------------------------------
  try({
    gsea_kegg <- clusterProfiler::gseKEGG(
      geneList      = ranked_vec,
      organism      = organism,
      keyType       = "kegg",
      minGSSize     = 10,
      maxGSSize     = 500,
      pvalueCutoff  = kegg_p_cutoff,
      pAdjustMethod = "BH",
      verbose       = FALSE
    )

    if (!is.null(gsea_kegg) && nrow(gsea_kegg) > 0) {
      gsea_kegg <- clusterProfiler::setReadable(gsea_kegg, OrgDb = orgdb_obj,
                                                 keyType = "ENTREZID")
      gsea_kegg_df <- as.data.frame(gsea_kegg)

      saveRDS(gsea_kegg, file.path(out_dir, "01.gsea_kegg_results.rds"))
      write.csv(gsea_kegg_df, file.path(out_dir, "01.gsea_kegg_results.csv"),
                row.names = FALSE)

      # -- Dotplot --
      p_dot_kegg <- enrichplot::dotplot(
        gsea_kegg, showCategory = cfg$show_n, split = ".sign") +
        ggplot2::facet_grid(. ~ .sign) +
        ggplot2::labs(title = paste0("KEGG GSEA: ", task_id)) +
        ggplot2::theme_classic(base_size = cfg$base_font_size)

      save_plot_pair(p_dot_kegg,
                     file.path(out_dir, "03.gsea_kegg_dotplot.pdf"),
                     file.path(out_dir, "03.gsea_kegg_dotplot.png"),
                     width = 12, height = 9)

      # -- Ridge Plot --
      try({
        p_ridge_kegg <- enrichplot::ridgeplot(gsea_kegg,
                                               showCategory = min(cfg$show_n, nrow(gsea_kegg))) +
          ggplot2::labs(x = "Enrichment Distribution",
                        title = paste0("KEGG GSEA Ridge: ", task_id)) +
          ggplot2::theme_classic(base_size = cfg$base_font_size)

        save_plot_pair(p_ridge_kegg,
                       file.path(out_dir, "04.gsea_kegg_ridge.pdf"),
                       file.path(out_dir, "04.gsea_kegg_ridge.png"),
                       width = 10, height = 8)
      }, silent = TRUE)

      # -- GSEA Enrichment Score Plot（Top 5 通路）--
      try({
        top_n <- min(5, nrow(gsea_kegg))
        p_gsea_curve <- enrichplot::gseaplot2(
          gsea_kegg,
          geneSetID = seq_len(top_n),
          title     = paste0("KEGG GSEA Score: ", task_id),
          pvalue_by = "p.adjust"
        )
        ggplot2::ggsave(
          filename = file.path(out_dir, "05.gsea_kegg_score.pdf"),
          plot = p_gsea_curve, width = 10, height = 3 + top_n * 1.5)
        ggplot2::ggsave(
          filename = file.path(out_dir, "05.gsea_kegg_score.png"),
          plot = p_gsea_curve, width = 10, height = 3 + top_n * 1.5, dpi = 300)
      }, silent = TRUE)

      kegg_stats$count <- nrow(gsea_kegg_df)
    }
  })

  # --- 4.2  GO GSEA ----------------------------------------------------------
  try({
    gsea_go <- clusterProfiler::gseGO(
      geneList      = ranked_vec,
      OrgDb         = orgdb_obj,
      keyType       = "ENTREZID",
      ont           = "ALL",
      minGSSize     = 10,
      maxGSSize     = 500,
      pvalueCutoff  = go_p_cutoff,
      pAdjustMethod = "BH",
      verbose       = FALSE
    )

    if (!is.null(gsea_go) && nrow(gsea_go) > 0) {
      gsea_go    <- clusterProfiler::setReadable(gsea_go, OrgDb = orgdb_obj,
                                                  keyType = "ENTREZID")
      gsea_go_df <- as.data.frame(gsea_go)

      saveRDS(gsea_go, file.path(out_dir, "02.gsea_go_results.rds"))
      write.csv(gsea_go_df, file.path(out_dir, "02.gsea_go_results.csv"),
                row.names = FALSE)

      # -- Dotplot --
      p_dot_go <- enrichplot::dotplot(
        gsea_go, showCategory = cfg$show_n, split = ".sign") +
        ggplot2::facet_grid(. ~ .sign) +
        ggplot2::labs(title = paste0("GO GSEA: ", task_id)) +
        ggplot2::theme_classic(base_size = cfg$base_font_size)

      save_plot_pair(p_dot_go,
                     file.path(out_dir, "06.gsea_go_dotplot.pdf"),
                     file.path(out_dir, "06.gsea_go_dotplot.png"),
                     width = 12, height = 9)

      # -- Ridge Plot --
      try({
        p_ridge_go <- enrichplot::ridgeplot(gsea_go,
                                             showCategory = min(cfg$show_n, nrow(gsea_go))) +
          ggplot2::labs(x = "Enrichment Distribution",
                        title = paste0("GO GSEA Ridge: ", task_id)) +
          ggplot2::theme_classic(base_size = cfg$base_font_size)

        save_plot_pair(p_ridge_go,
                       file.path(out_dir, "07.gsea_go_ridge.pdf"),
                       file.path(out_dir, "07.gsea_go_ridge.png"),
                       width = 10, height = 8)
      }, silent = TRUE)

      # -- GSEA Enrichment Score Plot（Top 5 通路）--
      try({
        top_n <- min(5, nrow(gsea_go))
        p_gsea_go_curve <- enrichplot::gseaplot2(
          gsea_go,
          geneSetID = seq_len(top_n),
          title     = paste0("GO GSEA Score: ", task_id),
          pvalue_by = "p.adjust"
        )
        ggplot2::ggsave(
          filename = file.path(out_dir, "08.gsea_go_score.pdf"),
          plot = p_gsea_go_curve, width = 10, height = 3 + top_n * 1.5)
        ggplot2::ggsave(
          filename = file.path(out_dir, "08.gsea_go_score.png"),
          plot = p_gsea_go_curve, width = 10, height = 3 + top_n * 1.5, dpi = 300)
      }, silent = TRUE)

      go_stats$count <- nrow(gsea_go_df)
    }
  })

  list(
    mode       = "GSEA",
    kegg_terms = kegg_stats$count,
    go_terms   = go_stats$count
  )
}


# ==============================================================================
# §5  CLI Adapter — 主程序入口
# ==============================================================================

message("\n", strrep("=", 70))
message("  智能基因富集分析 Skill  |  ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"))
message(strrep("=", 70))
message("  Input  : ", opt$input)
message("  OutDir : ", opt$outdir)
message("  Organism: ", opt$organism, "  |  OrgDb: ", opt$orgdb)

# --- 5.1  动态加载 OrgDb 包 --------------------------------------------------
if (!requireNamespace(opt$orgdb, quietly = TRUE)) {
  stop(sprintf("ERROR: OrgDb 包 '%s' 未安装，请先运行 BiocManager::install(\"%s\")",
               opt$orgdb, opt$orgdb))
}
suppressPackageStartupMessages(library(opt$orgdb, character.only = TRUE))
orgdb_obj <- get(opt$orgdb)

# --- 5.2  读取输入文件 -------------------------------------------------------
message("\n>>> [Adapter] Reading input file...")
df_raw <- tryCatch(
  data.table::fread(opt$input, header = "auto", stringsAsFactors = FALSE,
                    data.table = FALSE),
  error = function(e) stop("文件读取失败: ", e$message)
)
message(sprintf("  Loaded: %d rows x %d cols", nrow(df_raw), ncol(df_raw)))
if (ncol(df_raw) <= 5) {
  message("  Columns: ", paste(colnames(df_raw), collapse = ", "))
} else {
  message("  Columns (first 5): ", paste(head(colnames(df_raw), 5), collapse = ", "), " ...")
}

# --- 5.3  确定基因列 ---------------------------------------------------------
gene_col <- if (!is.null(opt$gene_col)) opt$gene_col else guess_gene_col(df_raw)

if (!gene_col %in% colnames(df_raw)) {
  stop(sprintf("ERROR: 指定的 gene_col='%s' 不存在于输入文件中！\n  可用列名: %s",
               gene_col, paste(colnames(df_raw), collapse = ", ")))
}
message("  Gene column  : ", gene_col)

# --- 5.4  分析模式路由 -------------------------------------------------------
# 确定模式
mode <- tolower(opt$mode)
if (mode == "auto") {
  if (!is.null(opt$metric_col) && opt$metric_col != "") {
    mode <- "gsea"
  } else {
    mode <- "ora"
  }
}

if (mode == "gsea") {
  if (is.null(opt$metric_col) || opt$metric_col == "") {
    stop("ERROR: GSEA 模式必须提供 --metric_col 指标列。")
  }
  message(">>> Mode: GSEA")
} else {
  message(">>> Mode: ORA")
}

# --- 执行路径选择 ---
if (mode == "gsea") {
  # -- GSEA 路径：需要 metric_col 且存在于文件中 --
  if (!opt$metric_col %in% colnames(df_raw)) {
    stop(sprintf("ERROR: --metric_col='%s' 不存在于输入文件中！", opt$metric_col))
  }
  message("  Metric column: ", opt$metric_col)

  # 过滤 NA 并构建 Ranked Named Vector
  df_gsea <- df_raw[!is.na(df_raw[[opt$metric_col]]) &
                      !is.na(df_raw[[gene_col]]) &
                      df_raw[[gene_col]] != "", ]
  metric_vec <- as.numeric(df_gsea[[opt$metric_col]])
  names(metric_vec) <- as.character(df_gsea[[gene_col]])

  # 去重（保留 metric 绝对值最大的那条）
  dup_genes  <- names(metric_vec)[duplicated(names(metric_vec))]
  if (length(dup_genes) > 0) {
    message(sprintf("  [Adapter] Removing %d duplicate gene entries (keeping max |metric|)",
                    length(dup_genes)))
    df_dedup    <- df_gsea[order(abs(metric_vec), decreasing = TRUE), ]
    df_dedup    <- df_dedup[!duplicated(df_dedup[[gene_col]]), ]
    metric_vec  <- as.numeric(df_dedup[[opt$metric_col]])
    names(metric_vec) <- as.character(df_dedup[[gene_col]])
  }

  metric_vec <- sort(metric_vec, decreasing = TRUE)
  message(sprintf("  Genes in ranked list: %d", length(metric_vec)))

  # ID 探测与转换
  id_type <- detect_id_type(names(metric_vec))
  message("  Detected ID type: ", id_type)

  convert_result <- convert_to_entrez(names(metric_vec), id_type, orgdb_obj)
  if (is.null(convert_result)) {
    stop("ERROR: ID 转换失败，没有有效的 Entrez ID 可供分析。")
  }

  # 将 ranked_vec 的 names 换成 ENTREZID
  gene_map       <- convert_result$gene_map
  metric_map     <- setNames(metric_vec[gene_map$INPUT], gene_map$ENTREZID)
  # 去重（同一个 ENTREZID 可能对应多个 Symbol）
  metric_map     <- metric_map[!duplicated(names(metric_map))]
  ranked_entrez  <- sort(metric_map, decreasing = TRUE)

  message(sprintf("  Ranked ENTREZID count: %d (mapping rate: %.1f%%)",
                  length(ranked_entrez),
                  length(ranked_entrez) / length(metric_vec) * 100))

  # 执行 GSEA
  final_stats <- run_gsea(
    ranked_vec    = ranked_entrez,
    orgdb_obj     = orgdb_obj,
    organism      = opt$organism,
    out_dir       = opt$outdir,
    task_id       = opt$task_id,
    cfg           = style_config,
    kegg_p_cutoff = opt$kegg_p_cutoff,
    go_p_cutoff   = opt$go_p_cutoff
  )

} else {
  # -- ORA 路径 --
  # (a) 如果提供了 pval_col，用它过滤显著基因
  if (!is.null(opt$pval_col)) {
    if (!opt$pval_col %in% colnames(df_raw)) {
      warning(sprintf("pval_col='%s' 不存在，跳过 p 值过滤，使用全部基因。", opt$pval_col))
      df_sig <- df_raw
    } else {
      df_sig <- df_raw[!is.na(df_raw[[opt$pval_col]]) &
                         as.numeric(df_raw[[opt$pval_col]]) < opt$pval_cutoff, ]
      message(sprintf("  Filtered by %s < %.3f: %d → %d genes",
                      opt$pval_col, opt$pval_cutoff, nrow(df_raw), nrow(df_sig)))
    }
  } else {
    df_sig <- df_raw
    message("  No pval_col specified, treating all rows as significant genes.")
  }

  df_sig  <- df_sig[!is.na(df_sig[[gene_col]]) & df_sig[[gene_col]] != "", ]
  gene_vec <- unique(as.character(df_sig[[gene_col]]))
  message(sprintf("  Target genes: %d", length(gene_vec)))

  if (length(gene_vec) < 5) {
    stop("ERROR: 显著基因数量不足 5 个，无法进行富集分析。")
  }

  # (b) 构建 fc_map（如果 metric_col 指定，即使是 ORA 也可以画山脊图）
  fc_map <- NULL
  if (!is.null(opt$metric_col) && opt$metric_col %in% colnames(df_sig)) {
    fc_vals <- as.numeric(df_sig[[opt$metric_col]])
    names(fc_vals) <- as.character(df_sig[[gene_col]])
    fc_map <- fc_vals[!is.na(fc_vals)]
  }

  # (c) 处理背景基因
  bg_entrez <- NULL
  if (!is.null(opt$background) && file.exists(opt$background)) {
    bg_raw    <- data.table::fread(opt$background, header = FALSE,
                                   data.table = FALSE)[[1]]
    bg_raw    <- unique(as.character(bg_raw[!is.na(bg_raw)]))
    message(sprintf("  Background genes loaded: %d", length(bg_raw)))
    bg_type   <- detect_id_type(bg_raw)
    bg_conv   <- convert_to_entrez(bg_raw, bg_type, orgdb_obj)
    if (!is.null(bg_conv)) bg_entrez <- bg_conv$entrez_ids
  }

  # (d) ID 探测与转换
  id_type <- detect_id_type(gene_vec)
  message("  Detected ID type: ", id_type)

  convert_result <- convert_to_entrez(gene_vec, id_type, orgdb_obj)
  if (is.null(convert_result)) {
    stop("ERROR: ID 转换失败，没有有效的 Entrez ID 可供分析。")
  }

  entrez_ids <- convert_result$entrez_ids
  gene_map   <- convert_result$gene_map
  message(sprintf("  Valid ENTREZID: %d (mapping rate: %.1f%%)",
                  length(entrez_ids), length(entrez_ids) / length(gene_vec) * 100))

  # 执行 ORA
  final_stats <- run_ora(
    entrez_ids    = entrez_ids,
    gene_map      = gene_map,
    fc_map        = fc_map,
    bg_entrez     = bg_entrez,
    orgdb_obj     = orgdb_obj,
    organism      = opt$organism,
    out_dir       = opt$outdir,
    task_id       = opt$task_id,
    cfg           = style_config,
    kegg_p_cutoff = opt$kegg_p_cutoff,
    go_p_cutoff   = opt$go_p_cutoff,
    go_q_cutoff   = opt$go_q_cutoff
  )
}


# ==============================================================================
# §6  汇总输出
# ==============================================================================

summary_info <- c(
  list(
    task_id    = opt$task_id,
    input_file = opt$input,
    organism   = opt$organism,
    timestamp  = format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  ),
  final_stats
)

summary_df <- as.data.frame(summary_info, stringsAsFactors = FALSE)
write.csv(summary_df,
          file.path(opt$outdir, "00.enrichment_summary.csv"),
          row.names = FALSE)

message("\n", strrep("=", 70))
message("  Analysis Completed!")
message("  Mode          : ", final_stats$mode)
message("  KEGG terms    : ", final_stats$kegg_terms)
message("  GO terms      : ", final_stats$go_terms)
message("  Results saved : ", opt$outdir)
message(strrep("=", 70), "\n")
