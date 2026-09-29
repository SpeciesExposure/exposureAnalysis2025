# Fig3_IUCN.r
# Figure 3: IUCN Red List status of the species exposed in 2025, and how their
# status differs from that of all study species in the same taxonomic group.
#
# Inputs:  singleMaxExposure_v1.qs, spAttributes_v8.qs
# Outputs: figure in plotDir

suppressPackageStartupMessages({
  library(qs2)
  library(dplyr)
  library(ggplot2)
  library(patchwork)
})

if (!exists("PROJECT_PATHS")) source(file.path("config", "paths.R"))
source(file.path(PROJECT_PATHS$src_r, "1_Setup.r"))

figure_version <- if (exists("RUN_VERSION") && is.character(RUN_VERSION) && length(RUN_VERSION) == 1 && nzchar(trimws(RUN_VERSION))) {
  toupper(trimws(RUN_VERSION))
} else {
  "V8"
}

intDir <- if (exists("intDir") && is.character(intDir) && length(intDir) == 1 && nzchar(trimws(intDir))) {
  path.expand(intDir)
} else {
  path.expand(PROJECT_PATHS$int_v8)
}

run_tag <- if (exists("RUN_NAME") && is.character(RUN_NAME) && length(RUN_NAME) == 1 && nzchar(trimws(RUN_NAME))) {
  trimws(RUN_NAME)
} else {
  figure_version
}

outDir <- if (exists("plotDir") && is.character(plotDir) && length(plotDir) == 1 && nzchar(trimws(plotDir))) {
  path.expand(plotDir)
} else {
  file.path(PROJECT_PATHS$output_base, tolower(run_tag), "figures", "main")
}
dir.create(outDir, recursive = TRUE, showWarnings = FALSE)

# ---- Locate files ----
sp_meta_file <- KEY_FILES$sp_attributes_v8
if (!file.exists(sp_meta_file)) stop("Cannot find spAttributes file.")

single_max_candidates <- unique(path.expand(c(
  file.path(intDir, "singleMaxExposure_v1.qs"),
  file.path(PROJECT_PATHS$int_v8, "singleMaxExposure_v1.qs"),
  file.path(PROJECT_PATHS$output_base, tolower(run_tag), paste0("Int_", figure_version), "singleMaxExposure_v1.qs")
)))
single_max_file <- single_max_candidates[file.exists(single_max_candidates)][1]
if (is.na(single_max_file)) stop("Cannot find singleMaxExposure_v1.qs.")

cat("Loading singleMaxExposure...\n")
ex <- tryCatch(qs_read(single_max_file), error = function(e) {
  cat("qs2 read failed, trying qs::qread...\n")
  qs::qread(single_max_file)
})

cat("Loading spAttributes...\n")
sp_meta <- qs_read(sp_meta_file) %>%
  distinct(spName, .keep_all = TRUE) %>%
  dplyr::select(spName, group, redlistCategory)

# ---- IUCN helpers ----
THREATENED_CATS <- c("Critically_Endangered", "Endangered", "Vulnerable")
LC_CATS <- c("Least_Concern", "Lower_Risk/least concern",
             "Lower_Risk/conservation dependent", "Lower_Risk/near threatened")
IUCN_FOCUS <- c("Threatened", "Near Threatened", "Least Concern", "Data Deficient")
iucn_order  <- c("Least Concern", "Near Threatened", "Data Deficient", "Threatened")

collapse_iucn <- function(x) {
  dplyr::case_when(
    x %in% THREATENED_CATS ~ "Threatened",
    x == "Near_Threatened"  ~ "Near Threatened",
    x %in% LC_CATS          ~ "Least Concern",
    x == "Data_Deficient"   ~ "Data Deficient",
    is.na(x)               ~ "Unknown",
    TRUE                   ~ "Other/Unknown"
  )
}

# ---- 2025 exposed species ----
ex_2025_sp <- ex %>%
  filter(year == 2025) %>%
  distinct(spName, group, redlistCategory) %>%
  mutate(iucn_simple = collapse_iucn(redlistCategory),
         group = if_else(is.na(group) | group == "", "Unknown", group))

# ---- Panel A: stacked bar by group ----
iucn_by_group <- ex_2025_sp %>%
  filter(!is.na(group), group != "Unknown", iucn_simple %in% IUCN_FOCUS) %>%
  count(group, iucn_simple)

group_order <- ex_2025_sp %>%
  filter(group != "Unknown") %>%
  count(group, sort = TRUE) %>%
  pull(group)

cm_cols1_safe <- if (exists("cm.cols1") && is.function(cm.cols1)) {
  cm.cols1
} else {
  function(x, bias = 1) {
    grDevices::colorRampPalette(
      c("grey90", "steelblue4", "steelblue1", "gold", "red1", "red4"),
      bias = bias
    )(x)
  }
}
iucn_cols <- stats::setNames(cm_cols1_safe(length(iucn_order), bias = 1), iucn_order)

p_bar <- ggplot(
  iucn_by_group %>% mutate(group = factor(group, levels = group_order)),
  aes(x = group, y = n, fill = factor(iucn_simple, levels = iucn_order))
) +
  geom_col(width = 0.7, colour = "white", linewidth = 0.3) +
  scale_y_continuous(name = "Number of species exposed in 2025") +
  scale_fill_manual(values = iucn_cols, name = "IUCN status") +
  labs(title = "A", x = NULL) +
  theme_bw(base_size = 14) +
  theme(legend.position = "right",
        legend.key.size  = unit(0.4, "cm"),
        legend.text      = element_text(size = 11),
        legend.title     = element_text(size = 11),
        axis.text.x      = element_text(angle = 30, hjust = 1))

# ---- Panel B: lollipop enrichment vs background null ----
bg_sp <- sp_meta %>%
  mutate(iucn_broad = collapse_iucn(redlistCategory),
         group = if_else(is.na(group) | group == "", "Unknown", group))

obs_counts <- ex_2025_sp %>%
  filter(group != "Unknown") %>%
  rename(iucn_broad = iucn_simple) %>%
  count(group, iucn_broad, name = "n_exposed") %>%
  group_by(group) %>%
  mutate(obs_pct = 100 * n_exposed / sum(n_exposed)) %>%
  ungroup()

null_counts <- bg_sp %>%
  filter(group != "Unknown") %>%
  count(group, iucn_broad, name = "n_all") %>%
  group_by(group) %>%
  mutate(null_pct = 100 * n_all / sum(n_all)) %>%
  ungroup()

enrich_df <- obs_counts %>%
  full_join(null_counts, by = c("group", "iucn_broad")) %>%
  mutate(
    n_exposed = replace_na(n_exposed, 0L),
    obs_pct   = replace_na(obs_pct, 0),
    null_pct  = replace_na(null_pct, 0),
    diff_ppt  = obs_pct - null_pct
  ) %>%
  filter(iucn_broad %in% IUCN_FOCUS, group %in% group_order)

iucn_plot_order <- c("Least Concern", "Near Threatened", "Threatened", "Data Deficient")
all_plot <- enrich_df %>%
  mutate(iucn_broad = factor(iucn_broad, levels = iucn_plot_order)) %>%
  arrange(iucn_broad, diff_ppt) %>%
  mutate(
    label_display = tools::toTitleCase(group),
    label = paste0(as.integer(iucn_broad), "__", label_display),
    label = factor(label, levels = label)
  )

y_labels <- setNames(all_plot$label_display, all_plot$label)

# compute midpoint y-position (integer) for each IUCN group to place group label
group_annot <- all_plot %>%
  mutate(y_int = as.integer(label)) %>%
  group_by(iucn_broad) %>%
  summarise(y_mid = mean(y_int), y_min = min(y_int), y_max = max(y_int), .groups = "drop")

n_rows <- nlevels(all_plot$label)

# narrow strip plot with just IUCN group labels, aligned to lollipop y-axis
p_iucn_strip <- ggplot(group_annot, aes(x = 0, y = y_mid, label = iucn_broad)) +
  geom_segment(aes(x = 0.6, xend = 0.6, y = y_min - 0.4, yend = y_max + 0.4),
               colour = "grey40", linewidth = 0.6) +
  geom_text(angle = 90, fontface = "bold", size = 4.2, colour = "grey20", hjust = 0.5) +
  scale_y_continuous(limits = c(0.5, n_rows + 0.5), expand = expansion(0)) +
  scale_x_continuous(limits = c(-1, 1)) +
  theme_void() +
  theme(plot.margin = margin(5, 0, 5, 5))

x_pad <- max(abs(all_plot$diff_ppt)) * 0.50

p_lollipop <- ggplot(all_plot, aes(x = diff_ppt, y = label, colour = diff_ppt > 0)) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey50") +
  geom_segment(aes(x = 0, xend = diff_ppt, y = label, yend = label),
               linewidth = 1.1, alpha = 0.6) +
  geom_point(size = 4) +
  geom_text(aes(label = sprintf("%+.1f%%", diff_ppt),
                hjust = if_else(diff_ppt > 0, -0.25, 1.25)),
            size = 3.84, colour = "black") +
  scale_colour_manual(values = c("TRUE" = "#d73027", "FALSE" = "#2166ac"), guide = "none") +
  scale_x_continuous(expand = expansion(add = c(x_pad, x_pad))) +
  scale_y_discrete(labels = y_labels, expand = expansion(add = c(0.6, 1.6))) +
  # label the two sides of the zero line so the direction is read at a glance
  annotate("text", x = -x_pad * 0.15, y = n_rows + 1.15, label = "under-represented",
           hjust = 1, size = 4.2, fontface = "italic", colour = "#2166ac") +
  annotate("text", x = x_pad * 0.15, y = n_rows + 1.15, label = "over-represented",
           hjust = 0, size = 4.2, fontface = "italic", colour = "#d73027") +
  labs(title = "B", x = "Observed % minus null %", y = NULL) +
  theme_bw(base_size = 14) +
  theme(axis.text.y = element_text(size = 12),
        axis.text.x = element_text(size = 12),
        plot.margin = margin(5, 10, 5, 2))

p_combined <- p_bar + p_iucn_strip + p_lollipop +
  plot_layout(widths = c(1, 0.08, 1.4)) &
  theme(plot.title = element_text(face = "bold"))

run_name_label <- tolower(gsub("[^A-Za-z0-9]+", "_", trimws(run_tag)))
plot_file <- file.path(outDir, paste0("iucn_2025exposed_", run_name_label, ".pdf"))
ggsave(plot_file, p_combined, width = 14, height = 8)
ggsave(sub("\\.pdf$", ".png", plot_file), p_combined, width = 14, height = 8, dpi = 300)
cat("Saved:", plot_file, "\n")
