# FigS4_IUCN_Chronic.r
# Figure S4: IUCN Red List status of chronically exposed species, those exposed
# in 6 or more years between 1990 and 2025, compared with all study species in
# the same taxonomic group.
#
# Inputs:  singleMaxExposure_v1.qs, spAttributes_v8.qs
# Outputs: figure in plotDir

suppressPackageStartupMessages({
  library(qs2)
  library(dplyr)
  library(ggplot2)
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

sp_meta_file <- KEY_FILES$sp_attributes_v8
if (!file.exists(sp_meta_file)) stop("Cannot find spAttributes file.")

CVAR_MAP <- c(
  temp__12_up = "temp_heat",
  temp__3__max_up = "temp_heat",
  temp__12_lo = "temp_cold",
  temp__3__min_lo = "temp_cold",
  precip__12_up = "precip_flood",
  precip__3__max_up = "precip_flood",
  precip__12_lo = "precip_dry",
  precip__3__min_lo = "precip_dry"
)

cat("Loading singleMaxExposure...\n")
single_max_candidates <- unique(path.expand(c(
  file.path(intDir, "singleMaxExposure_v1.qs"),
  if (exists("KEY_FILES") && is.list(KEY_FILES) && !is.null(KEY_FILES$single_max_exposure_v1)) {
    KEY_FILES$single_max_exposure_v1
  } else character(0),
  file.path(PROJECT_PATHS$int_v8, "singleMaxExposure_v1.qs"),
  file.path(PROJECT_PATHS$output_v8_intermediate, "singleMaxExposure_v1.qs"),
  file.path(PROJECT_PATHS$output_base, tolower(run_tag), paste0("Int_", figure_version), "singleMaxExposure_v1.qs"),
  file.path(PROJECT_PATHS$output_base, tolower(run_tag), "intermediate", "singleMaxExposure_v1.qs")
)))

single_max_file <- single_max_candidates[file.exists(single_max_candidates)][1]
if (is.na(single_max_file) || !nzchar(single_max_file)) {
  stop(
    "Could not locate singleMaxExposure_v1.qs. Checked:\n",
    paste0(" - ", single_max_candidates, collapse = "\n")
  )
}

cat("Using singleMaxExposure file:", single_max_file, "\n")
ex <- tryCatch(qs_read(single_max_file), error = function(e) {
  cat("qs2 read failed, trying qs::qread...\n")
  qs::qread(single_max_file)
})
ex <- ex %>% mutate(cvar = CVAR_MAP[var])

cat("Loading spAttributes...\n")
sp_meta <- qs_read(sp_meta_file) %>%
  distinct(spName, .keep_all = TRUE) %>%
  dplyr::select(spName, group, redlistCategory)

chronic_sp <- ex %>%
  distinct(spName, year) %>%
  count(spName, name = "n_yrs") %>%
  filter(n_yrs >= 6) %>%
  pull(spName)

dominant_cvar <- function(x) {
  x <- x[!is.na(x)]
  if (length(x) == 0) return(NA_character_)
  names(which.max(table(x)))[1]
}

sp_chronic_summary <- ex %>%
  filter(spName %in% chronic_sp) %>%
  group_by(spName) %>%
  summarise(
    n_yrs = n_distinct(year),
    dom_cvar = dominant_cvar(cvar),
    max_pExp = max(pExp, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  left_join(sp_meta, by = "spName")

THREATENED <- c("Critically_Endangered", "Endangered", "Vulnerable")
LC_CATS <- c(
  "Least_Concern",
  "Lower_Risk/least concern",
  "Lower_Risk/conservation dependent",
  "Lower_Risk/near threatened"
)

sp_chronic_summary <- sp_chronic_summary %>%
  mutate(
    iucn_simple = case_when(
      redlistCategory %in% THREATENED ~ "Threatened",
      redlistCategory == "Near_Threatened" ~ "Near Threatened",
      redlistCategory %in% LC_CATS ~ "Least Concern",
      redlistCategory == "Data_Deficient" ~ "Data Deficient",
      TRUE ~ "Other/NA"
    )
  )

iucn_by_group <- sp_chronic_summary %>%
  filter(!is.na(group), iucn_simple != "Other/NA") %>%
  count(group, iucn_simple)

iucn_order <- c("Least Concern", "Near Threatened", "Data Deficient", "Threatened")
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

# Group order: match Fig. 3 (most 2025-exposed species first)
group_order_bar <- ex %>%
  filter(year == 2025) %>%
  distinct(spName, group) %>%
  filter(!is.na(group), group != "") %>%
  count(group, sort = TRUE) %>%
  pull(group)

suppressPackageStartupMessages(library(patchwork))

p_iucn_bar <- ggplot(
  iucn_by_group %>% mutate(group = factor(group, levels = group_order_bar)),
  aes(x = group, y = n, fill = factor(iucn_simple, levels = iucn_order))
) +
  geom_col(width = 0.7, colour = "white", linewidth = 0.3) +
  scale_y_continuous(name = "Number of chronically exposed species") +
  scale_fill_manual(values = iucn_cols, name = "IUCN status") +
  labs(title = "A", x = NULL) +
  theme_bw(base_size = 14) +
  theme(legend.position = "right",
        legend.key.size  = unit(0.4, "cm"),
        legend.text      = element_text(size = 11),
        legend.title     = element_text(size = 11),
        axis.text.x      = element_text(angle = 30, hjust = 1))

# ---- Panel B: lollipop of % deviation from null (chronically exposed species, >=6 yrs) ----
IUCN_FOCUS <- c("Threatened", "Near Threatened", "Least Concern", "Data Deficient")

# sp_chronic_summary already has group + iucn_simple; rename for consistency
chr_sp_iucn <- sp_chronic_summary %>%
  filter(!is.na(group), iucn_simple %in% IUCN_FOCUS) %>%
  rename(iucn_broad = iucn_simple) %>%
  mutate(group = if_else(is.na(group) | group == "", "Unknown", group))

bg_sp <- sp_meta %>%
  distinct(spName, .keep_all = TRUE) %>%
  mutate(
    iucn_broad = case_when(
      redlistCategory %in% c("Critically_Endangered", "Endangered", "Vulnerable") ~ "Threatened",
      redlistCategory %in% c("Near_Threatened")                                   ~ "Near Threatened",
      redlistCategory %in% c("Least_Concern", "Lower_Risk/least concern",
                              "Lower_Risk/conservation dependent",
                              "Lower_Risk/near threatened")                        ~ "Least Concern",
      redlistCategory %in% c("Data_Deficient")                                    ~ "Data Deficient",
      TRUE ~ "Other/Unknown"
    ),
    group = if_else(is.na(group) | group == "", "Unknown", group)
  )

group_order_lollipop <- group_order_bar

obs_counts <- chr_sp_iucn %>%
  filter(group != "Unknown") %>%
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
  filter(iucn_broad %in% IUCN_FOCUS, group %in% group_order_lollipop)

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

group_annot <- all_plot %>%
  mutate(y_int = as.integer(label)) %>%
  group_by(iucn_broad) %>%
  summarise(y_mid = mean(y_int), y_min = min(y_int), y_max = max(y_int), .groups = "drop")

n_rows <- nlevels(all_plot$label)

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

p_combined <- p_iucn_bar + p_iucn_strip + p_lollipop +
  plot_layout(widths = c(1, 0.08, 1.4)) &
  theme(plot.title = element_text(face = "bold"))

run_name_label <- if (exists('RUN_NAME') && nzchar(trimws(RUN_NAME))) {
  RUN_NAME
} else {
  figure_version
}
run_name_label <- tolower(gsub('[^A-Za-z0-9]+', '_', trimws(run_name_label)))
plot_iucn <- file.path(outDir, paste0("iucn_gap_by_group_", run_name_label, ".pdf"))

ggsave(plot_iucn, p_combined, width = 14, height = 8)
ggsave(sub("\\.pdf$", ".png", plot_iucn), p_combined, width = 14, height = 8, dpi = 300)
cat("Saved:", plot_iucn, "\n")

