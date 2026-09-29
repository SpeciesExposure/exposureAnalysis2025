# enso_strip.R
# Colour strip of the NOAA Oceanic Niño Index (ONI v5, three-month running mean of Niño-3.4 SST
# anomaly) for the water years used in this analysis (October–September). Used beneath the
# time-series figures so readers can see El Niño (red) and La Niña (blue) conditions without a
# binary classification. Source file downloaded from
# https://www.cpc.ncep.noaa.gov/data/indices/oni.ascii.txt on 2026-09-15.
#
# oni_wateryear(years) -> data frame with one row per season, positioned within the water year
# oni_strip(years, ...) -> ggplot tile strip aligned to integer year positions (bars of width 1)

oni_file_default <- KEY_FILES$oni

oni_wateryear <- function(years, oni_file = oni_file_default) {
  oni <- utils::read.table(oni_file, header = TRUE, stringsAsFactors = FALSE)
  seas <- c("DJF","JFM","FMA","MAM","AMJ","MJJ","JJA","JAS","ASO","SON","OND","NDJ")
  # central calendar month of each season
  oni$month <- match(oni$SEAS, seas)
  # water year: October–September, named for the year in which it ends
  oni$wy <- ifelse(oni$month >= 10, oni$YR + 1, oni$YR)
  oni$wm <- ifelse(oni$month >= 10, oni$month - 9, oni$month + 3)   # Oct = 1 ... Sep = 12
  oni <- oni[oni$wy %in% years, ]
  # x position: bar for year y spans [y - 0.5, y + 0.5]; place each season at its month centre
  oni$x <- oni$wy - 0.5 + (oni$wm - 0.5) / 12
  oni[, c("wy", "wm", "SEAS", "YR", "ANOM", "x")]
}

oni_strip <- function(years, text_size = 8, show_legend = FALSE, strip_label = "ONI") {
  # legend off by default: the colour scale (blue -2.5 to red +2.5 degrees C) is stated in the figure captions
  d <- oni_wateryear(years)
  ggplot2::ggplot(d, ggplot2::aes(x = x, y = 1, fill = ANOM)) +
    ggplot2::geom_tile(width = 1 / 12, height = 1) +
    ggplot2::scale_fill_gradient2(low = "#2166ac", mid = "white", high = "#b2182b", midpoint = 0,
                                  limits = c(-2.5, 2.5), oob = scales::squish,
                                  name = "ONI (°C)", breaks = c(-2, -1, 0, 1, 2)) +
    ggplot2::scale_x_continuous(limits = c(min(years) - 0.5, max(years) + 0.5), expand = c(0, 0)) +
    ggplot2::scale_y_continuous(expand = c(0, 0)) +
    ggplot2::labs(x = NULL, y = strip_label) +
    ggplot2::theme_void() +
    ggplot2::theme(axis.title.y = ggplot2::element_text(size = text_size, angle = 0, vjust = 0.5, hjust = 1),
                   legend.position = if (show_legend) "right" else "none",
                   legend.key.height = grid::unit(0.3, "cm"), legend.key.width = grid::unit(0.25, "cm"),
                   legend.title = ggplot2::element_text(size = text_size), legend.text = ggplot2::element_text(size = text_size - 1),
                   plot.margin = ggplot2::margin(0, 0, 2, 0))
}

# ---- ENSO year lists: strictly from the NOAA CPC episode table (ONI v5, data/enso/) ----
# Each episode winter (a DJF season inside a >=5-consecutive-season run beyond +/-0.5 C) is assigned to the
# water year (Oct-Sep) containing that DJF. Strong El Nino = winter peak ONI (max of NDJ/DJF/JFM) >= 1.5 C.
# 2025 is included following NOAA's January 2025 La Nina advisory. Its winter minimum of -0.46 C sits
# just under the episode threshold in the final index.
EL_NINO_STRONG <- c(1992, 1998, 2010, 2016, 2024)
EL_NINO_OTHER  <- c(1995, 2003, 2005, 2007, 2015, 2019, 2020)
LA_NINA        <- c(1996, 1999, 2000, 2001, 2006, 2008, 2009, 2011, 2012, 2018, 2021, 2022, 2023, 2025)

enso_years <- function(years) {
  data.frame(year = years, el_nino = years %in% c(EL_NINO_STRONG, EL_NINO_OTHER),
             el_nino_strong = years %in% EL_NINO_STRONG, la_nina = years %in% LA_NINA)
}

# Add the symbols to a bar chart whose bars are at integer x = year. `bar_tops`: data frame with
# columns year and top (bar height). Red up-triangle = El Niño (larger = strong), blue down-triangle = La Niña.
enso_marks <- function(years, bar_tops, size = 2.4, gap_frac = 0.035) {
  # Each episode winter is drawn at the boundary between the two years it spans (x = water year - 0.5),
  # as in Fig. 1 of the 2023 report, raised above the taller of the two adjacent bars.
  e <- enso_years(years); t <- setNames(bar_tops$top, bar_tops$year)
  e$top <- pmax(t[as.character(e$year)], t[as.character(e$year - 1)], 0, na.rm = TRUE)
  e$x <- e$year - 0.5
  e$y <- e$top + gap_frac * max(bar_tops$top, na.rm = TRUE)
  e$sz <- ifelse(e$el_nino_strong, size * 1.6, size)
  list(
    ggplot2::geom_point(data = e[e$el_nino, ], ggplot2::aes(x = x, y = y, size = sz), inherit.aes = FALSE, shape = 24, fill = "#d73027", colour = "black", stroke = 0.3, show.legend = FALSE),
    ggplot2::geom_point(data = e[e$la_nina, ], ggplot2::aes(x = x, y = y), inherit.aes = FALSE, shape = 25, size = size, fill = "#2166ac", colour = "black", stroke = 0.3),
    ggplot2::scale_size_identity())
}

# In-panel key for the three symbols: place at (x, y) in data units; dy = row spacing
enso_key <- function(x, y, dy, size = 2.4, text_size = 3.4) {
  list(
    ggplot2::annotate("point", x = x, y = y, shape = 24, size = size * 1.6, fill = "#d73027", colour = "black", stroke = 0.3),
    ggplot2::annotate("text", x = x + 0.6, y = y, label = "Strong El Ni\u00f1o", hjust = 0, size = text_size),
    ggplot2::annotate("point", x = x, y = y - dy, shape = 24, size = size, fill = "#d73027", colour = "black", stroke = 0.3),
    ggplot2::annotate("text", x = x + 0.6, y = y - dy, label = "El Ni\u00f1o", hjust = 0, size = text_size),
    ggplot2::annotate("point", x = x, y = y - 2 * dy, shape = 25, size = size, fill = "#2166ac", colour = "black", stroke = 0.3),
    ggplot2::annotate("text", x = x + 0.6, y = y - 2 * dy, label = "La Ni\u00f1a", hjust = 0, size = text_size))
}
