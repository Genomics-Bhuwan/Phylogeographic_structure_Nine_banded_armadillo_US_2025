###set working directory.
setwd("E:/Armadillo_Analysis/Final_Armadillo_September_3/EEMS_Fourth_attempt/")
##install devtools
library("devtools")
install_github("dipetkov/reemsplots2")
library("reemsplots2")
library(terra)
library(tidyterra)
library(tidyverse)
library(ggspatial)
library(reemsplots2)

# Load US states shapefile
states <- vect("cb_2018_us_state_20m.shp")

# Load and filter occurrence data
occ <- read.csv("Dasy_original_GBIF+92_samples.csv") %>%
  filter(countryCode == "US") %>%
  filter(Latitude < 39.8)  # Filter to southern populations

# Create SpatVector from occurrence data
occ <- vect(occ, crs = "EPSG:4326", geom = c("Longitude", "Latitude"))

# Create study area boundary (without buffer)
occ_shp <- terra::convHull(occ) %>%
  terra::crop(states)

# Create state label data for all states
state_df <- as.data.frame(crds(centroids(states))) %>%
  mutate(NAME = states$NAME)

# Load EEMS results
mcmcpath <- "E:/Armadillo_Analysis/Final_Armadillo_September_3/EEMS_Fourth_attempt/aaaa/EEMS_Armadillo_not_FEEMS-output/EEMS_res_forBSB/res1"
plots <- make_eems_plots(mcmcpath, longlat = TRUE)

# Create final plot
final_plot <- plots$mrates01 +
  # Crop to study area
  coord_sf(
    xlim = ext(occ_shp)[1:2],
    ylim = ext(occ_shp)[3:4],
    expand = FALSE,
    crs = "EPSG:4326"
  ) +
  # Add color scale
  scale_fill_whitebox_c(
    palette = "muted",
    na.value = "gray10",
    guide = guide_colorbar(title = "log(m)")
  ) +
  # Add state outlines
  geom_spatvector(data = states, fill = NA, linewidth = 0.5, inherit.aes = FALSE) +
  geom_text(data = state_df, aes(x = x, y = y, label = NAME), size = 3) +
  # Add north arrow and scale
  annotation_north_arrow(
    location = "tr",
    which_north = "true",
    pad_x = unit(0.5, "cm"),
    pad_y = unit(0.5, "cm")
  ) +
  annotation_scale() +
  # Theme and axis customization
  theme_bw() +
  xlab("Longitude") +
  ylab("Latitude") +
  theme(
    axis.text.x = element_text(size = 14),
    axis.text.y = element_text(size = 14),
    axis.title.x = element_text(size = 18, face = "bold"),
    axis.title.y = element_text(size = 18, face = "bold"),
    # Adjusted legend position
    legend.position = c(0.75, 0.18),
    legend.direction = "horizontal",
    legend.box = "horizontal",
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 10),
    legend.background = element_rect(fill = "white", color = "black")
  )

# Save high-resolution image
ggsave(
  "EEMS_Final_Map.jpg",
  plot = final_plot,
  width = 12,
  height = 8,
  dpi = 600,
  device = "jpeg"
)
