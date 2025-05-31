#### I am assessing the population genetic structure with k=2 
#
## I am assessing pop.gen structure, individual and population (Ancestry barchart, pie-maps)

# Authors: Bhuwan Singh Bist
# Date: 11/07/2024

##Note: Save the file in jpeg and other four vector format.
# a. svg b. pdf c. wmf d)eps
# These files will be used in adobe illustrator for visualization.
# As I am editing, combining two images, SVG is best for me. 

## getwd
getwd()

### set working directory 
setwd("E:/Armadillo_Analysis/Final_Armadillo_September_3/mapMixture")

install.packages("mapmixture")
# install.packages("devtools")
devtools::install_github("Tom-Jenkins/mapmixture", force = TRUE)

# Load packages
library(maps)
library(mapmixture)
library(ggplot2)
library(gridExtra)
library(maps)
library(ggspatial)

# Set working directory
setwd("E:/Armadillo_Analysis/Final_Armadillo_September_3/mapMixture")

# Read in data
admixture1 <- read.csv("admixture2.csv")
coordinates <- read.csv("coordinates.csv")

# Get US state polygons
states <- map_data("state")

# Run mapmixture and build map
map5 <- mapmixture(
  admixture_df = admixture1,
  coords_df = coordinates,
  cluster_cols = c("lightgreen", "#6a3d9a"),
  cluster_names = c("Ancestry 1", "Ancestry 2"),
  crs = 4326,
  boundary = c(xmin = -103, xmax = -80.3, ymin = 25, ymax = 41),
  pie_size = 0.5
) +
  geom_polygon(data = states, aes(x = long, y = lat, group = group), fill = NA, color = "black", size = 0.5) +
  geom_text(data = subset(states, lat == min(lat)), aes(x = long, y = lat, label = region), size = 3, color = "black") +
  theme(
    legend.position = c(0.88, 0.08),
    legend.justification = c("right", "bottom"),
    legend.title = element_text(size = 14, face = "bold"),
    legend.text = element_text(size = 12),
    plot.margin = margin(l = 10, r = 10),
    axis.title.x = element_text(size = 16, face = "bold", margin = margin(t = 20)),  # Bold X-axis label
    axis.title.y = element_text(size = 16, face = "bold", margin = margin(r = 20)),  # Bold Y-axis label
    axis.text.x = element_text(size = 12, color = "black", hjust = 0.5, vjust = 0.5, margin = margin(t = -20)),
    axis.text.y = element_text(size = 12, color = "black", hjust = 0.5, vjust = 0.5, margin = margin(r = -35)),
    axis.ticks.x = element_line(size = 0.5, color = "black"),
    axis.ticks.y = element_line(size = 0.5, color = "black"),
    axis.ticks.length = unit(-0.2, "cm")
  ) +
  guides(fill = guide_legend(override.aes = list(size = 5, alpha = 1))) +
  labs(
    x = expression("Longitude"),
    y = expression("Latitude")
  )

# Show map
print(map5)

# Save plot
ggsave("E:/Armadillo_Analysis/Final_Armadillo_September_3/mapMixture/map_plot.jpeg",
       plot = map5, width = 14, height = 6, dpi = 600, device = "jpeg")

ggsave("E:/Armadillo_Analysis/Final_Armadillo_September_3/mapMixture/map_plot.tiff",
       plot = map5, width = 14, height = 6, dpi = 600, device = "tiff")



# Traditional structure barplot
##################
# Plot with structure_plot using the ordered data
structure_barplot <- structure_plot(
  admixture_df = admixture1,
  type = "structure",
  cluster_cols = c("lightgreen", "#6a3d9a"),
  site_dividers = TRUE,
  site_order = c(
    "FL1", "FL3", "FL4", "GA2", "MO5", "MO3", "TX5", "TX6", "MO4", 
    "TX7", "TX15", "TX10", "TX13", "TX14", "TX11", "TX4", "TX16", 
    "TX12", "TX9", "AL11", "IL2", "TX3", "AL2", "MO1", "TX18", 
    "TX1", "TX2", "MO2", "TX17", "AL12", "OK2", "OK4", "AR7", 
    "OK3", "AL13", "IL3", "OK1", "AR2", "MS3", "AL10", "TX19", 
    "MS4", "IL5", "MS5", "GA5", "FL5", "TX8", "AR5", "IL4", 
    "MS6", "IL6", "TX21", "AL3", "AL7", "TN1", "FL2", "AR1",
    "AR9", "GA4", "TN3", "AL5", "AR6", "AL9", "GA3", "TN2", 
    "AL1", "MO6", "AL14", "IL1", "GA1", "TX20", "AL4", "FL6",
    "MS1", "MS2", "AL6", "AL8", "AR3", "AR4", "AR8"
  ),
  divider_width = 0.4,
  labels = "site",
  flip_axis = FALSE,
  site_ticks_size = 0,
  site_labels_size = 0
) +
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_text(size = 12, angle = 90, hjust = 1, vjust = 0.5, face = "bold"),
    axis.text.y = element_text(size = 12, face = "bold"),  # Y-axis values bold and same size as x
    axis.title.y = element_text(size = 16, face = "bold", margin = margin(r = 15)),  # Larger y-label with space
    axis.ticks.x = element_line(),  # Show x-axis ticks
    plot.margin = margin(b = 20, l = 80, r = 20, t = 20),
    panel.grid = element_blank(),
    strip.text = element_blank()
  ) +
  ylab("Ancestry proportion")

# Save structure barplot
ggsave("E:/Armadillo_Analysis/Final_Armadillo_September_3/mapMixture/structure_barplot.jpeg", 
       plot = structure_barplot, 
       width = 20, 
       height = 4,
       dpi = 600,
       device = "jpeg")

# Save structure barplot as TIFF
ggsave("E:/Armadillo_Analysis/Final_Armadillo_September_3/mapMixture/structure_barplot.tiff", 
       plot = structure_barplot, 
       width = 20, 
       height = 4, 
       dpi = 600, 
       device = "tiff")

