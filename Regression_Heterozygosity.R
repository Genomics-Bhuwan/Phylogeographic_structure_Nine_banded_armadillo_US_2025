###Author: Bhuwan Singh Bist
###Date:11/11/2024
###Jezkova lab
###Landscape genomics: Obs. heterozygosity vs. Longitude and Latitude
###########################################################################
#############################################################################
#############################################################################
##Regression of the plots
# Load required libraries
# Load required libraries
# Load required libraries
library(ggplot2)
library(gridExtra)
library(broom)  # For tidy regression output

## Get and set working directory
getwd()

setwd("E:/Armadillo_Analysis/Final_Armadillo_September_3/Regression/")
# Read the CSV file
data <- read.csv("Armadillo_Regression_Heterozygosity.csv")

# Function to add direction to latitude and longitude
add_direction <- function(value, is_latitude = TRUE) {
  if (is_latitude) {
    direction <- ifelse(value >= 0, "N", "S")
    value <- abs(value)
  } else {
    direction <- ifelse(value >= 0, "E", "W")
    value <- abs(value)
  }
  
  return(paste(value, direction, sep = "° "))
}

# Apply function to the data to add direction
data$Latitude_label <- sapply(data$Latitude, add_direction, is_latitude = TRUE)
data$Longitude_label <- sapply(data$Longitude, add_direction, is_latitude = FALSE)

# Fit linear models
lat_model <- lm(Observed.Heterozygosity ~ Latitude, data = data)
long_model <- lm(Observed.Heterozygosity ~ Longitude, data = data)

# Get model summaries
lat_summary <- summary(lat_model)
long_summary <- summary(long_model)

# Extract R-squared, p-values, and correlation coefficient (R)
lat_r2 <- lat_summary$r.squared
lat_r2
long_r2 <- long_summary$r.squared
lat_p <- pf(lat_summary$fstatistic[1], lat_summary$fstatistic[2], lat_summary$fstatistic[3], lower.tail = FALSE)
lat_p
long_p <- pf(long_summary$fstatistic[1], long_summary$fstatistic[2], long_summary$fstatistic[3], lower.tail = FALSE)

# Calculate correlation coefficients (R) for Latitude and Longitude
lat_r <- cor(data$Latitude, data$Observed.Heterozygosity)
long_r <- cor(data$Longitude, data$Observed.Heterozygosity)

# Create plots with correlation coefficient and p-value in the bottom-left corner
lat_plot <- ggplot(data, aes(x = Latitude, y = Observed.Heterozygosity)) +
  geom_point(color = "blue", alpha = 0.8) +
  geom_smooth(method = "lm", color = "red", se = TRUE) +
  theme_bw() +
  labs(
    x = "Latitude",
    y = "Observed Heterozygosity"
  ) +
  scale_x_continuous(labels = function(x) sapply(x, add_direction, is_latitude = TRUE)) +
  annotate("text", x = min(data$Latitude), y = min(data$Observed.Heterozygosity), 
           label = paste("R = ", round(lat_r, 3), "\np = ", format(lat_p, scientific = FALSE)),
           hjust = 0, vjust = 0, size = 5, color = "red") +
  theme(
    plot.title = element_blank(),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 10, face = "bold"),  # Make axis values bold
    axis.ticks.length = unit(0.05, "inches"),  # Make ticks very small
    axis.ticks = element_line(size = 0.2),  # Make ticks smaller
    axis.text.x = element_text(margin = margin(t = 5)),  # Keep axis values same
    axis.text.y = element_text(margin = margin(r = 5))   # Keep axis values same
  )
lat_plot

## Observed heterozygosity vs. Longitude

# Define significance based on p-value
significance <- ifelse(long_p < 0.05, "**", "")

# Create the plot with asterisks for significance
long_plot <- ggplot(data, aes(x = Longitude, y = Observed.Heterozygosity)) +
  geom_point(color = "blue", alpha = 0.6) +
  geom_smooth(method = "lm", color = "red", se = TRUE) +
  theme_bw() +
  labs(
    x = "Longitude",
    y = "Observed Heterozygosity"
  ) +
  scale_x_continuous(labels = function(x) sapply(x, add_direction, is_latitude = FALSE)) +
  annotate("text", x = min(data$Longitude), y = min(data$Observed.Heterozygosity), 
           label = paste("R = ", round(long_r, 3), "\np = ", format(long_p, scientific = FALSE), significance),
           hjust = 0, vjust = 0, size = 5, color = "red") +
  theme(
    plot.title = element_blank(),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 10, face = "bold"),  # Make axis values bold
    axis.ticks.length = unit(0.05, "inches"),  # Make ticks very small
    axis.ticks = element_line(size = 0.2),  # Make ticks smaller
    axis.text.x = element_text(margin = margin(t = 5)),  # Keep axis values same
    axis.text.y = element_text(margin = margin(r = 5))   # Keep axis values same
  )
long_plot

# Arrange plots side by side
grid.arrange(long_plot, lat_plot, ncol = 2)

######### Save each of the plots ############

# JPEG format
ggsave("lat_plot.jpeg", plot = lat_plot, dpi = 600, width = 6, height = 4, 
       units = "in", device = "jpeg")
ggsave("long_plot.jpeg", plot = long_plot, dpi = 600, width = 6, height = 4, 
       units = "in", device = "jpeg")
ggsave("combined_plot.jpeg", plot = grid.arrange(long_plot, lat_plot, ncol = 2),
       dpi = 600, width = 16, height = 12, 
       units = "in", device = "jpeg")

# TIFF format
ggsave("lat_plot.tiff", plot = lat_plot, dpi = 600, width = 6, height = 4, 
       units = "in", device = "tiff")
ggsave("long_plot.tiff", plot = long_plot, dpi = 600, width = 6, height = 4, 
       units = "in", device = "tiff")
ggsave("combined_plot.tiff", plot = grid.arrange(long_plot, lat_plot, ncol = 2),
       dpi = 600, width = 16, height = 12, 
       units = "in", device = "tiff")


##*******************It was fun coding*************************
##*************************************************************
##*************************************************************
##*************************************************************
##*************************************************************
##*************************************************************


