###getwd
getwd()
##setwd
setwd("E:/Armadillo_Analysis/Final_Armadillo_September_3/Common_SNPs/")

##intersection of the SNPs or the outlier loci that are under selection in each methods.
# Set working directory (already set)
setwd("E:/Armadillo_Analysis/Final_Armadillo_September_3/Common_SNPs/")

########
library(VennDiagram)

# Set the working directory (already done in your code)
setwd("E:/Armadillo_Analysis/Final_Armadillo_September_3/Common_SNPs/")

# Step 1: Read the CSV file (Common_SNPs.csv)
data <- read.csv("Common_SNPs.csv")

# Step 2: Extract SNP lists from each method
lfmm_snps <- data$LFMM
rda_snps <- data$RDA
pcaadapt_snps <- data$pcaadapt

# Step 3: Convert to vectors for easier set operations (if not already)
lfmm_snps <- as.vector(lfmm_snps)
rda_snps <- as.vector(rda_snps)
pcaadapt_snps <- as.vector(pcaadapt_snps)

# Step 4: Find common SNPs
# Intersection of all three methods
common_all <- intersect(intersect(lfmm_snps, rda_snps), pcaadapt_snps)

# Intersection of any two methods
common_LFMM_RDA <- intersect(lfmm_snps, rda_snps)
common_LFMM_pcaadapt <- intersect(lfmm_snps, pcaadapt_snps)
common_RDA_pcaadapt <- intersect(rda_snps, pcaadapt_snps)

# Step 5: Find unique SNPs for each method
unique_LFMM <- setdiff(lfmm_snps, union(rda_snps, pcaadapt_snps))
unique_RDA <- setdiff(rda_snps, union(lfmm_snps, pcaadapt_snps))
unique_pcaadapt <- setdiff(pcaadapt_snps, union(lfmm_snps, rda_snps))

# Step 6: Customize and Create the Venn Diagram
venn.plot <- venn.diagram(
  x = list(
    LFMM = lfmm_snps,
    RDA = rda_snps,
    pcaadapt = pcaadapt_snps
  ),
  category.names = c("LFMM", "RDA", "pcadapt"),
  filename = NULL,  # Do not save the file yet
  output = TRUE,
  
  # Customizations for color-blind friendly color scheme
  fill = c("#E69F00", "#56B4E9", "#009E73"),  # Unique colors for each region
  #cat.col = c("#E69F00", "#56B4E9", "#009E73")  # Colors for the method labels
  cat.cex = 2,  # Increase label font size for categories (methods)
  cat.pos = c(90, 5, 0),  # Positioning of category labels inside the diagram
  cex = 1.5,  # Size of the numbers in the regions
  alpha = 0.5,  # Transparency of the regions
  euler.d = TRUE,  # Euler diagram (more accurate overlap sizes)
  scaled = TRUE,  # Scaled to reflect the actual size of overlap
  margin = 0.1,  # Padding around the diagram
  main = "SNP Overlap Across Methods",  # Title of the plot
  main.cex = 2,  # Title size
  cat.dist = 0.05,  # Distance between the categories and labels
  cat.fontface = "bold",  # Make category labels bold
  cex.number = 1.5,  # Adjust the size of the numbers in each section
)

# Display the Venn diagram
grid.draw(venn.plot)

# Step 7: Save the Venn diagram as a high-resolution JPEG
jpeg("armadillo_Venn_diagram_final.jpeg", width = 12, height = 12, units = "in", res = 600)
grid.draw(venn.plot)
dev.off()



################################################################
#############################################################
# Load required libraries
# Load required libraries
library(VennDiagram)
library(grid)

# Set working directory
setwd("E:/Armadillo_Analysis/Final_Armadillo_September_3/Common_SNPs/")

# Read the data and properly handle empty values
data <- read.csv("Common_SNPs.csv", stringsAsFactors = FALSE)

# Remove empty values and NA values explicitly
lfmm_snps <- data$LFMM[!is.na(data$LFMM) & data$LFMM != ""]
rda_snps <- data$RDA[!is.na(data$RDA) & data$RDA != ""]
pcaadapt_snps <- data$pcaadapt[!is.na(data$pcaadapt) & data$pcaadapt != ""]

# Create the Venn diagram with increased numbers inside the regions and no method names
venn.plot <- venn.diagram(
  x = list(
    LFMM = lfmm_snps,
    RDA = rda_snps,
    pcaadapt = pcaadapt_snps
  ),
  filename = NULL,  # Do not save the file yet
  output = TRUE,
  
  # Customizations for color-blind friendly color scheme
  fill = c("#E69F00", "#56B4E9", "#009E73"),  # Unique colors for each region
  
  # Numbers inside regions - Increase the size of the numbers in the sections
  cex = 4,  # Size of the numbers inside the regions
  cex.number = 5,  # Size of the numbers inside the sections (larger)
  
  # General appearance
  alpha = 0.5,  # Transparency of the regions
  euler.d = TRUE,  # Euler diagram (more accurate overlap sizes)
  scaled = TRUE,  # Scaled to reflect actual size of overlap
  margin = 0.1,  # Padding around the diagram
  
  # Completely remove category names (method names)
  cat.cex = 0,  # Remove category labels completely
  
  # Title
  main = "SNP Overlap Across Methods",
  main.cex = 2.5,  # Title size
  main.fontface = "bold"
)

# Save as high-resolution JPEG
jpeg(
  filename = "armadillo_Venn_diagram_final.jpeg",
  width = 12,
  height = 12,
  units = "in",
  res = 600,
  quality = 100
)

# Draw the diagram
grid.draw(venn.plot)

# Close the device
dev.off()

# Display the diagram in R
grid.draw(venn.plot)
