##Author: Bhuwan Singh Bist
##Date:08/26/2024
##Objectives: I want to use Mantel test to see the correlation between genetic distance and 
##geographic distance and see if they give a common result

##Install the packages
install.packages(c("vcfR","adegenet", "ade4","geosphere","dplyr"))
install.packages("fields")
# Load required libraries
library(vcfR)
library(adegenet)
library(ade4)
library(geosphere)
library(dplyr)
library(fields)
########## CODE FOR IBD ANALYSES ####################################
library(vcfR)
library(adegenet)
library(adegraphics)
library(pegas)
library(StAMPP)
library(lattice)
library(gplots)
library(ape)
library(ggmap)
library(spThin)

#Load in vcf file for %85 of data 
# Set working directory if needed
getwd()
###setworking directory
setwd("E:/Armadillo_Analysis/Final_Armadillo_September_3/IBD")

###read the vcf file
nf_85 <- read.vcfR("Armadillos_filtered.DP5g50maf05.recode.vcf")

#load in pop information 
nf_pop <- read.csv("E:/Armadillo_Analysis/Final_Armadillo_September_3/IBD/Armadillo_Pop.csv")
nf_pop
### convert to genlight
NF_85.genlight <- vcfR2genlight(nf_85, n.cores=2)
NF_85.genlight
#add pop names
pop(NF_85.genlight)<- nf_pop[,2]
pop(NF_85.genlight)

#generate list of genlight objects 
NF_SNPs <- list(NF_85.genlight)
NF_SNPs
#Calculation and visualization of Nei’s distances (using lapply)
### Calculate Nei's distances between individuals/pops
NF.ind <- lapply(NF_SNPs, stamppNeisD, pop = FALSE) # Nei's 1972 distance between indivs
NF.ind
NF.pop <- lapply(NF_SNPs, stamppNeisD, pop = TRUE) # Nei's 1972 distance between pops
NF.pop


NF_85.genlight@ploidy <- as.integer(ploidy(NF_85.genlight))
NF_85.genlight@ploidy

### Isolation by distance
coords <- read.csv ("E:/Armadillo_Analysis/Final_Armadillo_September_3/IBD/Armadillo_Pop.csv") # tab-separated file for all pops
coords
xy.coords.only<- subset(coords, select=c("Long","Lat"))
xy.coords.only
Dgeo <- dist(xy.coords.only)
Dgeo

###########################
## Calculate distance in km 
DistMat <- rdist.earth(x1 = coords[3:4], miles = FALSE)
DistMat <- as.dist(DistMat)


# create the dist objects used in analyses below
colnames(NF.ind[[1]]) <- rownames(NF.ind[[1]])
colnames(NF.ind[[1]])
NF.85.ind.dist<-as.dist(NF.ind[[1]], diag=T)
NF.85.ind.dist
attr(NF.85.ind.dist, "Labels")<-rownames(NF.ind[[1]]) # name the rows of a matrix


# Species-level Mantel test
species_mantel <- mantel.randtest(NF.85.ind.dist, DistMat, nrepet = 10000)
species_mantel


# Print species-level results
cat("Species-level Mantel test results:\n")
print(species_mantel)

#####################################
################################
#################################
# Plot species-level IBD as PDF
pdf("Species_IBD_plot.pdf")
plot(DistMat, NF.85.ind.dist, 
     xlab = "",  # Remove the default x-axis label
     ylab = "",  # Remove the default y-axis label
     main = "Species-level IBD",
     ylim = c(0.15, 0.35)) # Set ylim for the values between 0 and 0.4
abline(lm(as.vector(NF.85.ind.dist) ~ as.vector(DistMat)), col = "red", lwd = 2) # Increase line thickness

# Add custom x-axis and y-axis labels with controlled distance
mtext("Geographic Distance (km)", side = 1, line = 2)  # Adjust 'line' for distance from x-axis
mtext("Genetic Distance", side = 2, line = 2)          # Adjust 'line' for distance from y-axis

dev.off()

# Plot species-level IBD as JPEG
jpeg("Species_IBD_plot1.jpeg", width = 2000, height = 1600, res = 300) # Open a JPEG device with 300 DPI
plot(as.vector(DistMat), as.vector(NF.85.ind.dist), 
     xlab = "",  # Remove the default x-axis label
     ylab = "",  # Remove the default y-axis label
     ylim = c(0.15, 0.35)) # Set ylim for the values between 0 and 0.4
abline(lm(as.vector(NF.85.ind.dist) ~ as.vector(DistMat)), col = "red", lwd = 4) # Increase line thickness

# Add custom x-axis and y-axis labels with controlled distance
mtext("Geographic Distance (km)", side = 1, line = 2)  # Adjust 'line' as needed
mtext("Genetic Distance", side = 2, line = 2)          # Adjust 'line' as needed

dev.off() # Close the JPEG device

# Reset mgp to default for any future plots
par(mgp = c(3, 1, 0))

# Save all results to a text file
sink("Mantel_test_results.txt")
cat("Species-level Mantel test results:\n")
print(species_mantel)
cat("\nPopulation-level Mantel test results:\n")
sink()




############################
##########################
# Calculate linear model for p-value and R^2
# Calculate linear model for p-value and R^2
lm_model <- lm(as.vector(NF.85.ind.dist) ~ as.vector(DistMat))
summary_lm <- summary(lm_model)

# Extract R-squared and p-value
r_squared <- summary_lm$r.squared
p_value <- summary_lm$coefficients[2, 4]  # p-value for the slope

# Format R-squared to 5 decimal places
r_squared_text <- formatC(r_squared, format = "f", digits = 5)

# Define the p-value text based on its value
p_value_text <- ifelse(p_value < 0.05, "p < 0.05**", paste0("p = ", formatC(p_value, format = "f", digits = 5)))

# Plot species-level IBD as PDF
pdf("Species_IBD_plot.pdf")
plot(DistMat, NF.85.ind.dist, 
     xlab = "",  
     ylab = "",  
     main = "Species-level IBD",
     ylim = c(0.15, 0.35)) 
abline(lm_model, col = "red", lwd = 2)

# Add custom axis labels
mtext("Geographic Distance (km)", side = 1, line = 2)
mtext("Genetic Distance", side = 2, line = 2)

# Display R² and p-value in the bottom right corner
text(x = max(DistMat) * 0.95, y = 0.17, 
     labels = paste0("R² = ", r_squared_text, "\n", p_value_text),
     col = "red", adj = c(1, 0))

dev.off()

# Plot species-level IBD as JPEG
jpeg("Species_IBD_plot1.jpeg", width = 2000, height = 1600, res = 300)
plot(as.vector(DistMat), as.vector(NF.85.ind.dist), 
     xlab = "",  
     ylab = "",  
     ylim = c(0.15, 0.35)) 
abline(lm_model, col = "red", lwd = 4)

# Add custom axis labels
mtext("Geographic Distance (km)", side = 1, line = 2)
mtext("Genetic Distance", side = 2, line = 2)

# Display R² and p-value in the bottom right corner
text(x = max(DistMat) * 0.95, y = 0.17, 
     labels = paste0("R² = ", r_squared_text, "\n", p_value_text),
     col = "red", adj = c(1, 0))

dev.off()

