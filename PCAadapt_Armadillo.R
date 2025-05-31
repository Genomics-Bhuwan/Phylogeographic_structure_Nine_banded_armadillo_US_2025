##Author: Bhuwan Singh Bist
##Date: 08/29/2024
##pcadapt has been developed to detect the genetic markers involved in biological adaptation.
##it provides statistical tools for outlier detection based on PCA.


################
###install.packages
install.packages(c("vcfR", "pcadapt", "qvalue"))
install.packages("dplyr")

####library
library(vcfR)
library(pcadapt)
library(qvalue)
library(tidyverse)


## try http if https is not available
install.packages("devtools")
library("devtools")
install_github("jdstorey/qvalue")
browseVignettes(package = "qvalue")
library(qvalue)
### pcadapt is used to test for outliers and useful in terms of visualization.

##getwd
getwd()
##setwd
setwd("E:/Armadillo_Analysis/Final_Armadillo_September_3/Pcadapt")

# Assuming your population data is stored in a data frame 'pop_data' with 'Sample' and 'Population' columns
pop_data <- read.csv("E:/Armadillo_Analysis/Final_Armadillo_September_3/Pcadapt/Armadillo_Pop.csv")  # Load your population data
pop_data

#########################################################
####Used plink to convert vcf to .bed, .fam, .bim file.
#########################################################
###########################################################
##### Run pcadapt #########################################
###########################################################
###############################################################
# First step is to determine the number of principal components
###############################################################
Genotypes<- read.pcadapt("E:/Armadillo_Analysis/Final_Armadillo_September_3/Pcadapt/Armadillos_filtered_pcadapt.bed", type = "bed")
Genotypes  

###use plink to make the bedfile
###plink --vcf input_file.vcf --make-bed --out lenRef.pcadapt

x <- pcadapt(Genotypes, K = 20) 
x$pvalues


#####Screeplot
plot(x, option = "screeplot")
x <- pcadapt(Genotypes, K = 10)
x

# Open a JPEG device with specified dimensions and resolution
jpeg("Screeplot_PCA1.jpeg", width = 2000, height = 1600, res = 300)
plot(x, option = "screeplot")
dev.off()

####score plot
# Assuming pop_data is already loaded and contains the population names
# Create a vector with population names corresponding to each individual
poplist.names <- pop_data$Population
poplist.names
# Check that poplist.names aligns with the number of individuals in Genotypes
print(poplist.names)

# Plot the PCA scores with population names
plot(x, option = "scores", pop = poplist.names)
colors = c("Gold", "Limegreen")

# Open a JPEG device with specified dimensions and resolution
jpeg("Scoreplot_PCA1.jpeg", width = 2000, height = 1600, res = 300)
plot(x, option = "scores", i = , j = 2, pop = poplist.names)
plot(x, option = "scoreplot")
dev.off()


# Define a color palette that matches the number of unique populations
unique_pops <- length(unique(poplist.names))
unique_pops
colors <- rainbow(unique_pops)  # Generates a palette with as many colors as unique populations
# Interactive plotting with plotly
plot(x, option = "scores", pop = poplist.names, plt.pkg = "plotly", col = colors)

# project your PCA with corresponding colors
# use plotly code for an interactive plot

####here we plotted multiple principal component
##we plotted the principal component 1 and 2
plot(x, option = "scores", i = , j = 2, pop = poplist.names)

# Open a JPEG device with specified dimensions and resolution
jpeg("Scoreplot_PCA1.jpeg", width = 2000, height = 1600, res = 300)
plot(x, option = "scores", i = , j = 2, pop = poplist.names)
plot(x, option = "screeplot")
dev.off()

##############################################
# Second step is to compute our test statistic
##############################################
x <- pcadapt(Genotypes, K = 2)
x
#run based on your choice of K value

summary(x)
#look at your run


#####################################
# Third step is to graph your results
#####################################
jpeg("Manhattan_plot_1.jpeg", width = 2000, height = 1600, res = 300)
plot(x , option = "manhattan")
dev.off()

## creates a QQ-plot from your run
# creates a qqplot from your run
# looking for most p-values to be expected but some to stray on the side of significance 
jpeg("QQ_plot_1.jpeg", width = 2000, height = 1600, res = 300)
plot(x, option = "qqplot", threshold = 0.1)
dev.off()

# Plot the Mahalanobis distance
plot(x, option = "stat.distribution")

# plot a histogram of the pvalues, want to see a uniform distribution
hist(x$pvalues, xlab = "p-values", main = NULL, breaks = 50, col = "navy")




###################################
# Fourth step is to choose a cutoff
###################################

# view Storey 2010 for a good explanation of q-values
qval <- qvalue(x$pvalues)$qvalues
qval
alpha <- 0.1
outliersQ <- which(qval < alpha)
outliersQ
length(outliersQ)
outliers = which(qval < alpha)
length(outliersQ)

Qvalue<-as.data.frame(outliersQ)
Qvalue
candSNP = Qvalue
candSNP
candSNP = as.data.frame(unlist(candSNP))
candSNP
colnames(candSNP)<- "Position"
colnames(candSNP)
candSNP


write.csv(Qvalue, file = "QOutliers.csv")
# tests for SNPs with a qvalue less than 0.1
# gives us how many outliers were detected

###Let's test for Benjamini-Hochberg procedure
padj <- p.adjust(x$pvalues,method="BH")
alpha <- 0.1
outliers <- which(padj < alpha)
length(outliers)

##Bonferroni correction
padj <- p.adjust(x$pvalues,method="bonferroni")
alpha <- 0.1
outliers <- which(padj < alpha)
length(outliers)


######################################
##### Now let's make a cool plot #####
######################################
######################################
##### Now let's make a cool plot #####
######################################
####################
jpeg("Final_Manhattan_plot_1.jpeg", width = 2000, height = 1600, res = 300)
# Calculate -log10 of the p-values
log_pvalues <- -log10(x$pvalues)

# Define colors based on the condition
colors <- ifelse(log_pvalues > 5, "red", "skyblue")

# Set sizes based on the condition
sizes <- ifelse(log_pvalues > 5, 1.5, 1)  # Larger size for points above the threshold

# Plot the main data points with conditional colors and sizes
plot(log_pvalues, xlab = "SNP Id", ylab = "-log10(P-values)", 
     ylim = c(0, 20),
     cex = sizes, pch = 21,
     col = colors, bg = colors)  # Use the color vector here

# Add a horizontal abline at y = 5 to indicate the threshold
abline(h = 5, col = "blue", lty = 2, lwd = 2)  # Blue dashed line

# Save the plot
dev.off()







