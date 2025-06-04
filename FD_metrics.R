## R Code to Replicate FD Analysis ##
# Packages
library("knitr")
library("corrplot")
library("dplyr")
library("hypervolume")
library("BAT")
#################################################################
#################################################################
## Loading the trait and community databases and preparing the data
# 1- A data frame summarizing traits values for each species
# Load data
traits_birds <- read.delim("~/traits.txt", row.names=1)
as.data.frame(traits_birds)
# Evaluate trait correlation
corr_traits=cor(traits_birds)
head(round(corr_traits,2))
corrplot(corr_traits, method="number", type="lower")
# Drop correlated variables:
traits_birds <- subset(traits_birds, select = -c(tarsus))
traits_birds
# Standardize
traits_birds <- traits_birds %>% mutate_at(c("mass",
                                             "h.w",
                                             "wing"), ~(scale(.) %>% as.vector))
traits_birds
########
# 2- A matrix summarizing species assemblages.
# Excluding trees with zero species (26 trees)
commun_birds_26 <- read.delim("~/comm_26.txt", row.names=1)
commun_birds_26 <- as.matrix(commun_birds_26)
commun_birds_26
####
# Excluding trees with less than 2 species (22 trees)
commun_birds_22 <- read.delim("~/comm_22.txt", row.names=1)
commun_birds_22 <- as.matrix(commun_birds_22)
3
commun_birds_22
#################################################################
#################################################################
## Build hypervolumes ##
# Estimate bandwidth
BAND <- hypervolume::estimate_bandwidth(traits_birds)
BAND
# Hypervolumes - 22 trees
hv_22 <- BAT::kernel.build(comm = commun_birds_22, trait = traits_birds,
                           distance = "euclidean",
                           abund = TRUE,
                           method = "gaussian",
                           axes = 0,
                           kde.bandwidth = BAND,
                           cores = 0)
hv_22
# Hypervolumes - 26 trees
hv_26 <- BAT::kernel.build(comm = commun_birds_26, trait = traits_birds,
                           distance = "euclidean",
                           abund = TRUE,
                           method = "gaussian",
                           axes = 0,
                           kde.bandwidth = BAND,
                           cores = 0)
hv_26
#################################################################
## Calculating the functional richness of each plot using the kernel.alpha function
richness_26 <- kernel.alpha(comm = hv_26)
richness_22 <- kernel.alpha(comm = hv_22)
#################################################################
## Calculating the functional dispersion of each plot using the kernel.dispersion function
divergence_26 <- kernel.dispersion(comm = hv_26,frac=0.1, func= 'divergence')
divergence_22 <- kernel.dispersion(comm = hv_22,frac=0.1, func= 'divergence')
#################################################################
## Calculating the functional regularity of each plot using the kernel.evenness function
4
regularity_22 <- kernel.evenness(comm = hv_22)
#################################################################