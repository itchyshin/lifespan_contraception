# cause of death
# test

# loading packages

# loading
pacman::p_load(tidyverse,
               metafor,
               pander,
               knitr,
               stringr,
               ape,
               here,
               kableExtra,
               patchwork,
               lme4,
               readxl,
               #metaAidR,
               rotl,
               orchaRd,
               emmeans,
               clubSandwich,
               png,
               grid,
               here,
               cowplot,
               apextra,
               ggimage,
               ggstance,
               ggtree,
               phytools
)

# reading in the dataset


# main data
dat0 <- read_csv(here("data", "zoo2", "causeDeathAll.csv"), na = c("", "NA"))

# turning character strings into factors
dat0 <- dat0 %>% mutate(across(where(is.character), as.factor))

# length(unique(dat0$Species)) # the number of species 131

# phylogeny
tree <- read.tree(here("data", "zoo", "tree_zoo.tre"))

# taxonomy
tax <- read.csv(here("data", "zoo", "vertlife_taxonomy_translation_table.csv"))

dat0 %>% left_join(tax, by = c("Species" = "zims.species")) -> dat_full

dat <- dat_full

# talking out species with no data (Pseudocheirus peregrinus = likely to be mistaks in data)
# dat_full %>% filter(species != "Chrysocyon brachyurus" &
#                       species != "Crocuta crocuta" &
#                       species != "Neofelis nebulosa" &
#                       species != "Panthera uncia" &
#                       species != "Pseudocheirus peregrinus") %>% 
#   mutate(phylogeny = gsub(" ", "_", vertlife.species)) -> dat

# adding Cervus canadensis
#dat$vertlife.species[which(dat$Species == "Cervus canadensis")] <-"Cervus canadensis"
#dat$phylogeny[which(dat$species == "Cervus canadensis")] <-"Cervus_canadensis"

# fixing species name
dat$Species[dat$Species == "Equus asinus"] <- "Equus_africanus"
dat$Species[dat$Species == "Aonyx cinereus"] <- "Aonyx cinerea"
dat$Species[dat$Species == "Bubalus bubalis"] <- "Bubalus arnee"

dat$phylogeny <- gsub(" ", "_", dat$vertlife.species)

# life span data 
to_drop <-
  tree$tip.label[which(!(tree$tip.label %in% unique(dat$phylogeny)))]

tree <- drop.tip(tree, to_drop)

# checking the number of spp
#length(tree$tip.label)

tree <- force.ultrametric(tree)

#tree <- compute.brlen(tree)
cor_tree <- vcv(tree, corr = TRUE)

# meta-analysis plan
# we can use both log relative reisk and log response ratio (raito of means)

# function for raito of means

lnrrp <- function(m1, m2, n1, n2) {
  # arcsine transforamtion
  asin_trans <- function(p) { asin(sqrt(p)) }
  # SD for arcsine distribution (see Wiki - https://en.wikipedia.org/wiki/Arcsine_distribution)
  var1 <- 1/8
  var2 <- 1/8
  # lnRR - with 2nd order correction
  lnrr <- log(asin_trans(m1)/asin_trans(m2)) + 
    0.5 * ((var1 / (n1 * asin_trans(m1)^2)) - (var2 / (n2 * asin_trans(m2)^2)))	
  
  var <- var1 / (n1 * asin_trans(m1)^2) + var1^2 / (2 * n1^2 * asin_trans(m1)^4)  + 
    var2 / (n2 * asin_trans(m2)^2) + var2^2 / (2 * n2^2 * asin_trans(m2)^4) 
  
  invisible(data.frame(yi = lnrr , vi = var))
}

# function for relative risk

lnrrisk <- function(p1, p2, n1, n2) {
  e1 <- p1 * n1 # the nmber of dead in group 1
  e2 <- p2 * n2 # the number of dead in group 2
  
  a1 <- (1-p1)*n1 # the number of alive in group 1
  a2 <- (1-p2)*n2 # the number of alive in group 2
  
  # log relative risk
  
  lnrrisk <- log(p1/p2)
  
  vlnrrisk <- a1/(e1*(a1+e1)) + a2/(e2*(a2+e2))
  
  invisible(data.frame(yi2 = lnrrisk , vi2 = vlnrrisk))
}

