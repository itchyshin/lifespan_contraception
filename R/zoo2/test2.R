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

# ***************************************************************
#   *                          Note:                              *
#   *    force.ultrametric does not include a formal method to    *
#   *    ultrametricize a tree & should only be used to coerce    *
#   *   a phylogeny that fails is.ultrametric due to rounding --  *
#   *    not as a substitute for formal rate-smoothing methods.   *
#   ***************************************************************

#tree <- compute.brlen(tree)
cor_tree <- vcv(tree, corr = TRUE)




