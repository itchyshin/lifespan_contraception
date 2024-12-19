# make pre- and post purbety


rm(list = ls())

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
               ggtree
)


dat0 <- read_csv(here("data", "zoo", "timing_male.csv"), na = c("", "NA"))

# phylogeny
# read RData
#dford 

load(here("Rdata", "zoo", "maxCredTreeMammals.RData"))

tree <- maxCred #read.tree(here("data", "zoo", "tree_zoo.tre"))

# taxonomy
tax <- read.csv(here("data", "zoo", "species_merge_list.csv"))

dat0 %>% left_join(tax, by = c("Species" = "ZIMSSpecies")) -> dat_full

# talking out species with no data (Pseudocheirus peregrinus = likely to be mistaks in data)
# dat_full %>% filter(Species != "Chrysocyon brachyurus" &
#                     Species != "Crocuta crocuta" &
#                     Species != "Neofelis nebulosa" &
#                     Species != "Panthera uncia" &
#                     Species != "Pseudocheirus peregrinus") %>% 
#   mutate(phylogeny = gsub(" ", "_", vertlife.species)) -> dat

dat_full %>% 
  mutate(phylogeny = gsub(" ", "_", vertlife.treename)) -> dat

# adding Cervus canadensis
# dat$vertlife.species[which(dat$Species == "Cervus canadensis")] <-"Cervus canadensis"
# dat$phylogeny[which(dat$Species == "Cervus canadensis")] <-"Cervus_canadensis"

# fixing species names
#dat$Species[dat$Species == "Equus asinus"] <- "Equus_africanus"
dat$Species[dat$Species == "Aonyx cinereus"] <- "Aonyx cinerea"
#dat$Species[dat$Species == "Bubalus bubalis"] <- "Bubalus arnee"


# life span data 
to_drop <-
  tree$tip.label[which(!(tree$tip.label %in% unique(dat$phylogeny)))]

tree <- drop.tip(tree, to_drop)

# checking the number of spp
#length(tree$tip.label)

tree <- as.ultrametric(tree)

#tree <- compute.brlen(tree)
cor_tree <- vcv(tree, corr = TRUE)


missing <- which(is.na(dat$phylogeny))

dat$Species[missing]

# creating effect size 

dat <- escalc("ROM", 
                    m2i = LifeExpContraAft,
                    m1i = LifeExpContraBef,
                    sd2i = SEcontraAft*sqrt(NContraAft),
                    sd1i = SEcontraBef*sqrt(NContraBef), 
                    n2i = NContraAft,
                    n1i = NContraBef,
                    data = dat,
)

dat$species <- factor(dat$Species)

# meta-analysis

# variance-covariance matrix for sampling error assuming 0.5 correlation

mod_pp <- rma.mv(yi, V = vi,
                  random = list(
                    ~1|species,
                    ~1|phylogeny),
                  R = list(phylogeny = cor_tree),
                  data = dat)
summary(mod_pp)
round(i2_ml(mod_pp), 2)

#robust(mod_all, cluster = species)  

orchard_plot(mod_pp, xlab = "lnRR (all)", group = "species", g = FALSE)

# data set excluding primates

dat %>% filter(Primate == "No") -> dat_short


mod_pp2 <- rma.mv(yi, V = vi,
                 random = list(
                   ~1|species,
                   ~1|phylogeny),
                 R = list(phylogeny = cor_tree),
                 data = dat_short)
summary(mod_pp2)
round(i2_ml(mod_pp2), 2)

#robust(mod_all, cluster = species)  

orchard_plot(mod_pp2, xlab = "lnRR (all)", group = "species", g = FALSE)
