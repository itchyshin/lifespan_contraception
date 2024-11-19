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
               #apextra,
               ggimage,
               ggstance,
               ggtree,
               phytools,
               orchaRd
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

# there are missing species in the taxonomy table
dat0 %>% left_join(tax, by = c("Species" = "zims.species")) %>% 
  # make character strings into factors
  mutate(across(where(is.character), as.factor)) -> dat_full

dat <- dat_full


# talking out species with no data (Pseudocheirus peregrinus = likely to be mistakes in data)
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
# dat$Species[dat$Species == "Equus asinus"] <- "Equus africanus"
# dat$Species[dat$Species == "Aonyx cinereus"] <- "Aonyx cinerea"
# dat$Species[dat$Species == "Bubalus bubalis"] <- "Bubalus arnee"

# filling in NA in Species = Aonyx_cinerea

#dat$Species[which(is.na(dat$Species))] <- "Aonyx cinerea"


dat$Phylogeny <- gsub(" ", "_", dat$vertlife.species)

#TODO - we need to do phylogenetic tree matching

#life span data
to_drop <-
  tree$tip.label[which(!(tree$tip.label %in% unique(dat$Phylogeny)))]

tree <- drop.tip(tree, to_drop)

# checking the number of spp
#length(tree$tip.label)

tree <- compute.brlen(tree)

#tree <- compute.brlen(tree)
cor_tree <- vcv(tree, corr = TRUE)


# samller data set which does not have NA in dat$Phylogeny

dat <- dat %>% filter(!is.na(Phylogeny))

# Effect_ID is the unique identifier for the effect

dat$Effect_ID <- factor(1:nrow(dat))


# creating a varaible combining Sex and Type
dat$Sex_Type <- as.factor(paste0(dat$Sex, "_", dat$Type))


# meta-analysis plan
# we can use both log relative reisk and log response ratio (raito of means)

# function for raito of means

# lnrrp <- function(m1, m2, n1, n2) {
#   # if p1 or p2 = 0, turn into 0.025, if p1 or p2 = 1, turn into 0.975
#   m1[m1 == 0] <- 0.025
#   m1[m1 == 1] <- 0.975
#   m2[m2 == 0] <- 0.025
#   m2[m2 == 1] <- 0.975
#   # arcsine transforamtion
#   asin_trans <- function(p) { asin(sqrt(p)) }
#   # SD for arcsine distribution (see Wiki - https://en.wikipedia.org/wiki/Arcsine_distribution)
#   var1 <- 1/8
#   var2 <- 1/8
#   
#   # lnRR - with 2nd order correction
#   lnrr <- log(asin_trans(m1)/asin_trans(m2)) + 
#     0.5 * ((var1 / (n1 * asin_trans(m1)^2)) - (var2 / (n2 * asin_trans(m2)^2)))	
#   
#   var <- var1 / (n1 * asin_trans(m1)^2) + var1^2 / (2 * n1^2 * asin_trans(m1)^4)  + 
#     var2 / (n2 * asin_trans(m2)^2) + var2^2 / (2 * n2^2 * asin_trans(m2)^4) 
#   
#   invisible(data.frame(yi = lnrr , vi = var))
# }

# function for relative risk

# lnrrisk <- function(p1, p2, n1, n2) {
#   # if p1 or p2 = 0, turn into 0.025, if p1 or p2 = 1, turn into 0.975
#   if (p1 == 0) { p1 <- 0.025 }
#   
#   # if n1 or n2 = 0, lnrrisk and vlnrrisk should be NA
#   if (n1 == 0 | n2 == 0) {
#   
#     lnrrisk <- NA
#     vlnrrisk <- NA
#     
#   } else{
#   
#   e1 <- p1 * n1 # the nmber of dead in group 1
#   e2 <- p2 * n2 # the number of dead in group 2
#   
#   a1 <- (1-p1)*n1 # the number of alive in group 1
#   a2 <- (1-p2)*n2 # the number of alive in group 2
#   
#   # log relative risk
#   
#   lnrrisk <- log(p1/p2)
#   
#   vlnrrisk <- a1/(e1*(a1+e1)) + a2/(e2*(a2+e2))}
#   
#   invisible(data.frame(yi = lnrrisk , vi = vlnrrisk))
# }


# different kinds of causes of death
# Trauma
# Infectious Disease
# Non-infectious Disease
# Chronic Disease
# Death at birth
# Other

# calcuating 2 x effect sizes by applying 2 functions above to dat


# Trauma

# Low
dat_trauma <- escalc(measure = "RD", 
               ai = Contra_Trauma_Low*Contra_Trauma_N, 
               bi = (1-Contra_Trauma_Low)*Contra_Trauma_N, 
               ci = noContra_Trauma_Low*noContra_Trauma_N,
               di = (1-noContra_Trauma_Low)*noContra_Trauma_N,
               var.names = c("yi_trauma_low", "vi_trauma_low"),
               data = dat)
# Med

dat_trauma <- escalc(measure = "RD", 
               ai = Contra_Trauma_Med*Contra_Trauma_N, 
               bi = (1-Contra_Trauma_Med)*Contra_Trauma_N, 
               ci = noContra_Trauma_Med*noContra_Trauma_N,
               di = (1-noContra_Trauma_Med)*noContra_Trauma_N,
               var.names = c("yi_trauma_med", "vi_trauma_med"),
               data = dat_trauma)

# Upp

dat_trauma <- escalc(measure = "RD", 
               ai = Contra_Trauma_Upp*Contra_Trauma_N, 
               bi = (1-Contra_Trauma_Upp)*Contra_Trauma_N, 
               ci = noContra_Trauma_Upp*noContra_Trauma_N,
               di = (1-noContra_Trauma_Upp)*noContra_Trauma_N,
               var.names = c("yi_trauma_upp", "vi_trauma_upp"),
               data = dat_trauma)

dat_trauma %>% filter(Contra_Trauma_N > 0) -> dat_trauma

# create a long format of the data  using these 3 types of effect sizes (low, med, upp) yi and vi are the effect size and variance of the effect size

dat_long_trauma <- dat_trauma %>% select(Effect_ID, Species, Phylogeny, Sex_Type,
                           yi_trauma_low, vi_trauma_low, 
                           yi_trauma_med, vi_trauma_med, 
                           yi_trauma_upp, vi_trauma_upp) %>% 
  pivot_longer(cols = c(yi_trauma_low, yi_trauma_med, yi_trauma_upp, vi_trauma_low, vi_trauma_med, vi_trauma_upp), 
               names_to = c(".value", "type"), 
               names_pattern = "(yi|vi)_(.*)")

dat_long_trauma <- dat_long_trauma %>% filter(!is.na(yi))

str(dat_long_trauma)

dat_long_trauma$type <- factor(dat_long_trauma$type, levels = rev(c("trauma_low", "trauma_med", "trauma_upp")),
                        labels = rev(c("Lower", "Median", "Upper")))

# effect size level ID

dat_long_trauma$Effect_ID2 <- factor(1 : nrow(dat_long_trauma))


# VCV

VCV <- vcalc(dat_long_trauma$vi, 
            cluster = dat_long_trauma$Effect_ID, 
            obs = dat_long_trauma$Effect_ID2, 
            data = dat_long_trauma, rho = 0.5)


# meta-analysis using dat_long

mod_trauma <- rma.mv(yi = yi, V = VCV, 
                     random = list(
                       ~1|Species,
                       ~1|Phylogeny,
                       #~1|Effect_ID, 
                       ~1|Effect_ID2),
                     R = list(Phylogeny = cor_tree), 
                     data = dat_long_trauma,
                     method="REML", 
                     sparse=TRUE,
                     control=list(optimizer="optim", optmethod="Nelder-Mead")
)


summary(mod_trauma)

# meta-regression

mod_trauma_reg <- rma.mv(yi = yi, V = VCV, 
                         random = list(
                           ~1|Species,
                           ~1|Phylogeny,
                           ~1|Effect_ID2), 
                         #struct = "DIAG",
                         R = list(Phylogeny = cor_tree), 
                         data = dat_long_trauma,
                         method="REML", 
                         control=list(optimizer="optim", optmethod= "Nelder-Mead"),
                         mods = ~type - 1,
                         sparse=TRUE
)

summary(mod_trauma_reg)


p_trauma <- orchard_plot(mod_trauma_reg, mod = "type",
    xlab = "log risk ratio (Trauma)", group = "Species") + ylim(-0.7, 0.7)


p_trauma


# checking Sex_Type effect

# filter out Male_Hormonal
dat_long_trauma %>% filter(Sex_Type != "Male_Hormonal") -> dat_long_trauma2

# VCV2

VCV2 <- vcalc(dat_long_trauma2$vi, 
              cluster = dat_long_trauma2$Effect_ID, 
              obs = dat_long_trauma2$Effect_ID2, 
              data = dat_long_trauma2, rho = 0.5)

mod_trauma_reg2 <- rma.mv(yi = yi, V = VCV2, 
                         random = list(
                           ~1|Species,
                           ~1|Phylogeny,
                           ~1|Effect_ID2), 
                         R = list(Phylogeny = cor_tree), 
                         data = dat_long_trauma2,
                         method="REML", 
                         control=list(optimizer="optim", optmethod= "Nelder-Mead"),
                         mods = ~Sex_Type - 1,
                         sparse=TRUE
)

p_trauma2 <- orchard_plot(mod_trauma_reg2, mod = "Sex_Type",
    xlab = "log risk ratio (Trauma)", group = "Species", angle = 45) + ylim(-0.85, 0.7)

p_trauma2


# Infectious Disease

# Low

dat_infectious <- escalc(measure = "RD", 
               ai = Contra_InfectDisease_Low*Contra_InfectDisease_N,
               bi = (1-Contra_InfectDisease_Low)*Contra_InfectDisease_N,
               ci = noContra_InfectDisease_Low*noContra_InfectDisease_N,
               di = (1-noContra_InfectDisease_Low)*noContra_InfectDisease_N,
               var.names = c("yi_infectious_low", "vi_infectious_low"),
               data = dat)

# Med

dat_infectious <- escalc(measure = "RD",
               ai = Contra_InfectDisease_Med*Contra_InfectDisease_N, 
               bi = (1-Contra_InfectDisease_Med)*Contra_InfectDisease_N, 
               ci = noContra_InfectDisease_Med*noContra_InfectDisease_N,
               di = (1-noContra_InfectDisease_Med)*noContra_InfectDisease_N,
               var.names = c("yi_infectious_med", "vi_infectious_med"),
               data = dat_infectious)

# Upp

dat_infectious <- escalc(measure = "RD",
               ai = Contra_InfectDisease_Upp*Contra_InfectDisease_N, 
               bi = (1-Contra_InfectDisease_Upp)*Contra_InfectDisease_N, 
               ci = noContra_InfectDisease_Upp*noContra_InfectDisease_N,
               di = (1-noContra_InfectDisease_Upp)*noContra_InfectDisease_N,
               var.names = c("yi_infectious_upp", "vi_infectious_upp"),
               data = dat_infectious)

dat_infectious %>% filter(Contra_InfectDisease_N > 0) -> dat_infectious


# create a long format of the data using these 3 types of effect sizes (low, med, upp) yi and vi are the effect size and variance of the effect size

dat_long_infectious <- dat_infectious %>% select(Effect_ID, Species, Phylogeny, Sex_Type,
                           yi_infectious_low, vi_infectious_low, 
                           yi_infectious_med, vi_infectious_med, 
                           yi_infectious_upp, vi_infectious_upp) %>% 
  pivot_longer(cols = c(yi_infectious_low, yi_infectious_med, yi_infectious_upp, vi_infectious_low, vi_infectious_med, vi_infectious_upp), 
               names_to = c(".value", "type"), 
               names_pattern = "(yi|vi)_(.*)")

dat_long_infectious <- dat_long_infectious %>% filter(!is.na(yi))

str(dat_long_infectious)

dat_long_infectious$type <- factor(dat_long_infectious$type, 
                           levels = rev(c("infectious_low", "infectious_med", "infectious_upp")),
                           labels = rev(c("Lower", "Median", "Upper"))
)

# effect size level ID

dat_long_infectious$Effect_ID2 <- factor(1 : nrow(dat_long_infectious))

# VCV

VCV <- vcalc(dat_long_infectious$vi, 
            cluster = dat_long_infectious$Effect_ID, 
            obs = dat_long_infectious$Effect_ID2, 
            data = dat_long_infectious, rho = 0.5)


# meta-analysis using dat_long

mod_infectious <- rma.mv(yi = yi, V = VCV, 
                         random = list(
                           ~1|Species,
                           ~1|Phylogeny,
                           ~1|Effect_ID2), 
                         R = list(Phylogeny = cor_tree), 
                         data = dat_long_infectious,
                         test = "t",
                         method="REML", 
                         control=list(optimizer="optim", optmethod="Nelder-Mead")
)

summary(mod_infectious)

# meta-regression

mod_infectious_reg <- rma.mv(yi = yi, V = VCV, 
                             random = list(
                               ~1|Species,
                               ~1|Phylogeny,
                               ~1|Effect_ID2), 
                             #struct = "DIAG",
                             R = list(Phylogeny = cor_tree), 
                             data = dat_long_infectious,
                             method="REML", 
                             test = "t",
                             control=list(optimizer="optim", optmethod="Nelder-Mead"),
                             mods = ~type - 1,
                             sparse=TRUE
)

summary(mod_infectious_reg)

p_infectious <- orchard_plot(mod_infectious_reg, mod = "type",
    xlab = "log risk ratio (Infectious Disease)", group = "Species") + ylim(-0.85, 0.7)

p_infectious

# checking Sex_Type effect

mod_infectious_reg2 <- rma.mv(yi = yi, V = VCV, 
                             random = list(
                               ~1|Species,
                               ~1|Phylogeny,
                               ~1|Effect_ID2), 
                             R = list(Phylogeny = cor_tree), 
                             data = dat_long_infectious,
                             method="REML", 
                             control=list(optimizer="optim", optmethod="Nelder-Mead"),
                             mods = ~Sex_Type - 1,
                             test = "t",
                             sparse=TRUE
)

summary(mod_infectious_reg2)

p_infectious2 <- orchard_plot(mod_infectious_reg2, mod = "Sex_Type",
                              xlab = "log risk ratio (Infectious Diesase)", group = "Species", angle = 45) + ylim(-0.85, 0.7)

p_infectious2
    

# Non-infectious Disease

# Low

dat_noninfectious <- escalc(measure = "RD", 
               ai = Contra_NonInfectDisease_Low*Contra_NonInfectDisease_N,
               bi = (1-Contra_NonInfectDisease_Low)*Contra_NonInfectDisease_N,
               ci = noContra_NonInfectDisease_Low*noContra_NonInfectDisease_N,
               di = (1-noContra_NonInfectDisease_Low)*noContra_NonInfectDisease_N,
               var.names = c("yi_noninfectious_low", "vi_noninfectious_low"),
               data = dat)

# Med

dat_noninfectious <- escalc(measure = "RD",
               ai = Contra_NonInfectDisease_Med*Contra_NonInfectDisease_N, 
               bi = (1-Contra_NonInfectDisease_Med)*Contra_NonInfectDisease_N, 
               ci = noContra_NonInfectDisease_Med*noContra_NonInfectDisease_N,
               di = (1-noContra_NonInfectDisease_Med)*noContra_NonInfectDisease_N,
               var.names = c("yi_noninfectious_med", "vi_noninfectious_med"),
               data = dat_noninfectious)

# Upp

dat_noninfectious <- escalc(measure = "RD",
               ai = Contra_NonInfectDisease_Upp*Contra_NonInfectDisease_N, 
               bi = (1-Contra_NonInfectDisease_Upp)*Contra_NonInfectDisease_N, 
               ci = noContra_NonInfectDisease_Upp*noContra_NonInfectDisease_N,
               di = (1-noContra_NonInfectDisease_Upp)*noContra_NonInfectDisease_N,
               var.names = c("yi_noninfectious_upp", "vi_noninfectious_upp"),
               data = dat_noninfectious)

dat_noninfectious %>% filter(Contra_NonInfectDisease_N > 0) -> dat_noninfectious

# create a long format of the data using these 3 types of effect sizes (low, med, upp) yi and vi are the effect size and variance of the effect size

dat_long_noninfectious <- dat_noninfectious %>% select(Effect_ID, Species, Phylogeny, Sex_Type,
                           yi_noninfectious_low, vi_noninfectious_low, 
                           yi_noninfectious_med, vi_noninfectious_med, 
                           yi_noninfectious_upp, vi_noninfectious_upp) %>% 
  pivot_longer(cols = c(yi_noninfectious_low, yi_noninfectious_med, yi_noninfectious_upp, vi_noninfectious_low, vi_noninfectious_med, vi_noninfectious_upp), 
               names_to = c(".value", "type"), 
               names_pattern = "(yi|vi)_(.*)")

dat_long_noninfectious <- dat_long_noninfectious %>% filter(!is.na(yi))

str(dat_long_noninfectious)

dat_long_noninfectious$type <- factor(dat_long_noninfectious$type, 
                           levels = rev(c("noninfectious_low", "noninfectious_med", "noninfectious_upp")),
                           labels = rev(c("Lower", "Median", "Upper"))
                           )

# effect size level ID

dat_long_noninfectious$Effect_ID2 <- factor(1 : nrow(dat_long_noninfectious))

# VCV

VCV <- vcalc(dat_long_noninfectious$vi, 
            cluster = dat_long_noninfectious$Effect_ID, 
            obs = dat_long_noninfectious$Effect_ID2, 
            data = dat_long_noninfectious, rho = 0.5)

# meta-analysis using dat_long

mod_noninfectious <- rma.mv(yi = yi, V = VCV, 
                           random = list(
                             ~1|Species,
                             ~1|Phylogeny,
                             ~1|Effect_ID2), 
                           R = list(Phylogeny = cor_tree), 
                           data = dat_long_noninfectious,
                           method="REML", 
                           test = "t",
                           control=list(optimizer="optim", optmethod="Nelder-Mead")
)

summary(mod_noninfectious)

# meta-regression

mod_noninfectious_reg <- rma.mv(yi = yi, V = VCV, 
                               random = list(
                                 ~1|Species,
                                 ~1|Phylogeny,
                                 ~1|Effect_ID2), 
                               #struct = "DIAG",
                               R = list(Phylogeny = cor_tree), 
                               data = dat_long_noninfectious,
                               method="REML", 
                               control=list(optimizer="optim", optmethod="Nelder-Mead"),
                               mods = ~type - 1,
                               test = "t",
                               sparse=TRUE
)

summary(mod_noninfectious_reg)

p_noninfectious <- orchard_plot(mod_noninfectious_reg, mod = "type",
    xlab = "log risk ratio (Non-infectious Disease)", group = "Species") + ylim(-0.85, 0.7)

p_noninfectious

# Sex_Type effect

# filter out Male_Hormonal

dat_long_noninfectious %>% filter(Sex_Type != "Male_Hormonal") -> dat_long_noninfectious2

# VCV2

VCV2 <- vcalc(dat_long_noninfectious2$vi, 
            cluster = dat_long_noninfectious2$Effect_ID, 
            obs = dat_long_noninfectious2$Effect_ID2, 
            data = dat_long_noninfectious2, rho = 0.5)


mod_noninfectious_reg2 <- rma.mv(yi = yi, V = VCV2, 
                               random = list(
                                 ~1|Species,
                                 ~1|Phylogeny,
                                 ~1|Effect_ID2), 
                               R = list(Phylogeny = cor_tree), 
                               data = dat_long_noninfectious2,
                               method="REML", 
                               control=list(optimizer="optim", optmethod="Nelder-Mead"),
                               mods = ~Sex_Type - 1,
                               test = "t",
                               sparse=TRUE
)

summary(mod_noninfectious_reg2)

p_noninfectious2 <- orchard_plot(mod_noninfectious_reg2, mod = "Sex_Type",
    xlab = "log risk ratio (Non-infectious Disease)", group = "Species", angle = 45) + ylim(-0.85, 0.7)

p_noninfectious2


# Chronic Disease

# Low

dat_chronic <- escalc(measure = "RD", 
               ai = Contra_ChronDisease_Low*Contra_ChronDisease_N,
               bi = (1-Contra_ChronDisease_Low)*Contra_ChronDisease_N,
               ci = noContra_ChronDisease_Low*noContra_ChronDisease_N,
               di = (1-noContra_ChronDisease_Low)*noContra_ChronDisease_N,
               var.names = c("yi_chronic_low", "vi_chronic_low"),
               data = dat)

# Med

dat_chronic <- escalc(measure = "RD",
               ai = Contra_ChronDisease_Med*Contra_ChronDisease_N, 
               bi = (1-Contra_ChronDisease_Med)*Contra_ChronDisease_N, 
               ci = noContra_ChronDisease_Med*noContra_ChronDisease_N,
               di = (1-noContra_ChronDisease_Med)*noContra_ChronDisease_N,
               var.names = c("yi_chronic_med", "vi_chronic_med"),
               data = dat_chronic)

# Upp

dat_chronic <- escalc(measure = "RD",
               ai = Contra_ChronDisease_Upp*Contra_ChronDisease_N, 
               bi = (1-Contra_ChronDisease_Upp)*Contra_ChronDisease_N, 
               ci = noContra_ChronDisease_Upp*noContra_ChronDisease_N,
               di = (1-noContra_ChronDisease_Upp)*noContra_ChronDisease_N,
               var.names = c("yi_chronic_upp", "vi_chronic_upp"),
               data = dat_chronic)

dat_chronic %>% filter(Contra_ChronDisease_N > 0) -> dat_chronic

# create a long format of the data using these 3 types of effect sizes (low, med, upp) yi and vi are the effect size and variance of the effect size

dat_long_chronic <- dat_chronic %>% select(Effect_ID, Species, Phylogeny, Sex_Type,
                           yi_chronic_low, vi_chronic_low, 
                           yi_chronic_med, vi_chronic_med, 
                           yi_chronic_upp, vi_chronic_upp) %>% 
  pivot_longer(cols = c(yi_chronic_low, yi_chronic_med, yi_chronic_upp, vi_chronic_low, vi_chronic_med, vi_chronic_upp), 
               names_to = c(".value", "type"), 
               names_pattern = "(yi|vi)_(.*)")

dat_long_chronic <- dat_long_chronic %>% filter(!is.na(yi))

str(dat_long_chronic)

dat_long_chronic$type <- factor(dat_long_chronic$type, 
                           levels = rev(c("chronic_low", "chronic_med", "chronic_upp")),
                           labels = rev(c("Lower", "Median", "Upper"))
                           )

# effect size level ID

dat_long_chronic$Effect_ID2 <- factor(1 : nrow(dat_long_chronic))

# VCV

VCV <- vcalc(dat_long_chronic$vi, 
            cluster = dat_long_chronic$Effect_ID, 
            obs = dat_long_chronic$Effect_ID2, 
            data = dat_long_chronic, rho = 0.5)

# meta-analysis using dat_long

mod_chronic <- rma.mv(yi = yi, V = VCV, 
                     random = list(
                       ~1|Species,
                       ~1|Phylogeny,
                       ~1|Effect_ID2), 
                     R = list(Phylogeny = cor_tree), 
                     data = dat_long_chronic,
                     method="REML", 
                     test = "t",
                     dfs = "contain",
                     control=list(optimizer="optim", optmethod="Nelder-Mead")
)

summary(mod_chronic)

# meta-regression

mod_chronic_reg <- rma.mv(yi = yi, V = VCV, 
                         random = list(
                           ~1|Species,
                           ~1|Phylogeny,
                           ~1|Effect_ID2), 
                         struct = "DIAG",
                         R = list(Phylogeny = cor_tree), 
                         data = dat_long_chronic,
                         method="REML", 
                         control=list(optimizer="optim", optmethod="Nelder-Mead"),
                         mods = ~type - 1,
                         test = "t",
                         sparse=TRUE
)

summary(mod_chronic_reg)

p_chronic <- orchard_plot(mod_chronic_reg, mod = "type",
    xlab = "log risk ratio (Chronic Disease)", group = "Species") + ylim(-0.85, 0.7)

p_chronic

# Sex_Type effect

# filter out Male_Hormonal

dat_long_chronic2 <- dat_long_chronic %>% filter(Sex_Type != "Male_Hormonal") 

# VCV2
VCV2 <- vcalc(dat_long_chronic2$vi, 
             cluster = dat_long_chronic2$Effect_ID, 
             obs = dat_long_chronic2$Effect_ID2, 
             data = dat_long_chronic2, rho = 0.5)

mod_chronic_reg2 <- rma.mv(yi = yi, V = VCV2, 
                         random = list(
                           ~1|Species,
                           ~1|Phylogeny,
                           ~1|Effect_ID2), 
                         R = list(Phylogeny = cor_tree), 
                         data = dat_long_chronic2,
                         method="REML", 
                         control=list(optimizer="optim", optmethod="Nelder-Mead"),
                         mods = ~Sex_Type - 1,
                         test = "t",
                         sparse=TRUE
)

summary(mod_chronic_reg2)

p_chronic2 <- orchard_plot(mod_chronic_reg2, mod = "Sex_Type",
    xlab = "log risk ratio (Chronic Disease)", group = "Species" ,angle = 45) + ylim(-0.85, 0.7)

p_chronic2

# Death at birth

# Low

dat_deathAtBirth <- escalc(measure = "RD", 
               ai = Contra_deathAtBirth_Low*Contra_deathAtBirth_N,
               bi = (1-Contra_deathAtBirth_Low)*Contra_deathAtBirth_N,
               ci = noContra_deathAtBirth_Low*noContra_deathAtBirth_N,
               di = (1-noContra_deathAtBirth_Low)*noContra_deathAtBirth_N,
               var.names = c("yi_deathAtBirth_low", "vi_deathAtBirth_low"),
               data = dat)

# Med

dat_deathAtBirth <- escalc(measure = "RD",
               ai = Contra_deathAtBirth_Med*Contra_deathAtBirth_N, 
               bi = (1-Contra_deathAtBirth_Med)*Contra_deathAtBirth_N, 
               ci = noContra_deathAtBirth_Med*noContra_deathAtBirth_N,
               di = (1-noContra_deathAtBirth_Med)*noContra_deathAtBirth_N,
               var.names = c("yi_deathAtBirth_med", "vi_deathAtBirth_med"),
               data = dat_deathAtBirth)

# Upp

dat_deathAtBirth <- escalc(measure = "RD",
               ai = Contra_deathAtBirth_Upp*Contra_deathAtBirth_N, 
               bi = (1-Contra_deathAtBirth_Upp)*Contra_deathAtBirth_N, 
               ci = noContra_deathAtBirth_Upp*noContra_deathAtBirth_N,
               di = (1-noContra_deathAtBirth_Upp)*noContra_deathAtBirth_N,
               var.names = c("yi_deathAtBirth_upp", "vi_deathAtBirth_upp"),
               data = dat_deathAtBirth)

dat_deathAtBirth %>% filter(Contra_deathAtBirth_N > 0) -> dat_deathAtBirth


# create a long format of the data using these 3 types of effect sizes (low, med, upp) yi and vi are the effect size and variance of the effect size

dat_long_deathAtBirth <- dat_deathAtBirth %>% select(Effect_ID, Species, Phylogeny, Sex_Type,
                           yi_deathAtBirth_low, vi_deathAtBirth_low, 
                           yi_deathAtBirth_med, vi_deathAtBirth_med, 
                           yi_deathAtBirth_upp, vi_deathAtBirth_upp) %>% 
  pivot_longer(cols = c(yi_deathAtBirth_low, yi_deathAtBirth_med, yi_deathAtBirth_upp, vi_deathAtBirth_low, vi_deathAtBirth_med, vi_deathAtBirth_upp), 
               names_to = c(".value", "type"), 
               names_pattern = "(yi|vi)_(.*)")

dat_long_deathAtBirth <- dat_long_deathAtBirth %>% filter(!is.na(yi))

str(dat_long_deathAtBirth)

dat_long_deathAtBirth$type <- factor(dat_long_deathAtBirth$type, 
                           levels = rev(c("deathAtBirth_low", "deathAtBirth_med", "deathAtBirth_upp")),
                           labels = rev(c("Lower", "Median", "Upper"))
                           )

# effect size level ID

dat_long_deathAtBirth$Effect_ID2 <- factor(1 : nrow(dat_long_deathAtBirth))

# VCV

VCV <- vcalc(dat_long_deathAtBirth$vi, 
            cluster = dat_long_deathAtBirth$Effect_ID, 
            obs = dat_long_deathAtBirth$Effect_ID2, 
            data = dat_long_deathAtBirth, rho = 0.5)

# meta-analysis using dat_long

mod_deathAtBirth <- rma.mv(yi = yi, V = VCV, 
                         random = list(
                           ~1|Species,
                           ~1|Phylogeny,
                           ~1|Effect_ID2), 
                         R = list(Phylogeny = cor_tree), 
                         data = dat_long_deathAtBirth,
                         method="REML", 
                         test = "t",
                         control=list(optimizer="optim", optmethod="Nelder-Mead")
)

summary(mod_deathAtBirth)

# meta-regression

mod_deathAtBirth_reg <- rma.mv(yi = yi, V = VCV, 
                             random = list(
                               ~1|Species,
                               ~1|Phylogeny,
                               ~1|Effect_ID2), 
                             #struct = "DIAG",
                             R = list(Phylogeny = cor_tree), 
                             data = dat_long_deathAtBirth,
                             method="REML", 
                             control=list(optimizer="optim", optmethod="Nelder-Mead"),
                             mods = ~type - 1,
                             test = "t",
                             sparse=FALSE
)

summary(mod_deathAtBirth_reg)

p_deathAtBirth <- orchard_plot(mod_deathAtBirth_reg, mod = "type",
    xlab = "log risk ratio (Death at Birth)", group = "Species") + ylim(-0.85, 0.7)

p_deathAtBirth

# # Sex_Type effect
# 
# mod_deathAtBirth_reg2 <- rma.mv(yi = yi, V = VCV, 
#                              random = list(
#                                ~1|Species,
#                                ~1|Phylogeny,
#                                ~1|Effect_ID2), 
#                              R = list(Phylogeny = cor_tree), 
#                              data = dat_long_deathAtBirth,
#                              method="REML", 
#                              control=list(optimizer="optim", optmethod="Nelder-Mead"),
#                              mods = ~Sex_Type - 1,
#                              test = "t",
#                              sparse=TRUE
# )
# 
# summary(mod_deathAtBirth_reg2)
# 
# p_deathAtBirth2 <- orchard_plot(mod_deathAtBirth_reg2, mod = "Sex_Type",
#     xlab = "log risk ratio (Death at Birth)", group = "Species", angle = 45) + ylim(-0.85, 0.7)
# 
# p_deathAtBirth2


# Other

# Low

dat_other <- escalc(measure = "RD", 
               ai = Contra_Other_Low*Contra_Other_N,
               bi = (1-Contra_Other_Low)*Contra_Other_N,
               ci = noContra_Other_Low*noContra_Other_N,
               di = (1-noContra_Other_Low)*noContra_Other_N,
               var.names = c("yi_other_low", "vi_other_low"),
               data = dat)

# Med

dat_other <- escalc(measure = "RD",
               ai = Contra_Other_Med*Contra_Other_N, 
               bi = (1-Contra_Other_Med)*Contra_Other_N, 
               ci = noContra_Other_Med*noContra_Other_N,
               di = (1-noContra_Other_Med)*noContra_Other_N,
               var.names = c("yi_other_med", "vi_other_med"),
               data = dat_other)

# Upp

dat_other <- escalc(measure = "RD",
               ai = Contra_Other_Upp*Contra_Other_N, 
               bi = (1-Contra_Other_Upp)*Contra_Other_N, 
               ci = noContra_Other_Upp*noContra_Other_N,
               di = (1-noContra_Other_Upp)*noContra_Other_N,
               var.names = c("yi_other_upp", "vi_other_upp"),
               data = dat_other)

dat_other %>% filter(Contra_Other_N > 0) -> dat_other

# create a long format of the data using these 3 types of effect sizes (low, med, upp) yi and vi are the effect size and variance of the effect size

dat_long_other <- dat_other %>% select(Effect_ID, Species, Phylogeny, Sex_Type,
                           yi_other_low, vi_other_low, 
                           yi_other_med, vi_other_med, 
                           yi_other_upp, vi_other_upp) %>% 
  pivot_longer(cols = c(yi_other_low, yi_other_med, yi_other_upp, vi_other_low, vi_other_med, vi_other_upp), 
               names_to = c(".value", "type"), 
               names_pattern = "(yi|vi)_(.*)")

dat_long_other <- dat_long_other %>% filter(!is.na(yi))

str(dat_long_other)

dat_long_other$type <- factor(dat_long_other$type, 
                           levels = rev(c("other_low", "other_med", "other_upp")),
                           labels = rev(c("Lower", "Median", "Upper"))
                           )
# effect size level ID  

dat_long_other$Effect_ID2 <- factor(1 : nrow(dat_long_other))

# VCV

VCV <- vcalc(dat_long_other$vi, 
            cluster = dat_long_other$Effect_ID, 
            obs = dat_long_other$Effect_ID2, 
            data = dat_long_other, rho = 0.5)


# meta-analysis using dat_long

mod_other <- rma.mv(yi = yi, V = VCV, 
                     random = list(
                       ~1|Species,
                       ~1|Phylogeny,
                       ~1|Effect_ID2), 
                     R = list(Phylogeny = cor_tree), 
                     data = dat_long_other,
                     method="REML", 
                     test = "t",
                     control=list(optimizer="optim", optmethod="Nelder-Mead")
)

summary(mod_other)

# meta-regression

mod_other_reg <- rma.mv(yi = yi, V = VCV, 
                         random = list(
                           ~1|Species,
                           ~1|Phylogeny,
                           ~1|Effect_ID2), 
                         #struct = "DIAG",
                         R = list(Phylogeny = cor_tree), 
                         data = dat_long_other,
                         method="REML", 
                         control=list(optimizer="optim", optmethod="Nelder-Mead"),
                         mods = ~type - 1,
                         test = "t",
                         sparse=FALSE
)

summary(mod_other_reg)

p_other <- orchard_plot(mod_other_reg, mod = "type",
    xlab = "log risk ratio (Other)", group = "Species") + ylim(-0.85, 0.7)

p_other

# Sex_Type effect

# filter out Male_Hormonal

dat_long_other2 <- dat_long_other %>% 
  filter(Sex_Type != "Male_Hormonal")

# VCV2

VCV2 <- vcalc(dat_long_other2$vi, 
             cluster = dat_long_other2$Effect_ID, 
             obs = dat_long_other2$Effect_ID2, 
             data = dat_long_other2, rho = 0.5)

mod_other_reg2 <- rma.mv(yi = yi, V = VCV2, 
                         random = list(
                           ~1|Species,
                           ~1|Phylogeny,
                           ~1|Effect_ID2), 
                         R = list(Phylogeny = cor_tree), 
                         data = dat_long_other2,
                         method="REML", 
                         control=list(optimizer="optim", optmethod="Nelder-Mead"),
                         mods = ~Sex_Type - 1,
                         test = "t",
                         sparse=TRUE
                         )

summary(mod_other_reg2)

p_other2 <- orchard_plot(mod_other_reg2, mod = "Sex_Type", 
    xlab = "log risk ratio (Other)", group = "Species", angle = 45) + ylim(-0.85, 0.7)

p_other2


# combining all plots - use cowplot

p_all <- plot_grid(p_trauma, 
                   p_infectious, 
                   p_noninfectious, 
                   p_chronic, 
                   p_deathAtBirth, 
                   p_other,
                   ncol = 2)

p_all

p_all2 <- plot_grid(p_trauma2, 
                   p_infectious2, 
                   p_noninfectious2, 
                   p_chronic2, 
                   #p_deathAtBirth, 
                   p_other2,
                   ncol = 2)

p_all2
