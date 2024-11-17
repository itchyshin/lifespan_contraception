# test

# loading
pacman::p_load(tidyverse,
               metafor,
               pander,
               stringr,
               ape,
               kableExtra,
               patchwork,
               lme4,
               readxl,
               emmeans,
               rotl,
               orchaRd,
               clubSandwich,
               MuMIn,
               png,
               grid,
               here,
               formatR,
               naniar,
               GoodmanKruskal,
               ggalluvial
)

# need for metafor to understand MuMin 
eval(metafor:::.MuMIn)

# loading extra functions
# function for getting lnRR for proportional data (mortality)

lnrrp <- function(m1, m2, n1, n2) {
  # if p1 or p2 = 0, turn into 0.025, if p1 or p2 = 1, turn into 0.975
    m1[m1 == 0] <- 0.025
    m1[m1 == 1] <- 0.975
    m2[m2 == 0] <- 0.025
    m2[m2 == 1] <- 0.975
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


# function to get to lnRR for longevity data (CV required)
# The method proposed in Nakagawa et al (2022) - missing SD method

lnrrm <- function(m1, m2, n1, n2, cv21, cv22) {
  # lnRR - with 2nd order correction
  lnrr <- log(m1/m2) + 
    0.5 * ((cv21 /n1) - (cv22 / n2))	
  
  var <- (cv21 / n1) + ((cv21^2) / (2 * n1^2))  + 
    (cv22/ n2) + ((cv22^2) / (2 * n2^2) )
  
  invisible(data.frame(yi = lnrr , vi = var))
}

# for folded normal distribution see: https://en.wikipedia.org/wiki/Folded_normal_distribution

# folded mean
folded_mu <-function(mean, variance){
  mu <- mean
  sigma <- sqrt(variance)
  fold_mu <- sigma*sqrt(2/pi)*exp((-mu^2)/(2*sigma^2)) + mu*(1 - 2*pnorm(-mu/sigma))
  fold_mu
} 

# folded variance
folded_v <-function(mean, variance){
  mu <- mean
  sigma <- sqrt(variance)
  fold_mu <- sigma*sqrt(2/pi)*exp((-mu^2)/(2*sigma^2)) + mu*(1 - 2*pnorm(-mu/sigma))
  fold_se <- sqrt(mu^2 + sigma^2 - fold_mu^2)
  # adding se to make bigger mean
  fold_v <-fold_se^2
  fold_v
} 


# turning survival % into life expectancy
# note that if m1 = treatment lives 100% (happens twice in the data then we use the height observed in the data-set)

lnrre <- function(m1, m2, n1, n2) {
  
  # expected life span transformation
  
  ea1 <- -1/log(m1)
  ea2 <- -1/log(m2)    
  # var using the delta method
  # V is not var like above - it is not sampling SD but error
  V1 <- ((1-m1)*m1/n1)*(1/(m1*log(m1)^2))^2
  V2 <- ((1-m2)*m2/n2)*(1/(m2*log(m2)^2))^2
  
  lnrr <- log(ea1/ea2)
  
  var <- V1 + V2
  
  
  invisible(data.frame(yi = lnrr , vi = var))
}

#TODO  - we need do data checking and cleanring


# loading data

dat_full <- read.csv(here("data", "healthspan", "healthspan.csv"), na = c("", "NA"))

dat_full %>% 
  mutate_if(is.character, as.factor) -> dat

# Effect_ID is the unique identifier for the effect
dat$Effect_ID <- factor(1:nrow(dat))

effect_type <- ifelse(str_detect(dat$measurement_parameter, "rop"), "proportion", "other")



dat %>% group_by(Study) %>% 
  summarise(cv2_cont = mean((Error_control_SD/Control_value)^2, na.rm = T), 
            cv2_trt = mean((Error_experimental_SD/Experimental_value)^2, na.rm = T), 
            n_cont = mean(Sample_size_control, na.rm = T), 
            n_trt =  mean(Sample_size_experimental, na.rm = T)) %>% 
  ungroup() %>% 
  summarise(cv2_cont = weighted.mean(cv2_cont, n_cont, na.rm = T), 
            cv2_trt = weighted.mean(cv2_trt, n_trt, na.rm = T)) -> cvs

# lnRR

# survival proportion using arcsin transformation = lnrrp
# using CV
dat$yi <- ifelse(effect_type == "other", lnrrm(dat$Experimental_value,
                                               dat$Control_value, 
                                               dat$Sample_size_experimental, 
                                               dat$Sample_size_control, 
                                               cvs[["cv2_trt"]],
                                               cvs[["cv2_cont"]])[[1]],
                                         lnrrp(dat$Experimental_value, 
                                               dat$Control_value, 
                                               dat$Sample_size_experimental, 
                                               dat$Sample_size_control)[[1]])

dat$vi <- ifelse(effect_type == "other", lnrrm(dat$Experimental_value,
                                               dat$Control_value, 
                                               dat$Sample_size_experimental, 
                                               dat$Sample_size_control, 
                                               cvs[["cv2_trt"]],
                                               cvs[["cv2_cont"]])[[2]],
                                         lnrrp(dat$Experimental_value, 
                                               dat$Control_value, 
                                               dat$Sample_size_experimental, 
                                               dat$Sample_size_control)[[2]])

# meta-analysis

VCV <- vcalc(vi = vi, cluster = Study, obs = Effect_ID, #subgroup = Measure,
             data = dat, rho = 0.5)

mod <-  rma.mv(yi = yi, 
               V = vi, 
               #mod = ~ 1, 
               random = list(~1|Strain, 
                             ~ 1|Study, 
                             ~ 1|Effect_ID), 
               #struct = "DIAG",
               data = dat, 
               test = "t",
               sparse = TRUE,
               control=list(optimizer="optim", optmethod="Nelder-Mead")
)
summary(mod) 

round(i2_ml(mod),2) # almost no phylogenetic effect

# visualizing the result
orchard_plot(mod, xlab = "log response ratio (lnRR)", group = "Study")


# measure 

VCV <- vcalc(vi = vi, cluster = Study, obs = Effect_ID, subgroup = Measure,
             data = dat, rho = 0.5)

mod1 <-  rma.mv(yi = yi, 
               V = vi, 
               mod = ~ Measure - 1, 
               random = list(~ 1|Strain, 
                             ~ 1|Study, 
                             ~ Measure|Effect_ID), 
               struct = "DIAG",
               data = dat, 
               test = "t",
               sparse = TRUE,
               control=list(optimizer="optim", optmethod="BFGS")
)
summary(mod1) 

# visualizing the result
orchard_plot(mod1, mod = "Measure", 
             xlab = "log response ratio (lnRR)", group = "Study", angle = 45)

# Measure x Sex

dat$MesSex <- paste0(dat$Measure, "_", dat$Sex)

VCV <- vcalc(vi = vi, cluster = Study, obs = Effect_ID, subgroup = MesSex,
             data = dat, rho = 0.5)

mod2 <-  rma.mv(yi = yi, 
               V = vi, 
               mod = ~ MesSex - 1, 
               random = list(~ 1|Strain, 
                             ~ 1|Study, 
                             ~ MesSex|Effect_ID), 
               struct = "DIAG",
               data = dat, 
               test = "t",
               sparse = TRUE,
               control=list(optimizer="optim", optmethod="BFGS")
)

summary(mod2)

# visualizing the result

orchard_plot(mod2, mod = "MesSex", 
             xlab = "log response ratio (lnRR)", group = "Study", angle = 45)


