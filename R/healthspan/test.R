# test - health span data

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

# flipping directions

dat$direction <- ifelse(dat$direction.of.improved.health == "Increased", 1, -1)

dat$yi <- dat$yi * dat$direction


# meta-analysis

VCV <- vcalc(vi = vi, cluster = Study, obs = Effect_ID, #subgroup = Measure,
             data = dat, rho = 0.5)

mod <-  rma.mv(yi = yi, 
               V = VCV, 
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


# male and female difference

# VCV1 <- vcalc(vi = vi, cluster = Study, obs = Effect_ID, subgroup = Sex, data = dat, rho = 0.5)
# 
# mod1b <-  rma.mv(yi = yi, 
#                 V = VCV1, 
#                 mod = ~ Sex - 1, 
#                 random = list(~ 1|Strain, 
#                               ~ 1|Study, 
#                               ~ Sex|Effect_ID), 
#                 struct = "DIAG",
#                 data = dat, 
#                 test = "t",
#                 sparse = TRUE,
#                 control=list(optimizer="optim", optmethod="BFGS")
# )
# summary(mod1b) 


mod1 <-  rma.mv(yi = yi, 
                V = VCV, 
                mod = ~ Sex - 1, 
                random = list(~ 1|Strain, 
                              ~ 1|Study, 
                              ~ 1|Effect_ID), 
                #struct = "DIAG",
                data = dat, 
                test = "t",
                sparse = TRUE,
                control=list(optimizer="optim", optmethod="BFGS")
)
summary(mod1) 

# visualizing the result
orchard_plot(mod1, mod = "Sex", 
             xlab = "log response ratio (lnRR)", group = "Study", angle = 45, cb = F) # colour = T)


# attr(lm_result, "class") <- NULL

# orchard plot

main <- mod_results(mod, group = "Study")
attr(main, "class") <- NULL
main$mod_table$name <- gsub("Intrcpt", "Overall", main$mod_table$name)
main$mod_table$name <- factor(main$mod_table$name)
main$data$moderator <- gsub("Intrcpt", "Overall", main$data$moderator)
main$data$moderator <- factor(main$data$moderator)

class(main) <- c("orchard", "data.frame")

p_overall <- orchard_plot(main, xlab = "log response ratio (lnRR)", angle = 0, group = "Study") + ylim(-2.2,2)


sex_diff <- mod_results(mod1, mod = "Sex", group = "Study")


p_sex_diff <-orchard_plot(sex_diff, mod = "Sex", 
                          xlab = "log response ratio (lnRR)", group = "Study", angle = 0, cb = F) + ylim(-2.2,2)

combined <- submerge(sex_diff, main)

# changing the name (intercept) to:
# TODO - this should be done in the orchard_plot function
combined$mod_table$name <- gsub("Intrcpt", "Overall", combined$mod_table$name)
combined$mod_table$name <- factor(combined$mod_table$name)

combined$data$moderator <- gsub("Intrcpt", "Overall", combined$data$moderator)


orchard_plot(combined,
             xlab = "log response ratio (lnRR)", group = "Study", angle = 0, cb = F) + 
  scale_colour_manual(values = rev(c("#999999", "#88CCEE", "#CC6677"))) +
  scale_fill_manual(values = rev(c("#999999", "#88CCEE", "#CC6677")))




# mod1 better
#anova(mod1, mod1b)

# measure 

#VCV2 <- vcalc(vi = vi, cluster = Study, obs = Effect_ID, subgroup = Sub.measure,
#            data = dat, rho = 0.5)

# mod2b <-  rma.mv(yi = yi, 
#                V = VCV2, 
#                mod = ~ Sub.measure - 1, 
#                random = list(~ 1|Strain, 
#                              ~ 1|Study, 
#                              ~ Sub.measure|Effect_ID), 
#                struct = "DIAG",
#                data = dat, 
#                test = "t",
#                sparse = TRUE,
#                control=list(optimizer="optim", optmethod="BFGS")
# )
# summary(mod2) 

mod2 <-  rma.mv(yi = yi, 
                V = VCV, 
                mod = ~ Sub.measure - 1, 
                random = list(~ 1|Strain, 
                              ~ 1|Study, 
                              ~ 1|Effect_ID), 
                #struct = "DIAG",
                data = dat, 
                test = "t",
                sparse = TRUE,
                control=list(optimizer="optim", optmethod="BFGS")
)
summary(mod2) 


# visualizing the result
orchard_plot(mod2, mod = "Sub.measure", 
             xlab = "log response ratio (lnRR)", group = "Study", angle = 45)



# Measure x Sex
# need to clearn up Sub.measure

dat$Sub.measure <- factor(dat$Sub.measure, levels = levels(dat$Sub.measure),
                          labels = c("Cardiac\nfunction/\npathology", "Cardiac size", "Cognition", "Frailty",
                                     "Immune\nfunction", "Metabolism", "Muscle size", "Non-tumor\npathology",
                                     "Sensory\nfunction", "Strength/\nbalance", "Tumor\nmammory", "Tumor\nnonmammory",
                                     "Voluntary\nactivity"))


dat$MesSex <- paste0(dat$Sub.measure, "_", dat$Sex)

#VCV <- vcalc(vi = vi, cluster = Study, obs = Effect_ID, subgroup = MesSex,
#             data = dat, rho = 0.5)

# mod3b <-  rma.mv(yi = yi, 
#                V = vi, 
#                mod = ~ MesSex - 1, 
#                random = list(~ 1|Strain, 
#                              ~ 1|Study, 
#                              ~ MesSex|Effect_ID), 
#                struct = "DIAG",
#                data = dat, 
#                test = "t",
#                sparse = TRUE,
#                control=list(optimizer="optim", optmethod="BFGS")
# )
# 
# summary(mod3b)


mod3 <-  rma.mv(yi = yi, 
                 V = vi, 
                 mod = ~ MesSex - 1, 
                 random = list(~ 1|Strain, 
                               ~ 1|Study, 
                               ~ 1|Effect_ID), 
                 #struct = "DIAG",
                 data = dat, 
                 test = "t",
                 sparse = TRUE,
                 control=list(optimizer="optim", optmethod="BFGS")
)

summary(mod3)


# visualizing the result

orchard_plot(mod3, mod = "MesSex", 
             xlab = "log response ratio (lnRR)", group = "Study", angle = 0, cb = F)


# taking out 



# some experiments 

res3 <- mod_results(mod3, mod = "MesSex", group = "Study")

mod_table_m <- res3$mod_table[c(2,4,6, 7, 10, 12, 14, 16, 18, 21, 23), ]
mod_table_f <- res3$mod_table[c(1,3,5, 8, 9, 11, 13, 15, 17, 19, 20, 22), ]

data_m <- res3$data[ !is.na(match(res3$data$moderator, mod_table_m$name)), ]
data_f <- res3$data[ !is.na(match(res3$data$moderator, mod_table_f$name)), ]
  
res3_male <- list(mod_table = mod_table_m, data = data_m)
res3_female <- list(mod_table = mod_table_f, data = data_f)

# changing names
# male
res3_male$mod_table$name <- gsub("_Male", "", res3_male$mod_table$name)
res3_male$mod_table$name <- factor(res3_male$mod_table$name)

res3_male$data$moderator <- gsub("_Male", "", res3_male$data$moderator)

# female
res3_female$mod_table$name <- gsub("_Female", "", res3_female$mod_table$name)
res3_female$mod_table$name <- factor(res3_female$mod_table$name)

res3_female$data$moderator <- gsub("_Female", "", res3_female$data$moderator)

class(res3_male) <- c("orchard", "data.frame")
class(res3_female) <- c("orchard", "data.frame")





p_male <- orchard_plot(res3_male, mod = "MesSex", 
             xlab = "log response ratio (lnRR)", group = "Study", angle = 0, cb = F) + labs(title = "Male") + ylim(-2.2,2)


# p_male <- p_male + plot_annotation(title = "Male",
#                          theme = theme(plot.title = element_text(size = 20)))

p_female <- orchard_plot(res3_female, mod = "MesSex", 
             xlab = "log response ratio (lnRR)", group = "Study", angle = 0, cb = F) + labs(title = "Female") + ylim(-2.2,2)

# p_female <- p_female + plot_annotation(title = "Female",
#                          theme = theme(plot.title = element_text(size = 20)))

p_overall + p_sex_diff + p_female + p_male +  plot_layout(widths = c(2, 2), heights = c(1,2.5)) 


# creating a forest plot using dat

str(dat)

# selecting

fdat <- dat %>% select(Study, Measure, Effect_ID, Measurement.type, Sub.measure, yi, vi)

# add lower.ci and upper.ci

fdat <- dat %>% mutate(lower.ci = yi - 1.96*sqrt(vi), upper.ci = yi + 1.96*sqrt(vi))

# p_forest <- ggplot(data = fdat, aes(x = yi, y = Measurement.type)) +
#   geom_errorbarh(aes(xmin = lower.ci, xmax = upper.ci, colour = Sub.measure), 
#                  height = 1, show.legend = TRUE, size = 1, alpha = 0.8, position = position_dodge2(width = 1)) +
#   geom_point(aes(col = Sub.measure), fill = "white", size = 2, shape = 21, position =position_dodge2(width = 1)) +
#   geom_vline(xintercept = 0, linetype = 2, colour = "black", alpha = 0.3) +
#   geom_vline(xintercept = mod$b, linetype = 2, colour = "red", alpha = 0.3) +
#   labs(x = "lnRR (effect size)", y = "")

fdat$Measurement.type <- factor(
  fdat$Measurement.type,
  levels = unique(fdat$Measurement.type)
)

dodge <- position_dodge2(width = 0.6, preserve = "single")

p_forest <- ggplot(data = fdat, aes(x = yi, y = Measurement.type)) +
  geom_errorbarh(
    aes(xmin = lower.ci, xmax = upper.ci, colour = Sub.measure), 
    height = 0.3,
    position = dodge,
    size = 1, alpha = 0.8
  ) +
  geom_point(
    aes(colour = Sub.measure),
    fill = "white",
    shape = 21,
    size = 2,
    position = dodge
  ) +
  geom_vline(xintercept = 0, linetype = 2, colour = "black", alpha = 0.3) +
  geom_vline(xintercept = mod$b, linetype = 1, linewidth = 2, colour = "red", alpha = 0.3) +
  scale_y_discrete(expand = c(0,0)) +
  # Rename legend title
  scale_colour_discrete(name = "Sub-measure") +
  labs(x = "lnRR (effect size)", y = "") +
  theme_minimal() +
  # Place legend below plot
  theme(
    legend.position = "bottom",
    axis.text.y = element_text(size = 9)
  )

p_forest

########################
# different drawing .... 
########################

# example....






