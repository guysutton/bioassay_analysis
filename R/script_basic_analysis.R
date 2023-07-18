# Script: Basic workthrough for bioassay analyses (proportion survival/mortality)

# -----------------------------------------------------------------------------
# Session setup
# -----------------------------------------------------------------------------

# Load required packages
if (!require("pacman"))
  install.packages("pacman")
pacman::p_load(
  tidyverse,
  raster,
  tidyr,
  DHARMa
)

# Change ggplot theme
theme_set(
  theme_classic() +
    theme(
      panel.border = element_rect(colour = "black",
                                  fill = NA),
      axis.text = element_text(colour = "black"),
      axis.title.x = element_text(margin = unit(c(2, 0, 0, 0),
                                                "mm")),
      axis.title.y = element_text(margin = unit(c(0, 4, 0, 0),
                                                "mm")),
      legend.position = "none"
    )
)

# Function to plot doses as X ^10x (still needs work!)
scientific_10 <- function(x) {
  parse(text=gsub("e\\+*", " %*% 10^", scales::scientific_format()(x)))
}

# Write function to calculate the lower and upper 95% confidence intervals for 
# LC50, LC90 and LC99 from a probit GLM
calc_lc_ci <- function(model){
  xp <- MASS::dose.p(model, p=c(0.50, 0.90, 0.99))  
  xp.ci <- xp + attr(xp, "SE") %*% matrix(qnorm(1 - 0.05/2)*c(-1,1), nrow=1)
  zp.est <- cbind(xp, attr(xp, "SE"), xp.ci[,1], xp.ci[,2])
  dimnames(zp.est)[[2]] <- c("LD", "SE", "LCL","UCL")
  zp.est
}

# -----------------------------------------------------------------------------
# Import data
# -----------------------------------------------------------------------------

# Import raw dataset 
df <- readxl::read_xlsx(
  here::here("./data/data_raw/bioassay.xlsx")
)

# Check data imported correctly
head(df)

# Check the columns are the correct format 
dplyr::glimpse(df)

# - Does it make sense to include 'dose' as a numeric or categorical predictor? 
#   - I am going to assume that dose should be a numeric... 

# -----------------------------------------------------------------------------
# Exploratory data analysis
# -----------------------------------------------------------------------------

# Plot relationship between dose and proportion mortality, and how it differs
# between insecticide treatments? 
# - You should be able to guess what the statistics say from this graph 
df %>%
  # Calculate a column containing proportion insects that died (2 dec places)
  dplyr::mutate(
    prop_dead = round(dead / total, 2)
  ) %>%
  # Calculate standard error 
  dplyr::group_by(insecticide, dose) %>%
  dplyr::summarise(
    mean_p = mean(prop_dead),
    sd_p = sd(prop_dead),
    n_p = n(), 
    se = sd_p / sqrt(n_p),
    ymin = mean_p - se,
    ymax = mean_p + se
  ) %>% 
  ggplot(
    data = ., 
    aes(
      x = dose, 
      y = mean_p,
      colour = insecticide
    )
  ) +
  geom_point() + 
  geom_line() +
  theme(
    legend.position = "right"
  )

# - I think that it is quite clear that the 'GV75-NPV25' treatment produces 
#   the highest mortality rate, higher than the other three treatments. 


# -----------------------------------------------------------------------------
# Model building
# -----------------------------------------------------------------------------

# In this section, we are going to fit a probit logistic regression model 
# whereby we compare the proportion of insects that died by treatment
# - This model will test the following hypotheses:
#   - 1. If there is a signficant difference in the proportion of insects that 
#        died between the four insecticide treatments (indicated by 'treatment'),
#        averaged across all the different doses.
#   - 2. If there is a signficant difference in the proportion of insects that 
#        died between the seven dose treatments (indicated by 'dose'), 
#        averaged across the four insecticide treatments.
#   - 3. Lastly, and most importantly, if there is a signficant interaction term
#        between 'treatment' and 'dose', which biologically translates into a 
#        hypothesis something like, "does the treatment effect vary between doses",
#        or "does the dose effect depend on the treatment"... 

# Fit model 
m1 <- glmmTMB::glmmTMB(
  cbind(dead, total - dead) ~ 
    insecticide + dose + insecticide:dose,
  data = df,
  family = binomial(link = "probit")
)

# Perform Likelihood Ratio Test 
car::Anova(
  mod = m1, 
  test = "Chisq",
  type = "III"
)

# - When there is a significant interaction term, the interpretation of the 
#   fixed effects (e.g. insecticide, dose), usually are meaningless.
# - So, there is support for a statistically significant interaction term 
#   between treatment and dose (X2 = 86.24, df = 3, P < 0.001). 
#   - This means that the effect of 'dose' on proportion mortality differed 
#     between the treatments... 


# -----------------------------------------------------------------------------
# Marginal predictions 
# -----------------------------------------------------------------------------

# Extract marginal predicted means from the fitted model 
preds <- ggeffects::ggeffect(
  model = m1, 
  terms = c("insecticide", "dose [0:1530000 by = 1000]"),
  type = "fe"
  ) %>%
  as.data.frame(.) %>%
  # Change column names
  dplyr::rename(
    insecticide = x,
    dose = group
  ) %>%
  # Convert dose into a numeric
  dplyr::mutate(
    dose = readr::parse_number(as.character(dose))
  )
head(preds)

# Make plot
p1 <- preds %>%
  ggplot(
    data = ., 
    aes(
      x = dose, 
      y = predicted, 
      colour = insecticide, 
      fill = insecticide
      )
    ) +
  geom_line() +
  scale_y_continuous(
    limits = c(0, 1),
    breaks = seq(0, 1, 0.2)
  ) +
 # scale_x_continuous(
 #   label = scientific_10
 # ) + 
  labs(
    x = "Dose",
    y = "Proportion mortality",
    colour = "Insecticide"
  ) +
  theme(
    legend.position = c(0.85, 0.3)
  )
p1

# Make plot (each treatment in separate panel for ci)
p2 <- preds %>%
  ggplot(
    data = ., 
    aes(
      x = dose, 
      y = predicted, 
      colour = insecticide, 
      fill = insecticide
    )
  ) +
  geom_ribbon(
    aes(
      ymin = conf.low,
      ymax = conf.high
    ),
    alpha = 0.2
  ) + 
  geom_line() +
  facet_wrap(~ insecticide, nrow = 1) + 
  scale_y_continuous(
    limits = c(0, 1),
    breaks = seq(0, 1, 0.2)
  ) +
  # scale_x_continuous(
  #   label = scientific_10
  # ) + 
  labs(
    x = "Dose",
    y = "Proportion mortality",
    fill = "Insecticide"
  ) 
p2

# Combine plots (just an idea)
plots <- cowplot::plot_grid(
  p1,
  p2,
  nrow = 2,
  labels = c("(A)", "(B)")
)
plots

# Save figure 
ggsave(
  here::here("./figures/fig_shared_legend.png"),
  dpi = 600,
  height = 6,
  width = 6
)

# -----------------------------------------------------------------------------
# Calculating LC values per treatment 
# -----------------------------------------------------------------------------

# Setup a loop to iterate through the different treatments, and fit 
# a probit model to each treatment subset, and extract LC values, as desired 

# Split data into different treatments 
by_treatment <- split(df, df$insecticide)

# Loop over each treatment and fit probit model and extract LC and associated
# 95% confidence intervals 
mod_list <- by_treatment %>% 
  # Iterate model fitting 
  purrr::map(
    ~ glm(
      cbind(dead, total - dead) ~ 
        dose,
      data = .x,
      family = binomial(link = "probit"))
  ) %>%
  # Calculate 95% confidence intervals around the LC values 
  purrr::map(
    ~ calc_lc_ci(.x)
  )
mod_list

# Convert results into a dataframe 
z <- mod_list %>% 
  tibble::enframe(.) %>%
  tidyr::unnest_longer(value) 
z

# Value_1 = LC estimate
# Value_2 = LC standard error
# Value_3 = LC lower 95% CI
# Value_4 = LC upper 95% CI

# -----------------------------------------------------------------------------
# Post-hoc comparisons 
# -----------------------------------------------------------------------------

# What doses would we want to compare? 
df %>%
  dplyr::distinct(dose)

# Post-hoc comparisons ( compare treatments at 1.6 x 10^4 )
emm <- emmeans::emmeans(
  m1,
  specs = pairwise ~ insecticide,
  regrid = "response",
  at = list(dose = c(16000)),
  adjust = "bonferroni"
)
emm

# Post-hoc comparisons ( compare treatments at 4.0 x 10^4 )
emm <- emmeans::emmeans(
  m1,
  specs = pairwise ~ insecticide,
  regrid = "response",
  at = list(dose = c(40000)),
  adjust = "bonferroni"
)
emm

# Post-hoc comparisons ( compare treatments at 1.0 x 10^5 )
emm <- emmeans::emmeans(
  m1,
  specs = pairwise ~ insecticide,
  regrid = "response",
  at = list(dose = c(100000)),
  adjust = "bonferroni"
)
emm

# Post-hoc comparisons ( compare treatments at 2.5 x 10^5 )
emm <- emmeans::emmeans(
  m1,
  specs = pairwise ~ insecticide,
  regrid = "response",
  at = list(dose = c(250000)),
  adjust = "bonferroni"
)
emm

# Post-hoc comparisons ( compare treatments at 6.25 x 10^5 )
emm <- emmeans::emmeans(
  m1,
  specs = pairwise ~ insecticide,
  regrid = "response",
  at = list(dose = c(625000)),
  adjust = "bonferroni"
)
emm

# Post-hoc comparisons ( compare treatments at 1.5 x 10^6 )
emm <- emmeans::emmeans(
  m1,
  specs = pairwise ~ insecticide,
  regrid = "response",
  at = list(dose = c(1500000)),
  adjust = "bonferroni"
)
emm






















