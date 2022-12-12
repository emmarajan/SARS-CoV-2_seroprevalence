library(ggplot2)
library(tidyverse)

lc_logistic_regression <- read.csv(file = 'lc_logistic_regression.csv')

ggplot(lc_logistic_regression, aes(y = estimate, reorder(x = term, estimate))) +
  geom_pointrange(aes(ymin = conf.low, ymax = conf.high),
                  size = 0.5) +
  geom_hline(yintercept = 1.0, linetype = "dotted", size = 1) +
  scale_y_log10(breaks = c(0.01, 0.1, 1.0, 10, 100, 1000)) +
  labs(y = "Odds ratio (95% CI)", x = "post COVID-19 symptoms") +
  coord_flip(ylim = c(0.01, 1000)) +
  theme_bw()

