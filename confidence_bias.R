library(tidyverse)
library(hmetad)
library(tidybayes)
library(paletteer)
library(patchwork)

PALETTE <- paletteer_d("fishualize::Phractocephalus_hemioliopterus")[c(1, 3)]

## Number of simulated datasets and trials/dataset
N_trials <- 1000000

## fit sample models
m.absolute <- fit_metad(bf(N ~ 0 + Intercept),
                        data=sim_metad(N=N_trials, c2_0_diff=.5, c2_1_diff=.5),
                        backend='cmdstanr',
                        algorithm='laplace',
                        prior=prior(normal(0, 1)) +
                          prior(normal(0, 1), class=dprime) +
                          prior(normal(0, 1), class=c) +
                          prior(lognormal(0, 1), class=metac2zero1diff) +
                          prior(lognormal(0, 1), class=metac2one1diff))
m.absolute

m.relative <- fit_metad(bf(N ~ 0 + Intercept),
                        data=sim_metad(N=N_trials, c2_0_diff=.5, c2_1_diff=.5),
                        backend='cmdstanr', metac_absolute=FALSE,
                algorithm='laplace',
                prior=prior(normal(0, 1)) +
                  prior(normal(0, 1), class=dprime) +
                  prior(normal(0, 1), class=c) +
                  prior(lognormal(0, 1), class=metac2zero1diff) +
                  prior(lognormal(0, 1), class=metac2one1diff))
m.relative

#' Fit a model repeatedly until convergence
#' @param model the model to fit
#' @param newdata the new data to fit the model to
#' @param times the maximum number of times to attempt model fitting
fitUntilConvergence <- function(model, newdata, times=10) {
  fit <- tryCatch(update(model, newdata=newdata, init=0),
                  error=function(e) NULL)
  i <- 1

  while (is.null(fit) && i < times) {
    fit <- tryCatch(update(model, newdata=newdata, init=2),
                    error=function(e) NULL)
    i <- i + 1
  }

  fit
}

#' For a range of parameters, fit models to simulated data
#' and calculate measures of metacognitive bias
#' @param parameters a tibble of model parameters for data simulation
#' @param metac_absolute if TRUE, use absolute parameterization in data simulation (metac=c).
#' if FALSE, use relative parameterization (metac = M*c)
#' @param model model to re-fit to simulated data
simulate_conf_bias <- function(parameters, metac_absolute=TRUE, model=m.absolute) {
  parameters |>
    mutate(data=pmap(list(N_trials, dprime, c, log_M, c2_0_diff, c2_1_diff),
                     function(...) aggregate_metad(sim_metad(..., metac_absolute=metac_absolute))),
           model=map(data, partial(fitUntilConvergence, model=model)),
           mean_conf0=map_dbl(model,
                              ~ tibble(.row=1) |>
                                add_mean_confidence_draws(., by_stimulus=FALSE) |>
                                filter(response==0) |>
                                median_qi(.epred) |>
                                pull(.epred)),
           mean_conf1=map_dbl(model,
                              ~ tibble(.row=1) |>
                                add_mean_confidence_draws(., by_stimulus=FALSE) |>
                                filter(response==1) |>
                                median_qi(.epred) |>
                                pull(.epred)),
           mdistance0=map_dbl(model,
                              ~ tibble(.row=1) |>
                                add_linpred_draws(., dpar=TRUE, transform=TRUE) |>
                                mutate(mdistance=rowMeans(pick(starts_with('metac2zero'))) /
                                         (mu*dprime)) |>
                                median_qi(mdistance) |>
                                pull(mdistance)),
           mdistance1=map_dbl(model,
                              ~ tibble(.row=1) |>
                                add_linpred_draws(., dpar=TRUE, transform=TRUE) |>
                                mutate(mdistance=rowMeans(pick(starts_with('metac2one'))) /
                                         (mu*dprime)) |>
                                median_qi(mdistance) |>
                                pull(mdistance)),
           metadelta0=map_dbl(model,
                              ~ tibble(.row=1) |>
                                add_metacognitive_bias_draws(.) |>
                                filter(response==0) |>
                                median_qi(metacognitive_bias) |>
                                pull(metacognitive_bias)),
           metadelta1=map_dbl(model,
                              ~ tibble(.row=1) |>
                                add_metacognitive_bias_draws(.) |>
                                filter(response==1) |>
                                median_qi(metacognitive_bias) |>
                                pull(metacognitive_bias)))
}


## varying type 2 criteria
fits.c2_0 <- expand_grid(N=N_trials, dprime=2, c=.5, log_M=log(.5),
                         c2_0_diff=seq(.1, 3, by=.1), c2_1_diff=.5) |>
  simulate_conf_bias(metac_absolute=TRUE, model=m.absolute)

p.c2_0.conf <- fits.c2_0 |>
  pivot_longer(mean_conf0:mean_conf1, names_to='response') |>
  mutate(response=str_match(response, '[[:digit:]]')) |>
  ggplot(aes(x=c2_0_diff, y=value, color=response)) +
  geom_line() +
  scale_y_continuous('Mean Confidence', expand=expansion()) +
  scale_x_continuous(expression('dmeta-c'[2]^0), expand=expansion()) +
  coord_cartesian(xlim=c(0, 3), ylim=c(1, 2)) +
  theme_bw(18)
p.c2_0.conf

p.c2_0.mdistance <- fits.c2_0 |>
  pivot_longer(mdistance0:mdistance1, names_to='response') |>
  mutate(response=str_match(response, '[[:digit:]]')) |>
  ggplot(aes(x=c2_0_diff, y=value, color=response)) +
  geom_line() +
  scale_y_continuous('mdistance', expand=expansion()) +
  scale_x_continuous(expression('dmeta-c'[2]^0), expand=expansion()) +
  coord_cartesian(xlim=c(0, 3), ylim=c(0, 3)) +
  theme_bw(18)
p.c2_0.mdistance

p.c2_0.metadelta <- fits.c2_0 |>
  pivot_longer(metadelta0:metadelta1, names_to='response') |>
  mutate(response=str_match(response, '[[:digit:]]')) |>
  ggplot(aes(x=c2_0_diff, y=value, color=response)) +
  geom_line() +
  scale_y_continuous(expression('meta-'*Delta), expand=expansion()) +
  scale_x_continuous(expression('dmeta-c'[2]^0), expand=expansion()) +
  coord_cartesian(xlim=c(0, 3), ylim=c(0, 3)) +
  theme_bw(18)
p.c2_0.metadelta


fits.c2_1 <- expand_grid(N=N_trials, dprime=2, c=.5, log_M=log(.5),
                         c2_0_diff=.5, c2_1_diff=seq(.1, 3, by=.1)) |>
  simulate_conf_bias(metac_absolute=TRUE, model=m.absolute)

p.c2_1.conf <- fits.c2_1 |>
  pivot_longer(mean_conf0:mean_conf1, names_to='response') |>
  mutate(response=str_match(response, '[[:digit:]]')) |>
  ggplot(aes(x=c2_1_diff, y=value, color=response)) +
  geom_line() +
  scale_y_continuous('Mean Confidence', expand=expansion()) +
  scale_x_continuous(expression('dmeta-c'[2]^1), expand=expansion()) +
  coord_cartesian(xlim=c(0, 3), ylim=c(1, 2)) +
  theme_bw(18)
p.c2_1.conf

p.c2_1.mdistance <- fits.c2_1 |>
  pivot_longer(mdistance0:mdistance1, names_to='response') |>
  mutate(response=str_match(response, '[[:digit:]]')) |>
  ggplot(aes(x=c2_1_diff, y=value, color=response)) +
  geom_line() +
  scale_y_continuous('mdistance', expand=expansion()) +
  scale_x_continuous(expression('dmeta-c'[2]^1), expand=expansion()) +
  coord_cartesian(xlim=c(0, 3), ylim=c(0, 3)) +
  theme_bw(18)
p.c2_1.mdistance

p.c2_1.metadelta <- fits.c2_1 |>
  pivot_longer(metadelta0:metadelta1, names_to='response') |>
  mutate(response=str_match(response, '[[:digit:]]')) |>
  ggplot(aes(x=c2_1_diff, y=value, color=response)) +
  geom_line() +
  scale_y_continuous(expression('meta-'*Delta), expand=expansion()) +
  scale_x_continuous(expression('dmeta-c'[2]^1), expand=expansion()) +
  coord_cartesian(xlim=c(0, 3), ylim=c(0, 3)) +
  theme_bw(18)
p.c2_1.metadelta


((p.c2_0.conf | p.c2_1.conf ) /
   (p.c2_0.metadelta | p.c2_1.metadelta)) +
  plot_layout(guides='collect') +
  plot_annotation(tag_levels='A') &
  scale_color_manual('Response', values=PALETTE)
ggsave('plots/estimation_sensitivity.pdf', width=8, height=6)




## varying type 1 criterion
fits.c <- expand_grid(N=N_trials, dprime=2, c=seq(-1.5, 1.5, by=.1), log_M=log(.5),
                      c2_0_diff=.5, c2_1_diff=.5) |>
  simulate_conf_bias(metac_absolute=TRUE, model=m.absolute)

p.c.conf <- fits.c |>
  pivot_longer(mean_conf0:mean_conf1, names_to='response') |>
  mutate(response=str_match(response, '[[:digit:]]')) |>
  ggplot(aes(x=c, y=value, color=response)) +
  geom_line() +
  scale_y_continuous('Mean Confidence', expand=expansion()) +
  scale_x_continuous(breaks=-1:1, expand=expansion()) +
  coord_cartesian(xlim=c(-1.5, 1.5), ylim=c(1, 2)) +
  theme_bw(18)
p.c.conf

p.c.mdistance <- fits.c |>
  pivot_longer(mdistance0:mdistance1, names_to='response') |>
  mutate(response=str_match(response, '[[:digit:]]')) |>
  ggplot(aes(x=c, y=value, color=response)) +
  geom_line() +
  scale_y_continuous('mdistance', expand=expansion()) +
  scale_x_continuous(breaks=-1:1,expand=expansion()) +
  coord_cartesian(xlim=c(-1.5, 1.5), ylim=c(0, 1)) +
  theme_bw(18)
p.c.mdistance

p.c.metadelta <- fits.c |>
  pivot_longer(metadelta0:metadelta1, names_to='response') |>
  mutate(response=str_match(response, '[[:digit:]]')) |>
  ggplot(aes(x=c, y=value, color=response)) +
  geom_line() +
  scale_y_continuous(expression('meta-'*Delta), expand=expansion()) +
  scale_x_continuous(breaks=-1:1,expand=expansion()) +
  coord_cartesian(xlim=c(-1.5, 1.5), ylim=c(0, 1)) +
  theme_bw(18)
p.c.metadelta






## varying type 1 sensitivity (fixing meta-dprime)
fits.dprime <- expand_grid(N=N_trials,
                           dprime=seq(.5, 3, by=.1), c=0.5,
                           meta_dprime=1,
                           c2_0_diff=.5, c2_1_diff=.5) |>
  mutate(log_M=log(meta_dprime/dprime)) |>
  simulate_conf_bias(metac_absolute=TRUE, model=m.absolute)

p.dprime.conf <- fits.dprime |>
  pivot_longer(starts_with('mean_conf'), names_to='response') |>
  mutate(response=str_match(response, '[[:digit:]]')) |>
  ggplot(aes(x=dprime, y=value, color=response)) +
  geom_line() +
  scale_y_continuous('Mean Confidence', expand=expansion()) +
  scale_x_continuous("d'", expand=expansion()) +
  coord_cartesian(xlim=c(0, 3), ylim=c(1, 2)) +
  theme_bw(18)
p.dprime.conf

p.dprime.mdistance <- fits.dprime |>
  pivot_longer(mdistance0:mdistance1, names_to='response') |>
  mutate(response=str_match(response, '[[:digit:]]')) |>
  ggplot(aes(x=dprime, y=value, color=response)) +
  geom_line() +
  scale_y_continuous('mdistance', expand=expansion()) +
  scale_x_continuous("d'", expand=expansion()) +
  coord_cartesian(xlim=c(0, 3), ylim=c(0, 1)) +
  theme_bw(18)
p.dprime.mdistance

p.dprime.metadelta <- fits.dprime |>
  pivot_longer(metadelta0:metadelta1, names_to='response') |>
  mutate(response=str_match(response, '[[:digit:]]')) |>
  ggplot(aes(x=dprime, y=value, color=response)) +
  geom_line() +
  scale_y_continuous(expression('meta-'*Delta), expand=expansion()) +
  scale_x_continuous("d'", expand=expansion()) +
  coord_cartesian(xlim=c(0, 3), ylim=c(0, 1)) +
  theme_bw(18)
p.dprime.metadelta



## varying type 2 sensitivity
fits.meta_dprime <- expand_grid(N=N_trials, dprime=2, c=.5,
                                meta_dprime=seq(.5, 3, by=.1),
                                c2_0_diff=.5, c2_1_diff=.5) |>
  mutate(log_M=log(meta_dprime/dprime)) |>
  simulate_conf_bias(metac_absolute=TRUE, model=m.absolute)

p.meta_dprime.conf <- fits.meta_dprime |>
  pivot_longer(starts_with('mean_conf'), names_to='response') |>
  mutate(response=str_match(response, '[[:digit:]]')) |>
  ggplot(aes(x=meta_dprime, y=value, color=response)) +
  geom_line() +
  scale_y_continuous('Mean Confidence', expand=expansion()) +
  scale_x_continuous("meta-d'", expand=expansion()) +
  coord_cartesian(xlim=c(0, 3), ylim=c(1, 2)) +
  theme_bw(18)
p.meta_dprime.conf

p.meta_dprime.mdistance <- fits.meta_dprime |>
  pivot_longer(mdistance0:mdistance1, names_to='response') |>
  mutate(response=str_match(response, '[[:digit:]]')) |>
  ggplot(aes(x=meta_dprime, y=value, color=response)) +
  geom_line() +
  scale_y_continuous('mdistance', expand=expansion()) +
  scale_x_continuous("meta-d'",  expand=expansion()) +
  coord_cartesian(xlim=c(0, 3), ylim=c(0, 1)) +
  theme_bw(18)
p.meta_dprime.mdistance

p.meta_dprime.metadelta <- fits.meta_dprime |>
  pivot_longer(metadelta0:metadelta1, names_to='response') |>
  mutate(response=str_match(response, '[[:digit:]]')) |>
  ggplot(aes(x=meta_dprime, y=value, color=response)) +
  geom_line() +
  scale_y_continuous(expression('meta-'*Delta),  expand=expansion()) +
  scale_x_continuous("meta-d'", expand=expansion()) +
  coord_cartesian(xlim=c(0, 3), ylim=c(0, 1)) +
  theme_bw(18)
p.meta_dprime.metadelta



((p.c.conf | p.dprime.conf | p.meta_dprime.conf) /
   (p.c.metadelta | p.dprime.metadelta | p.meta_dprime.metadelta)) +
  plot_layout(guides='collect') +
  plot_annotation(tag_levels='A') &
  scale_color_manual('Response', values=PALETTE)
ggsave('plots/estimation.pdf', width=10, height=6)








#################################################################################################
#                                  Fit under incorrect generative model
#################################################################################################
## varying type 1 criterion
fits.c <- expand_grid(N=N_trials, dprime=2, c=seq(-1.5, 1.5, by=.1), log_M=log(.5),
                      c2_0_diff=.5, c2_1_diff=.5) |>
  simulate_conf_bias(metac_absolute=FALSE, model=m.absolute)

p.c.conf <- fits.c |>
  pivot_longer(mean_conf0:mean_conf1, names_to='response') |>
  mutate(response=str_match(response, '[[:digit:]]')) |>
  ggplot(aes(x=c, y=value, color=response)) +
  geom_line() +
  scale_y_continuous('Mean Confidence', expand=expansion()) +
  scale_x_continuous(breaks=-1:1, expand=expansion()) +
  coord_cartesian(xlim=c(-1.5, 1.5), ylim=c(1, 2)) +
  theme_bw(18)
p.c.conf

p.c.mdistance <- fits.c |>
  pivot_longer(mdistance0:mdistance1, names_to='response') |>
  mutate(response=str_match(response, '[[:digit:]]')) |>
  ggplot(aes(x=c, y=value, color=response)) +
  geom_line() +
  scale_y_continuous('mdistance', expand=expansion()) +
  scale_x_continuous(breaks=-1:1,expand=expansion()) +
  coord_cartesian(xlim=c(-1.5, 1.5), ylim=c(0, 1)) +
  theme_bw(18)
p.c.mdistance

p.c.metadelta <- fits.c |>
  pivot_longer(metadelta0:metadelta1, names_to='response') |>
  mutate(response=str_match(response, '[[:digit:]]')) |>
  ggplot(aes(x=c, y=value, color=response)) +
  geom_hline(yintercept=.5, linetype='dashed') +
  geom_line() +
  scale_y_continuous(expression('meta-'*Delta), expand=expansion()) +
  scale_x_continuous(breaks=-1:1,expand=expansion()) +
  coord_cartesian(xlim=c(-1.5, 1.5), ylim=c(0, 1)) +
  theme_bw(18)
p.c.metadelta


## varying type 1 sensitivity (fixing meta-dprime)
fits.dprime <- expand_grid(N=N_trials,
                           dprime=seq(.5, 3, by=.1), c=0.5,
                           meta_dprime=1,
                           c2_0_diff=.5, c2_1_diff=.5) |>
  mutate(log_M=log(meta_dprime/dprime)) |>
  simulate_conf_bias(metac_absolute=FALSE, model=m.absolute)

p.dprime.conf <- fits.dprime |>
  pivot_longer(starts_with('mean_conf'), names_to='response') |>
  mutate(response=str_match(response, '[[:digit:]]')) |>
  ggplot(aes(x=dprime, y=value, color=response)) +
  geom_line() +
  scale_y_continuous('Mean Confidence', expand=expansion()) +
  scale_x_continuous("d'", expand=expansion()) +
  coord_cartesian(xlim=c(0, 3), ylim=c(1, 2)) +
  theme_bw(18)
p.dprime.conf

p.dprime.mdistance <- fits.dprime |>
  pivot_longer(mdistance0:mdistance1, names_to='response') |>
  mutate(response=str_match(response, '[[:digit:]]')) |>
  ggplot(aes(x=dprime, y=value, color=response)) +
  geom_line() +
  scale_y_continuous('mdistance', expand=expansion()) +
  scale_x_continuous("d'", expand=expansion()) +
  coord_cartesian(xlim=c(0, 3), ylim=c(0, 1)) +
  theme_bw(18)
p.dprime.mdistance

p.dprime.metadelta <- fits.dprime |>
  pivot_longer(metadelta0:metadelta1, names_to='response') |>
  mutate(response=str_match(response, '[[:digit:]]')) |>
  ggplot(aes(x=dprime, y=value, color=response)) +
  geom_hline(yintercept=.5, linetype='dashed') +
  geom_line() +
  scale_y_continuous(expression('meta-'*Delta), expand=expansion()) +
  scale_x_continuous("d'", expand=expansion()) +
  coord_cartesian(xlim=c(0, 3), ylim=c(0, 1)) +
  theme_bw(18)
p.dprime.metadelta


## varying type 2 sensitivity
fits.meta_dprime <- expand_grid(N=N_trials, dprime=2, c=.5,
                                meta_dprime=seq(.5, 3, by=.1),
                                c2_0_diff=.5, c2_1_diff=.5) |>
  mutate(log_M=log(meta_dprime/dprime)) |>
  simulate_conf_bias(metac_absolute=FALSE, model=m.absolute)

p.meta_dprime.conf <- fits.meta_dprime |>
  pivot_longer(starts_with('mean_conf'), names_to='response') |>
  mutate(response=str_match(response, '[[:digit:]]')) |>
  ggplot(aes(x=meta_dprime, y=value, color=response)) +
  geom_line() +
  scale_y_continuous('Mean Confidence', expand=expansion()) +
  scale_x_continuous("meta-d'", expand=expansion()) +
  coord_cartesian(xlim=c(0, 3), ylim=c(1, 2)) +
  theme_bw(18)
p.meta_dprime.conf

p.meta_dprime.mdistance <- fits.meta_dprime |>
  pivot_longer(mdistance0:mdistance1, names_to='response') |>
  mutate(response=str_match(response, '[[:digit:]]')) |>
  ggplot(aes(x=meta_dprime, y=value, color=response)) +
  geom_line() +
  scale_y_continuous('mdistance', expand=expansion()) +
  scale_x_continuous("meta-d'",  expand=expansion()) +
  coord_cartesian(xlim=c(0, 3), ylim=c(0, 1)) +
  theme_bw(18)
p.meta_dprime.mdistance

p.meta_dprime.metadelta <- fits.meta_dprime |>
  pivot_longer(metadelta0:metadelta1, names_to='response') |>
  mutate(response=str_match(response, '[[:digit:]]')) |>
  ggplot(aes(x=meta_dprime, y=value, color=response)) +
  geom_hline(yintercept=.5, linetype='dashed') +
  geom_line() +
  scale_y_continuous(expression('meta-'*Delta),  expand=expansion()) +
  scale_x_continuous("meta-d'", expand=expansion()) +
  coord_cartesian(xlim=c(0, 3), ylim=c(0, 1)) +
  theme_bw(18)
p.meta_dprime.metadelta



((p.c.conf | p.dprime.conf | p.meta_dprime.conf) /
   (p.c.metadelta | p.dprime.metadelta | p.meta_dprime.metadelta)) +
  plot_layout(guides='collect') +
  plot_annotation(tag_levels='A') &
  scale_color_manual('Response', values=PALETTE)
ggsave('plots/misspecification.pdf', width=10, height=6)
