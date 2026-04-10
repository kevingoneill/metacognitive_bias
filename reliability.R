library(tidyverse)
library(hmetad)
library(tidybayes)
library(patchwork)

#################################################################################################
##                                    Data processing
#################################################################################################
## read in and format data
data.reliability <- read_csv('20260330_124909_main_post_exclusion.csv') |>
  rename_with(tolower) |>
  rename(participant=id,
         correct=correct1st,
         confidence=confidence1st,
         stimulus=target) |>
  mutate(participant=factor(participant),
         session=factor(session),
         day=as.integer(day),
         block=factor(block),
         trial=factor(trial),
         stimulus=to_unsigned(stimulus),
         response=ifelse(correct, stimulus, 1-stimulus),
         confidence=as.integer(factor(confidence, levels=c(67, 83, 100),
                                      labels=1:3)),
         dotdiff=scale(dotdiff)[,1]) |>
  select(participant, session, day, block, trial, infocost, dotdiff,
         stimulus, response, correct, confidence, seeagain)
data.reliability

## count number of missing cells per participant (ignoring session)
counts_participant <- data.reliability |>
  mutate(stimulus=factor(stimulus),
         joint_response=factor(joint_response(response, confidence, K=3))) |>
  count(participant, stimulus, joint_response, .drop=FALSE) |>
  group_by(participant) |>
  summarize(n_participant=sum(n==0))

## count number of missing cells per participant/session
counts_session <- data.reliability |>
  mutate(stimulus=factor(stimulus),
         joint_response=factor(joint_response(response, confidence, K=3))) |>
  count(participant, session, stimulus, joint_response, .drop=FALSE) |>
  group_by(participant, session) |>
  summarize(n_session=sum(n==0))

## exclude sessions with 3 or more empty stimulus/response/confidence cells
## and pre-aggregate data
data.aggregated <- data.reliability |>
  left_join(counts_participant) |>
  left_join(counts_session) |>
  filter(n_participant < 3,
         n_session < 3) |>
  select(-n_participant, -n_session) |>
  mutate(participant=factor(participant)) |>
  aggregate_metad(participant, session) |>
  filter(N_0 > 0, N_1 > 0)
data.aggregated

#################################################################################################
##                                    Model fitting
#################################################################################################
## show parameters for 3 confidence levels
metad(n_distinct(data.reliability$confidence))

# sample from the prior
prior.reliability <- fit_metad(bf(N ~ (0 + session | participant),
                                  dprime+c ~ (0 + session | participant),
                                  metac2zero1diff+metac2zero2diff ~ (0 + session | participant),
                                  metac2one1diff+metac2one2diff ~ (0 + session | participant)),
                               data.aggregated, aggregate=FALSE,
                               init=0, cores=4, sample_prior='only', 
                               prior=prior(normal(0, .25), class=Intercept) +
                                 prior(normal(1, .25), class=Intercept, dpar=dprime) +
                                 prior(normal(0, .25), class=Intercept, dpar=c) +
                                 set_prior('normal(-0.5, .5)', class='Intercept',
                                           dpar=c('metac2zero1diff',
                                                  'metac2zero2diff',
                                                  'metac2one1diff',
                                                  'metac2one2diff')) +
                                 prior(normal(0, 1), class=sd) +
                                 set_prior('normal(0, 1)', class='sd',
                                           dpar=c('dprime', 'c',
                                                  'metac2zero1diff',
                                                  'metac2zero2diff',
                                                  'metac2one1diff',
                                                  'metac2one2diff')),
                               backend='cmdstanr')
prior.reliability

# prior predictive check on joint response probabilities
data.reliability |>
  distinct(session) |>
  add_epred_rvars_metad(prior.reliability, re_formula=NA) |>
  ggplot(aes(x=factor(joint_response), group=session)) +
  stat_pointinterval(aes(ydist=.epred)) +
  facet_grid( ~ stimulus) +
  theme_classic()

# sample from posterior
m.reliability <- update(prior.reliability, sample_prior='no', 
                        cores=4, file='reliability',
                        control=list(adapt_delta=.99),
                        backend='cmdstanr')
print(m.reliability, prior=TRUE)



#################################################################################################
##                                  Compute reliability correlations
#################################################################################################
## obtain model estimates
draws.linpred <- data.aggregated |>
  expand(participant, session) |>
  add_linpred_draws(m.reliability, value='M', dpar=c('dprime', 'c'), transform=TRUE) |>
  left_join(data.aggregated |>
              expand(participant, session) |>
              add_metacognitive_bias_draws(m.reliability) |>
              pivot_wider(names_from=response, values_from=metacognitive_bias,
                          names_prefix='meta_delta') |>
              mutate(meta_delta=map2_dbl(meta_delta0, meta_delta1, ~ mean(c(.x, .y))))) |>
  pivot_longer(M:meta_delta, names_to='.variable', values_to='.value') |>
  group_by(participant, session, .variable) |>
  select(-.row)
draws.linpred

#' plot participant-level estimates across three sessions
#' along with correlations between sessions
session_plot <- function(draws, variable='M', buffer=.25, label=NULL,
                         text_x=NULL, text_y=NULL, text_size=4,
                         size=1, linewidth=.1, alpha=.4) {
  d <- draws |>
    filter(.variable==variable) |>
    median_qi(.value)

  if (is.null(label))
    label <- variable
  
  lims <- range(d$.value) + c(-buffer, buffer)
  d <- d |>
    pivot_wider(names_from=session, values_from=c(.value, .lower, .upper))

  if (is.null(text_x))
    text_x <- lims[1] + (lims[2]-lims[1])*3/4
  if (is.null(text_y))
    text_y <- lims[1] + (lims[2]-lims[1])/8
  
  
  d.cor <- draws |>
    filter(.variable==variable) |>
    pivot_wider(names_from=session, values_from=.value, names_prefix='session') |>
    group_by(.variable, .draw) |>
    summarize(r12=cor(session1, session2, use='complete.obs'),
              r13=cor(session1, session3, use='complete.obs'),
              r23=cor(session2, session3, use='complete.obs')) |>
    median_qi(r12, r13, r23) |>
    mutate(r_1_2=sprintf('r = %.2f [%.2f, %.2f]', r12, r12.lower, r12.upper),
           r_1_3=sprintf('r = %.2f [%.2f, %.2f]', r13, r13.lower, r13.upper),
           r_2_3=sprintf('r = %.2f [%.2f, %.2f]', r23, r23.lower, r23.upper),
           x=text_x,
           y=text_y) |>
    select(.variable, r_1_2, r_1_3, r_2_3, x, y)
  
  
  p1 <- ggplot(d, aes(x=.value_1, y=.value_2)) +
    geom_errorbar(aes(xmin=.lower_1, xmax=.upper_1), linewidth=linewidth, alpha=alpha) +
    geom_errorbar(aes(ymin=.lower_2, ymax=.upper_2), linewidth=linewidth, alpha=alpha) +
    geom_abline(slope=1, intercept=0, linetype='dashed') +
    geom_point(size=size) +
    xlab(parse(text=paste0(label, '~(Session~1)'))) +
    ylab(parse(text=paste0(label, '~(Session~2)'))) +
    geom_label(aes(x=x, y=y, label=r_1_2), data=d.cor,
               hjust=.5, vjust=.5, size=text_size)

  p2 <- ggplot(d, aes(x=.value_1, y=.value_3)) +
    geom_errorbar(aes(xmin=.lower_1, xmax=.upper_1), linewidth=linewidth, alpha=alpha) +
    geom_errorbar(aes(ymin=.lower_3, ymax=.upper_3), linewidth=linewidth, alpha=alpha) +
    geom_abline(slope=1, intercept=0, linetype='dashed') +
    geom_point(size=size) +
    xlab(parse(text=paste0(label, '~(Session~1)'))) +
    ylab(parse(text=paste0(label, '~(Session~3)'))) +
    geom_label(aes(x=x, y=y, label=r_1_3), data=d.cor,
               hjust=.5, vjust=.5, size=text_size)

  p3 <- ggplot(d, aes(x=.value_2, y=.value_3)) +
    geom_errorbar(aes(xmin=.lower_2, xmax=.upper_2), linewidth=linewidth, alpha=alpha) +
    geom_errorbar(aes(ymin=.lower_3, ymax=.upper_3), linewidth=linewidth, alpha=alpha) +
    geom_abline(slope=1, intercept=0, linetype='dashed') +
    geom_point(size=size) +
    xlab(parse(text=paste0(label, '~(Session~2)'))) +
    ylab(parse(text=paste0(label, '~(Session~3)'))) +
    geom_label(aes(x=x, y=y, label=r_2_3), data=d.cor,
               hjust=.5, vjust=.5, size=text_size)
  
  (p1 | p3 | p2) &
    coord_fixed(xlim=lims, ylim=lims) &
    theme_classic(18) &
    theme(panel.border=element_rect(linewidth=1.5),
          axis.line=element_blank())
}

session_plot(draws.linpred, 'M', text_y=0.25) &
  coord_cartesian(xlim=c(0, 2), ylim=c(0, 2))
ggsave('plots/reliability_mratio.pdf', width=12, height=4)

session_plot(draws.linpred, 'dprime', text_x=1.25, text_y=0.15) &
  coord_cartesian(xlim=c(0, 2.5), ylim=c(0, 2.5))
ggsave('plots/reliability_dprime.pdf', width=12, height=4)

session_plot(draws.linpred, 'c', text_y=-.95)
ggsave('plots/reliability_c.pdf', width=12, height=4)

((session_plot(draws.linpred, 'meta_delta0', label='"meta-"*Delta[0]', text_x=2, text_y=1/3) /
  session_plot(draws.linpred, 'meta_delta1', label='"meta-"*Delta[1]', text_x=2, text_y=1/3) /
  session_plot(draws.linpred, 'meta_delta', label='"meta-"*Delta', text_x=2, text_y=1/3)) &
  coord_fixed(xlim=c(0, 3), ylim=c(0, 3))) +
  plot_annotation(tag_levels='A')
ggsave('plots/reliability.pdf', width=12, height=12)



## compute posteriors over test-retest correlations
draws.cor <- draws.linpred |>
  pivot_wider(names_from=session, values_from=.value, names_prefix='session') |>
  group_by(.variable, .draw) |>
  summarize(r12=cor(session1, session2),
            r13=cor(session1, session3),
            r23=cor(session2, session3))

draws.cor |>
  median_qi() |>
  mutate(r_1_2=sprintf('%.2f [%.2f, %.2f]', r12, r12.lower, r12.upper),
         r_2_3=sprintf('%.2f [%.2f, %.2f]', r23, r23.lower, r23.upper),
         r_1_3=sprintf('%.2f [%.2f, %.2f]', r13, r13.lower, r13.upper)) |>
  select(.variable, r_1_2, r_1_3, r_2_3)
