library(tidyverse)
library(hmetad)
library(tidybayes)
library(patchwork)

#################################################################################################
##                                    Data processing
#################################################################################################
## read in and format data
data.reliability <- read_csv('20260421_182620_main_post_exclusion.csv') |>
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

(ggplot(counts_participant, aes(x=n_participant, fill=n_participant<=3)) /
   ggplot(counts_session, aes(x=n_session, fill=n_session<3))) &
  geom_histogram(binwidth=1) &
  theme_classic(18)


## exclude sessions with 3 or more empty stimulus/response/confidence cells
count(data.reliability, participant, session)
data.reliability <- data.reliability |>
  left_join(counts_participant) |>
  left_join(counts_session) |>
  filter(n_participant <= 3,
         n_session < 3) |>
  select(-n_participant, -n_session) |>
  mutate(participant=factor(participant))
count(data.reliability, participant, session)

#################################################################################################
##                                    Model fitting
#################################################################################################
## show parameters for 3 confidence levels
metad(n_distinct(data.reliability$confidence))
dpar_metac2 <- c('metac2zero1diff', 'metac2zero2diff', 'metac2one1diff', 'metac2one2diff')

# sample from the prior
prior.reliability <- fit_metad(
  bf(N ~ (0 + session | participant),
     dprime+c ~ (0 + session | participant),
     metac2zero1diff+metac2zero2diff ~ (0 + session |p0| participant),
     metac2one1diff+metac2one2diff ~ (0 + session |p1| participant)),
  data.reliability,
  init=0, cores=4, sample_prior='only', backend='cmdstanr',
  prior=prior(normal(0, .25), class=Intercept) +
    prior(normal(1, .25), class=Intercept, dpar=dprime) +
    prior(normal(0, .1), class=Intercept, dpar=c) +
    set_prior('normal(-0.5, 1)', class='Intercept', dpar=dpar_metac2) +
    prior(normal(0, 0.5), class=sd) +
    set_prior('normal(0, 0.33)', class='sd', dpar=c('dprime', 'c')) +
    set_prior('normal(0, 0.5)', class='sd', dpar=dpar_metac2) +
    prior(lkj(2), class=cor))
prior.reliability

# prior predictive check on joint response probabilities
tibble(session=1) |>
  add_epred_rvars_metad(prior.reliability, re_formula=NA) |>
  ggplot(aes(x=factor(joint_response),
             color=factor(response), alpha=factor(confidence),
             group=session)) +
  stat_pointinterval(aes(ydist=.epred)) +
  scale_alpha_discrete() +
  facet_grid( ~ stimulus, labeller=label_both) +
  theme_classic()

expand_grid(session=1,
            participant=first(levels(data.reliability$participant))) |>
  add_epred_rvars_metad(prior.reliability) |>
  ggplot(aes(x=factor(joint_response),
             color=factor(response), alpha=factor(confidence),
             group=session)) +
  stat_pointinterval(aes(ydist=.epred)) +
  scale_alpha_discrete() +
  facet_grid( ~ stimulus, labeller=label_both) +
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
draws.linpred <- data.reliability |>
  expand(participant, session) |>
  add_linpred_draws(m.reliability, value='M', dpar=c('dprime', 'c'), transform=TRUE) |>
  left_join(data.reliability |>
              expand(participant, session) |>
              add_metacognitive_bias_draws(m.reliability) |>
              pivot_wider(names_from=response, values_from=metacognitive_bias,
                          names_prefix='meta_delta') |>
              mutate(meta_delta=map2_dbl(meta_delta0, meta_delta1, ~ mean(c(.x, .y))))) |>
  pivot_longer(M:meta_delta, names_to='.variable', values_to='.value') |>
  group_by(participant, session, .variable) |>
  select(-.row)
draws.linpred


## compute posteriors over test-retest correlations
draws.cor <- draws.linpred |>
  pivot_wider(names_from=session, values_from=.value, names_prefix='session') |>
  group_by(.variable, .draw) |>
  summarize(r12=cor(session1, session2, use='complete.obs'),
            r13=cor(session1, session3, use='complete.obs'),
            r14=cor(session1, session4, use='complete.obs'),
            r23=cor(session2, session3, use='complete.obs'),
            r24=cor(session2, session4, use='complete.obs'),
            r34=cor(session3, session4, use='complete.obs')) |>
  median_qi(r12, r13, r14, r23, r24, r34) |>
  mutate(r_1_2=sprintf('r = %.2f [%.2f, %.2f]', r12, r12.lower, r12.upper),
         r_1_3=sprintf('r = %.2f [%.2f, %.2f]', r13, r13.lower, r13.upper),
         r_1_4=sprintf('r = %.2f [%.2f, %.2f]', r14, r14.lower, r14.upper),
         r_2_3=sprintf('r = %.2f [%.2f, %.2f]', r23, r23.lower, r23.upper),
         r_2_4=sprintf('r = %.2f [%.2f, %.2f]', r24, r24.lower, r24.upper),
         r_3_4=sprintf('r = %.2f [%.2f, %.2f]', r34, r34.lower, r34.upper)) |>
  select(.variable, r_1_2, r_1_3, r_1_4, r_2_3, r_2_4, r_3_4)
draws.cor


#' plot participant-level estimates across three sessions
#' along with correlations between sessions
session_plot <- function(draws, draws.cor, variable='M', buffer=.25, label=NULL,
                         text_x=NULL, text_y=NULL, text_size=4,
                         size=1, linewidth=.1, alpha=.4,
                         layout=`+`, all_pairs=TRUE) {
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
  
  
  d.cor <- draws.cor |>
    filter(.variable==variable) |>
    mutate(x=text_x, y=text_y)
  
  pairs_plot <- function(session1, session2) {
    ggplot(d) +
      aes_string(x=paste0('.value_', session1), y=paste0('.value_', session2)) +
      geom_errorbar(aes_string(xmin=paste0('.lower_', session1),
                               xmax=paste0('.upper_', session1)),
                    linewidth=linewidth, alpha=alpha) +
      geom_errorbar(aes_string(ymin=paste0('.lower_', session2),
                               ymax=paste0('.upper_', session2)),
                    linewidth=linewidth, alpha=alpha) +
      geom_abline(slope=1, intercept=0, linetype='dashed') +
      geom_point(size=size) +
      xlab(parse(text=paste0(label, '~(Session~', session1, ')'))) +
      ylab(parse(text=paste0(label, '~(Session~', session2, ')'))) +
      geom_label(aes_string(x='x', y='y', label=paste0('r_', session1, '_', session2)),
                 data=d.cor, hjust=.5, vjust=.5, size=text_size) +
      coord_fixed(xlim=lims, ylim=lims) +
      theme_classic(18) +
      theme(panel.border=element_rect(linewidth=1.5),
            axis.line=element_blank())
  }

  sessions <- draws |> pull(session) |> unique() |> as.integer()

  if (all_pairs) {
    sessions <- expand_grid(session1=sessions, session2=sessions)
  } else {
    sessions <- expand_grid(session1=first(sessions), session2=sessions)
  }

  sessions |>
    filter(session1 < session2) |>
    arrange(session2 - session1) |>
    mutate(plots=map2(session1, session2, pairs_plot)) |>
    pull(plots) |>
    reduce(layout)
}

session_plot(draws.linpred, draws.cor, 'M', text_y=0.25) &
  coord_cartesian(xlim=c(0, 2), ylim=c(0, 2))
ggsave('plots/reliability_mratio.pdf', width=12, height=8)

session_plot(draws.linpred, draws.cor, 'dprime', text_x=1.25, text_y=0.15) &
  coord_cartesian(xlim=c(0, 2.5), ylim=c(0, 2.5))
ggsave('plots/reliability_dprime.pdf', width=12, height=8)

session_plot(draws.linpred, draws.cor, 'c', text_y=-.95)
ggsave('plots/reliability_c.pdf', width=12, height=8)

((session_plot(draws.linpred, draws.cor,
               'meta_delta0', label='"meta-"*Delta[0]',
               text_x=2, text_y=1/3, layout=`/`) |
    session_plot(draws.linpred, draws.cor,
                 'meta_delta1', label='"meta-"*Delta[1]',
                 text_x=2, text_y=1/3, layout=`/`) |
    session_plot(draws.linpred, draws.cor,
                 'meta_delta', label='"meta-"*Delta',
                 text_x=2, text_y=1/3, layout=`/`)) &
  coord_fixed(xlim=c(0, 3), ylim=c(0, 3))) +
  plot_annotation(tag_levels='A')
ggsave('plots/reliability_4.pdf', width=12, height=24)



((session_plot(draws.linpred, draws.cor,
               'meta_delta0', label='"meta-"*Delta[0]',
               text_x=2, text_y=1/3, layout=`|`, all_pairs=FALSE) /
    session_plot(draws.linpred, draws.cor,
                 'meta_delta1', label='"meta-"*Delta[1]',
                 text_x=2, text_y=1/3, layout=`|`, all_pairs=FALSE) /
    session_plot(draws.linpred, draws.cor,
                 'meta_delta', label='"meta-"*Delta',
                 text_x=2, text_y=1/3, layout=`|`, all_pairs=FALSE)) &
  coord_fixed(xlim=c(0, 3), ylim=c(0, 3))) +
  plot_annotation(tag_levels='A')
ggsave('plots/reliability.pdf', width=12, height=12)
