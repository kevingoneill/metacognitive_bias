library(tidyverse)
library(distributional)
library(tidybayes)
library(paletteer)
library(patchwork)

## set SDT distributions
d <- tibble(stimulus=factor(0:1),
            d_prime=2, c=.5, meta_d_prime=1,
            mu1=c(-1, 1)*d_prime/2, sd=1,
            mu2=c(-1, 1)*meta_d_prime/2)


PALETTE <- paletteer_d("fishualize::Phractocephalus_hemioliopterus")[c(1, 3)]
DIST_COLOR <- paletteer_d("fishualize::Phractocephalus_hemioliopterus")[4]
DELTA_COLOR <- paletteer_d("fishualize::Phractocephalus_hemioliopterus")[2]

THEME_SDT <- function(xlab='Evidence', limits=c(-2.5, 2.5), base_size=18,
                      expand=expansion(mult=c(.1, .1))) {
  list(scale_x_continuous(xlab, limits=limits, breaks=seq(-4, 4, by=2), expand=expand),
       scale_color_manual('Stimulus', values=PALETTE),
       scale_fill_manual('Stimulus', values=PALETTE),
       scale_y_continuous(expand=c(0, 0)),
       theme_classic(base_size=base_size),
       theme(axis.title.y=element_blank(),
             axis.text.y=element_blank(),
             axis.ticks.y=element_blank(),
             axis.line.y=element_blank()))
}

geom_metadelta <- function(c=.5, c2_diff=c(.5, .5), response=0,
                           y=.75, size=6, gap=.025, ygap=.125, width=.025,
                           rel_y=2/3, rel_size=1, hjust=1/2,
                           dist_color=DIST_COLOR, delta_color=DELTA_COLOR) {
  if (response) {
    c2 <- c + cumsum(c2_diff)
  } else {
    c2 <- c - cumsum(c2_diff)
    gap <- -gap
    
  }
  
  i <- seq_along(c2)
  prev <- lag(c2, default=c)
  avg_c2 <- mean(c2)
  
  c(map(c2, \(c2) geom_vline(aes(xintercept=c2), linetype='dashed')),
    pmap(list(c2, prev, y), \(c2, prev, y)
         geom_errorbar(aes(y=y, xmin=prev+gap, xmax=c2-gap), width=width, color=dist_color)),
    pmap(list(c2, prev, y, i), \(c2, prev, y, i)
         geom_text(aes(y=y+ygap, x=(c2+prev)/2),
                   label=parse(text=paste0("'dmeta-c'[2*','*", i, "]^", response)),
                   size=size, color=dist_color)),
    geom_segment(aes(y=0, yend=y*rel_y, x=avg_c2, xend=avg_c2),
                 linetype='dotted', linewidth=.5, color=delta_color),
    geom_errorbar(aes(y=y*rel_y, xmin=c+gap, xmax=avg_c2-gap), width=width, color=delta_color),
    geom_text(aes(y=y*rel_y+ygap*2/3, x=c+(avg_c2-c)*hjust),
              label=parse(text=paste0('"meta-"*Delta[', response, ']')),
              color=delta_color, size=size*rel_size)
    )
}

p.type1 <- ggplot(d, aes(xdist=dist_normal(mu1, sd))) +
  stat_slab(aes(fill=stimulus), color=NA, alpha=.2, scale=.8, show.legend=FALSE) +
  geom_vline(aes(xintercept=c)) +
  THEME_SDT(xlab='Type 1 Evidence', expand=expansion())
p.type2.0 <- ggplot(d, aes(xdist=dist_truncated(dist_normal(mu2, sd), upper=c))) +
  stat_slab(aes(fill=stimulus), color=NA, alpha=.2, scale=.8, show.legend=FALSE) +
  geom_vline(aes(xintercept=c)) +
  geom_metadelta(c=.5, c2_diff=c(1.25, .75), response=0) +
  THEME_SDT(xlab='Type 2 Evidence\n("0" Response)', limits=c(-2.5, .5),
            expand=expansion(add=c(0, 0.01)))
p.type2.1 <- ggplot(d, aes(xdist=dist_truncated(dist_normal(mu2, sd), lower=c))) +
  stat_slab(aes(fill=stimulus), color=NA, alpha=.2, scale=.8, show.legend=FALSE) +
  geom_vline(aes(xintercept=c)) +
  geom_metadelta(c=.5, response=1, c2_diff=c(.75, 1), hjust=1) +
  THEME_SDT(xlab='Type 2 Evidence\n("1" Response)', limits=c(.5, 2.5),
            expand=expansion(add=c(.01, 0)))

(p.type1 /
   ((p.type2.0 | p.type2.1) +
      plot_layout(widths=c(3, 2)))) +
  plot_annotation(tag_levels='A')
ggsave('plots/meta_delta.pdf', width=10, height=6)
