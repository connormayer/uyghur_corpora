require(lme4)
require(tidyverse)
require(archive)

setwd("C:/Users/conno/ling/uyghur_corpora")
rfa_file <- archive_read("corpora/rfa/output/conservative_parses.zip")
rfa_data <- read_csv(rfa_file, col_types = cols())
rfa_data$source <- 'rfa'
awazi_file <- archive_read("corpora/awazi/output/conservative_parses.zip")
awazi_data <- read_csv(awazi_file, col_types = cols())
awazi_data$source <- 'awazi'
akademiye_file <- archive_read("corpora/akademiye/output/conservative_parses.zip")
akademiye_data <- read_file(akademiye_file, col_types = cols())
akademiye_data$source <- 'akademiye'

full_data <- rbind(rfa_data, awazi_data, akademiye_data)
n = sum(full_data$count)

counts <- full_data %>% 
  group_by(root) %>%
  summarize(total_count = sum(count),
            # Convert counts to counts per million
            norm_count = 1000000 * (total_count / n),
            raised_form_prop = sum(raised) / total_count)

full_data <- inner_join(full_data, counts, by="root")

possible_opaque_harmonizers <- full_data %>%
  filter(last_two == 'BF' | last_two == 'FB') %>%
  filter(raising_candidate) %>%
  filter(back_count + front_count > 0)

possible_opaque_harmonizers <- possible_opaque_harmonizers %>% 
  mutate(
    has_blocker = grepl("gpr_rsub3|ger_past2|prog|px2pl", possible_opaque_harmonizers$tags),
    has_ye = grepl("ye$", possible_opaque_harmonizers$root),
    has_ane = grepl("ane$", possible_opaque_harmonizers$root),
    has_che = grepl("che$", possible_opaque_harmonizers$root),
    has_name = grepl('name$', possible_opaque_harmonizers$root),
  )

# Remove tokens with suffixes that behave oddly
possible_opaque_harmonizers[possible_opaque_harmonizers$root == 'anche',]$has_che <- FALSE
possible_opaque_harmonizers[possible_opaque_harmonizers$root == 'bunche',]$has_che <- FALSE
possible_opaque_harmonizers[possible_opaque_harmonizers$root == 'qoshumche',]$has_che <- FALSE
possible_opaque_harmonizers[possible_opaque_harmonizers$root == 'qanche',]$has_che <- FALSE
possible_opaque_harmonizers[possible_opaque_harmonizers$root == 'birqanche',]$has_che <- FALSE

possible_opaque_harmonizers[possible_opaque_harmonizers$root == 'bahane',]$has_ane <- FALSE
possible_opaque_harmonizers[possible_opaque_harmonizers$root == 'durdane',]$has_ane <- FALSE
possible_opaque_harmonizers[possible_opaque_harmonizers$root == 'epsane',]$has_ane <- FALSE
possible_opaque_harmonizers[possible_opaque_harmonizers$root == "i'ane",]$has_ane <- FALSE
possible_opaque_harmonizers[possible_opaque_harmonizers$root == 'jerimane',]$has_ane <- FALSE
possible_opaque_harmonizers[possible_opaque_harmonizers$root == 'perwane',]$has_ane <- FALSE
possible_opaque_harmonizers[possible_opaque_harmonizers$root == 'perghane',]$has_ane <- FALSE
possible_opaque_harmonizers[possible_opaque_harmonizers$root == 'péshane',]$has_ane <- FALSE

# Filter out a few bogus roots
possible_opaque_harmonizers <- possible_opaque_harmonizers %>%
  filter(
     # Not a real word
     root != 'qure' & 
     # Correct Uyghur word is komitét, but seems to often be written
     # as kométit
     root != 'kométet' &
     # Correct root is teqqasla
     root != 'teqqasle' &
     # Correct root is etiwarla
     root != 'etiwarle' &
     # Correct root is mejburla
     root != 'mejburle' &
     # Eli is the name 'Ali', not a raised form of 'ela'
     root != 'ela' &
     # Alem is a real word, but these all look to be cases of alim 'scholar'
     root != 'alem' &
     # Nur'ela is a name, but these stories refer to a man named Nur'eli
     root != "nur'ela" &
     # These correspond to esme 'memory', not esma
     root != 'esma'
     # Is this right?
     # root != 'xarabe'
  )

tokens <- possible_opaque_harmonizers %>% filter(
  !has_che & !has_ye & !has_ane & !has_blocker, !has_name
) %>%
  select(-has_blocker, -has_che, -has_ye, -has_ane, -has_name)

# Does opacity vary as a result of overall frequency

harmonic_raisers <- tokens %>%
  filter(raised) %>%
  mutate(
    opaque = case_when(back_count > 0 & last_two == 'BF' ~ 1,
                       front_count > 0 & last_two == 'FB' ~ 1,
                       TRUE ~ 0),
    log_total_count = log(total_count)
  )

root_agg <- harmonic_raisers %>% 
  group_by(root, last_two) %>%
  summarise(n=sum(count),
            percent_back = sum(back_count) / n,
            total_count = mean(total_count),
            raised_form_prop = mean(raised_form_prop))

m1 <- glmer(opaque ~ log_total_count + raised_form_prop + last_two + (1|root), data=harmonic_raisers, family="binomial")
summary(m1)


# Get all raisers

raisers_no_suffixes_full <- stem_agg[(stem_agg$LastTwo == "BF" | stem_agg$LastTwo == "FB" | stem_agg$LastTwo == 'BB' | stem_agg$LastTwo == 'FF') &
                                  !stem_agg$HasChe &
                                  !stem_agg$HasAne &
                                  !stem_agg$HasYe &
                                  !stem_agg$HasName, ]


raised_forms <-  full_data %>% group_by(Stem) %>%
  summarise(raised_form_prop = sum(Count[Raises %in% TRUE]) / mean(total_count))
raised_forms$raised_form_prop[is.na(raised_forms$raised_form_prop)] <- 0


raising_rates <- possible_raisers %>% group_by(Stem) %>%
  summarise(
    rate=sum(Count[Raises %in% TRUE]) / sum(Count)
  )

raisers_no_suffixes_full <- merge(raisers_no_suffixes_full, raised_forms, by="Stem")
raisers_no_suffixes_full <- merge(raisers_no_suffixes_full, raising_rates, by="Stem")
raisers_no_suffixes_full <- distinct(raisers_no_suffixes_full)

raisers_no_suffixes <- stem_agg[(stem_agg$LastTwo == "BF" | stem_agg$LastTwo == "FB") &
  !stem_agg$HasChe &
  !stem_agg$HasAne &
  !stem_agg$HasYe &
  !stem_agg$HasName, ]

raisers_no_suffixes$percent_opaque <- raisers_no_suffixes$percent_back
raisers_no_suffixes[raisers_no_suffixes$LastTwo == "BF",]$percent_opaque <- 1 - raisers_no_suffixes[raisers_no_suffixes$LastTwo == "BF",]$percent_back

raised_forms <-  full_data %>% group_by(Stem) %>%
  summarise(raised_form_prop = sum(Count[Raises %in% TRUE]) / mean(total_count))
raised_forms$raised_form_prop[is.na(raised_forms$raised_form_prop)] <- 0

raising_rates <- possible_raisers %>% group_by(Stem) %>%
  summarise(
    rate=sum(Count[Raises %in% TRUE]) / sum(Count)
  )

raisers_no_suffixes <- merge(raisers_no_suffixes, raised_forms, by="Stem")
raisers_no_suffixes <- merge(raisers_no_suffixes, raising_rates, by="Stem")
raisers_no_suffixes <- distinct(raisers_no_suffixes)

raisers_no_suffixes$raised_tokens <- raisers_no_suffixes$total_count * raisers_no_suffixes$raised_form_prop
raisers_no_suffixes$unraised_tokens <- raisers_no_suffixes$total_count - raisers_no_suffixes$raised_tokens

raisers_no_suffixes <- droplevels(raisers_no_suffixes)

foo <- glm(percent_opaque ~ raised_form_prop + log(total_count) + LastTwo, data=raisers_no_suffixes)
summary(foo)

bar <- gam(percent_opaque ~ s(raised_form_prop) + s(log(total_count)) + LastTwo, data=raisers_no_suffixes)

foo <- glm(percent_opaque ~ raised_form_prop + log(total_count), data=raisers_no_suffixes[raisers_no_suffixes$LastTwo == 'BF',])
summary(foo)

foo <- glm(percent_opaque ~ log(total_count) + LastTwo, data=raisers_no_suffixes)
summary(foo)

cor.test(raisers_no_suffixes$percent_opaque, log(raisers_no_suffixes$total_count), method="pearson")
cor.test(raisers_no_suffixes$percent_opaque, raisers_no_suffixes$raised_form_prop, method='pearson')
cor.test(raisers_no_suffixes$percent_opaque, raisers_no_suffixes$rate, method='pearson')

learned_weights <- fread('learned_spu_weights.csv', encoding="UTF-8")
comb <- cbind(raisers_no_suffixes, learned_weights)
comb$unraised_form_prop <- 1-comb$raised_form_prop

cor.test(comb$Weight, log(comb$total_count), method="pearson")
cor.test(comb$Weight, comb$raised_form_prop, method='pearson')
cor.test(comb$Weight, comb$rate, method='pearson')
cor.test(comb$Weight, log(comb$raised_tokens), method='pearson')

comb[comb$unraised_tokens < 1,]$unraised_tokens <- 1
cor.test(comb[comb$unraised_tokens >= 1,]$Weight, log(comb[comb$unraised_tokens >= 1,]$unraised_tokens), method='pearson')

foo <- glm(Weight ~ raised_form_prop + log(total_count) + LastTwo, data=comb)
foo <- glm(Weight ~ raised_form_prop + log(total_count) + LastTwo, data=comb)
summary(foo)

plot(comb$Weight, comb$raised_form_prop)
plot(comb$Weight, log(comb$total_count))

ggplot(data=comb, aes(x=log(total_count), y=Weight)) +
  geom_point(size=4) +
  xlab("Log root token count") +
  ylab("Indexed constraint weight") +
  theme(axis.text.x = element_text(size=20),
        axis.text.y = element_text(size=20),
        axis.title.x = element_text(size=25),
        axis.title.y = element_text(size=25),
        plot.title = element_text(size=30, hjust=0.5)) +
 geom_smooth(method="glm",
              se=FALSE) +
ggsave("freq_plot_weights.png")

ggplot(data=comb, aes(x=raised_form_prop, y=Weight)) +
  geom_point(size=4) +
  xlab("Proportion of raised forms") +
  ylab("Indexed constraint weight") +
  theme(axis.text.x = element_text(size=20),
        axis.text.y = element_text(size=20),
        axis.title.x = element_text(size=25),
        axis.title.y = element_text(size=25),
        plot.title = element_text(size=30, hjust=0.5)) +
  geom_smooth(method="glm",
              se=FALSE) +
ggsave("raising_prop_plot_weights.png")

ggplot(data=comb, aes(x=log(unraised_tokens), y=Weight)) +
  geom_point(size=4) +
  xlab("Log number of unraised tokens") +
  ylab("Indexed constraint weight") +
  theme(axis.text.x = element_text(size=20),
        axis.text.y = element_text(size=20),
        axis.title.x = element_text(size=25),
        axis.title.y = element_text(size=25),
        plot.title = element_text(size=30, hjust=0.5)) +
  geom_smooth(method="glm",
              se=FALSE) +
ggsave("unraised_tokens_plot_weights.png")

ggplot(data=comb, aes(x=log(raised_tokens), y=Weight)) +
  geom_point(size=4) +
  xlab("Log number of raised tokens") +
  ylab("Indexed constraint weight") +
  theme(axis.text.x = element_text(size=20),
        axis.text.y = element_text(size=20),
        axis.title.x = element_text(size=25),
        axis.title.y = element_text(size=25),
        plot.title = element_text(size=30, hjust=0.5)) +
  geom_smooth(method="glm",
              se=FALSE) +
ggsave("raised_tokens_plot_weights.png")



shit <- possible_raisers[possible_raisers$IsNeutral == FALSE,]
shit <- shit[shit$HasBlocker == FALSE,]
shit <- shit[shit$Ambiguous == FALSE,]
shit <- shit[shit$LastTwo == "BF" | shit$LastTwo == "FB",]

stem_agg_bar_plot <- shit %>% group_by(Stem, Raises, LastTwo, RaiserClass, HasChe, HasYe, HasAne, HasName) %>%
  summarise(n=sum(Back.Count) + sum(Front.Count),
            percent_back = sum(Back.Count) / n,
            percent_back_sd = sd(Back.Count))

stem_agg_bar_plot <- stem_agg_bar_plot[(stem_agg_bar_plot$LastTwo == "BF" | stem_agg_bar_plot$LastTwo == "FB") &
                                         !stem_agg_bar_plot$HasChe &
                                         !stem_agg_bar_plot$HasAne &
                                         !stem_agg_bar_plot$HasYe &
                                         !stem_agg_bar_plot$HasName, ]

nrow(stem_agg_bar_plot)
nrow(stem_agg_bar_plot[stem_agg_bar_plot$LastTwo == 'BF',])
nrow(stem_agg_bar_plot[stem_agg_bar_plot$LastTwo == 'FB',])

percent_graph <- raisers_no_suffixes %>% group_by(LastTwo) %>%
  summarise(opaque_harmony = mean(percent_opaque),
            transparent_harmony = 1 - opaque_harmony)

token_graph <- raisers_no_suffixes %>% group_by(LastTwo) %>%
  summarise(opaque_harmony = round(mean(percent_opaque) * sum(n)),
            transparent_harmony = round((1 - mean(percent_opaque)) * sum(n))
  )


percent_graph <- melt(as.data.frame(percent_graph), idvar="LastTwo")
token_graph <- melt(as.data.frame(token_graph), idvar="LastTwo")
graph_data <- merge(percent_graph, token_graph, by=c("LastTwo", "variable"))

ggplot(data=graph_data, aes(x=LastTwo, y=value.x, fill=factor(variable), label=value.y)) +
  geom_bar(stat="identity") +
  geom_text(size=6, position=position_stack(vjust=0.6), aes(label = paste("n =", value.y))) +
  theme(axis.text.x = element_text(angle = -60, hjust = 0, size=20),
        axis.text.y = element_text(size=20),
        axis.title.x = element_text(size=25),
        axis.title.y = element_text(size=25),
        #plot.title = element_text(size=2, hjust=0.5),
        legend.title = element_blank(),
        legend.text=element_text(size=20)) +
  xlab("Underlying root template") +
  ylab("Proportion of tokens") +
  scale_y_continuous(breaks=seq(0,1,0.25), limits=c(0,1)) +
  scale_fill_discrete(labels=c("Opaque harmony", "Transparent harmony")) +
  #ggtitle("Harmonizing behavior in raised tokens with harmonizing suffixes")
  ggsave("full_harmonic_raisers.png")

ggplot(raisers_no_suffixes[raisers_no_suffixes$LastTwo == "BF",]) +
  geom_histogram(aes(percent_opaque), binwidth=0.1) +
  xlab("Proportion of opaque responses for BF roots") +
  ylab("Root type count") +
  theme(axis.text.x = element_text(hjust = 0, size=20),
        axis.text.y = element_text(size=20),
        axis.title.x = element_text(size=25),
        axis.title.y = element_text(size=25))
ggsave('bf_histogram.png')

ggplot(raisers_no_suffixes[raisers_no_suffixes$LastTwo == "FB",]) +
  geom_histogram(aes(percent_opaque), binwidth=0.1) +
  xlab("Proportion of opaque responses for FB roots") +
  ylab("Root type count") +
  theme(axis.text.x = element_text(hjust = 0, size=20),
        axis.text.y = element_text(size=20),
        axis.title.x = element_text(size=25),
        axis.title.y = element_text(size=25))
ggsave('fb_histogram.png')


ggplot(data=raisers_no_suffixes, aes(x=log(total_count), y=percent_opaque)) +
  #geom_point(size=4) +
  geom_jitter(size=3, height=0.01) +
  xlab("Log root token count") +
  ylab("Proportion of opaque tokens") +
  theme(axis.text.x = element_text(size=20),
        axis.text.y = element_text(size=20),
        axis.title.x = element_text(size=25),
        axis.title.y = element_text(size=25),
        plot.title = element_text(size=30, hjust=0.5)) +
 geom_smooth(method="glm",
              se=FALSE) +
ggsave("freq_plot.png")


ggplot(data=raisers_no_suffixes, aes(x=raised_form_prop, y=percent_opaque)) +
  #geom_point(size=4) +
  geom_jitter(size=3, height=0.01) +
  xlab("Proportion of raised tokens") +
  ylab("Proportion of opaque tokens") +
  theme(axis.text.x = element_text(size=20),
        axis.text.y = element_text(size=20),
        axis.title.x = element_text(size=25),
        axis.title.y = element_text(size=25),
        plot.title = element_text(size=30, hjust=0.5)) +
  geom_smooth(method="glm",
              se=FALSE) +
ggsave("raised_form_plot.png")
#
# ggplot(raisers_no_suffixes[raisers_no_suffixes$LastTwo == "BF",]) +
#   geom_histogram(aes(percent_opaque), binwidth=0.1) +
#   xlab("Proportion of opaque responses for BF stems") +
#   ylab("Stem count") +
#   theme(axis.text.x = element_text(hjust = 0, size=20),
#         axis.text.y = element_text(size=20),
#         axis.title.x = element_text(size=25),
#         axis.title.y = element_text(size=25))
# ggsave('bf_histogram.png')
#
# ggplot(raisers_no_suffixes[raisers_no_suffixes$LastTwo == "FB",]) +
#   geom_histogram(aes(percent_opaque), binwidth=0.1) +
#   xlab("Proportion of opaque responses for FB stems") +
#   ylab("Stem count") +
#   theme(axis.text.x = element_text(hjust = 0, size=20),
#         axis.text.y = element_text(size=20),
#         axis.title.x = element_text(size=25),
#         axis.title.y = element_text(size=25))
# ggsave('fb_histogram.png')
#
#
# ggplot(data=raisers_no_suffixes, aes(x=rate, y=percent_opaque)) +
#   geom_point(size=4) +
#   xlab("Raising rate") +
#   ylab("Proportion of opaque tokens") +
#   theme(axis.text.x = element_text(size=20),
#         axis.text.y = element_text(size=20),
#         axis.title.x = element_text(size=25),
#         axis.title.y = element_text(size=25),
#         plot.title = element_text(size=30, hjust=0.5)) +
#   geom_smooth(method="glm",
#               se=FALSE)
# ggsave("raised_form_plot.png")
#
#
# raisers_no_suffixes <- droplevels(raisers_no_suffixes)
# foo <- aov(percent_opaque ~ log(total_count) + raised_form_prop, data=raisers_no_suffixes)
# summary(foo)
#
#
#
# shit <- raisers_no_suffixes[raisers_no_suffixes$unraised_tokens > 0,]
#
# cor.test(log(shit$unraised_tokens), shit$percent_opaque, method="pearson", na.action='na.omit')
#
# shit <- droplevels(shit)
# foo <- glm(percent_opaque ~ log(unraised_tokens) + raised_form_prop + LastTwo, data=shit)
# summary(foo)
#
# foo <- aov(percent_opaque ~  raised_form_prop * log(unraised_tokens), data=shit)
# summary(foo)
#
# ggplot(data=shit, aes(x=log(shit$unraised_tokens), y=percent_opaque)) +
#   geom_point(size=4) +
#   xlab("Log unraised token count") +
#   ylab("Proportion of opaque tokens") +
#   theme(axis.text.x = element_text(size=20),
#         axis.text.y = element_text(size=20),
#         axis.title.x = element_text(size=25),
#         axis.title.y = element_text(size=25),
#         plot.title = element_text(size=30, hjust=0.5)) +
#   geom_smooth(method="glm",
#               se=FALSE)
# ggsave("freq_plot.png")
#
#
#


vascillators <- raisers_no_suffixes[raisers_no_suffixes$percent_opaque < 1,]

cor.test(vascillators$percent_opaque, vascillators$raised_form_prop, method='pearson')

foo <- glm(percent_opaque ~ log(total_count) + raised_form_prop, data=vascillators)
summary(foo)

ggplot(data=vascillators, aes(x=raised_form_prop, y=percent_opaque)) +
  geom_point(size=4) +
  xlab("Proportion of raised observations") +
  ylab("Proportion of opaque tokens") +
  theme(axis.text.x = element_text(size=20),
        axis.text.y = element_text(size=20),
        axis.title.x = element_text(size=25),
        axis.title.y = element_text(size=25),
        plot.title = element_text(size=30, hjust=0.5)) +
  geom_smooth(method="glm",
              se=FALSE)
ggsave("raised_form_plot.png")

cor.test(vascillators$percent_opaque, vascillators$total_count, method='pearson')

ggplot(data=vascillators, aes(x=log(total_count), y=percent_opaque)) +
  geom_point(size=4) +
  xlab("Log root frequency") +
  ylab("Proportion of opaque tokens") +
  theme(axis.text.x = element_text(size=20),
        axis.text.y = element_text(size=20),
        axis.title.x = element_text(size=25),
        axis.title.y = element_text(size=25),
        plot.title = element_text(size=30, hjust=0.5)) +
  geom_smooth(method="glm",
              se=FALSE) +
ggsave("vascillators_freq_plot.png")

cor.test(vascillators$percent_opaque, vascillators$raised_form_prop, method='pearson')

# ggplot(data=vascillators, aes(x=raised_form_prop, y=percent_opaque)) +
#   geom_point(size=4) +
#   xlab("Proportion of raiased forms") +
#   ylab("Proportion of opaque tokens") +
#   theme(axis.text.x = element_text(size=20),
#         axis.text.y = element_text(size=20),
#         axis.title.x = element_text(size=25),
#         axis.title.y = element_text(size=25),
#         plot.title = element_text(size=30, hjust=0.5)) +
#   geom_smooth(method="glm",
#               se=FALSE) +
#   ggsave("vascillators_freq_plot.png")

# Fully indexed output
filename <- 'indexed_output.csv'
n_general_constraints <- 4
n_total_constraints <- n_general_constraints + nrow(raisers_no_suffixes_full)

headers <- c('', '', '', 'VAgreeBack', 'VAgreeFront', 'DoRaising', '*Sf')
for (i in 1:nrow(raisers_no_suffixes_full)) {
  headers <- c(headers, paste('SPU(VAgree)-', raisers_no_suffixes_full[i,]$Stem, sep=''))
}
write.table(matrix(headers, nrow=1), file=filename, row.names=FALSE, col.names = FALSE, sep=',')
write.table(matrix(headers, nrow=1), file=filename, row.names=FALSE, col.names = FALSE, sep=',', append=TRUE)
write.table(matrix(c('', '', '', rep(0, n_total_constraints)), nrow=1), file=filename, row.names=FALSE, col.names = FALSE, sep=',', append=TRUE)

for (i in 1:nrow(raisers_no_suffixes_full)) {
  row <- raisers_no_suffixes_full[i,]
  no_raise_f <- c(row$Stem, paste(row$Stem, '-F', sep=''))
  no_raise_b <- c(row$Stem, paste(row$Stem, '-B', sep=''))
  raise_f <- c(row$Stem, paste(row$Stem, '-RAISE-F', sep=''))
  raise_b <- c(row$Stem, paste(row$Stem, '-RAISE-B', sep=''))

  # Add frequencies
  no_raise_f <- c(no_raise_f, 0)
  no_raise_b <- c(no_raise_b, 0)
  raise_f <- c(raise_f, (row$n * (1 - row$percent_back)))
  raise_b <- c(raise_b, (row$n * row$percent_back))

  # Calculate VAgreeBack violations
  if (row$LastTwo == 'BB') {
    no_raise_f <- c(no_raise_f, 1)
    no_raise_b <- c(no_raise_b, '')
    raise_f <- c(raise_f, 1)
    raise_b <- c(raise_b, '')
  }
  if (row$LastTwo == 'FF') {
    no_raise_f <- c(no_raise_f, '')
    no_raise_b <- c(no_raise_b, '')
    raise_f <- c(raise_f, '')
    raise_b <- c(raise_b, '')
  }
  if (row$LastTwo == 'FB') {
    no_raise_f <- c(no_raise_f, 1)
    no_raise_b <- c(no_raise_b, '')
    raise_f <- c(raise_f, '')
    raise_b <- c(raise_b, '')
  }
  if (row$LastTwo == 'BF') {
    no_raise_f <- c(no_raise_f, '')
    no_raise_b <- c(no_raise_b, '')
    raise_f <- c(raise_f, 1)
    raise_b <- c(raise_b, '')
  }

  # Calculate VAgreeFront violations
  if (row$LastTwo == 'BB') {
    no_raise_f <- c(no_raise_f, '')
    no_raise_b <- c(no_raise_b, '')
    raise_f <- c(raise_f, '')
    raise_b <- c(raise_b, '')
  }
  if (row$LastTwo == 'FF') {
    no_raise_f <- c(no_raise_f, '')
    no_raise_b <- c(no_raise_b, 1)
    raise_f <- c(raise_f, '')
    raise_b <- c(raise_b, 1)
  }
  if (row$LastTwo == 'FB') {
    no_raise_f <- c(no_raise_f, '')
    no_raise_b <- c(no_raise_b, '')
    raise_f <- c(raise_f, '')
    raise_b <- c(raise_b, 1)
  }
  if (row$LastTwo == 'BF') {
    no_raise_f <- c(no_raise_f, '')
    no_raise_b <- c(no_raise_b, 1)
    raise_f <- c(raise_f, '')
    raise_b <- c(raise_b, '')
  }

  # Do raising violations
  no_raise_f <- c(no_raise_f, 1)
  raise_f <- c(raise_f, '')
  no_raise_b <- c(no_raise_b, 1)
  raise_b <- c(raise_b, '')

  # Do *Sf violations
  no_raise_f <- c(no_raise_f, 1)
  raise_f <- c(raise_f, 1)
  no_raise_b <- c(no_raise_b, 0)
  raise_b <- c(raise_b, 0)

  # Add empty spaces for previous constraints
  n_prev_constraints <- i - 1
  prev_vals <- rep('', n_prev_constraints)
  no_raise_f <- c(no_raise_f, prev_vals)
  no_raise_b <- c(no_raise_b, prev_vals)
  raise_f <- c(raise_f, prev_vals)
  raise_b <- c(raise_b, prev_vals)

  # Calculate SPU violation
  if (row$LastTwo == 'BB') {
    no_raise_f <- c(no_raise_f, '')
    no_raise_b <- c(no_raise_b, '')
    raise_f <- c(raise_f, 1)
    raise_b <- c(raise_b, '')
  }
  if (row$LastTwo == 'FF') {
    no_raise_f <- c(no_raise_f, '')
    no_raise_b <- c(no_raise_b, '')
    raise_f <- c(raise_f, '')
    raise_b <- c(raise_b, 1)
  }
  if (row$LastTwo == 'FB') {
    no_raise_f <- c(no_raise_f, '')
    no_raise_b <- c(no_raise_b, '')
    raise_f <- c(raise_f, 1)
    raise_b <- c(raise_b, '')
  }
  if (row$LastTwo == 'BF') {
    no_raise_f <- c(no_raise_f, '')
    no_raise_b <- c(no_raise_b, '')
    raise_f <- c(raise_f, '')
    raise_b <- c(raise_b, 1)
  }

  # Add empty spaces for following constraints
  n_following_constraints <- n_total_constraints - i
  following_vals <- rep('', n_following_constraints)
  no_raise_f <- c(no_raise_f, following_vals)
  no_raise_b <- c(no_raise_b, following_vals)
  raise_f <- c(raise_f, following_vals)
  raise_b <- c(raise_b, following_vals)

  write.table(matrix(no_raise_f, nrow=1), file=filename, row.names=FALSE, col.names = FALSE, append=TRUE, sep=',')
  write.table(matrix(no_raise_b, nrow=1), file=filename, row.names=FALSE, col.names = FALSE, append=TRUE, sep=',')
  write.table(matrix(raise_f, nrow=1), file=filename, row.names=FALSE, col.names = FALSE, append=TRUE, sep=',')
  write.table(matrix(raise_b, nrow=1), file=filename, row.names=FALSE, col.names = FALSE, append=TRUE, sep=',')
}

# Only conflicting forms indexed output
filename <- 'partially_indexed_output.csv'
n_general_constraints <- 4
n_indexed_constraints <- nrow(raisers_no_suffixes_full[raisers_no_suffixes_full$LastTwo == 'BF' | raisers_no_suffixes_full$LastTwo == 'FB',])
n_total_constraints <- n_general_constraints + n_indexed_constraints

headers <- c('', '', '', 'VAgreeBack', 'VAgreeFront', 'DoRaising', '*Sf')
for (i in 1:nrow(raisers_no_suffixes_full)) {
  if (raisers_no_suffixes_full[i,]$LastTwo %in% c('BF', 'FB')) {
    headers <- c(headers, paste('SPU(VAgree)-', raisers_no_suffixes_full[i,]$Stem, sep=''))
  }
}
write.table(matrix(headers, nrow=1), file=filename, row.names=FALSE, col.names = FALSE, sep=',')
write.table(matrix(headers, nrow=1), file=filename, row.names=FALSE, col.names = FALSE, sep=',', append=TRUE)
write.table(matrix(c('', '', '', rep(0, n_total_constraints)), nrow=1), file=filename, row.names=FALSE, col.names = FALSE, sep=',', append=TRUE)

normal_forms <- raisers_no_suffixes_full[raisers_no_suffixes_full$LastTwo %in% c('BB', 'FF'),]
conflict_forms <- raisers_no_suffixes_full[raisers_no_suffixes_full$LastTwo %in% c('BF', 'FB'),]

for (i in 1:nrow(normal_forms)) {
  row <- normal_forms[i,]
  no_raise_f <- c(row$Stem, paste(row$Stem, '-F', sep=''))
  no_raise_b <- c(row$Stem, paste(row$Stem, '-B', sep=''))
  raise_f <- c(row$Stem, paste(row$Stem, '-RAISE-F', sep=''))
  raise_b <- c(row$Stem, paste(row$Stem, '-RAISE-B', sep=''))

  # Add frequencies
  no_raise_f <- c(no_raise_f, 0)
  no_raise_b <- c(no_raise_b, 0)
  raise_f <- c(raise_f, (row$n * (1 - row$percent_back)))
  raise_b <- c(raise_b, (row$n * row$percent_back))

  # Calculate VAgreeBack violations
  if (row$LastTwo == 'BB') {
    no_raise_f <- c(no_raise_f, 1)
    no_raise_b <- c(no_raise_b, '')
    raise_f <- c(raise_f, 1)
    raise_b <- c(raise_b, '')
  }
  if (row$LastTwo == 'FF') {
    no_raise_f <- c(no_raise_f, '')
    no_raise_b <- c(no_raise_b, '')
    raise_f <- c(raise_f, '')
    raise_b <- c(raise_b, '')
  }
  if (row$LastTwo == 'FB') {
    no_raise_f <- c(no_raise_f, 1)
    no_raise_b <- c(no_raise_b, '')
    raise_f <- c(raise_f, '')
    raise_b <- c(raise_b, '')
  }
  if (row$LastTwo == 'BF') {
    no_raise_f <- c(no_raise_f, '')
    no_raise_b <- c(no_raise_b, '')
    raise_f <- c(raise_f, 1)
    raise_b <- c(raise_b, '')
  }

  # Calculate VAgreeFront violations
  if (row$LastTwo == 'BB') {
    no_raise_f <- c(no_raise_f, '')
    no_raise_b <- c(no_raise_b, '')
    raise_f <- c(raise_f, '')
    raise_b <- c(raise_b, '')
  }
  if (row$LastTwo == 'FF') {
    no_raise_f <- c(no_raise_f, '')
    no_raise_b <- c(no_raise_b, 1)
    raise_f <- c(raise_f, '')
    raise_b <- c(raise_b, 1)
  }
  if (row$LastTwo == 'FB') {
    no_raise_f <- c(no_raise_f, '')
    no_raise_b <- c(no_raise_b, '')
    raise_f <- c(raise_f, '')
    raise_b <- c(raise_b, 1)
  }
  if (row$LastTwo == 'BF') {
    no_raise_f <- c(no_raise_f, '')
    no_raise_b <- c(no_raise_b, 1)
    raise_f <- c(raise_f, '')
    raise_b <- c(raise_b, '')
  }

  # Do raising violations
  no_raise_f <- c(no_raise_f, 1)
  raise_f <- c(raise_f, '')
  no_raise_b <- c(no_raise_b, 1)
  raise_b <- c(raise_b, '')

  # Do *Sf violations
  no_raise_f <- c(no_raise_f, 1)
  raise_f <- c(raise_f, 1)
  no_raise_b <- c(no_raise_b, 0)
  raise_b <- c(raise_b, 0)

  # Add indexed constraint violations (blank)
  blank_vals <- rep('', n_indexed_constraints)
  no_raise_f <- c(no_raise_f, blank_vals)
  no_raise_b <- c(no_raise_b, blank_vals)
  raise_f <- c(raise_f, blank_vals)
  raise_b <- c(raise_b, blank_vals)

  write.table(matrix(no_raise_f, nrow=1), file=filename, row.names=FALSE, col.names = FALSE, append=TRUE, sep=',')
  write.table(matrix(no_raise_b, nrow=1), file=filename, row.names=FALSE, col.names = FALSE, append=TRUE, sep=',')
  write.table(matrix(raise_f, nrow=1), file=filename, row.names=FALSE, col.names = FALSE, append=TRUE, sep=',')
  write.table(matrix(raise_b, nrow=1), file=filename, row.names=FALSE, col.names = FALSE, append=TRUE, sep=',')
}


for (i in 1:nrow(conflict_forms)) {
  row <- conflict_forms[i,]
  no_raise_f <- c(row$Stem, paste(row$Stem, '-F', sep=''))
  no_raise_b <- c(row$Stem, paste(row$Stem, '-B', sep=''))
  raise_f <- c(row$Stem, paste(row$Stem, '-RAISE-F', sep=''))
  raise_b <- c(row$Stem, paste(row$Stem, '-RAISE-B', sep=''))

  # Add frequencies
  no_raise_f <- c(no_raise_f, 0)
  no_raise_b <- c(no_raise_b, 0)
  raise_f <- c(raise_f, (row$n * (1 - row$percent_back)))
  raise_b <- c(raise_b, (row$n * row$percent_back))

  # Calculate VAgreeBack violations
  if (row$LastTwo == 'BB') {
    no_raise_f <- c(no_raise_f, 1)
    no_raise_b <- c(no_raise_b, '')
    raise_f <- c(raise_f, 1)
    raise_b <- c(raise_b, '')
  }
  if (row$LastTwo == 'FF') {
    no_raise_f <- c(no_raise_f, '')
    no_raise_b <- c(no_raise_b, '')
    raise_f <- c(raise_f, '')
    raise_b <- c(raise_b, '')
  }
  if (row$LastTwo == 'FB') {
    no_raise_f <- c(no_raise_f, 1)
    no_raise_b <- c(no_raise_b, '')
    raise_f <- c(raise_f, '')
    raise_b <- c(raise_b, '')
  }
  if (row$LastTwo == 'BF') {
    no_raise_f <- c(no_raise_f, '')
    no_raise_b <- c(no_raise_b, '')
    raise_f <- c(raise_f, 1)
    raise_b <- c(raise_b, '')
  }

  # Calculate VAgreeFront violations
  if (row$LastTwo == 'BB') {
    no_raise_f <- c(no_raise_f, '')
    no_raise_b <- c(no_raise_b, '')
    raise_f <- c(raise_f, '')
    raise_b <- c(raise_b, '')
  }
  if (row$LastTwo == 'FF') {
    no_raise_f <- c(no_raise_f, '')
    no_raise_b <- c(no_raise_b, 1)
    raise_f <- c(raise_f, '')
    raise_b <- c(raise_b, 1)
  }
  if (row$LastTwo == 'FB') {
    no_raise_f <- c(no_raise_f, '')
    no_raise_b <- c(no_raise_b, '')
    raise_f <- c(raise_f, '')
    raise_b <- c(raise_b, 1)
  }
  if (row$LastTwo == 'BF') {
    no_raise_f <- c(no_raise_f, '')
    no_raise_b <- c(no_raise_b, 1)
    raise_f <- c(raise_f, '')
    raise_b <- c(raise_b, '')
  }

  # Do raising violations
  no_raise_f <- c(no_raise_f, 1)
  raise_f <- c(raise_f, '')
  no_raise_b <- c(no_raise_b, 1)
  raise_b <- c(raise_b, '')

  # Do *Sf violations
  no_raise_f <- c(no_raise_f, 1)
  raise_f <- c(raise_f, 1)
  no_raise_b <- c(no_raise_b, 0)
  raise_b <- c(raise_b, 0)

  # Add empty spaces for previous constraints
  n_prev_constraints <- i - 1
  prev_vals <- rep('', n_prev_constraints)
  no_raise_f <- c(no_raise_f, prev_vals)
  no_raise_b <- c(no_raise_b, prev_vals)
  raise_f <- c(raise_f, prev_vals)
  raise_b <- c(raise_b, prev_vals)

  # Calculate SPU violation
  if (row$LastTwo == 'FB') {
    no_raise_f <- c(no_raise_f, '')
    no_raise_b <- c(no_raise_b, '')
    raise_f <- c(raise_f, 1)
    raise_b <- c(raise_b, '')
  }
  if (row$LastTwo == 'BF') {
    no_raise_f <- c(no_raise_f, '')
    no_raise_b <- c(no_raise_b, '')
    raise_f <- c(raise_f, '')
    raise_b <- c(raise_b, 1)
  }

  # Add empty spaces for following constraints
  n_following_constraints <- n_indexed_constraints - i
  following_vals <- rep('', n_following_constraints)
  no_raise_f <- c(no_raise_f, following_vals)
  no_raise_b <- c(no_raise_b, following_vals)
  raise_f <- c(raise_f, following_vals)
  raise_b <- c(raise_b, following_vals)

  write.table(matrix(no_raise_f, nrow=1), file=filename, row.names=FALSE, col.names = FALSE, append=TRUE, sep=',')
  write.table(matrix(no_raise_b, nrow=1), file=filename, row.names=FALSE, col.names = FALSE, append=TRUE, sep=',')
  write.table(matrix(raise_f, nrow=1), file=filename, row.names=FALSE, col.names = FALSE, append=TRUE, sep=',')
  write.table(matrix(raise_b, nrow=1), file=filename, row.names=FALSE, col.names = FALSE, append=TRUE, sep=',')
}


# Non-indexed version (surface-true)
filename <- 'nonindexed_output.csv'

headers <- c('', '', '', 'VAgreeBack', 'VAgreeFront', 'DoRaising', '*Sf', 'SPU(VAgree)')
write.table(matrix(headers, nrow=1), file=filename, row.names=FALSE, col.names = FALSE, sep=',')
write.table(matrix(headers, nrow=1), file=filename, row.names=FALSE, col.names = FALSE, sep=',', append=TRUE)
write.table(matrix(c('', '', '', rep(0, 5)), nrow=1), file=filename, row.names=FALSE, col.names = FALSE, sep=',', append=TRUE)

for (i in 1:nrow(raisers_no_suffixes_full)) {
  row <- raisers_no_suffixes_full[i,]
  no_raise_f <- c(row$Stem, paste(row$Stem, '-F', sep=''))
  no_raise_b <- c(row$Stem, paste(row$Stem, '-B', sep=''))
  raise_f <- c(row$Stem, paste(row$Stem, '-RAISE-F', sep=''))
  raise_b <- c(row$Stem, paste(row$Stem, '-RAISE-B', sep=''))

  # Add frequencies
  no_raise_f <- c(no_raise_f, 0)
  no_raise_b <- c(no_raise_b, 0)
  raise_f <- c(raise_f, (row$n * (1 - row$percent_back)))
  raise_b <- c(raise_b, (row$n * row$percent_back))

  # Calculate VAgreeBack violations
  if (row$LastTwo == 'BB') {
    no_raise_f <- c(no_raise_f, 1)
    no_raise_b <- c(no_raise_b, '')
    raise_f <- c(raise_f, 1)
    raise_b <- c(raise_b, '')
  }
  if (row$LastTwo == 'FF') {
    no_raise_f <- c(no_raise_f, '')
    no_raise_b <- c(no_raise_b, '')
    raise_f <- c(raise_f, '')
    raise_b <- c(raise_b, '')
  }
  if (row$LastTwo == 'FB') {
    no_raise_f <- c(no_raise_f, 1)
    no_raise_b <- c(no_raise_b, '')
    raise_f <- c(raise_f, '')
    raise_b <- c(raise_b, '')
  }
  if (row$LastTwo == 'BF') {
    no_raise_f <- c(no_raise_f, '')
    no_raise_b <- c(no_raise_b, '')
    raise_f <- c(raise_f, 1)
    raise_b <- c(raise_b, '')
  }

  # Calculate VAgreeFront violations
  if (row$LastTwo == 'BB') {
    no_raise_f <- c(no_raise_f, '')
    no_raise_b <- c(no_raise_b, '')
    raise_f <- c(raise_f, '')
    raise_b <- c(raise_b, '')
  }
  if (row$LastTwo == 'FF') {
    no_raise_f <- c(no_raise_f, '')
    no_raise_b <- c(no_raise_b, 1)
    raise_f <- c(raise_f, '')
    raise_b <- c(raise_b, 1)
  }
  if (row$LastTwo == 'FB') {
    no_raise_f <- c(no_raise_f, '')
    no_raise_b <- c(no_raise_b, '')
    raise_f <- c(raise_f, '')
    raise_b <- c(raise_b, 1)
  }
  if (row$LastTwo == 'BF') {
    no_raise_f <- c(no_raise_f, '')
    no_raise_b <- c(no_raise_b, 1)
    raise_f <- c(raise_f, '')
    raise_b <- c(raise_b, '')
  }

  # Do raising violations
  no_raise_f <- c(no_raise_f, 1)
  raise_f <- c(raise_f, '')
  no_raise_b <- c(no_raise_b, 1)
  raise_b <- c(raise_b, '')

  # Do *Sf violations
  no_raise_f <- c(no_raise_f, 1)
  raise_f <- c(raise_f, 1)
  no_raise_b <- c(no_raise_b, 0)
  raise_b <- c(raise_b, 0)

  # Calculate SPU violation
  if (row$LastTwo == 'BB') {
    no_raise_f <- c(no_raise_f, '')
    no_raise_b <- c(no_raise_b, '')
    raise_f <- c(raise_f, 1)
    raise_b <- c(raise_b, '')
  }
  if (row$LastTwo == 'FF') {
    no_raise_f <- c(no_raise_f, '')
    no_raise_b <- c(no_raise_b, '')
    raise_f <- c(raise_f, '')
    raise_b <- c(raise_b, 1)
  }
  if (row$LastTwo == 'FB') {
    no_raise_f <- c(no_raise_f, '')
    no_raise_b <- c(no_raise_b, '')
    raise_f <- c(raise_f, 1)
    raise_b <- c(raise_b, '')
  }
  if (row$LastTwo == 'BF') {
    no_raise_f <- c(no_raise_f, '')
    no_raise_b <- c(no_raise_b, '')
    raise_f <- c(raise_f, '')
    raise_b <- c(raise_b, 1)
  }

  write.table(matrix(no_raise_f, nrow=1), file=filename, row.names=FALSE, col.names = FALSE, append=TRUE, sep=',')
  write.table(matrix(no_raise_b, nrow=1), file=filename, row.names=FALSE, col.names = FALSE, append=TRUE, sep=',')
  write.table(matrix(raise_f, nrow=1), file=filename, row.names=FALSE, col.names = FALSE, append=TRUE, sep=',')
  write.table(matrix(raise_b, nrow=1), file=filename, row.names=FALSE, col.names = FALSE, append=TRUE, sep=',')
}

# write.table(matrix(headers, nrow=1), file='test_output.csv', row.names=FALSE, col.names = FALSE, sep=',')
# write.table(matrix(headers, nrow=1), file='test_output.csv', row.names=FALSE, col.names = FALSE, sep=',', append=TRUE)
# write.table(matrix(c('', '', '', 0, 0, 0, 0, 0), nrow=1), file='test_output.csv', row.names=FALSE, col.names = FALSE, sep=',', append=TRUE)
#
# for (i in 1:nrow(raisers_no_suffixes_full)) {
#   row <- raisers_no_suffixes_full[i,]
#   no_raise_f <- c(row$Stem, paste(row$Stem, '-F', sep=''))
#   no_raise_b <- c(row$Stem, paste(row$Stem, '-B', sep=''))
#   raise_f <- c(row$Stem, paste(row$Stem, '-RAISE-F', sep=''))
#   raise_b <- c(row$Stem, paste(row$Stem, '-RAISE-B', sep=''))
#
#   # Add frequencies
#   no_raise_f <- c(no_raise_f, 0)
#   no_raise_b <- c(no_raise_b, 0)
#   raise_f <- c(raise_f, (row$n * (1 - row$percent_back)))
#   raise_b <- c(raise_b, (row$n * row$percent_back))
#
#   # Calculate VAgreeBack violations
#   if (row$LastTwo == 'BB') {
#     no_raise_f <- c(no_raise_f, 1)
#     no_raise_b <- c(no_raise_b, '')
#     raise_f <- c(raise_f, 1)
#     raise_b <- c(raise_b, '')
#   }
#   if (row$LastTwo == 'FF') {
#     no_raise_f <- c(no_raise_f, '')
#     no_raise_b <- c(no_raise_b, '')
#     raise_f <- c(raise_f, '')
#     raise_b <- c(raise_b, '')
#   }
#   if (row$LastTwo == 'FB') {
#     no_raise_f <- c(no_raise_f, 1)
#     no_raise_b <- c(no_raise_b, '')
#     raise_f <- c(raise_f, '')
#     raise_b <- c(raise_b, '')
#   }
#   if (row$LastTwo == 'BF') {
#     no_raise_f <- c(no_raise_f, '')
#     no_raise_b <- c(no_raise_b, '')
#     raise_f <- c(raise_f, 1)
#     raise_b <- c(raise_b, '')
#   }
#
#   # Calculate VAgreeFront violations
#   if (row$LastTwo == 'BB') {
#     no_raise_f <- c(no_raise_f, '')
#     no_raise_b <- c(no_raise_b, '')
#     raise_f <- c(raise_f, '')
#     raise_b <- c(raise_b, '')
#   }
#   if (row$LastTwo == 'FF') {
#     no_raise_f <- c(no_raise_f, '')
#     no_raise_b <- c(no_raise_b, 1)
#     raise_f <- c(raise_f, '')
#     raise_b <- c(raise_b, 1)
#   }
#   if (row$LastTwo == 'FB') {
#     no_raise_f <- c(no_raise_f, '')
#     no_raise_b <- c(no_raise_b, '')
#     raise_f <- c(raise_f, '')
#     raise_b <- c(raise_b, 1)
#   }
#   if (row$LastTwo == 'BF') {
#     no_raise_f <- c(no_raise_f, '')
#     no_raise_b <- c(no_raise_b, 1)
#     raise_f <- c(raise_f, '')
#     raise_b <- c(raise_b, '')
#   }
#
#   # Calculate SPU(back) violations
#   # unraised_tokens <- max(log(row$total_count * (1 - row$raised_form_prop)), 0)
#   # if (is.na(unraised_tokens)) {
#   #   unraised_tokens <- 0
#   # }
#   #unraised_tokens <- max(row$total_count * (1 - row$raised_form_prop), 0) / 10000
#   unraised_tokens <- 1 - row$raised_form_prop
#   if (row$LastTwo == 'BB') {
#     no_raise_f <- c(no_raise_f, '')
#     no_raise_b <- c(no_raise_b, '')
#     raise_f <- c(raise_f, unraised_tokens)
#     raise_b <- c(raise_b, '')
#   }
#   if (row$LastTwo == 'FF') {
#     no_raise_f <- c(no_raise_f, '')
#     no_raise_b <- c(no_raise_b, '')
#     raise_f <- c(raise_f, '')
#     raise_b <- c(raise_b, unraised_tokens)
#   }
#   if (row$LastTwo == 'FB') {
#     no_raise_f <- c(no_raise_f, '')
#     no_raise_b <- c(no_raise_b, '')
#     raise_f <- c(raise_f, unraised_tokens)
#     raise_b <- c(raise_b, '')
#   }
#   if (row$LastTwo == 'BF') {
#     no_raise_f <- c(no_raise_f, '')
#     no_raise_b <- c(no_raise_b, '')
#     raise_f <- c(raise_f, '')
#     raise_b <- c(raise_b, unraised_tokens)
#   }
#
#   # Do raising violations
#   no_raise_f <- c(no_raise_f, 1)
#   raise_f <- c(raise_f, '')
#   no_raise_b <- c(no_raise_b, 1)
#   raise_b <- c(raise_b, '')
#
#   # Do *Sf violations
#   no_raise_f <- c(no_raise_f, 1)
#   raise_f <- c(raise_f, 1)
#   no_raise_b <- c(no_raise_b, 0)
#   raise_b <- c(raise_b, 0)
#
#   if (row$raised_form_prop != 1) {
#
#     write.table(matrix(no_raise_f, nrow=1), file='test_output.csv', row.names=FALSE, col.names = FALSE, append=TRUE, sep=',')
#     write.table(matrix(no_raise_b, nrow=1), file='test_output.csv', row.names=FALSE, col.names = FALSE, append=TRUE, sep=',')
#     write.table(matrix(raise_f, nrow=1), file='test_output.csv', row.names=FALSE, col.names = FALSE, append=TRUE, sep=',')
#     write.table(matrix(raise_b, nrow=1), file='test_output.csv', row.names=FALSE, col.names = FALSE, append=TRUE, sep=',')
#   }
# }


# OLD VERSION
# headers <- c('', '', '', 'VAgreeBack', 'VAgreeFront', 'SPU(back)', 'DoRaising', '*Sf')
# write.table(matrix(headers, nrow=1), file='test_output.csv', row.names=FALSE, col.names = FALSE, sep=',')
# write.table(matrix(headers, nrow=1), file='test_output.csv', row.names=FALSE, col.names = FALSE, sep=',', append=TRUE)
# write.table(matrix(c('', '', '', 0, 0, 0, 0, 0), nrow=1), file='test_output.csv', row.names=FALSE, col.names = FALSE, sep=',', append=TRUE)
#
# for (i in 1:nrow(raisers_no_suffixes_full)) {
#   row <- raisers_no_suffixes_full[i,]
#   no_raise_f <- c(row$Stem, paste(row$Stem, '-F', sep=''))
#   no_raise_b <- c(row$Stem, paste(row$Stem, '-B', sep=''))
#   raise_f <- c(row$Stem, paste(row$Stem, '-RAISE-F', sep=''))
#   raise_b <- c(row$Stem, paste(row$Stem, '-RAISE-B', sep=''))
#
#   # Add frequencies
#   no_raise_f <- c(no_raise_f, 0)
#   no_raise_b <- c(no_raise_b, 0)
#   raise_f <- c(raise_f, (row$n * (1 - row$percent_back)))
#   raise_b <- c(raise_b, (row$n * row$percent_back))
#
#   # Calculate VAgreeBack violations
#   if (row$LastTwo == 'BB') {
#     no_raise_f <- c(no_raise_f, 1)
#     no_raise_b <- c(no_raise_b, '')
#     raise_f <- c(raise_f, 1)
#     raise_b <- c(raise_b, '')
#   }
#   if (row$LastTwo == 'FF') {
#     no_raise_f <- c(no_raise_f, '')
#     no_raise_b <- c(no_raise_b, '')
#     raise_f <- c(raise_f, '')
#     raise_b <- c(raise_b, '')
#   }
#   if (row$LastTwo == 'FB') {
#     no_raise_f <- c(no_raise_f, 1)
#     no_raise_b <- c(no_raise_b, '')
#     raise_f <- c(raise_f, '')
#     raise_b <- c(raise_b, '')
#   }
#   if (row$LastTwo == 'BF') {
#     no_raise_f <- c(no_raise_f, '')
#     no_raise_b <- c(no_raise_b, '')
#     raise_f <- c(raise_f, 1)
#     raise_b <- c(raise_b, '')
#   }
#
#   # Calculate VAgreeFront violations
#   if (row$LastTwo == 'BB') {
#     no_raise_f <- c(no_raise_f, '')
#     no_raise_b <- c(no_raise_b, '')
#     raise_f <- c(raise_f, '')
#     raise_b <- c(raise_b, '')
#   }
#   if (row$LastTwo == 'FF') {
#     no_raise_f <- c(no_raise_f, '')
#     no_raise_b <- c(no_raise_b, 1)
#     raise_f <- c(raise_f, '')
#     raise_b <- c(raise_b, 1)
#   }
#   if (row$LastTwo == 'FB') {
#     no_raise_f <- c(no_raise_f, '')
#     no_raise_b <- c(no_raise_b, '')
#     raise_f <- c(raise_f, '')
#     raise_b <- c(raise_b, 1)
#   }
#   if (row$LastTwo == 'BF') {
#     no_raise_f <- c(no_raise_f, '')
#     no_raise_b <- c(no_raise_b, 1)
#     raise_f <- c(raise_f, '')
#     raise_b <- c(raise_b, '')
#   }
#
#   # Calculate SPU(back) violations
#   # unraised_tokens <- max(log(row$total_count * (1 - row$raised_form_prop)), 0)
#   # if (is.na(unraised_tokens)) {
#   #   unraised_tokens <- 0
#   # }
#   #unraised_tokens <- max(row$total_count * (1 - row$raised_form_prop), 0) / 10000
#   unraised_tokens <- 1 - row$raised_form_prop
#   if (row$LastTwo == 'BB') {
#     no_raise_f <- c(no_raise_f, '')
#     no_raise_b <- c(no_raise_b, '')
#     raise_f <- c(raise_f, unraised_tokens)
#     raise_b <- c(raise_b, '')
#   }
#   if (row$LastTwo == 'FF') {
#     no_raise_f <- c(no_raise_f, '')
#     no_raise_b <- c(no_raise_b, '')
#     raise_f <- c(raise_f, '')
#     raise_b <- c(raise_b, unraised_tokens)
#   }
#   if (row$LastTwo == 'FB') {
#     no_raise_f <- c(no_raise_f, '')
#     no_raise_b <- c(no_raise_b, '')
#     raise_f <- c(raise_f, unraised_tokens)
#     raise_b <- c(raise_b, '')
#   }
#   if (row$LastTwo == 'BF') {
#     no_raise_f <- c(no_raise_f, '')
#     no_raise_b <- c(no_raise_b, '')
#     raise_f <- c(raise_f, '')
#     raise_b <- c(raise_b, unraised_tokens)
#   }
#
#   # Do raising violations
#   no_raise_f <- c(no_raise_f, 1)
#   raise_f <- c(raise_f, '')
#   no_raise_b <- c(no_raise_b, 1)
#   raise_b <- c(raise_b, '')
#
#   # Do *Sf violations
#   no_raise_f <- c(no_raise_f, 1)
#   raise_f <- c(raise_f, 1)
#   no_raise_b <- c(no_raise_b, 0)
#   raise_b <- c(raise_b, 0)
#
#   if (row$raised_form_prop != 1) {
#
#     write.table(matrix(no_raise_f, nrow=1), file='test_output.csv', row.names=FALSE, col.names = FALSE, append=TRUE, sep=',')
#     write.table(matrix(no_raise_b, nrow=1), file='test_output.csv', row.names=FALSE, col.names = FALSE, append=TRUE, sep=',')
#     write.table(matrix(raise_f, nrow=1), file='test_output.csv', row.names=FALSE, col.names = FALSE, append=TRUE, sep=',')
#     write.table(matrix(raise_b, nrow=1), file='test_output.csv', row.names=FALSE, col.names = FALSE, append=TRUE, sep=',')
#   }
# }
