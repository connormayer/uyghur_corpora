library(data.table)
library(lme4)
library(tidyverse)

load_corpus <- function(path) {
  corpus_data <- fread(path, encoding="UTF-8")
  corpus_data$Front.Count <- as.numeric(as.character(corpus_data$Front.Count))
  corpus_data$Back.Count <- as.numeric(as.character(corpus_data$Back.Count))
  return(corpus_data)
}

#data <- load_corpus('C://Users/conno/Dropbox/ling/Dissertation/code/processed_corpora/harmonic_raisers.csv')
data <- read_csv('C://Users/conno/Dropbox/ling/Dissertation/code/processed_corpora/harmonic_raisers.csv')
# data <- data[, V1 := NULL]
# data <- data[, V1 := NULL]
# data <- data[data$LastTwo == 'BF' | data$LastTwo == 'FB',]
# data$opaque <- 0
# data[Front.Count > 0 &  LastOne == 'F', opaque := 1]
# data[Back.Count > 0 & LastOne == 'B', opaque := 1]
# data <- data[data$HasChe == FALSE,]
# data <- data[data$HasYe == FALSE,]
# data <- data[data$HasAne == FALSE,]
# data <- data[data$HasName == FALSE,]
#
# data <- merge(data, raised_forms, by="Stem")
#
# write.table(data, file="harmonic_raisers.csv", row.names=FALSE, sep=",")
# Hacky way to transform frequencies to per million
data <- data %>%
  mutate(log_total_count=log(total_count) - log(15.7))

foo <- glm(opaque ~ log_total_count + raised_form_prop + LastTwo, data=data, family="binomial")
foo2 <- glmer(opaque ~ log_total_count + raised_form_prop + LastTwo + (1|Stem) + (1|Source), data=data, family="binomial")
summary(foo2)

data %>%
  group_by(LastTwo, log_total_count, raised_form_prop) %>% 
  summarize(opaque=mean(opaque)) %>%
  filter(opaque < 0.99) %>%
  ggplot(aes(x=raised_form_prop, y=opaque)) +
  geom_point()

percent_data <- data %>% group_by(LastTwo) %>%
  summarise(opaque_harmony = mean(opaque),
            transparent_harmony = 1 - opaque_harmony)

token_data <- data %>% group_by(LastTwo) %>%
  summarise(opaque_harmony = round(mean(opaque) * n()),
            transparent_harmony = round((1 - mean(opaque)) * n())
  )
percent_data <- melt(as.data.frame(percent_data), idvar="LastTwo")
token_data <- melt(as.data.frame(token_data), idvar="LastTwo")
graph_data <- merge(percent_data, token_data, by=c("LastTwo", "variable"))

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
  scale_fill_discrete(labels=c("Opaque harmony", "Transparent harmony"))
  #ggtitle("Harmonizing behavior in raised tokens with harmonizing suffixes")
  ggsave("full_harmonic_raisers.png")
