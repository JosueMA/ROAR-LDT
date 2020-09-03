library(mirt)
library(dplyr)
library(ggplot2)
library(gridExtra)
library(cowplot)
library(corrplot)
library(RColorBrewer)
library(stringr)
library(psych)
library(viridis)
library(wesanderson)

setwd('~/git/ROAR-LDT/Study2/')

# Read data
dfv1 <- read.csv('~/git/ROAR-LDT/data_allsubs/LDT_alldata_wide.csv')
dfv2 <- read.csv('~/git/ROAR-LDT/Study2/data/LDT_alldata_wide_v2.csv')
df <- full_join(dfv1,dfv2)
metadata <- read.csv('~/git/ROAR-LDT/Study2/data/metadata_all_roundeddates.csv')
sub.summary <- read.csv('~/git/ROAR-LDT/Study2/data/LDT_summarymeasures_wide_v1v2.csv')
metadata <- rename(metadata,subj = record_id)

# Fit IRT model and compute ability
m1 <- mirt(select(df,-subj), model = 1, itemtype = 'Rasch', guess=0.5, technical=list(NCYCLES=5000))
# Just for version 2
m2 <- mirt(select(dfv2,-subj), model = 1, itemtype = 'Rasch', guess=0.5, technical=list(NCYCLES=5000))

irt.estimates <- select(df,subj)
irt.estimates$theta1pl <- fscores(m1)
# Note whether they are V1 or V2
irt.estimates$version <- factor(is.na(df$throomba),levels = c(TRUE,FALSE), labels = c('v1','v2'))

# Compute percent correct
irt.estimates$ncor <- rowSums(select(df,-subj),na.rm=TRUE)

# Filter out V1 for the time being
irt.estimates <- filter(irt.estimates,version=='v2')

# Join metadata and irt estiates
sub.data <- left_join(irt.estimates, metadata)
# Get data without nan
sub.data <- filter(sub.data,!is.na(wj_lwid_raw))
# Remove subjects with very old scores
sub.data <- filter(sub.data,MonthsSinceTesting<11)


## Figure
# Set color limits
clims = c(80, 120)
sub.data$wj_clip <- sub.data$wj_lwid_ss
sub.data$wj_clip[sub.data$wj_clip<clims[1]]<-clims[1]
sub.data$wj_clip[sub.data$wj_clip>clims[2]]<-clims[2]
agerange <- c(min(sub.data$visit_age/12),max(sub.data$visit_age/12))
p1 <- ggplot(sub.data, aes(x=ncor, y=wj_lwid_raw)) +
  geom_point(aes(colour=agenow/12),size=3,alpha=1) + stat_smooth(method="lm", se=TRUE, color='gray30') + xlab('Number correct (out of 252)') +
  scale_color_gradientn(colours = c( 'dodgerblue1','firebrick1','goldenrod1')) + 
  labs(colour = "Age") + theme(legend.position = "right")+
  geom_label(x=min(sub.data$ncor),y=max(sub.data$wj_lwid_raw)-1, label=sprintf('r = %.2f',cor(select(sub.data, ncor,wj_lwid_raw))[1,2]),hjust=0, vjust=0,size=3)
p1
p2 <- ggplot(sub.data, aes(x=theta1pl, y=wj_lwid_raw)) +
  geom_point(aes(colour=agenow/12),size=3,alpha=1) + stat_smooth(method="lm", se=TRUE, color='gray30') + xlab('Ability estimate (1PL model)') +
  # scale_color_viridis(option='viridis',direction=-1)+
  scale_color_gradientn(colours = wes_palette(n=5, name="Zissou1"),limits = agerange) + 
  #scale_color_gradientn(colours = c( 'dodgerblue1','firebrick1','goldenrod1')) + 
  labs(colour = "Age") + theme(legend.position = "right")+
  geom_label(x=min(sub.data$theta1pl),y=max(sub.data$wj_lwid_raw)-1, label=sprintf('r = %.2f',cor(select(sub.data, theta1pl,wj_lwid_raw))[1,2]),hjust=0, vjust=0,size=3)
p2
grid.arrange(p1,p2,nrow=1)
g = arrangeGrob(p1,p2,nrow=1)
ggsave('ROAR-LDT_v2.pdf',g,width=6, height=3)

# Now make figures just for the young children
sub.data <- filter(sub.data,agenow<=84)
p3 <- ggplot(sub.data, aes(x=ncor, y=wj_lwid_raw)) +
  geom_point(aes(colour=agenow/12),size=3,alpha=1) + stat_smooth(method="lm", se=TRUE, color='gray30') + xlab('Number correct (out of 252)') +
  scale_color_gradientn(colours = c( 'dodgerblue1','firebrick1','goldenrod1'),limits = agerange) + 
  labs(colour = "Age") + theme(legend.position = "right")+
  geom_label(x=min(sub.data$ncor),y=max(sub.data$wj_lwid_raw)-1, label=sprintf('r = %.2f',cor(select(sub.data, ncor,wj_lwid_raw))[1,2]),hjust=0, vjust=0,size=3)
p3
p4 <- ggplot(sub.data, aes(x=theta1pl, y=wj_lwid_raw)) +
  geom_point(aes(colour=agenow/12),size=3,alpha=1) + stat_smooth(method="lm", se=TRUE, color='gray30') + xlab('Ability estimate (1PL model)') +
  #scale_color_gradientn(colours = c( 'dodgerblue1','firebrick1','goldenrod1'),limits = agerange) + 
  scale_color_gradientn(colours = wes_palette(n=5, name="Zissou1"),limits = agerange) + 
  labs(colour = "Age") + theme(legend.position = "right")+
  geom_label(x=min(sub.data$theta1pl),y=max(sub.data$wj_lwid_raw)-1, label=sprintf('r = %.2f',cor(select(sub.data, theta1pl,wj_lwid_raw))[1,2]),hjust=0, vjust=0,size=3)
p4
grid.arrange(p1,p2,p3,p4,nrow=2)
g = arrangeGrob(p1,p2,p3,p4,nrow=2)
ggsave('ROAR-LDT_v2_2rows.pdf',g,width=6, height=5)

## Analyze performance on Study 2 words

v2words = c('sit','listen', 'lunch', 'night', 'tis', 'stenil', 'nulch', 'ginth',
            'fun' ,'cold', 'bathroom','teacher','nuf', 'dolc', 'throomba', 'chareet',
            'hello', 'name', 'good', 'hungry', 'loleh', 'eamn', 'dogo', 'gurynh')
df.v2words <- select(df, all_of(c('subj',v2words)))
df.v2words.1st <- filter(df.v2words, subj %in% sub.data$subj)
print(colSums(df.v2words.1st)/dim(df.v2words.1st)[1])

# Make figre looking at item difficulty and item correlation with overal performance
item.cors = as.data.frame(cor(select(dfv2,-subj),rowSums(select(dfv2,-subj),dim=1)/(dim(dfv2)[2]-1)))
item.cors <- rename(item.cors,r.performance = V1)
item.cors$pcor <- colSums(select(dfv2,-subj))/(dim(dfv2)[1])
item.cors$items <- row.names(item.cors)
item.cors$v2newwords <- item.cors$items %in% v2words
# Scatter plot of corr values
gt1 = ggplot(item.cors,aes(x=r.performance,y=1-pcor,label=items,colour=v2newwords)) +
  geom_text(size=2,fontface='bold',alpha=0.7)+
  scale_color_manual(values = c('deepskyblue1','firebrick1'))+
  
  theme(legend.position='none') + 
  xlab('Item correlation with overall proportion correct') + ylab('Item difficulty (1 - proportion correct)')
gt1
ggsave('itemcors.pdf',gt1,width=5,height=5)
