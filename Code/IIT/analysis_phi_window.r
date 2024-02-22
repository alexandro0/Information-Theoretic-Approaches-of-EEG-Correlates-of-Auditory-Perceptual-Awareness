################################################################################
# SCRIPT FOR PLOTTING AND ANALYSES WINDOWED PHI DATA
################################################################################

# CTRL + ALT + R

################################################################################
# LIBRARY LOAD
################################################################################

library(ggplot2)
library(ggsignif)
library(ggpubr)
library(doBy)
library(pracma)
library(lme4)
library(nlme)
library(lmerTest)
library(emmeans)
library(multcomp)
library(doBy)
library(stringr)
library(report)
library(xtable)

################################################################################
# FUNCTION
################################################################################

summarySE <- function(data=NULL, measurevar, groupvars=NULL, na.rm=FALSE, 
                      conf.interval=.95) {
  library(doBy)
  # New version of length which can handle NA's: if na.rm==T, don't count them
  length2 <- function (x, na.rm=FALSE) {
    if (na.rm) sum(!is.na(x))
    else       length(x)
  }
  # Collapse the data
  formula <- as.formula(paste(measurevar, paste(groupvars, collapse=" + "), 
                              sep=" ~ "))
  datac <- summaryBy(formula, data=data, FUN=c(length2,mean,sd), na.rm=na.rm)
  # Rename columns
  names(datac)[ names(datac) == paste(measurevar, ".mean", sep="")] = measurevar
  names(datac)[ names(datac) == paste(measurevar, ".sd", sep="")] = "sd"
  names(datac)[ names(datac) == paste(measurevar, ".length2", sep="")] = "N"
  datac$se <- datac$sd / sqrt(datac$N)  # Calculate standard error of the mean
  # Confidence interval multiplier for standard error
  # Calculate t-statistic for confidence interval: 
  # e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
  ciMult <- qt(conf.interval/2 + .5, datac$N-1)
  datac$ci <- datac$se * ciMult
  return(datac)
}

################################################################################
# DATA LOAD
################################################################################

setwd("/home/link/Documents/thèse_onera/experimentation/eeg/manipeeg/results/")
setwd("EEG/data_for_stats/hits_miss_around_detection/125Hz/")
df = read.table("phi_sagittal_measures.csv", header=TRUE, sep=',', dec=".")

df$Condition = df$Time
df$Condition[which(df$Condition > 0.)] = "After"
df$Condition[which(df$Condition < 0.)] = "Before"
df$Condition[which(df$Condition == 0.)] = "Before"

################################################################################
# REFERENCING FOR PLOTTING
################################################################################

root_path = "/home/link/Documents/thèse_onera"
these_path = "thèse_alex/Figures/illustrations/Exp_EEG/Inf_Int"
article_path = "/home/link/Documents/thèse_onera/articles_alex/"
slides_path = "/home/link/Documents/thèse_onera/diapos_phd_thesis/"
slides_path_acc = "images/EEG/info_int/"
slides_path_final = str_c(slides_path, slides_path_acc, sep="")

################################################################################
# PARAMETERS FOR PLOTTING
################################################################################

theme_set(theme_bw())

pd_eb = position_dodge(0.)
sizepoint_eb = 6.
ebwidth_eb = .6
cexx_eb = 0.6

pd_w = position_dodge(0.)
sizepoint_w = 6.
ebwidth_w = .6
cexx_w = 0.6

tt = element_text(size=50, family="Times New Roman", color="black", face="bold")
atx = element_text(color="black", size=15, angle=90, hjust=1, vjust=0)
atx_diff = element_text(color="black", size=35, angle=90, hjust=1, vjust=0)

################################################################################
# FACTOR
################################################################################

df$Sujet = as.factor(df$Sujet)
df$Detection = as.factor(df$Detection)
df$Condition = as.factor(df$Condition)
df$Fenetre = as.factor(df$Fenetre)
df$Tau = NULL
df$Time = NULL

################################################################################
################################################################################
################################################################################
# PHI MULTI INFORMATION
################################################################################
################################################################################
################################################################################

# Plot Detection * Fenetre
n=length(unique(df$Fenetre))
fenetre_label = round(linspace(-2.15, 3, n),2)
detection_pos = 28

phi_MI_eb = summarySE(df, measurevar="Phi_MI", 
                      groupvars=c("Detection","Fenetre"))
phi_MI_stats = summarySE(df, measurevar="Phi_MI", 
                         groupvars=c("Sujet", "Detection", "Fenetre"))

file_name = "Phi_MI_windows_sagittal.jpeg"
figure_to_save = str_c(root_path, these_path, file_name, sep="/")
jpeg(file = figure_to_save, width=1209, height=648, units = "px")
par(mar=c(4,4,0,0)+0.2)
phi_MI_time_sagittal = 
  ggplot(phi_MI_eb, aes(x=rev(Fenetre), y=Phi_MI, group=Detection, 
                        colour=Detection)) + 
  geom_errorbar(aes(ymin=Phi_MI-ci, ymax=Phi_MI+ci), width=ebwidth_w) + 
  geom_point(size=sizepoint_w) + 
  scale_x_discrete(labels=fenetre_label) + 
  geom_vline(xintercept=detection_pos, colour="red", size=1.) +
  xlab("Time (sec)") + ylab(expression(phi^MI)) +
  theme_bw() + theme(legend.position = c(0.88,0.86)) + 
  theme(text=tt, axis.text.x=atx)
phi_MI_time_sagittal
dev.off()

file_name = "Phi_MI_windows_sagittal_by_subject.jpeg"
figure_to_save = str_c(root_path, these_path, file_name, sep="/")
jpeg(file = figure_to_save, width=1209, height=648, units = "px")
par(mar=c(4,4,0,0)+0.2)
phi_MI_time_sagittal_by_subject = 
  ggplot(phi_MI_stats, aes(x=rev(Fenetre), y=Phi_MI, group=Detection, 
                           fill=Sujet, colour=Detection)) + 
  geom_errorbar(aes(ymin=Phi_MI-ci, ymax=Phi_MI+ci), width=.1) + 
  geom_point(size=1) + xlab("Time (sec)") +
  scale_x_discrete(labels=fenetre_label) + ylab(expression(phi^MI)) +
  geom_vline(xintercept=detection_pos, colour="red") +
  facet_wrap(~ Sujet) + guides(fill=FALSE) +
  theme_bw() + theme(legend.position = c(0.95, 0.03)) + 
  theme(text=tt, axis.text.x = element_text(color="black", size=4, angle=90, 
                                            hjust=1, vjust=0))
phi_MI_time_sagittal_by_subject
dev.off()

# Statistical Analysis Detection * Fenetre
phi_MI_lmm = lme(Phi_MI ~ Detection * Fenetre, data=phi_MI_stats, 
                 random = ~1|Sujet, method= "ML")

residuals = resid(phi_MI_lmm)
plot(fitted(phi_MI_lmm), residuals)
abline(0,0)
plot(fitted(phi_MI_lmm), phi_MI_lmm$data$Phi_MI)
qqnorm(residuals)
qqline(residuals)

summary(phi_MI_lmm)
anova(phi_MI_lmm)

emmip(phi_MI_lmm, Detection ~ Fenetre)

phi_MI_lmm.emm = emmeans(phi_MI_lmm, ~ Detection * Fenetre)
contrast(phi_MI_lmm.emm, "eff", by="Detection")
contrast(phi_MI_lmm.emm, "eff", by="Fenetre")
contrast(phi_MI_lmm.emm, interaction = c("eff", "pairwise"), 
         adjust="bonferroni")

phi_MI_lmm = lmer(Phi_MI ~ Detection * Fenetre + (1|Sujet), data=phi_MI_stats)
anova(phi_MI_lmm)
phi_MI.emm = emmeans(phi_MI_lmm, ~ Detection*Fenetre)

contrast(phi_MI_lmm.emm, "eff", by="Detection")
contrast(phi_MI_lmm.emm, "eff", by="Fenetre")
contrast(phi_MI.emm, interaction = c("eff", "pairwise"), adjust = "bonferroni")

# post-hocs Detection 
phi_MI_stats$res = interaction(phi_MI_stats$Detection)
phi_MI_stats$res = factor(phi_MI_stats$res)
phi_MI_mixmodel_det = lme(Phi_MI ~ res, random=~1|Sujet, data=phi_MI_stats, 
                          method = "ML")
mcp_det = glht(phi_MI_mixmodel_det, linfct=mcp(res="Tukey"))
summary(mcp_det, test=adjusted(type="bonferroni"))
pq_det = summary(mcp_det, test=adjusted(type="bonferroni"))$test
mtests_det = cbind(pq_det$coefficients, pq_det$sigma, pq_det$tstat, 
                   pq_det$pvalues)
error_det = attr(pq_det$pvalues, "error")
pname_det = switch(
  mcp_det$alternativ,
  less = paste("Pr(<", ifelse(mcp_det$df ==0, "z", "t"), ")", sep=""),
  greater = paste("Pr(>", ifelse(mcp_det$df == 0, "z", "t"), ")", sep=""),
  two.sided = paste("Pr(>|", ifelse(mcp_det$df == 0, "z", "t"), "|)", sep=""))
colnames(mtests_det) = c("Estimate", "Std. Error", 
                         ifelse(mcp_det$df==0,"z value","t value"),pname_det)

# post-hocs Fenetre 
phi_MI_stats$res = interaction(phi_MI_stats$Fenetre)
phi_MI_stats$res = factor(phi_MI_stats$res)
phi_MI_mixmodel_fen = lme(Phi_MI ~ res, random=~1|Sujet, data=phi_MI_stats, 
                          method = "ML")
mcp_fen = glht(phi_MI_mixmodel_fen, linfct=mcp(res="Tukey"))
summary(mcp_fen, test=adjusted(type="bonferroni"))
pq_fen = summary(mcp_fen, test=adjusted(type="bonferroni"))$test
mtests_fen = cbind(pq_fen$coefficients, pq_fen$sigma, pq_fen$tstat, 
                   pq_fen$pvalues)
error_fen = attr(pq_fen$pvalues, "error")
pname_fen = switch(
  mcp_fen$alternativ,
  less = paste("Pr(<", ifelse(mcp_fen$df ==0, "z", "t"), ")", sep=""),
  greater = paste("Pr(>", ifelse(mcp_fen$df == 0, "z", "t"), ")", sep=""),
  two.sided = paste("Pr(>|", ifelse(mcp_fen$df == 0, "z", "t"), "|)", sep=""))
colnames(mtests_fen) = c("Estimate", "Std. Error", 
                         ifelse(mcp_fen$df==0,"z value","t value"),pname_fen)

# post-hocs Detection * Fenetre 
phi_MI_stats$res = interaction(phi_MI_stats$Detection, phi_MI_stats$Fenetre)
phi_MI_stats$res = factor(phi_MI_stats$res)
phi_MI_mixmodel_det_fen = lme(Phi_MI ~ res, random=~1|Sujet, 
                              data=phi_MI_stats, method = "ML")
mcp_det_fen = glht(phi_MI_mixmodel_det_fen, linfct=mcp(res="Tukey"))
summary(mcp_det_fen, test=adjusted(type="bonferroni"))
pq_det_fen = summary(mcp_det_fen, test=adjusted(type="bonferroni"))$test
mtests_det_fen = cbind(pq_det_fen$coefficients, pq_det_fen$sigma, 
                       pq_det_fen$tstat, pq_det_fen$pvalues)
error_det_fen = attr(pq_det_fen$pvalues, "error")
pname_det_fen = switch(
  mcp_det_fen$alternativ,
  less = paste("Pr(<", ifelse(mcp_det_fen$df ==0, "z", "t"), ")", sep=""),
  greater = paste("Pr(>", ifelse(mcp_det_fen$df == 0, "z", "t"), ")", sep=""),
  two.sided = paste("Pr(>|", ifelse(mcp_det_fen$df == 0,"z","t"),"|)", sep=""))
colnames(mtests_det_fen) = c("Estimate", "Std. Error", 
                             ifelse(mcp_det_fen$df==0,"z value","t value"),
                             pname_det_fen)

# RESUME
report(phi_MI_lmm)
report(anova(phi_MI_lmm))

xtable(anova(phi_MI_lmm))
xtable(mtests_det)
xtable(mtests_fen)
xtable(mtests_det_fen)

# xtable(contrast(phi_MI_lmm.emm, "eff", by="Detection"))
# xtable(contrast(phi_MI_lmm.emm, "eff", by="Fenetre"))
# xtable(contrast(phi_MI.emm, interaction = c("eff", "pairwise"), 
#                 adjust = "bonferroni"))
# xtable(report_table(phi_MI_lmm))
# xtable(report_table(anova(phi_MI_lmm)))

# Plot Phi_MI for significant windows Hit VS Miss
y_min = min(phi_MI_stats$Phi_MI)-min(phi_MI_stats$se)
y_max = max(phi_MI_stats$Phi_MI)+max(phi_MI_stats$se)

win = c(1:n)
# sign_win = c(1:14,34:65) 
# 65-34=31 / 65-14=51
# sign_win = c(1:30,51:65)
# non_sign_win = win[-sign_win]

file_name = "Phi_MI_windows_sagittal_stats.jpeg"
figure_to_save = str_c(root_path, these_path, file_name, sep="/")
jpeg(file = figure_to_save, width=1209, height=648, units = "px")
par(mar=c(4,4,0,0)+0.2)
phi_MI_time_sagittal_stats = 
  ggplot(phi_MI_eb, aes(x=rev(Fenetre), y=Phi_MI, 
                        colour=Detection, group=Detection)) + 
  geom_line(position=pd_w) + geom_point(position=pd_w, size=sizepoint_w) +
  geom_errorbar(aes(ymin=Phi_MI-ci, ymax=Phi_MI+ci), 
                position=pd_w, width=ebwidth_w, cex=cexx_w) + 
  # geom_vline(xintercept=sign_win, colour="black", size=1.) +
  geom_vline(xintercept=detection_pos, colour="red", size=1.) +
  xlab("Time (sec)") + ylab(expression(phi^MI)) +
  scale_x_discrete(labels=fenetre_label) + 
  theme_bw() + theme(legend.position = c(0.88,0.86)) + 
  theme(text=tt, axis.text.x=atx)
phi_MI_time_sagittal_stats
dev.off()

# Plot Phi_MI HIT-MISS DIFFERENCE vs window For FC Cluster 
hit = phi_MI_eb[phi_MI_eb$Detection == "hits",] 
miss = phi_MI_eb[phi_MI_eb$Detection == "miss",] 
diff = hit$Phi_MI-miss$Phi_MI
se = hit$se-miss$se
fen = as.vector(hit$Fenetre)
Phi_MI_diff = data.frame(fenetre_label, diff, se)

y_min = min(Phi_MI_diff$diff)+min(Phi_MI_diff$diff)/10
y_max = max(Phi_MI_diff$diff)+max(Phi_MI_diff$diff)/10

file_name = "Phi_MI_windows_sagittal_diff_hitmiss.jpeg"
figure_to_save = str_c(root_path, these_path, file_name, sep="/")
jpeg(file = figure_to_save, width=1209, height=648, units = "px")
par(mar=c(4,4,0,0)+0.2)
phi_MI_time_sagittal_diff = 
  ggplot(Phi_MI_diff, aes(x=rev(fenetre_label), y=diff)) + 
  geom_line(position=pd_w) + geom_point(position=pd_w, size=sizepoint_w) +
  xlab("Time (sec)") + ylab("Diff Hit - Miss") +
  geom_vline(xintercept=0., colour="red", size=1.) +
  geom_hline(yintercept=0., colour="darkblue", size=1.) +
  theme_bw() + ylim(y_min, y_max) + 
  theme(text=tt, axis.text.x=atx_diff)
phi_MI_time_sagittal_diff
dev.off()

diff = miss$Phi_MI-hit$Phi_MI
se = miss$se-hit$se
fen = as.vector(hit$Fenetre)
Phi_MI_diff = data.frame(fenetre_label, diff, se)

y_min = min(Phi_MI_diff$diff)+min(Phi_MI_diff$diff)/10
y_max = max(Phi_MI_diff$diff)+max(Phi_MI_diff$diff)/10

file_name = "Phi_MI_windows_sagittal_diff_misshit.jpeg"
figure_to_save = str_c(root_path, these_path, file_name, sep="/")
jpeg(file = figure_to_save, width=1209, height=648, units = "px")
par(mar=c(4,4,0,0)+0.2)
phi_MI_time_sagittal_diff = 
  ggplot(Phi_MI_diff, aes(x=rev(fenetre_label), y=diff)) + 
  geom_line(position=pd_w) + geom_point(position=pd_w, size=sizepoint_w) +
  xlab("Time (sec)") + ylab("Diff Miss - Hit") +
  geom_vline(xintercept=0., colour="red", size=1.) +
  geom_hline(yintercept=0., colour="darkblue", size=1.) +
  theme_bw() + ylim(y_min, y_max) + 
  theme(text=tt, axis.text.x=atx_diff)
phi_MI_time_sagittal_diff
dev.off()

################################################################################
################################################################################
################################################################################
# PHI H
################################################################################
################################################################################
################################################################################

# Plot Detection * Fenetre
n=length(unique(df$Fenetre))
fenetre_label = round(linspace(-2.15, 3, n),2)
detection_pos = 28

phi_H_eb = summarySE(df, measurevar="Phi_H", 
                     groupvars=c("Detection","Fenetre"))
phi_H_stats = summarySE(df, measurevar="Phi_H", 
                        groupvars=c("Sujet", "Detection", "Fenetre"))

file_name = "Phi_H_windows_sagittal.jpeg"
figure_to_save = str_c(root_path, these_path, file_name, sep="/")
jpeg(file = figure_to_save, width=1209, height=648, units = "px")
par(mar=c(4,4,0,0)+0.2)
phi_H_time_sagittal = 
  ggplot(phi_H_eb, aes(x=rev(Fenetre), y=Phi_H, group=Detection, 
                       colour=Detection)) + 
  geom_errorbar(aes(ymin=Phi_H-ci, ymax=Phi_H+ci), width=ebwidth_w) + 
  geom_point(size=sizepoint_w) + 
  scale_x_discrete(labels=fenetre_label) + 
  geom_vline(xintercept=detection_pos, colour="red", size=1.) +
  xlab("Time (sec)") + ylab(expression(phi^H)) +
  theme_bw() + theme(legend.position = c(0.88,0.86)) + 
  theme(text=tt, axis.text.x=atx)
phi_H_time_sagittal
dev.off()

file_name = "Phi_H_windows_sagittal_by_subject.jpeg"
figure_to_save = str_c(root_path, these_path, file_name, sep="/")
jpeg(file = figure_to_save, width=1209, height=648, units = "px")
par(mar=c(4,4,0,0)+0.2)
phi_H_time_sagittal_by_subject = 
  ggplot(phi_H_stats, aes(x=rev(Fenetre), y=Phi_H, group=Detection, 
                          fill=Sujet, colour=Detection)) + 
  geom_errorbar(aes(ymin=Phi_H-ci, ymax=Phi_H+ci), width=.1) + 
  geom_point(size=1) + xlab("Time (sec)") +
  scale_x_discrete(labels=fenetre_label) + ylab(expression(phi^H)) +
  geom_vline(xintercept=detection_pos, colour="red") +
  facet_wrap(~ Sujet) + guides(fill=FALSE) +
  theme_bw() + theme(legend.position = c(0.95, 0.03)) + 
  theme(text=tt, axis.text.x = element_text(color="black", size=4, angle=90, 
                                            hjust=1, vjust=0))
phi_H_time_sagittal_by_subject
dev.off()

# Statistical Analysis Detection * Fenetre
phi_H_lmm = lme(Phi_H ~ Detection * Fenetre, data=phi_H_stats, 
                random = ~1|Sujet, method= "ML")

residuals = resid(phi_H_lmm)
plot(fitted(phi_H_lmm), residuals)
abline(0,0)
plot(fitted(phi_H_lmm), phi_H_lmm$data$Phi_H)
qqnorm(residuals)
qqline(residuals)

summary(phi_H_lmm)
anova(phi_H_lmm)

emmip(phi_H_lmm, Detection ~ Fenetre)

phi_H_lmm.emm = emmeans(phi_H_lmm, ~ Detection * Fenetre)
contrast(phi_H_lmm.emm, "eff", by="Detection")
contrast(phi_H_lmm.emm, "eff", by="Fenetre")
contrast(phi_H_lmm.emm, interaction = c("eff", "pairwise"), 
         adjust="bonferroni")

phi_H_lmm = lmer(Phi_H ~ Detection * Fenetre + (1|Sujet), data=phi_H_stats)
anova(phi_H_lmm)
phi_H.emm = emmeans(phi_H_lmm, ~ Detection*Fenetre)

contrast(phi_H.emm, "eff", by="Detection")
contrast(phi_H.emm, "eff", by="Fenetre")
contrast(phi_H.emm, interaction = c("eff", "pairwise"), adjust="bonferroni")

# post-hocs Detection 
phi_H_stats$res = interaction(phi_H_stats$Detection)
phi_H_stats$res = factor(phi_H_stats$res)
phi_H_mixmodel_det = lme(Phi_H ~ res, random=~1|Sujet, data=phi_H_stats, 
                         method = "ML")
mcp_det = glht(phi_H_mixmodel_det, linfct=mcp(res="Tukey"))
summary(mcp_det, test=adjusted(type="bonferroni"))
pq_det = summary(mcp_det, test=adjusted(type="bonferroni"))$test
mtests_det = cbind(pq_det$coefficients, pq_det$sigma, pq_det$tstat, 
                   pq_det$pvalues)
error_det = attr(pq_det$pvalues, "error")
pname_det = switch(
  mcp_det$alternativ,
  less = paste("Pr(<", ifelse(mcp_det$df ==0, "z", "t"), ")", sep=""),
  greater = paste("Pr(>", ifelse(mcp_det$df == 0, "z", "t"), ")", sep=""),
  two.sided = paste("Pr(>|", ifelse(mcp_det$df == 0, "z", "t"), "|)", sep=""))
colnames(mtests_det) = c("Estimate", "Std. Error", 
                         ifelse(mcp_det$df==0,"z value","t value"),pname_det)

# post-hocs Fenetre 
phi_H_stats$res = interaction(phi_H_stats$Fenetre)
phi_H_stats$res = factor(phi_H_stats$res)
phi_H_mixmodel_fen = lme(Phi_H ~ res, random=~1|Sujet, data=phi_H_stats, 
                         method = "ML")
mcp_fen = glht(phi_H_mixmodel_fen, linfct=mcp(res="Tukey"))
summary(mcp_fen, test=adjusted(type="bonferroni"))
pq_fen = summary(mcp_fen, test=adjusted(type="bonferroni"))$test
mtests_fen = cbind(pq_fen$coefficients, pq_fen$sigma, pq_fen$tstat, 
                   pq_fen$pvalues)
error_fen = attr(pq_fen$pvalues, "error")
pname_fen = switch(
  mcp_fen$alternativ,
  less = paste("Pr(<", ifelse(mcp_fen$df ==0, "z", "t"), ")", sep=""),
  greater = paste("Pr(>", ifelse(mcp_fen$df == 0, "z", "t"), ")", sep=""),
  two.sided = paste("Pr(>|", ifelse(mcp_fen$df == 0, "z", "t"), "|)", sep=""))
colnames(mtests_fen) = c("Estimate", "Std. Error", 
                         ifelse(mcp_fen$df==0,"z value","t value"),pname_fen)

# post-hocs Detection * Fenetre 
phi_H_stats$res = interaction(phi_H_stats$Detection, phi_H_stats$Fenetre)
phi_H_stats$res = factor(phi_H_stats$res)
phi_H_mixmodel_det_fen = lme(Phi_H ~ res, random=~1|Sujet, 
                             data=phi_H_stats, method = "ML")
mcp_det_fen = glht(phi_H_mixmodel_det_fen, linfct=mcp(res="Tukey"))
summary(mcp_det_fen, test=adjusted(type="bonferroni"))
pq_det_fen = summary(mcp_det_fen, test=adjusted(type="bonferroni"))$test
mtests_det_fen = cbind(pq_det_fen$coefficients, pq_det_fen$sigma, 
                       pq_det_fen$tstat, pq_det_fen$pvalues)
error_det_fen = attr(pq_det_fen$pvalues, "error")
pname_det_fen = switch(
  mcp_det_fen$alternativ,
  less = paste("Pr(<", ifelse(mcp_det_fen$df ==0, "z", "t"), ")", sep=""),
  greater = paste("Pr(>", ifelse(mcp_det_fen$df == 0, "z", "t"), ")", sep=""),
  two.sided = paste("Pr(>|", ifelse(mcp_det_fen$df == 0,"z","t"),"|)", sep=""))
colnames(mtests_det_fen) = c("Estimate", "Std. Error", 
                             ifelse(mcp_det_fen$df==0,"z value","t value"),
                             pname_det_fen)

# RESUME
report(phi_H_lmm)
report(anova(phi_H_lmm))

xtable(anova(phi_H_lmm))
xtable(mtests_det)
xtable(mtests_fen)
xtable(mtests_det_fen)

# xtable(contrast(phi_H_lmm.emm, "eff", by="Detection"))
# xtable(contrast(phi_H_lmm.emm, "eff", by="Fenetre"))
# xtable(contrast(phi_H.emm, interaction = c("eff", "pairwise"), 
#                 adjust = "bonferroni"))
# xtable(report_table(phi_H_lmm))
# xtable(report_table(anova(phi_H_lmm)))

# Plot Phi_H HIT-MISS DIFFERENCE vs window For FC Cluster 
hit = phi_H_eb[phi_H_eb$Detection == "hits",] 
miss = phi_H_eb[phi_H_eb$Detection == "miss",] 
diff = hit$Phi_H-miss$Phi_H
se = hit$se-miss$se
fen = as.vector(hit$Fenetre)
phi_H_diff = data.frame(fenetre_label, diff, se)

y_min = min(phi_H_diff$diff)+min(phi_H_diff$diff)/10
y_max = max(phi_H_diff$diff)+max(phi_H_diff$diff)/10

file_name = "Phi_H_windows_sagittal_diff_hitmiss.jpeg"
figure_to_save = str_c(root_path, these_path, file_name, sep="/")
jpeg(file = figure_to_save, width=1209, height=648, units = "px")
par(mar=c(4,4,0,0)+0.2)
phi_H_time_sagittal_diff = 
  ggplot(phi_H_diff, aes(x=rev(fenetre_label), y=diff)) + 
  geom_line(position=pd_w) + geom_point(position=pd_w, size=sizepoint_w) +
  xlab("Time (sec)") + ylab("Diff Hit - Miss") +
  geom_vline(xintercept=0., colour="red", size=1.) +
  geom_hline(yintercept=0., colour="darkblue", size=1.) +
  theme_bw() + ylim(y_min, y_max) + 
  theme(text=tt, axis.text.x=atx_diff)
phi_H_time_sagittal_diff
dev.off()

diff = miss$Phi_H-hit$Phi_H
se = miss$se-hit$se
fen = as.vector(hit$Fenetre)
phi_H_diff = data.frame(fenetre_label, diff, se)

y_min = min(phi_H_diff$diff)+min(phi_H_diff$diff)/10
y_max = max(phi_H_diff$diff)+max(phi_H_diff$diff)/10

file_name = "Phi_H_windows_sagittal_diff_misshit.jpeg"
figure_to_save = str_c(root_path, these_path, file_name, sep="/")
jpeg(file = figure_to_save, width=1209, height=648, units = "px")
par(mar=c(4,4,0,0)+0.2)
phi_H_time_sagittal_diff = 
  ggplot(phi_H_diff, aes(x=rev(fenetre_label), y=diff)) + 
  geom_line(position=pd_w) + geom_point(position=pd_w, size=sizepoint_w) +
  xlab("Time (sec)") + ylab("Diff Miss - Hit") +
  geom_vline(xintercept=0., colour="red", size=1.) +
  geom_hline(yintercept=0., colour="darkblue", size=1.) +
  theme_bw() + ylim(y_min, y_max) + 
  theme(text=tt, axis.text.x=atx_diff)
phi_H_time_sagittal_diff
dev.off()

################################################################################
################################################################################
################################################################################
# PHI STAR
################################################################################
################################################################################
################################################################################

# Plot Detection * Fenetre
n=length(unique(df$Fenetre))
fenetre_label = round(linspace(-2.15, 3, n),2)
detection_pos = 28

phi_star_eb = summarySE(df, measurevar="Phi_Star", 
                        groupvars=c("Detection","Fenetre"))
phi_star_stats = summarySE(df, measurevar="Phi_Star", 
                           groupvars=c("Sujet", "Detection", "Fenetre"))

file_name = "Phi_Star_windows_sagittal.jpeg"
figure_to_save = str_c(root_path, these_path, file_name, sep="/")
jpeg(file = figure_to_save, width=1209, height=648, units = "px")
par(mar=c(4,4,0,0)+0.2)
phi_star_time_sagittal = 
  ggplot(phi_star_eb, aes(x=rev(Fenetre), y=Phi_Star, group=Detection, 
                          colour=Detection)) + 
  geom_errorbar(aes(ymin=Phi_Star-ci, ymax=Phi_Star+ci), width=ebwidth_w) + 
  geom_point(size=sizepoint_w) + 
  scale_x_discrete(labels=fenetre_label) + 
  geom_vline(xintercept=detection_pos, colour="red", size=1.) +
  xlab("Time (sec)") + ylab(expression(phi^x)) +
  theme_bw() + theme(legend.position = c(0.88,0.86)) + 
  theme(text=tt, axis.text.x=atx)
phi_star_time_sagittal
dev.off()

file_name = "Phi_Star_windows_sagittal_by_subject.jpeg"
figure_to_save = str_c(root_path, these_path, file_name, sep="/")
jpeg(file = figure_to_save, width=1209, height=648, units = "px")
par(mar=c(4,4,0,0)+0.2)
phi_star_time_sagittal_by_subject = 
  ggplot(phi_star_stats, aes(x=rev(Fenetre), y=Phi_Star, group=Detection, 
                             fill=Sujet, colour=Detection)) + 
  geom_errorbar(aes(ymin=Phi_Star-ci, ymax=Phi_Star+ci), width=.1) + 
  geom_point(size=1) + xlab("Time (sec)") +
  scale_x_discrete(labels=fenetre_label) + ylab(expression(phi^x)) +
  geom_vline(xintercept=detection_pos, colour="red") +
  facet_wrap(~ Sujet) + guides(fill=FALSE) +
  theme_bw() + theme(legend.position = c(0.95, 0.03)) + 
  theme(text=tt, axis.text.x = element_text(color="black", size=4, angle=90, 
                                            hjust=1, vjust=0))
phi_star_time_sagittal_by_subject
dev.off()

# Statistical Analysis Detection * Fenetre
phi_star_lmm = lme(Phi_Star ~ Detection * Fenetre, data=phi_star_stats, 
                   random = ~1|Sujet, method= "ML")

residuals = resid(phi_star_lmm)
plot(fitted(phi_star_lmm), residuals)
abline(0,0)
plot(fitted(phi_star_lmm), phi_star_lmm$data$Phi_star)
qqnorm(residuals)
qqline(residuals)

summary(phi_star_lmm)
anova(phi_star_lmm)

emmip(phi_star_lmm, Detection ~ Fenetre)

phi_star_lmm.emm = emmeans(phi_star_lmm, ~ Detection * Fenetre)
contrast(phi_star_lmm.emm, "eff", by="Detection")
contrast(phi_star_lmm.emm, "eff", by="Fenetre")
contrast(phi_star_lmm.emm, interaction = c("eff", "pairwise"), 
         adjust="bonferroni")

phi_star_lmm = lmer(Phi_Star ~ Detection * Fenetre + (1|Sujet), 
                    data=phi_star_stats)
anova(phi_star_lmm)
phi_star.emm = emmeans(phi_star_lmm, ~ Detection*Fenetre)

contrast(phi_star.emm, "eff", by="Detection")
contrast(phi_star.emm, "eff", by="Fenetre")
contrast(phi_star.emm, interaction = c("eff", "pairwise"), adjust="bonferroni")

# post-hocs Detection 
phi_star_stats$res = interaction(phi_star_stats$Detection)
phi_star_stats$res = factor(phi_star_stats$res)
phi_star_mixmodel_det = lme(Phi_Star ~ res, random=~1|Sujet, data=phi_star_stats, 
                            method = "ML")
mcp_det = glht(phi_star_mixmodel_det, linfct=mcp(res="Tukey"))
summary(mcp_det, test=adjusted(type="bonferroni"))
pq_det = summary(mcp_det, test=adjusted(type="bonferroni"))$test
mtests_det = cbind(pq_det$coefficients, pq_det$sigma, pq_det$tstat, 
                   pq_det$pvalues)
error_det = attr(pq_det$pvalues, "error")
pname_det = switch(
  mcp_det$alternativ,
  less = paste("Pr(<", ifelse(mcp_det$df ==0, "z", "t"), ")", sep=""),
  greater = paste("Pr(>", ifelse(mcp_det$df == 0, "z", "t"), ")", sep=""),
  two.sided = paste("Pr(>|", ifelse(mcp_det$df == 0, "z", "t"), "|)", sep=""))
colnames(mtests_det) = c("Estimate", "Std. Error", 
                         ifelse(mcp_det$df==0,"z value","t value"),pname_det)

# post-hocs Fenetre 
phi_star_stats$res = interaction(phi_star_stats$Fenetre)
phi_star_stats$res = factor(phi_star_stats$res)
phi_star_mixmodel_fen = lme(Phi_Star ~ res, random=~1|Sujet, data=phi_star_stats, 
                            method = "ML")
mcp_fen = glht(phi_star_mixmodel_fen, linfct=mcp(res="Tukey"))
summary(mcp_fen, test=adjusted(type="bonferroni"))
pq_fen = summary(mcp_fen, test=adjusted(type="bonferroni"))$test
mtests_fen = cbind(pq_fen$coefficients, pq_fen$sigma, pq_fen$tstat, 
                   pq_fen$pvalues)
error_fen = attr(pq_fen$pvalues, "error")
pname_fen = switch(
  mcp_fen$alternativ,
  less = paste("Pr(<", ifelse(mcp_fen$df ==0, "z", "t"), ")", sep=""),
  greater = paste("Pr(>", ifelse(mcp_fen$df == 0, "z", "t"), ")", sep=""),
  two.sided = paste("Pr(>|", ifelse(mcp_fen$df == 0, "z", "t"), "|)", sep=""))
colnames(mtests_fen) = c("Estimate", "Std. Error", 
                         ifelse(mcp_fen$df==0,"z value","t value"),pname_fen)

# post-hocs Detection * Fenetre 
phi_star_stats$res = interaction(phi_star_stats$Detection, phi_star_stats$Fenetre)
phi_star_stats$res = factor(phi_star_stats$res)
phi_star_mixmodel_det_fen = lme(Phi_Star ~ res, random=~1|Sujet, 
                                data=phi_star_stats, method = "ML")
mcp_det_fen = glht(phi_star_mixmodel_det_fen, linfct=mcp(res="Tukey"))
summary(mcp_det_fen, test=adjusted(type="bonferroni"))
pq_det_fen = summary(mcp_det_fen, test=adjusted(type="bonferroni"))$test
mtests_det_fen = cbind(pq_det_fen$coefficients, pq_det_fen$sigma, 
                       pq_det_fen$tstat, pq_det_fen$pvalues)
error_det_fen = attr(pq_det_fen$pvalues, "error")
pname_det_fen = switch(
  mcp_det_fen$alternativ,
  less = paste("Pr(<", ifelse(mcp_det_fen$df ==0, "z", "t"), ")", sep=""),
  greater = paste("Pr(>", ifelse(mcp_det_fen$df == 0, "z", "t"), ")", sep=""),
  two.sided = paste("Pr(>|", ifelse(mcp_det_fen$df == 0,"z","t"),"|)", sep=""))
colnames(mtests_det_fen) = c("Estimate", "Std. Error", 
                             ifelse(mcp_det_fen$df==0,"z value","t value"),
                             pname_det_fen)

# RESUME
report(phi_star_lmm)
report(anova(phi_star_lmm))

xtable(anova(phi_star_lmm))
xtable(mtests_det)
xtable(mtests_fen)
xtable(mtests_det_fen)

# xtable(contrast(phi_star_lmm.emm, "eff", by="Detection"))
# xtable(contrast(phi_star_lmm.emm, "eff", by="Fenetre"))
# xtable(contrast(phi_star.emm, interaction = c("eff", "pairwise"), 
#                 adjust = "bonferroni"))
# xtable(report_table(phi_star_lmm))
# xtable(report_table(anova(phi_star_lmm)))

# Plot Phi_Star HIT-MISS DIFFERENCE vs window For FC Cluster 
hit = phi_star_eb[phi_star_eb$Detection == "hits",] 
miss = phi_star_eb[phi_star_eb$Detection == "miss",] 
diff = hit$Phi_Star-miss$Phi_Star
se = hit$se-miss$se
fen = as.vector(hit$Fenetre)
phi_star_diff = data.frame(fenetre_label, diff, se)

y_min = min(phi_star_diff$diff)+min(phi_star_diff$diff)/10
y_max = max(phi_star_diff$diff)+max(phi_star_diff$diff)/10

file_name = "Phi_Star_windows_sagittal_diff_hitmiss.jpeg"
figure_to_save = str_c(root_path, these_path, file_name, sep="/")
jpeg(file = figure_to_save, width=1209, height=648, units = "px")
par(mar=c(4,4,0,0)+0.2)
phi_star_time_sagittal_diff = 
  ggplot(phi_star_diff, aes(x=rev(fenetre_label), y=diff)) + 
  geom_line(position=pd_w) + geom_point(position=pd_w, size=sizepoint_w) +
  xlab("Time (sec)") + ylab("Diff Hit - Miss") +
  geom_vline(xintercept=0., colour="red", size=1.) +
  geom_hline(yintercept=0., colour="darkblue", size=1.) +
  theme_bw() + ylim(y_min, y_max) + 
  theme(text=tt, axis.text.x=atx_diff)
phi_star_time_sagittal_diff
dev.off()

diff = miss$Phi_Star-hit$Phi_Star
se = miss$se-hit$se
fen = as.vector(hit$Fenetre)
phi_star_diff = data.frame(fenetre_label, diff, se)

y_min = min(phi_star_diff$diff)+min(phi_star_diff$diff)/10
y_max = max(phi_star_diff$diff)+max(phi_star_diff$diff)/10

file_name = "Phi_Star_windows_sagittal_diff_misshit.jpeg"
figure_to_save = str_c(root_path, these_path, file_name, sep="/")
jpeg(file = figure_to_save, width=1209, height=648, units = "px")
par(mar=c(4,4,0,0)+0.2)
phi_star_time_sagittal_diff = 
  ggplot(phi_star_diff, aes(x=rev(fenetre_label), y=diff)) + 
  geom_line(position=pd_w) + geom_point(position=pd_w, size=sizepoint_w) +
  xlab("Time (sec)") + ylab("Diff Miss - Hit") +
  geom_vline(xintercept=0., colour="red", size=1.) +
  geom_hline(yintercept=0., colour="darkblue", size=1.) +
  theme_bw() + ylim(y_min, y_max) + 
  theme(text=tt, axis.text.x=atx_diff)
phi_star_time_sagittal_diff
dev.off()

################################################################################
################################################################################
################################################################################
# PHI GEO
################################################################################
################################################################################
################################################################################

# Plot Detection * Fenetre
n=length(unique(df$Fenetre))
fenetre_label = round(linspace(-2.15, 3, n),2)
detection_pos = 28

phi_geo_eb = summarySE(df, measurevar="Phi_Geo", 
                       groupvars=c("Detection","Fenetre"))
phi_geo_stats = summarySE(df, measurevar="Phi_Geo", 
                          groupvars=c("Sujet", "Detection", "Fenetre"))

file_name = "Phi_Geo_windows_sagittal.jpeg"
figure_to_save = str_c(root_path, these_path, file_name, sep="/")
jpeg(file = figure_to_save, width=1209, height=648, units = "px")
par(mar=c(4,4,0,0)+0.2)
phi_geo_time_sagittal = 
  ggplot(phi_geo_eb, aes(x=rev(Fenetre), y=Phi_Geo, group=Detection, 
                         colour=Detection)) + 
  geom_errorbar(aes(ymin=Phi_Geo-ci, ymax=Phi_Geo+ci), width=ebwidth_w) + 
  geom_point(size=sizepoint_w) + 
  scale_x_discrete(labels=fenetre_label) + 
  geom_vline(xintercept=detection_pos, colour="red", size=1.) +
  xlab("Time (sec)") + ylab(expression(phi^G)) +
  theme_bw() + theme(legend.position = c(0.88,0.86)) + 
  theme(text=tt, axis.text.x=atx)
phi_geo_time_sagittal
dev.off()

file_name = "Phi_Geo_windows_sagittal_by_subject.jpeg"
figure_to_save = str_c(root_path, these_path, file_name, sep="/")
jpeg(file = figure_to_save, width=1209, height=648, units = "px")
par(mar=c(4,4,0,0)+0.2)
phi_geo_time_sagittal_by_subject = 
  ggplot(phi_geo_stats, aes(x=rev(Fenetre), y=Phi_Geo, group=Detection, 
                            fill=Sujet, colour=Detection)) + 
  geom_errorbar(aes(ymin=Phi_Geo-ci, ymax=Phi_Geo+ci), width=.1) + 
  geom_point(size=1) + xlab("Time (sec)") +
  scale_x_discrete(labels=fenetre_label) + ylab(expression(phi^G)) +
  geom_vline(xintercept=detection_pos, colour="red") +
  facet_wrap(~ Sujet) + guides(fill=FALSE) +
  theme_bw() + theme(legend.position = c(0.95, 0.03)) + 
  theme(text=tt, axis.text.x = element_text(color="black", size=4, angle=90, 
                                            hjust=1, vjust=0))
phi_geo_time_sagittal_by_subject
dev.off()

# Statistical Analysis Detection * Fenetre
phi_geo_lmm = lme(Phi_Geo ~ Detection * Fenetre, data=phi_geo_stats, 
                  random = ~1|Sujet, method= "ML")

residuals = resid(phi_geo_lmm)
plot(fitted(phi_geo_lmm), residuals)
abline(0,0)
plot(fitted(phi_geo_lmm), phi_geo_lmm$data$Phi_geo)
qqnorm(residuals)
qqline(residuals)

summary(phi_geo_lmm)
anova(phi_geo_lmm)

emmip(phi_geo_lmm, Detection ~ Fenetre)

phi_geo_lmm.emm = emmeans(phi_geo_lmm, ~ Detection * Fenetre)
contrast(phi_geo_lmm.emm, "eff", by="Detection")
contrast(phi_geo_lmm.emm, "eff", by="Fenetre")
contrast(phi_geo_lmm.emm, interaction = c("eff", "pairwise"), 
         adjust="bonferroni")

phi_geo_lmm = lmer(Phi_Geo ~ Detection * Fenetre + (1|Sujet), 
                   data=phi_geo_stats)
anova(phi_geo_lmm)
phi_geo.emm = emmeans(phi_geo_lmm, ~ Detection*Fenetre)

contrast(phi_geo.emm, "eff", by="Detection")
contrast(phi_geo.emm, "eff", by="Fenetre")
contrast(phi_geo.emm, interaction = c("eff", "pairwise"), adjust="bonferroni")

# post-hocs Detection 
phi_geo_stats$res = interaction(phi_geo_stats$Detection)
phi_geo_stats$res = factor(phi_geo_stats$res)
phi_geo_mixmodel_det = lme(Phi_Geo ~ res, random=~1|Sujet, data=phi_geo_stats, 
                           method = "ML")
mcp_det = glht(phi_geo_mixmodel_det, linfct=mcp(res="Tukey"))
summary(mcp_det, test=adjusted(type="bonferroni"))
pq_det = summary(mcp_det, test=adjusted(type="bonferroni"))$test
mtests_det = cbind(pq_det$coefficients, pq_det$sigma, pq_det$tstat, 
                   pq_det$pvalues)
error_det = attr(pq_det$pvalues, "error")
pname_det = switch(
  mcp_det$alternativ,
  less = paste("Pr(<", ifelse(mcp_det$df ==0, "z", "t"), ")", sep=""),
  greater = paste("Pr(>", ifelse(mcp_det$df == 0, "z", "t"), ")", sep=""),
  two.sided = paste("Pr(>|", ifelse(mcp_det$df == 0, "z", "t"), "|)", sep=""))
colnames(mtests_det) = c("Estimate", "Std. Error", 
                         ifelse(mcp_det$df==0,"z value","t value"),pname_det)

# post-hocs Fenetre 
phi_geo_stats$res = interaction(phi_geo_stats$Fenetre)
phi_geo_stats$res = factor(phi_geo_stats$res)
phi_geo_mixmodel_fen = lme(Phi_Geo ~ res, random=~1|Sujet, data=phi_geo_stats, 
                           method = "ML")
mcp_fen = glht(phi_geo_mixmodel_fen, linfct=mcp(res="Tukey"))
summary(mcp_fen, test=adjusted(type="bonferroni"))
pq_fen = summary(mcp_fen, test=adjusted(type="bonferroni"))$test
mtests_fen = cbind(pq_fen$coefficients, pq_fen$sigma, pq_fen$tstat, 
                   pq_fen$pvalues)
error_fen = attr(pq_fen$pvalues, "error")
pname_fen = switch(
  mcp_fen$alternativ,
  less = paste("Pr(<", ifelse(mcp_fen$df ==0, "z", "t"), ")", sep=""),
  greater = paste("Pr(>", ifelse(mcp_fen$df == 0, "z", "t"), ")", sep=""),
  two.sided = paste("Pr(>|", ifelse(mcp_fen$df == 0, "z", "t"), "|)", sep=""))
colnames(mtests_fen) = c("Estimate", "Std. Error", 
                         ifelse(mcp_fen$df==0,"z value","t value"),pname_fen)

# post-hocs Detection * Fenetre 
phi_geo_stats$res = interaction(phi_geo_stats$Detection, phi_geo_stats$Fenetre)
phi_geo_stats$res = factor(phi_geo_stats$res)
phi_geo_mixmodel_det_fen = lme(Phi_Geo ~ res, random=~1|Sujet, 
                               data=phi_geo_stats, method = "ML")
mcp_det_fen = glht(phi_geo_mixmodel_det_fen, linfct=mcp(res="Tukey"))
summary(mcp_det_fen, test=adjusted(type="bonferroni"))
pq_det_fen = summary(mcp_det_fen, test=adjusted(type="bonferroni"))$test
mtests_det_fen = cbind(pq_det_fen$coefficients, pq_det_fen$sigma, 
                       pq_det_fen$tstat, pq_det_fen$pvalues)
error_det_fen = attr(pq_det_fen$pvalues, "error")
pname_det_fen = switch(
  mcp_det_fen$alternativ,
  less = paste("Pr(<", ifelse(mcp_det_fen$df ==0, "z", "t"), ")", sep=""),
  greater = paste("Pr(>", ifelse(mcp_det_fen$df == 0, "z", "t"), ")", sep=""),
  two.sided = paste("Pr(>|", ifelse(mcp_det_fen$df == 0,"z","t"),"|)", sep=""))
colnames(mtests_det_fen) = c("Estimate", "Std. Error", 
                             ifelse(mcp_det_fen$df==0,"z value","t value"),
                             pname_det_fen)

# RESUME
report(phi_geo_lmm)
report(anova(phi_geo_lmm))

xtable(anova(phi_geo_lmm))
xtable(mtests_det)
xtable(mtests_fen)
xtable(mtests_det_fen)

# xtable(contrast(phi_geo_lmm.emm, "eff", by="Detection"))
# xtable(contrast(phi_geo_lmm.emm, "eff", by="Fenetre"))
# xtable(contrast(phi_geo.emm, interaction = c("eff", "pairwise"), 
#                 adjust = "bonferroni"))
# xtable(report_table(phi_geo_lmm))
# xtable(report_table(anova(phi_geo_lmm)))

# Plot Phi_Geo HIT-MISS DIFFERENCE vs window For FC Cluster 
hit = phi_geo_eb[phi_geo_eb$Detection == "hits",] 
miss = phi_geo_eb[phi_geo_eb$Detection == "miss",] 
diff = hit$Phi_Geo-miss$Phi_Geo
se = hit$se-miss$se
fen = as.vector(hit$Fenetre)
phi_geo_diff = data.frame(fenetre_label, diff, se)

y_min = min(phi_geo_diff$diff)+min(phi_geo_diff$diff)/10
y_max = max(phi_geo_diff$diff)+max(phi_geo_diff$diff)/10

file_name = "Phi_Geo_windows_sagittal_diff_hitmiss.jpeg"
figure_to_save = str_c(root_path, these_path, file_name, sep="/")
jpeg(file = figure_to_save, width=1209, height=648, units = "px")
par(mar=c(4,4,0,0)+0.2)
phi_geo_time_sagittal_diff = 
  ggplot(phi_geo_diff, aes(x=rev(fenetre_label), y=diff)) + 
  geom_line(position=pd_w) + geom_point(position=pd_w, size=sizepoint_w) +
  xlab("Time (sec)") + ylab("Diff Hit - Miss") +
  geom_vline(xintercept=0., colour="red", size=1.) +
  geom_hline(yintercept=0., colour="darkblue", size=1.) +
  theme_bw() + ylim(y_min, y_max) + 
  theme(text=tt, axis.text.x=atx_diff)
phi_geo_time_sagittal_diff
dev.off()

diff = miss$Phi_Geo-hit$Phi_Geo
se = miss$se-hit$se
fen = as.vector(hit$Fenetre)
phi_geo_diff = data.frame(fenetre_label, diff, se)

y_min = min(phi_geo_diff$diff)+min(phi_geo_diff$diff)/10
y_max = max(phi_geo_diff$diff)+max(phi_geo_diff$diff)/10

file_name = "Phi_Geo_windows_sagittal_diff_misshit.jpeg"
figure_to_save = str_c(root_path, these_path, file_name, sep="/")
jpeg(file = figure_to_save, width=1209, height=648, units = "px")
par(mar=c(4,4,0,0)+0.2)
phi_geo_time_sagittal_diff = 
  ggplot(phi_geo_diff, aes(x=rev(fenetre_label), y=diff)) + 
  geom_line(position=pd_w) + geom_point(position=pd_w, size=sizepoint_w) +
  xlab("Time (sec)") + ylab("Diff Miss - Hit") +
  geom_vline(xintercept=0., colour="red", size=1.) +
  geom_hline(yintercept=0., colour="darkblue", size=1.) +
  theme_bw() + ylim(y_min, y_max) + 
  theme(text=tt, axis.text.x=atx_diff)
phi_geo_time_sagittal_diff
dev.off()

