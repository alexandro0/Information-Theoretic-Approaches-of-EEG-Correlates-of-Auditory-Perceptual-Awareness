################################################################################
# SCRIPT FOR PLOTTING AND ANALYSES PHI DATA
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

phi_MI_eb = summarySE(df, measurevar="Phi_MI", 
                      groupvars=c("Detection","Condition"))
phi_MI_stats = summarySE(df, measurevar="Phi_MI", 
                         groupvars=c("Sujet", "Detection", "Condition"))

cond = c("Before", "After")
phi_MI_stats$Condition = factor(phi_MI_stats$Condition, levels = cond)

# Plot Detection * Condition
file_name = "Phi_MI_sagittal.jpeg"
figure_to_save = str_c(root_path, these_path, file_name, sep="/")
jpeg(file = figure_to_save, width=1209, height=648, units = "px")
par(mar=c(4,4,0,0)+0.2)
phi_MI_sagittal = 
  ggplot(phi_MI_stats, aes(x=Condition, y=Phi_MI, fill=Detection)) +  
  geom_violin() +
  geom_dotplot(binaxis='y', stackdir='center', position=position_dodge(0.9)) +
  ylab(expression(phi^MI)) + 
  theme_bw() + theme(legend.position = c(0.88,0.86)) + theme(text=tt)
phi_MI_sagittal
dev.off()

# Statistical Analysis Detection * Condition
phi_MI_lmm = lme(Phi_MI ~ Detection * Condition, data=phi_MI_stats, 
                 random = ~1|Sujet, method= "ML")

residuals = resid(phi_MI_lmm)
plot(fitted(phi_MI_lmm), residuals)
abline(0,0)
plot(fitted(phi_MI_lmm), phi_MI_lmm$data$Phi_MI)
qqnorm(residuals)
qqline(residuals)

summary(phi_MI_lmm)
anova(phi_MI_lmm)

emmip(phi_MI_lmm, Detection ~ Condition)

phi_MI_lmm.emm = emmeans(phi_MI_lmm, ~ Detection * Condition)
contrast(phi_MI_lmm.emm, "eff", by="Detection")
contrast(phi_MI_lmm.emm, "eff", by="Condition")
contrast(phi_MI_lmm.emm, interaction = c("eff", "pairwise"), 
         adjust="bonferroni")

phi_MI_lmm = lmer(Phi_MI ~ Detection * Condition + (1|Sujet), data=phi_MI_stats)
anova(phi_MI_lmm)
phi_MI.emm = emmeans(phi_MI_lmm, ~ Detection*Condition)

contrast(phi_MI_lmm.emm, "eff", by="Detection")
contrast(phi_MI_lmm.emm, "eff", by="Condition")
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

# post-hocs Condition 
phi_MI_stats$res = interaction(phi_MI_stats$Condition)
phi_MI_stats$res = factor(phi_MI_stats$res)
phi_MI_mixmodel_cond = lme(Phi_MI ~ res, random=~1|Sujet, data=phi_MI_stats, 
                           method = "ML")
mcp_cond = glht(phi_MI_mixmodel_cond, linfct=mcp(res="Tukey"))
summary(mcp_cond, test=adjusted(type="bonferroni"))
pq_cond = summary(mcp_cond, test=adjusted(type="bonferroni"))$test
mtests_cond = cbind(pq_cond$coefficients, pq_cond$sigma, pq_cond$tstat, 
                    pq_cond$pvalues)
error_cond = attr(pq_cond$pvalues, "error")
pname_cond = switch(
  mcp_cond$alternativ,
  less = paste("Pr(<", ifelse(mcp_cond$df ==0, "z", "t"), ")", sep=""),
  greater = paste("Pr(>", ifelse(mcp_cond$df == 0, "z", "t"), ")", sep=""),
  two.sided = paste("Pr(>|", ifelse(mcp_cond$df == 0, "z", "t"), "|)", sep=""))
colnames(mtests_cond) = c("Estimate", "Std. Error", 
                          ifelse(mcp_cond$df==0,"z value","t value"),pname_cond)

# post-hocs Detection * Condition 
phi_MI_stats$res = interaction(phi_MI_stats$Detection, phi_MI_stats$Condition)
phi_MI_stats$res = factor(phi_MI_stats$res)
phi_MI_mixmodel_det_cond = lme(Phi_MI ~ res, random=~1|Sujet, 
                               data=phi_MI_stats, method = "ML")
mcp_det_cond = glht(phi_MI_mixmodel_det_cond, linfct=mcp(res="Tukey"))
summary(mcp_det_cond, test=adjusted(type="bonferroni"))
pq_det_cond = summary(mcp_det_cond, test=adjusted(type="bonferroni"))$test
mtests_det_cond = cbind(pq_det_cond$coefficients, pq_det_cond$sigma, 
                        pq_det_cond$tstat, pq_det_cond$pvalues)
error_det_cond = attr(pq_det_cond$pvalues, "error")
pname_det_cond = switch(
  mcp_det_cond$alternativ,
  less = paste("Pr(<", ifelse(mcp_det_cond$df ==0, "z", "t"), ")", sep=""),
  greater = paste("Pr(>", ifelse(mcp_det_cond$df == 0, "z", "t"), ")", sep=""),
  two.sided = paste("Pr(>|", ifelse(mcp_det_cond$df == 0,"z","t"),"|)", sep=""))
colnames(mtests_det_cond) = c("Estimate", "Std. Error", 
                              ifelse(mcp_det_cond$df==0,"z value","t value"),
                              pname_det_cond)

# RESUME
report(phi_MI_lmm)
report(anova(phi_MI_lmm))

xtable(anova(phi_MI_lmm))
xtable(mtests_det)
xtable(mtests_cond)
xtable(mtests_det_cond)

# xtable(contrast(phi_MI_lmm.emm, "eff", by="Detection"))
# xtable(contrast(phi_MI_lmm.emm, "eff", by="Condition"))
# xtable(contrast(phi_MI.emm, interaction = c("eff", "pairwise"), 
#                 adjust = "bonferroni"))
# xtable(report_table(phi_MI_lmm))
# xtable(report_table(anova(phi_MI_lmm)))

################################################################################
################################################################################
################################################################################
# PHI H
################################################################################
################################################################################
################################################################################

phi_H_eb = summarySE(df, measurevar="Phi_H", 
                     groupvars=c("Detection","Condition"))
phi_H_stats = summarySE(df, measurevar="Phi_H", 
                        groupvars=c("Sujet", "Detection", "Condition"))

cond = c("Before", "After")
phi_H_stats$Condition = factor(phi_H_stats$Condition, levels = cond)

# Plot Detection * Condition
file_name = "Phi_H_sagittal.jpeg"
figure_to_save = str_c(root_path, these_path, file_name, sep="/")
jpeg(file = figure_to_save, width=1209, height=648, units = "px")
par(mar=c(4,4,0,0)+0.2)
phi_H_sagittal = 
  ggplot(phi_H_stats, aes(x=Condition, y=Phi_H, fill=Detection)) +  
  geom_violin() +
  geom_dotplot(binaxis='y', stackdir='center', position=position_dodge(0.9)) +
  ylab(expression(phi^H)) + 
  theme_bw() + theme(legend.position = c(0.88,0.86)) + theme(text=tt)
phi_H_sagittal
dev.off()

# Statistical Analysis Detection * Condition
phi_H_lmm = lme(Phi_H ~ Detection * Condition, data=phi_H_stats, 
                random = ~1|Sujet, method= "ML")

residuals = resid(phi_H_lmm)
plot(fitted(phi_H_lmm), residuals)
abline(0,0)
plot(fitted(phi_H_lmm), phi_H_lmm$data$Phi_H)
qqnorm(residuals)
qqline(residuals)

summary(phi_H_lmm)
anova(phi_H_lmm)

emmip(phi_H_lmm, Detection ~ Condition)

phi_H_lmm.emm = emmeans(phi_H_lmm, ~ Detection * Condition)
contrast(phi_H_lmm.emm, "eff", by="Detection")
contrast(phi_H_lmm.emm, "eff", by="Condition")
contrast(phi_H_lmm.emm, interaction = c("eff", "pairwise"), adjust="bonferroni")

phi_H_lmm = lmer(Phi_H ~ Detection * Condition + (1|Sujet), data=phi_H_stats)
anova(phi_H_lmm)

phi_H.emm = emmeans(phi_H_lmm, ~ Detection*Condition)

contrast(phi_H_lmm.emm, "eff", by="Detection")
contrast(phi_H_lmm.emm, "eff", by="Condition")
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

# post-hocs Condition 
phi_H_stats$res = interaction(phi_H_stats$Condition)
phi_H_stats$res = factor(phi_H_stats$res)
phi_H_mixmodel_cond = lme(Phi_H ~ res, random=~1|Sujet, data=phi_H_stats, 
                          method = "ML")
mcp_cond = glht(phi_H_mixmodel_cond, linfct=mcp(res="Tukey"))
summary(mcp_cond, test=adjusted(type="bonferroni"))
pq_cond = summary(mcp_cond, test=adjusted(type="bonferroni"))$test
mtests_cond = cbind(pq_cond$coefficients, pq_cond$sigma, pq_cond$tstat, 
                    pq_cond$pvalues)
error_cond = attr(pq_cond$pvalues, "error")
pname_cond = switch(
  mcp_cond$alternativ,
  less = paste("Pr(<", ifelse(mcp_cond$df ==0, "z", "t"), ")", sep=""),
  greater = paste("Pr(>", ifelse(mcp_cond$df == 0, "z", "t"), ")", sep=""),
  two.sided = paste("Pr(>|", ifelse(mcp_cond$df == 0, "z", "t"), "|)", sep=""))
colnames(mtests_cond) = c("Estimate", "Std. Error", 
                          ifelse(mcp_cond$df==0,"z value","t value"),pname_cond)

# post-hocs Detection * Condition 
phi_H_stats$res = interaction(phi_H_stats$Detection, phi_H_stats$Condition)
phi_H_stats$res = factor(phi_H_stats$res)
phi_H_mixmodel_det_cond = lme(Phi_H ~ res, random=~1|Sujet, 
                              data=phi_H_stats, method = "ML")
mcp_det_cond = glht(phi_H_mixmodel_det_cond, linfct=mcp(res="Tukey"))
summary(mcp_det_cond, test=adjusted(type="bonferroni"))
pq_det_cond = summary(mcp_det_cond, test=adjusted(type="bonferroni"))$test
mtests_det_cond = cbind(pq_det_cond$coefficients, pq_det_cond$sigma, 
                        pq_det_cond$tstat, pq_det_cond$pvalues)
error_det_cond = attr(pq_det_cond$pvalues, "error")
pname_det_cond = switch(
  mcp_det_cond$alternativ,
  less = paste("Pr(<", ifelse(mcp_det_cond$df ==0, "z", "t"), ")", sep=""),
  greater = paste("Pr(>", ifelse(mcp_det_cond$df == 0, "z", "t"), ")", sep=""),
  two.sided = paste("Pr(>|", ifelse(mcp_det_cond$df == 0,"z","t"),"|)", sep=""))
colnames(mtests_det_cond) = c("Estimate", "Std. Error", 
                              ifelse(mcp_det_cond$df==0,"z value","t value"),
                              pname_det_cond)

# RESUME
report(phi_H_lmm)
report(anova(phi_H_lmm))

xtable(anova(phi_H_lmm))
xtable(mtests_det)
xtable(mtests_cond)
xtable(mtests_det_cond)

# xtable(contrast(phi_H_lmm.emm, "eff", by="Detection"))
# xtable(contrast(phi_H_lmm.emm, "eff", by="Condition"))
# xtable(contrast(phi_H.emm, interaction = c("eff", "pairwise"), 
#                 adjust = "bonferroni"))
# xtable(report_table(phi_H_lmm))
# xtable(report_table(anova(phi_H_lmm)))

################################################################################
################################################################################
################################################################################
# PHI STAR
################################################################################
################################################################################
################################################################################

phi_star_eb = summarySE(df, measurevar="Phi_Star", 
                        groupvars=c("Detection","Condition"))
phi_star_stats = summarySE(df, measurevar="Phi_Star", 
                           groupvars=c("Sujet", "Detection", "Condition"))

cond = c("Before", "After")
phi_star_stats$Condition = factor(phi_star_stats$Condition, levels = cond)

# Plot Detection * Condition
file_name = "Phi_Star_sagittal.jpeg"
figure_to_save = str_c(root_path, these_path, file_name, sep="/")
jpeg(file = figure_to_save, width=1209, height=648, units = "px")
par(mar=c(4,4,0,0)+0.2)
phi_star_sagittal = 
  ggplot(phi_star_stats, aes(x=Condition, y=Phi_Star, fill=Detection)) +  
  geom_violin() +
  geom_dotplot(binaxis='y', stackdir='center', position=position_dodge(0.9)) +
  ylab(expression(phi^x)) + 
  theme_bw() + theme(legend.position = c(0.88,0.86)) + theme(text=tt)
phi_star_sagittal
dev.off()

# Statistical Analysis Detection * Condition
phi_star_lmm = lme(Phi_Star ~ Detection * Condition, data=phi_star_stats, 
                   random = ~1|Sujet, method= "ML")

residuals = resid(phi_star_lmm)
plot(fitted(phi_star_lmm), residuals)
abline(0,0)
plot(fitted(phi_star_lmm), phi_star_lmm$data$Phi_star)
qqnorm(residuals)
qqline(residuals)

summary(phi_star_lmm)
anova(phi_star_lmm)

emmip(phi_star_lmm, Detection ~ Condition)

phi_star_lmm.emm = emmeans(phi_star_lmm, ~ Detection * Condition)
contrast(phi_star_lmm.emm, "eff", by="Detection")
contrast(phi_star_lmm.emm, "eff", by="Condition")
contrast(phi_star_lmm.emm, interaction = c("eff", "pairwise"), 
         adjust="bonferroni")

phi_star_lmm = lmer(Phi_Star ~ Detection * Condition + (1|Sujet), 
                    data=phi_star_stats)
anova(phi_star_lmm)

phi_star.emm = emmeans(phi_star_lmm, ~ Detection*Condition)

contrast(phi_star.emm, "eff", by="Detection")
contrast(phi_star.emm, "eff", by="Condition")
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

# post-hocs Condition 
phi_star_stats$res = interaction(phi_star_stats$Condition)
phi_star_stats$res = factor(phi_star_stats$res)
phi_star_mixmodel_cond = lme(Phi_Star ~ res, random=~1|Sujet, data=phi_star_stats, 
                             method = "ML")
mcp_cond = glht(phi_star_mixmodel_cond, linfct=mcp(res="Tukey"))
summary(mcp_cond, test=adjusted(type="bonferroni"))
pq_cond = summary(mcp_cond, test=adjusted(type="bonferroni"))$test
mtests_cond = cbind(pq_cond$coefficients, pq_cond$sigma, pq_cond$tstat, 
                    pq_cond$pvalues)
error_cond = attr(pq_cond$pvalues, "error")
pname_cond = switch(
  mcp_cond$alternativ,
  less = paste("Pr(<", ifelse(mcp_cond$df ==0, "z", "t"), ")", sep=""),
  greater = paste("Pr(>", ifelse(mcp_cond$df == 0, "z", "t"), ")", sep=""),
  two.sided = paste("Pr(>|", ifelse(mcp_cond$df == 0, "z", "t"), "|)", sep=""))
colnames(mtests_cond) = c("Estimate", "Std. Error", 
                          ifelse(mcp_cond$df==0,"z value","t value"),pname_cond)

# post-hocs Detection * Condition 
phi_star_stats$res = interaction(phi_star_stats$Detection, phi_star_stats$Condition)
phi_star_stats$res = factor(phi_star_stats$res)
phi_star_mixmodel_det_cond = lme(Phi_Star ~ res, random=~1|Sujet, 
                                 data=phi_star_stats, method = "ML")
mcp_det_cond = glht(phi_star_mixmodel_det_cond, linfct=mcp(res="Tukey"))
summary(mcp_det_cond, test=adjusted(type="bonferroni"))
pq_det_cond = summary(mcp_det_cond, test=adjusted(type="bonferroni"))$test
mtests_det_cond = cbind(pq_det_cond$coefficients, pq_det_cond$sigma, 
                        pq_det_cond$tstat, pq_det_cond$pvalues)
error_det_cond = attr(pq_det_cond$pvalues, "error")
pname_det_cond = switch(
  mcp_det_cond$alternativ,
  less = paste("Pr(<", ifelse(mcp_det_cond$df ==0, "z", "t"), ")", sep=""),
  greater = paste("Pr(>", ifelse(mcp_det_cond$df == 0, "z", "t"), ")", sep=""),
  two.sided = paste("Pr(>|", ifelse(mcp_det_cond$df == 0,"z","t"),"|)", sep=""))
colnames(mtests_det_cond) = c("Estimate", "Std. Error", 
                              ifelse(mcp_det_cond$df==0,"z value","t value"),
                              pname_det_cond)

# RESUME
report(phi_star_lmm)
report(anova(phi_star_lmm))

xtable(anova(phi_star_lmm))
xtable(mtests_det)
xtable(mtests_cond)
xtable(mtests_det_cond)

# xtable(report_table(phi_star_lmm))
# xtable(report_table(anova(phi_star_lmm)))
# xtable(contrast(phi_star_lmm.emm, "eff", by="Detection"))
# xtable(contrast(phi_star_lmm.emm, "eff", by="Condition"))
# xtable(contrast(phi_star.emm, interaction = c("eff", "pairwise"), 
#                 adjust="bonferroni"))

################################################################################
################################################################################
################################################################################
# PHI GEO
################################################################################
################################################################################
################################################################################

phi_geo_eb = summarySE(df, measurevar="Phi_Geo", 
                       groupvars=c("Detection","Condition"))
phi_geo_stats = summarySE(df, measurevar="Phi_Geo", 
                          groupvars=c("Sujet", "Detection", "Condition"))

cond = c("Before", "After")
phi_geo_stats$Condition = factor(phi_geo_stats$Condition, levels = cond)

# Plot Detection * Condition
file_name = "Phi_Geo_sagittal.jpeg"
figure_to_save = str_c(root_path, these_path, file_name, sep="/")
jpeg(file = figure_to_save, width=1209, height=648, units = "px")
par(mar=c(4,4,0,0)+0.2)
phi_geo_sagittal = 
  ggplot(phi_geo_stats, aes(x=Condition, y=Phi_Geo, fill=Detection)) +  
  geom_violin() +
  geom_dotplot(binaxis='y', stackdir='center', position=position_dodge(0.9)) +
  ylab(expression(phi^G)) + 
  theme_bw() + theme(legend.position = c(0.88,0.86)) + theme(text=tt)
phi_geo_sagittal
dev.off()

# Statistical Analysis Detection * Condition
phi_geo_lmm = lme(Phi_Geo ~ Detection * Condition, data=phi_geo_stats, 
                  random = ~1|Sujet, method= "ML")

residuals = resid(phi_geo_lmm)
plot(fitted(phi_geo_lmm), residuals)
abline(0,0)
plot(fitted(phi_geo_lmm), phi_geo_lmm$data$Phi_geo)
qqnorm(residuals)
qqline(residuals)

summary(phi_geo_lmm)
anova(phi_geo_lmm)

emmip(phi_geo_lmm, Detection ~ Condition)

phi_geo_lmm.emm = emmeans(phi_geo_lmm, ~ Detection * Condition)
contrast(phi_geo_lmm.emm, "eff", by="Detection")
contrast(phi_geo_lmm.emm, "eff", by="Condition")
contrast(phi_geo_lmm.emm, interaction=c("eff", "pairwise"), adjust="bonferroni")

phi_geo_lmm = lmer(Phi_Geo ~ Detection * Condition + (1|Sujet), 
                   data=phi_geo_stats)
anova(phi_geo_lmm)

phi_geo.emm = emmeans(phi_geo_lmm, ~ Detection*Condition)

contrast(phi_geo.emm, "eff", by="Detection")
contrast(phi_geo.emm, "eff", by="Condition")
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

# post-hocs Condition 
phi_geo_stats$res = interaction(phi_geo_stats$Condition)
phi_geo_stats$res = factor(phi_geo_stats$res)
phi_geo_mixmodel_cond = lme(Phi_Geo ~ res, random=~1|Sujet, data=phi_geo_stats, 
                            method = "ML")
mcp_cond = glht(phi_geo_mixmodel_cond, linfct=mcp(res="Tukey"))
summary(mcp_cond, test=adjusted(type="bonferroni"))
pq_cond = summary(mcp_cond, test=adjusted(type="bonferroni"))$test
mtests_cond = cbind(pq_cond$coefficients, pq_cond$sigma, pq_cond$tstat, 
                    pq_cond$pvalues)
error_cond = attr(pq_cond$pvalues, "error")
pname_cond = switch(
  mcp_cond$alternativ,
  less = paste("Pr(<", ifelse(mcp_cond$df ==0, "z", "t"), ")", sep=""),
  greater = paste("Pr(>", ifelse(mcp_cond$df == 0, "z", "t"), ")", sep=""),
  two.sided = paste("Pr(>|", ifelse(mcp_cond$df == 0, "z", "t"), "|)", sep=""))
colnames(mtests_cond) = c("Estimate", "Std. Error", 
                          ifelse(mcp_cond$df==0,"z value","t value"),pname_cond)

# post-hocs Detection * Condition 
phi_geo_stats$res = interaction(phi_geo_stats$Detection, phi_geo_stats$Condition)
phi_geo_stats$res = factor(phi_geo_stats$res)
phi_geo_mixmodel_det_cond = lme(Phi_Geo ~ res, random=~1|Sujet, 
                                data=phi_geo_stats, method = "ML")
mcp_det_cond = glht(phi_geo_mixmodel_det_cond, linfct=mcp(res="Tukey"))
summary(mcp_det_cond, test=adjusted(type="bonferroni"))
pq_det_cond = summary(mcp_det_cond, test=adjusted(type="bonferroni"))$test
mtests_det_cond = cbind(pq_det_cond$coefficients, pq_det_cond$sigma, 
                        pq_det_cond$tstat, pq_det_cond$pvalues)
error_det_cond = attr(pq_det_cond$pvalues, "error")
pname_det_cond = switch(
  mcp_det_cond$alternativ,
  less = paste("Pr(<", ifelse(mcp_det_cond$df ==0, "z", "t"), ")", sep=""),
  greater = paste("Pr(>", ifelse(mcp_det_cond$df == 0, "z", "t"), ")", sep=""),
  two.sided = paste("Pr(>|", ifelse(mcp_det_cond$df == 0,"z","t"),"|)", sep=""))
colnames(mtests_det_cond) = c("Estimate", "Std. Error", 
                              ifelse(mcp_det_cond$df==0,"z value","t value"),
                              pname_det_cond)

# RESUME
report(phi_geo_lmm)
report(anova(phi_geo_lmm))

xtable(anova(phi_geo_lmm))
xtable(mtests_det)
xtable(mtests_cond)
xtable(mtests_det_cond)

# xtable(report_table(phi_geo_lmm))
# xtable(report_table(anova(phi_geo_lmm)))
# xtable(contrast(phi_geo_lmm.emm, "eff", by="Detection"))
# xtable(contrast(phi_geo_lmm.emm, "eff", by="Condition"))
# xtable(contrast(phi_geo.emm, interaction = c("eff", "pairwise"), 
#                 adjust="bonferroni"))

################################################################################