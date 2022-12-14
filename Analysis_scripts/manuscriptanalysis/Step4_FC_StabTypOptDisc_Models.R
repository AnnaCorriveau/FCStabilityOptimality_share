# Construct linear models to predict performance using connectome features
# include motion variables and random intercept of dataset

# AC June 2022


# clear 
rm(list=ls())

library(apastats)
library(ggplot2)
library(data.table)
library(lmerTest)
# ========================================================================
# ========================================================================
# Sustained attention runs
df <- read.table('/Users/annacorriveau/Documents/GitHub/FCStability_Typicality/outputs/FC_DOTS_GLMinput_SA_NoNorm.txt')
names <- c('behav', 'stab', 'typ', 'disc', 'opt', 'dataset','saCPM_strength','mot_avg','mot_diff')
setnames(df, new = names)
dataset <- factor(df$dataset)
df5 <- scale(df[1:5])
dfmot <- scale(df[8:9])
df <- cbind(df5, dataset, dfmot)

df <- data.frame(df)

# Delete any NAs so we can compare across all variables
df <- df[!is.na(df$behav),]
df <- df[!is.na(df$stab),]
df <- df[!is.na(df$typ),]
df <- df[!is.na(df$opt),]
df <- df[!is.na(df$disc),]


# does stability predict SA performance across datasets?
model_SA_stab <- lmer(behav ~ stab + mot_avg + mot_diff + (1 | dataset), data=df)
summary(model_SA_stab)

# does typicality predict SA performance across datasets?
model_SA_typ <- lmer(behav ~ typ + mot_avg + mot_diff + (1 | dataset), data=df)
summary(model_SA_typ)

# does optimality predict SA performance across datasets?
model_SA_opt <- lmer(behav ~ opt + mot_avg + mot_diff +  (1 | dataset), data=df)
summary(model_SA_opt)

# does discriminability predict SA performance across datasets?
model_SA_disc <- lmer(behav ~ disc + mot_avg + mot_diff + (1 | dataset), data=df)
summary(model_SA_disc)



model_SA_stab_typ <- lmer(behav ~ stab + typ + mot_avg + mot_diff + (1 | dataset), data=df)
summary(model_SA_stab_typ)

model_SA_stabtyp <- lmer(behav ~ stab * typ + mot_avg + mot_diff + (1 | dataset), data=df)
summary(model_SA_stabtyp)

model_SA_stab_opt <- lmer(behav ~ stab + opt + mot_avg + mot_diff + (1 | dataset), data=df)
summary(model_SA_stab_opt)

model_SA_stabopt <- lmer(behav ~ stab * opt + mot_avg + mot_diff + (1 | dataset), data=df)
summary(model_SA_stabopt)

AIC(model_SA_stab, model_SA_typ, model_SA_opt, model_SA_disc, model_SA_stab_typ, model_SA_stabtyp, model_SA_stab_opt, model_SA_stabopt)

# Are any models significantly better than others?
anova(model_SA_opt, model_SA_stab_opt)
anova(model_SA_opt, model_SA_stabopt)
anova(model_SA_opt, model_SA_stab)


# ========================================================================
# Working memory
# ========================================================================
dw <- read.table('/Users/annacorriveau/Documents/GitHub/FCStability_Typicality/outputs/FC_DOTS_GLMinput_WM_NoNorm.txt')
names <- c('behav', 'stab', 'typ', 'disc', 'opt', 'dataset','wmCPM_strength','mot_avg','mot_diff')
setnames(dw, new = names)
dataset <- factor(dw$dataset)
dw5 <- scale(dw[1:5])
dwmot <- scale(dw[8:9])
dw <- cbind(dw5, dataset, dwmot)
dw <- data.frame(dw)

# Delete any NAs so we can compare across all variables
dw <- dw[!is.na(dw$behav),]
dw <- dw[!is.na(dw$stab),]
dw <- dw[!is.na(dw$typ),]
dw <- dw[!is.na(dw$opt),]
dw <- dw[!is.na(dw$disc),]


# does stability predict WM performance across datasets?
model_WM_stab <- lmer(behav ~ stab + mot_avg + mot_diff + (1 | dataset), data=dw)
summary(model_WM_stab)

# does typicality predict WM performance across datasets?
model_WM_typ <- lmer(behav ~ typ + mot_avg + mot_diff + (1 | dataset), data=dw)
summary(model_WM_typ)


# does optimality predict WM performance across datasets?
model_WM_opt <- lmer(behav ~ opt + mot_avg + mot_diff + (1 | dataset), data=dw)
summary(model_WM_opt)


# does discriminability predict WM performance across datasets?
model_WM_disc <- lmer(behav ~ disc + mot_avg + mot_diff + (1 | dataset), data=dw)
summary(model_WM_disc)



model_WM_stab_typ <- lmer(behav ~ stab + typ + mot_avg + mot_diff + (1 | dataset), data=dw)
summary(model_WM_stab_typ)

model_WM_stabtyp <- lmer(behav ~ stab * typ + mot_avg + mot_diff + (1 | dataset), data=dw)
summary(model_WM_stabtyp)

model_WM_stab_opt <- lmer(behav ~ stab + opt + mot_avg + mot_diff + (1 | dataset), data=dw)
summary(model_WM_stab_opt)

model_WM_stabopt <- lmer(behav ~ stab * opt + mot_avg + mot_diff +  (1 | dataset), data=dw)
summary(model_WM_stabopt)


# Model comparison
AIC(model_WM_stab, model_WM_typ, model_WM_opt, model_WM_disc, model_WM_stab_typ, model_WM_stabtyp, model_WM_stab_opt, model_WM_stabopt)

# Are models significantly better than others?
anova(model_WM_stab, model_WM_stab_opt)
anova(model_WM_stab, model_WM_stabopt)


# ==============================================================================
# =======# How do saCPM strength and FC features compare?=======================
# ==============================================================================
# Reload datasets -- do not need to scale across datasets since we will be modeling within
df <- read.table('/Users/annacorriveau/Documents/GitHub/FCStability_Typicality/outputs/FC_DOTS_GLMinput_SA_NormWithinDataset.txt')
names <- c('behav', 'stab', 'typ', 'disc', 'opt', 'dataset','saCPM_strength','mot_avg','mot_diff')
setnames(df, new = names)

# Delete any NAs so we can compare across all variables
df <- df[!is.na(df$behav),]
df <- df[!is.na(df$stab),]
df <- df[!is.na(df$typ),]
df <- df[!is.na(df$opt),]
df <- df[!is.na(df$disc),]

df1 <- df[df$dataset == 1,]
df2 <- df[df$dataset == 2,]
df3 <- df[df$dataset == 3,]

dw <- read.table('/Users/annacorriveau/Documents/GitHub/FCStability_Typicality/outputs/FC_DOTS_GLMinput_WM_NormWithinDataset.txt')
names <- c('behav', 'stab', 'typ', 'disc', 'opt', 'dataset','wmCPM_strength','mot_avg','mot_diff')
setnames(dw, new = names)
# Delete any NAs so we can compare across all variables
dw <- dw[!is.na(dw$behav),]
dw <- dw[!is.na(dw$stab),]
dw <- dw[!is.na(dw$typ),]
dw <- dw[!is.na(dw$opt),]
dw <- dw[!is.na(dw$disc),]


dw2 <- dw[dw$dataset == 2,]
dw3 <- dw[dw$dataset == 3,]

# Sustained attention DS2
model_SA2_saCPM <- glm(behav~saCPM_strength + mot_avg + mot_diff, data = df2)
summary(model_SA2_saCPM)
model_SA2_stab <- glm(behav~stab + mot_avg + mot_diff, data = df2)
summary(model_SA2_stab)
model_SA2_typ <- glm(behav~typ + mot_avg + mot_diff, data = df2)
summary(model_SA2_typ)
model_SA2_opt <- glm(behav~opt + mot_avg + mot_diff, data = df2)
summary(model_SA2_opt)

model_SA2_saCPM_stab <- glm(behav~saCPM_strength*stab + mot_avg + mot_diff, data = df2)
summary(model_SA2_saCPM_stab)
model_SA2_saCPM_typ <- glm(behav~saCPM_strength*typ + mot_avg + mot_diff, data = df2)
summary(model_SA2_saCPM_typ)
model_SA2_saCPM_opt <- glm(behav~saCPM_strength*opt + mot_avg + mot_diff, data = df2)
summary(model_SA2_saCPM_opt)



# Sustained attention DS3
model_SA3_saCPM <- glm(behav~saCPM_strength + mot_avg + mot_diff, data = df3)
summary(model_SA3_saCPM)
model_SA3_stab <- glm(behav~stab + mot_avg + mot_diff, data = df3)
summary(model_SA3_stab)
model_SA3_typ <- glm(behav~typ + mot_avg + mot_diff, data = df3)
summary(model_SA3_typ)
model_SA3_opt <- glm(behav~opt + mot_avg + mot_diff, data = df3)
summary(model_SA3_opt)

model_SA3_saCPM_stab <- glm(behav~saCPM_strength*stab + mot_avg + mot_diff, data = df3)
summary(model_SA3_saCPM_stab)
model_SA3_saCPM_typ <- glm(behav~saCPM_strength*typ + mot_avg + mot_diff, data = df3)
summary(model_SA3_saCPM_typ)
model_SA3_saCPM_opt <- glm(behav~saCPM_strength*opt + mot_avg + mot_diff, data = df3)
summary(model_SA3_saCPM_opt)


# Working memory DS2
model_WM2_wmCPM <- glm(behav~wmCPM_strength + mot_avg + mot_diff, data = dw2)
summary(model_WM2_wmCPM)
model_WM2_stab <- glm(behav~stab + mot_avg + mot_diff, data = dw2)
summary(model_WM2_stab)
model_WM2_typ <- glm(behav~typ + mot_avg + mot_diff, data = dw2)
summary(model_WM2_typ)
model_WM2_opt <- glm(behav~opt + mot_avg + mot_diff, data = dw2)
summary(model_WM2_opt)

model_WM2_wmCPM_stab <- glm(behav~wmCPM_strength*stab + mot_avg + mot_diff, data = dw2)
summary(model_WM2_wmCPM_stab)
model_WM2_wmCPM_typ <- glm(behav~wmCPM_strength*typ + mot_avg + mot_diff, data = dw2)
summary(model_WM2_wmCPM_typ)
model_WM2_wmCPM_opt <- glm(behav~wmCPM_strength*opt + mot_avg + mot_diff, data = dw2)
summary(model_WM2_wmCPM_opt)



 
#===============================================================================
# Model Identifiability with stability and typicality 
#
df <- read.table('/Users/annacorriveau/Documents/GitHub/FCStability_Typicality/outputs/FC_DOTS_GLMinput_SA_NoNorm.txt')
ID_ATTN <- read.table('/Users/annacorriveau/Documents/GitHub/FCStability_Typicality/outputs/FC_DOTS_ID_ATTN.txt')
ID_unscale <- ID_ATTN$V1
dataset <- factor(df$V6)
df5 <- scale(df[1:5])
dfmot <- scale(df[8:9])
ID <- scale(ID_unscale)
df <- cbind(df5, dataset, dfmot, ID)
df <- data.frame(df)
names <- c('behav', 'stab', 'typ', 'disc', 'opt', 'dataset','mot_avg','mot_diff','ID')
setnames(df, new = names)


model_ID_attn <- lmer(ID ~ stab + typ + (1|dataset), data=df)
summary (model_ID_attn)
 

dw <- read.table('/Users/annacorriveau/Documents/GitHub/FCStability_Typicality/outputs/FC_DOTS_GLMinput_WM_NoNorm.txt')
ID_WM <- read.table('/Users/annacorriveau/Documents/GitHub/FCStability_Typicality/outputs/FC_DOTS_ID_WM.txt')
ID_unscale <- ID_WM$V1
dataset <- factor(dw$V6)
dw5 <- scale(dw[1:5])
dwmot <- scale(dw[8:9])
ID <- scale(ID_unscale)
dw <- cbind(dw5, dataset, dwmot, ID)
dw <- data.frame(dw)
names <- c('behav', 'stab', 'typ', 'disc', 'opt', 'dataset','mot_avg','mot_diff','ID')
setnames(dw, new = names)


model_ID_wm <- lmer(ID ~ stab + typ + (1|dataset), data=dw)
summary(model_ID_wm)
