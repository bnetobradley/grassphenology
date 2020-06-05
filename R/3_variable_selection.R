library(dplyr)
library(lme4)

flr_dat <- read.csv("data/merged_herbarium_climate_na.csv")
flr_dat$binomial <- sub("Deschampsia_caespitosa", "Deschampsia_cespitosa", flr_dat$binomial)
flr_dat <- flr_dat %>% filter(anthers == "anthers")
flr_dat <- flr_dat %>% dplyr::select(-anthers)
mean_flr <- flr_dat %>% group_by(binomial) %>% summarise(sp_mean = mean(doy), sp_sd = sd(doy), sp_sterr = sd(doy)/sqrt(length(doy)))
#write.csv(mean_flr, "data/slopes/mean_flowering.csv", row.names=FALSE)
flr_dat <- full_join(flr_dat, mean_flr, by = c("binomial" = "binomial"))

# create new dataframe with climate variables (prior to flowering) scaled relative to species' mean flowering times ####
## ie. if a species' mean flowering date is in May, 'tave_month_prior' is the average temperature in April
new <- flr_dat %>% mutate(tave_month_prior = case_when(
  sp_mean > 121 & sp_mean < 152 ~ Tave04,
  sp_mean > 151 & sp_mean < 182 ~ Tave05,
  sp_mean > 181 & sp_mean < 213 ~ Tave06,
  sp_mean > 212 & sp_mean < 244 ~ Tave07,
  sp_mean > 243 & sp_mean < 274 ~ Tave08),
  tave_2month_prior = case_when(
    sp_mean > 121 & sp_mean < 152 ~ Tave03,
    sp_mean > 151 & sp_mean < 182 ~ Tave04,
    sp_mean > 181 & sp_mean < 213 ~ Tave05,
    sp_mean > 212 & sp_mean < 244 ~ Tave06,
    sp_mean > 243 & sp_mean < 274 ~ Tave07),
  tave_3month_prior = case_when(
    sp_mean > 121 & sp_mean < 152 ~ Tave02,
    sp_mean > 151 & sp_mean < 182 ~ Tave03,
    sp_mean > 181 & sp_mean < 213 ~ Tave04,
    sp_mean > 212 & sp_mean < 244 ~ Tave05,
    sp_mean > 243 & sp_mean < 274 ~ Tave06),
  tave_4month_prior = case_when(
    sp_mean > 121 & sp_mean < 152 ~ Tave01,
    sp_mean > 151 & sp_mean < 182 ~ Tave02,
    sp_mean > 181 & sp_mean < 213 ~ Tave03,
    sp_mean > 212 & sp_mean < 244 ~ Tave04,
    sp_mean > 243 & sp_mean < 274 ~ Tave05),
  ppt_1month_prior = case_when(
    sp_mean > 121 & sp_mean < 152 ~ PPT04,
    sp_mean > 151 & sp_mean < 182 ~ PPT05,
    sp_mean > 181 & sp_mean < 213 ~ PPT06,
    sp_mean > 212 & sp_mean < 244 ~ PPT07,
    sp_mean > 243 & sp_mean < 274 ~ PPT08),
  ppt_2month_prior = case_when(
    sp_mean > 121 & sp_mean < 152 ~ PPT03,
    sp_mean > 151 & sp_mean < 182 ~ PPT04,
    sp_mean > 181 & sp_mean < 213 ~ PPT05,
    sp_mean > 212 & sp_mean < 244 ~ PPT06,
    sp_mean > 243 & sp_mean < 274 ~ PPT07),
  ppt_3month_prior = case_when(
    sp_mean > 121 & sp_mean < 152 ~ PPT02,
    sp_mean > 151 & sp_mean < 182 ~ PPT03,
    sp_mean > 181 & sp_mean < 213 ~ PPT04,
    sp_mean > 212 & sp_mean < 244 ~ PPT05,
    sp_mean > 243 & sp_mean < 274 ~ PPT06),
  ppt_4month_prior = case_when(
    sp_mean > 121 & sp_mean < 152 ~ PPT01,
    sp_mean > 151 & sp_mean < 182 ~ PPT02,
    sp_mean > 181 & sp_mean < 213 ~ PPT03,
    sp_mean > 212 & sp_mean < 244 ~ PPT04,
    sp_mean > 243 & sp_mean < 274 ~ PPT05),
)


# test which climate variable (relative to species' mean flowering date) has the lowest AIC value ####
aic <- data.frame(value=1:8, variable=c("tavemonth", "tave2month", "tave3month", "tave4month", "pptmonth", "ppt2month", "ppt3month", "ppt4month"))
for (i in 265:273) {
  test <- lmer(data = new, doy ~ new[,i] + (1|binomial))
  aic[i-264,1] <- (AIC(test))
  aic[i-264,3:7] <- (anova(test))
}

tave <- lmer(data = new, doy ~ new$tave_month_prior + (1|binomial))
tave2 <- lmer(data = new, doy ~ new$tave_2month_prior + (1|binomial))
tave3 <- lmer(data = new, doy ~ new$tave_3month_prior + (1|binomial))
ppt <- lmer(data = new, doy ~ new$ppt_1month_prior + (1|binomial))
ppt2 <- lmer(data = new, doy ~ new$ppt_2month_prior + (1|binomial))
ppt3 <- lmer(data = new, doy ~ new$ppt_3month_prior + (1|binomial))
anova(tave, tave2, tave3, ppt, ppt2, ppt3)

# get slope estimates for relationship between flowering and average temperature in the month and two months prior (as well as ppt 2 months prior)
lm_results <- data.frame()
for (i in unique(new$binomial)) {
  speciesx<-filter(new, binomial == i)
  lm_speciesx<-lm(speciesx$doy~speciesx$ppt_2month_prior)
  tmp_res <- broom::tidy(lm_speciesx)
  tmp_res$binomial <- rep(i, 2)
  tmp_res$rsquare <- summary(lm_speciesx)$r.squared
  lm_results <- rbind(lm_results, tmp_res)
}
ppt2month <- filter(lm_results, term == "speciesx$ppt_2month_prior")
write.csv(ppt2month, "data/slopes/ppt_2month_prior_slopes.csv",row.names=FALSE)
