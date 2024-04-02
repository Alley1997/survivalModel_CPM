# Setting Working Directory -----------------------------------------------
setwd("/Users/alleythawani/Documents/MSc_HDS/Stats_Modelling/Assessment")


# Loading Libraries -------------------------------------------------------

library(tidyverse)
library(janitor)
library(survival)
library(gtsummary)
library(grid)
library(gridExtra)
library(rms)
library(stats)
library(survival)
library(conflicted)
library(DescTools)
library(survminer)

# Data Loading, Cleaning and Preparation ----------------------------------

cvd.df <- read_csv("cvddata.csv") |> 
  mutate(
        ID = as.character(ID),
        Sex = as.factor(Sex),
        Smoking_Status = as.factor(Smoking_Status),
        Diabetic = as.factor(Diabetes),
        CKD = as.factor(CKD),
        atrial_fibrillation = as.factor(AF),
        TEVENT = TEVENT
        )|> 
          select(-c(AF,ID,Diabetes)) |> 
  clean_names() #use globally accepted valiable naming standards

# Summary characteristics  --------------------------------------------------------

cvd.df |> select(
  c(sex,smoking_status,diabetic,status,ckd,atrial_fibrillation)
) |> 
  mutate(
    sex = ifelse(sex == 0,"female","male"),
    smoking_status  = ifelse(smoking_status == 0,"no-smoker","smoker/ever_smoked"),
    diabetic = ifelse(diabetic == 0,"no","yes"),
    chronic_kidney_disease = ifelse( ckd == 0,"no","yes"),
    atrial_fibrillation = ifelse(atrial_fibrillation == 0,"no","yes"),
    status = ifelse(status == 0,"no event/ censored","event")
  ) |> 
  tbl_summary(by = status,
              type = all_dichotomous() ~ "categorical",
              statistic = list(
                all_categorical() ~ "{n} ({p})%"
              )) |> 
  add_overall()
              

#Distribution of Continuous variables

grid.arrange(cvd.df |> ggplot( aes(x = age)) +
               geom_histogram(binwidth = 1, fill = "black", color = "white") +
               labs(title = "Age Distribution") +
               theme(plot.title = element_text(hjust = 0.5)),
             
             cvd.df |> ggplot( aes(x = bmi)) +
               geom_histogram(binwidth = 1, fill = "black", color = "white") +
               labs(title = "BMI Distribution",
                    x = "Body Mass Index") +
               theme(plot.title = element_text(hjust = 0.5)) +
               theme(plot.title = element_text(hjust = 0.5)),ncol = 2, heights = c(1, 0.2, 0))

# Cox PH Modelling with RMS package ------------------------------------------------------------
units(cvd.df$tevent) <- "years"
dd = datadist(cvd.df)
options(datadist='dd')

surv.obj <- Surv(cvd.df$tevent,cvd.df$status)

model_simple <- cph(surv.obj ~ age + sex + smoking_status + bmi + diabetic + ckd + atrial_fibrillation,
             data =cvd.df ,
             surv = TRUE,
             x = TRUE, 
             y = TRUE,
             method = 'breslow',
             time.inc = 5)

summary(model_simple)
# Proportional Hazard Test Simple Model ------------------------------------------------
ph.test.simple <- cox.zph(model_simple)

ph.test.simple

# Linearity Test for age --------------------------------------------------

normal_age_model <- cph(surv.obj ~ age,
                        data = cvd.df ,
                        surv = TRUE,
                        x = TRUE, 
                        y = TRUE,
                        method = 'breslow')

rcs_age_model <- cph(surv.obj ~ rcs(age,3),
                     data = cvd.df ,
                     surv = TRUE,
                     x = TRUE, 
                     y = TRUE,
                     method = 'breslow')

lr_test  = lrtest(normal_age_model,rcs_age_model)

# Fitting Complex Model ---------------------------------------------------
model_complex <- cph(surv.obj ~ rcs(age,3) + strat(sex) + smoking_status + bmi + diabetic + ckd + atrial_fibrillation,
             data =cvd.df ,
             surv = TRUE,
             x = TRUE, 
             y = TRUE,
             method = 'breslow',
             time.inc = 5)

# Proportional Hazard Test Simple Model -----------------------------------

ph.test.complex <- cox.zph(model_complex)

ph.test.complex


# AIC TESTS ---------------------------------------------------------------

AIC(model_simple)
AIC(model_complex)
#AIC(final_model_complex)


# STEP WISE MODEL SELLECTION ----------------------------------------------

final <- fastbw(model_complex, rule = "aic")


# FINAL OPTIMAL MODEL -----------------------------------------------------

final_model_complex <- cph(surv.obj ~ rcs(age,3) + strat(sex) + bmi + diabetic + ckd,
                     data = cvd.df ,
                     surv = TRUE,
                     x = TRUE, 
                     y = TRUE,
                     method = 'breslow',
                     time.inc = 5)

adjusted_hazard_ratios <- exp(final_model_complex$coefficients)
ahr_conf_int <- exp(confint(final_model_complex))
f <- update(final_model_complex,x=TRUE,y=TRUE)

# INTERNAL VALIDATION WITH BOOTSTRAP METHOD -------------------------------

set.seed(123)  # Set seed for reproducibility

validation_results_boot <- validate(final_model_complex,method = "boot", B = 200,u = 5)

reds <- predict(f, type = "lp", times = 5)

c_index <- rcorr.cens(preds,surv.obj)

c_indexx <- c_index['C Index']
se <- c_index['S.D.']/2

low <- c_indexx-1.96*se
hi <- c_indexx +1.96*se
# INTERNAL CALIBRATION WITH BOOTSTRAP METHOD ------------------------------

calibrate_results_boot <- calibrate(final_model_complex,method = "boot",dxy = T, B = 200,u = 5)

plot(calibrate_results_boot,las = 1)


# MODEL PRESENTATION WITH NOMOGRAM ----------------------------------------
nom <- nomogram(final_model_complex)

plot(nom)



#####Nomogram

surv <- Survival(final_model_complex)

surv2y <- function(x) surv(5, lp = x)

ss <- c(-2, -1,0.05, 0.2, 0.4, 0.6, 0.7, 0.8, 0.9, 1,2)

nomo_gram <- nomogram(final_model_complex, fun = surv2y, fun.at = ss, lp = TRUE, funlabel = "5 year survival")

plot(nomo_gram)



survdiff(surv.obj~sex,data = cvd.df) |> plot()

cox_model <- survfit(surv.obj ~ sex, data = cvd.df)

ggsurvplot(cox_model, data = cvd.df, risk.table = TRUE)

# Create a survival plot using survplot
survplot(cox_model, data = cvd.df, col = c("blue", "red"), lty = c(1, 2), lwd = c(2, 2), xlab = "Time", ylab = "Survival Probability", main = "Survival Curve by Sex")

# Linearity Test for age --------------------------------------------------

normal_age_model <- cph(surv.obj ~ bmi,
                        data = cvd.df ,
                        surv = TRUE,
                        x = TRUE, 
                        y = TRUE,
                        method = 'breslow')

rcs_age_model <- cph(surv.obj ~ rcs(bmi,3),
                     data = cvd.df ,
                     surv = TRUE,
                     x = TRUE, 
                     y = TRUE,
                     method = 'breslow')

lr_test  = lrtest(normal_age_model,rcs_age_model)
