#= = = = = = = = = = = = = = = =
#- - - - - - - - - - - - - - - -
#       PROJET MITIGATE +
#- - - - - - - - - - - - - - - -
#= = = = = = = = = = = = = = = =

#Code created by Thierry Brunelle on 16/11/2023 under behalf of the Pretag Project to produce the paper Economic Mechanisms of Pesticide Use in Cameroon


### -----------------------------------
###            Import packages
### -----------------------------------


library(ggplot2)
library(ivreg)
library(dplyr)
library(readxl)
library(corrplot)
library(car)
library(caTools)
library(lmtest)
library(tseries)
library(ARDL)
library(aTSA)
library(urca)
library(Rmisc)
library(AER)
library(stargazer)
library(xtable)
library(scales)
library(descr)



### -----------------------------------
###            Input data
### -----------------------------------

data_elasticity_estimation <- read_excel("Data-econometric-estim-pest-elasticity.xlsx", 
                                         sheet = "data")
lend=29 #2018
lbeg=6 #1995

PestUseCmr_init<-data_elasticity_estimation$`Pesticides Utilisation (tons)`[lbeg:lend]
PestUseperHaCmr_init<-data_elasticity_estimation$`Pesticide Use Per Ha`[lbeg:lend]
HerbicideUseCmr_init<-data_elasticity_estimation$`Herbicide use tons`[lbeg:lend]
HerbicideUseperHaCmr_init<-data_elasticity_estimation$`Herbicide use per ha (kg/ha)`[lbeg:lend]
InsecticideUseperHaCmr_init<-data_elasticity_estimation$`Insecticide use per ha (kg/ha)`[lbeg:lend]
FongicideUseperHaCmr_init<-data_elasticity_estimation$`Fongicide use per ha (kg/ha)`[lbeg:lend]
NonHerbicideUseperHaCmr_init<-data_elasticity_estimation$NonHerbUseperHa[lbeg:lend]
PestPriceCmr_init<-data_elasticity_estimation$`Pesticide price`[lbeg:lend]
wAgriExportPrice_init<-data_elasticity_estimation$`Weighted Agri Export Price`[lbeg:lend]
ShareChinaPestExp_init<-data_elasticity_estimation$`Share China PesticideExport`[lbeg:lend]
GDPperWorkerCmr_init<-data_elasticity_estimation$`Agriculture, forestry, and fishing, value added per worker (constant 2015 US$)`[lbeg:lend]
ShareLaborAgri_init<-data_elasticity_estimation$`Share labor force in Agriculture, World Bank`[lbeg:lend]
UnemploymentRate_init<-data_elasticity_estimation$`Unemployment % labor force (World Bank)`[lbeg:lend]/100
dummy_legislation <- data_elasticity_estimation$Dummy_legislation[lbeg:lend]
Year <-  data_elasticity_estimation$Year[lbeg:lend]

PestUseCmr<-scale(log(PestUseCmr_init))
PestUseperHaCmr<-scale(log(PestUseperHaCmr_init))
HerbicideUseCmr<-scale(log(HerbicideUseCmr_init))
HerbicideUseperHaCmr<-scale(log(HerbicideUseperHaCmr_init))
InsecticideUseperHaCmr<-scale(log(InsecticideUseperHaCmr_init))
FongicideUseperHaCmr<-scale(log(FongicideUseperHaCmr_init))
NonHerbicideUseperHaCmr<-scale(log(NonHerbicideUseperHaCmr_init))
PestPriceCmr<-scale(log(PestPriceCmr_init))
wAgriExportPrice<-scale(wAgriExportPrice_init)
ShareChinaPestExp<-scale(log(ShareChinaPestExp_init))
GDPperWorkerCmr<-log(GDPperWorkerCmr_init)
ShareLaborAgri<-log(ShareLaborAgri_init)
UnemploymentRate<-log(UnemploymentRate_init)


TAB<-data.frame(PestUseCmr,PestUseperHaCmr,HerbicideUseCmr, HerbicideUseperHaCmr, NonHerbicideUseperHaCmr, PestPriceCmr, wAgriExportPrice, ShareChinaPestExp, GDPperWorkerCmr, ShareLaborAgri, UnemploymentRate, dummy_legislation, Year)
TAB_init<-data.frame(PestUseCmr_init, PestUseperHaCmr_init*1000, HerbicideUseperHaCmr_init*1000, InsecticideUseperHaCmr_init*1000, FongicideUseperHaCmr_init*1000, PestPriceCmr_init, wAgriExportPrice_init, ShareChinaPestExp_init, UnemploymentRate_init)
Xvar<-cbind(PestPriceCmr,wAgriExportPrice,UnemploymentRate,dummy_legislation)


### -----------------------------------
###       Data description
### -----------------------------------

summary(TAB)
cor(Xvar) -> cor_matrice
corrplot(cor_matrice,type="lower") 
descr(TAB_init)

# Test stationarité : Kwiatkowski-Phillips-Schmidt-Shin test
kpss.test(PestUseperHaCmr) 
kpss.test(HerbicideUseperHaCmr)
kpss.test(NonHerbicideUseperHaCmr)
kpss.test(PestPriceCmr)
kpss.test(wAgriExportPrice)
kpss.test(UnemploymentRate)

# Test stationarité : Augmented Dickey Fuller test
adf.test(PestUseperHaCmr) 
adf.test(HerbicideUseperHaCmr)
adf.test(NonHerbicideUseperHaCmr)
adf.test(PestPriceCmr)
adf.test(wAgriExportPrice)
adf.test(UnemploymentRate)



### -----------------------------------
###       Estimation élasticités 
### -----------------------------------


#######################################
### Model 1: Pesticide total use per Ha
#######################################
reg_ols_mod1 <- lm(PestUseperHaCmr~PestPriceCmr+wAgriExportPrice+UnemploymentRate+dummy_legislation, data=TAB)
summary(reg_ols_mod1)
AIC(reg_ols_mod1)

# Alternative specification
reg_ols_mod1b <- lm(PestUseperHaCmr~PestPriceCmr, data=TAB)
summary(reg_ols_mod1b)
AIC(reg_ols_mod1b)

# Mod 1: Test multicollinéarité
vif(reg_ols_mod1)
## ---> ok

# Mod 1: Test homoscedasticité
bptest(reg_ols_mod1)
## ---> p-value = 0.03547 >>> presence d'heteroscedaticité

# Weighted least square regression 
#wt <- 1 / lm(abs(reg_ols_mod1$residuals) ~ reg_ols_mod1$fitted.values)$fitted.values^2
#wt <- 1 / reg_ols_mod1$residuals^2
#reg_wls_mod1 <- lm(PestUseperHaCmr~PestPriceCmr+wAgriExportPrice+UnemploymentRate+dummy_legislation, weights = wt, data=TAB)
#summary(reg_wls_mod1)
#AIC(reg_wls_mod1)


# Cointegration test avec 1 lag pour données annuelles
adf.test(reg_ols_mod1$residuals, nlag=1)
## --> pvalue <= 0.01 --> les résidus sont stationnaires, relation de cointégration attestée

# Mod 1: Régression avec variable instrumentale 
# Relation de court terme
reg_iv_mod1_ST <- ivreg(PestUseperHaCmr~PestPriceCmr+lag(PestUseperHaCmr)+wAgriExportPrice+UnemploymentRate+dummy_legislation | ShareChinaPestExp+lag(PestUseperHaCmr)+wAgriExportPrice+UnemploymentRate+dummy_legislation, data=TAB )
summary(reg_iv_mod1_ST, diagnostics=TRUE)
# Relation de long terme
reg_iv_mod1_LT <- ivreg(PestUseperHaCmr~PestPriceCmr+wAgriExportPrice+UnemploymentRate+dummy_legislation| ShareChinaPestExp+wAgriExportPrice+UnemploymentRate+dummy_legislation, data=TAB )
summary(reg_iv_mod1_LT, diagnostics=TRUE)
confint(reg_iv_mod1_LT)


# Hausman test rejeté (p-value > 0.05). L'estimateur ols est plus consistent que IV 
reg_ols_mod1_ST <- lm(PestUseperHaCmr~lag(PestUseperHaCmr)+PestPriceCmr+wAgriExportPrice+UnemploymentRate+dummy_legislation, data=TAB)
summary(reg_ols_mod1_ST)

# Newey West correction
reg_ols_mod1_LT<-reg_ols_mod1
coef_mod1_ST<-coeftest(reg_ols_mod1_ST, vcov=NeweyWest(reg_ols_mod1_ST, verbose=T))
coef_mod1_LT<-coeftest(reg_ols_mod1_LT, vcov=NeweyWest(reg_ols_mod1_LT, verbose=T))
# Même coefficient que OLS, je peux utiliser le même R2

#newvar=cbind(x1 = PestPriceCmr[1], x2 = c(24:1))
#newd1 <- data.frame(PestPriceCmr,wAgriExportPrice,UnemploymentRate,dummy_legislation)
#newd2 <- data.frame(newd2 <- data.frame(log(data_elasticity_estimation$`Pesticide price`[lbeg])*rep(1,24),wAgriExportPrice,UnemploymentRate,dummy_legislation),wAgriExportPrice,UnemploymentRate,dummy_legislation)
#exp(predict(reg_ols_mod1 , newd1))
#exp(predict(reg_ols_mod1 , newd2))

        

#######################################
### Model 2: Herbicide Use per Hectare
#######################################
reg_ols_mod2 <- lm(HerbicideUseperHaCmr~PestPriceCmr+wAgriExportPrice+UnemploymentRate+dummy_legislation, data=TAB)
summary(reg_ols_mod2)
AIC(reg_ols_mod2)

# Alternative specification 
reg_ols_mod2b <- lm(HerbicideUseCmr~UnemploymentRate, data=TAB)
summary(reg_ols_mod2b)
AIC(reg_ols_mod2b)

# Mod2: Test multicollinéarité
vif(reg_ols_mod2)
## ---> ok

# Mod 2: Test homoscedasticité
bptest(reg_ols_mod2)
## ---> ok

# Cointegration test
kpss.test(reg_ols_mod2$residuals)
adf.test(reg_ols_mod2$residuals)
##  pvalue <= 0.01 --> les résidus sont stationnaires, relation de cointégration attestée

# Mod2: Régression avec variable instrumental
# Relation de court terme
reg_iv_mod2_ST <- ivreg(HerbicideUseperHaCmr~lag(HerbicideUseperHaCmr)+PestPriceCmr+wAgriExportPrice+UnemploymentRate+dummy_legislation | lag(HerbicideUseperHaCmr)+ShareChinaPestExp+wAgriExportPrice+UnemploymentRate+dummy_legislation, data=TAB)
summary(reg_iv_mod2_ST, diagnostics=TRUE)
# Relation de long terme
reg_iv_mod2_LT <- ivreg(HerbicideUseperHaCmr~PestPriceCmr+wAgriExportPrice+UnemploymentRate+dummy_legislation  | ShareChinaPestExp+wAgriExportPrice+UnemploymentRate+dummy_legislation, data=TAB)
summary(reg_iv_mod2_LT, diagnostics=TRUE)
# Hausman test n'est pas rejeté (p-value < 0.05). L'estimateur IV est plus consistent que ols. 

# Regression de court terme pour comparaison 
reg_ols_mod2_ST <- lm(HerbicideUseperHaCmr~lag(HerbicideUseperHaCmr)+PestPriceCmr+wAgriExportPrice+UnemploymentRate+dummy_legislation, data=TAB)
summary(reg_ols_mod2_ST)
vif(reg_ols_mod2_ST)
reg_ols_mod2_LT<-reg_ols_mod2

#########################################
## Model 3: Insecticide Use per Hectare
#######################################
reg_ols_mod3 <- lm(InsecticideUseperHaCmr~PestPriceCmr+wAgriExportPrice+UnemploymentRate+dummy_legislation, data=TAB)
summary(reg_ols_mod3)
AIC(reg_ols_mod3)

# Alternative specification
reg_ols_mod3b <- lm(InsecticideUseperHaCmr~PestPriceCmr, data=TAB)
summary(reg_ols_mod3b)
AIC(reg_ols_mod3b)

# Mod 3: Test multicollinéarité
vif(reg_ols_mod3)

# Mod 3: Test homoscedasticité
bptest(reg_ols_mod3)
# p-value = 0.01378 >> présence d'hétéroscedasticité

# # Weighted least square regression 
# wt <- 1 / lm(abs(reg_ols_mod3$residuals) ~ reg_ols_mod3$fitted.values)$fitted.values^2
# #wt <- 1 / reg_ols_mod3$residuals^2
# reg_wls_mod3 <- lm(InsecticideUseperHaCmr~PestPriceCmr+wAgriExportPrice+UnemploymentRate+dummy_legislation, weights = wt, data=TAB)
# summary(reg_wls_mod3)
# AIC(reg_wls_mod3)

# Cointegration test
adf.test(reg_ols_mod3$residuals, nlag =1)
##  pvalue <= 0.01 --> les résidus sont stationnaires, relation de cointégration attestée

# Mod 3: Régression avec variable instrumentale 
# Relation de court terme
reg_iv_mod3_ST <- ivreg(InsecticideUseperHaCmr~PestPriceCmr+lag(InsecticideUseperHaCmr)+wAgriExportPrice+UnemploymentRate+dummy_legislation | ShareChinaPestExp+lag(InsecticideUseperHaCmr)+wAgriExportPrice+UnemploymentRate+dummy_legislation, data=TAB)
summary(reg_iv_mod3_ST, diagnostics=TRUE)
# Hausman ok (p-value < 0.05). L'estimateur IV est plus consistent que OLS

# Relation de long terme
reg_iv_mod3_LT <- ivreg(InsecticideUseperHaCmr~PestPriceCmr+wAgriExportPrice+UnemploymentRate+dummy_legislation | ShareChinaPestExp+wAgriExportPrice+UnemploymentRate+dummy_legislation, data=TAB)
summary(reg_iv_mod3_LT, diagnostics=TRUE)
# Hausman test rejeté (p-value > 0.05). L'estimateur ols est plus consistent que IV 

reg_ols_mod3_ST <- lm(InsecticideUseperHaCmr~lag(InsecticideUseperHaCmr)+PestPriceCmr+wAgriExportPrice+UnemploymentRate+dummy_legislation, data=TAB)
summary(reg_ols_mod3_ST)

# Newey West correction
reg_ols_mod3_LT<-reg_ols_mod3
coef_mod3_ST<- coeftest(reg_iv_mod3_ST, vcov=NeweyWest(reg_iv_mod3_ST, verbose=T))
coef_mod3_LT <- coeftest(reg_ols_mod3_LT, vcov=NeweyWest(reg_ols_mod3_LT, verbose=T))
# Même coefficient que OLS, je peux utiliser le même R2


#########################################
## Model 4: Fongicide Use per Hectare
#######################################
reg_ols_mod4 <- lm(FongicideUseperHaCmr~PestPriceCmr+wAgriExportPrice+UnemploymentRate+dummy_legislation , data=TAB)
summary(reg_ols_mod4)
AIC(reg_ols_mod4)

# Alternative specification
reg_ols_mod4b <- lm(FongicideUseperHaCmr~PestPriceCmr, data=TAB)
summary(reg_ols_mod4b)
AIC(reg_ols_mod4b)

# Mod 3: Test multicollinéarité
vif(reg_ols_mod4)

# Mod 3: Test homoscedasticité
bptest(reg_ols_mod4)
# pvalue = 0.09017 >> on ne rejette pas H0 (homoscedasticité)

# Cointegration test
adf.test(reg_ols_mod4$residuals, nlag =1)
##  pvalue <= 0.01 --> les résidus sont stationnaires, relation de cointégration attestée

# Mod 3: Régression avec variable instrumentale 
# Relation de court terme
reg_iv_mod4_ST <- ivreg(FongicideUseperHaCmr~PestPriceCmr+lag(FongicideUseperHaCmr)+wAgriExportPrice+UnemploymentRate+dummy_legislation | ShareChinaPestExp+lag(FongicideUseperHaCmr)+wAgriExportPrice+UnemploymentRate+dummy_legislation, data=TAB)
summary(reg_iv_mod4_ST, diagnostics=TRUE)
# Relation de long terme
reg_iv_mod4_LT <- ivreg(FongicideUseperHaCmr~PestPriceCmr+wAgriExportPrice+UnemploymentRate+dummy_legislation | ShareChinaPestExp+wAgriExportPrice+UnemploymentRate+dummy_legislation, data=TAB)
summary(reg_iv_mod4_LT, diagnostics=TRUE)

# Hausman test rejeté (p-value > 0.05). L'estimateur ols est plus consistent que IV 
reg_ols_mod4_ST <- lm(FongicideUseperHaCmr~lag(FongicideUseperHaCmr)+PestPriceCmr+wAgriExportPrice+UnemploymentRate+dummy_legislation, data=TAB)
summary(reg_ols_mod4_ST)
reg_ols_mod4_LT<-reg_ols_mod4



#########################################
## Model 4: Robustness test on lag agricultural price
#######################################
test_iv_mod1_LT <- ivreg(PestUseperHaCmr~PestPriceCmr+lag(wAgriExportPrice)+UnemploymentRate+dummy_legislation| ShareChinaPestExp+lag(wAgriExportPrice)+UnemploymentRate+dummy_legislation, data=TAB )
summary(test_iv_mod1_LT, diagnostics=TRUE)
test_iv_mod2_LT <- ivreg(HerbicideUseperHaCmr~PestPriceCmr+lag(wAgriExportPrice)+UnemploymentRate+dummy_legislation  | ShareChinaPestExp+lag(wAgriExportPrice)+UnemploymentRate+dummy_legislation, data=TAB)
summary(test_iv_mod2_LT, diagnostics=TRUE)
test_ols_mod3 <- lm(InsecticideUseperHaCmr~PestPriceCmr+lag(wAgriExportPrice)+UnemploymentRate+dummy_legislation, data=TAB)
summary(test_ols_mod3, diagnostics=TRUE)
test_ols_mod4 <- lm(FongicideUseperHaCmr~PestPriceCmr+lag(wAgriExportPrice)+UnemploymentRate+dummy_legislation , data=TAB)
summary(test_ols_mod4, diagnostics=TRUE)


#########################################
## Export des résultats
#######################################

olsm2ST<- reg_ols_mod2_ST
ivm2ST<-reg_iv_mod2_ST
coefm3ST<-coef_mod3_ST
ivm3ST<-reg_iv_mod3_ST
olsm4ST<-reg_ols_mod4_ST
ivm4ST <- reg_iv_mod4_ST
table_ST<-stargazer(olsm2ST, ivm2ST , coefm3ST, ivm3ST, olsm4ST , ivm4ST) 
olsm2LT<- reg_ols_mod2_LT
ivm2LT<-reg_iv_mod2_LT
coefm3LT<-coef_mod3_LT
ivm3LT<-reg_iv_mod3_LT
olsm4LT<-reg_ols_mod4_LT
ivm4LT <- reg_iv_mod4_LT
table_LT<-stargazer(olsm2LT, ivm2LT , coef_mod3_LT, ivm3LT, olsm4LT , ivm4LT) 
table_TotPest<-stargazer(coef_mod1_ST,reg_iv_mod1_ST,coef_mod1_LT,reg_iv_mod1_LT)


fs <- lm(FongicideUseperHaCmr~PestPriceCmr+wAgriExportPrice+UnemploymentRate+dummy_legislation+lag(FongicideUseperHaCmr), data=TAB)
fn <- lm(FongicideUseperHaCmr~lag(FongicideUseperHaCmr)+wAgriExportPrice+UnemploymentRate+dummy_legislation, data=TAB)
waldtest(fs, fn, vcov = vcovHC(fs, type="HC0"))$F[2]
fs <- lm(HerbicideUseperHaCmr~PestPriceCmr+wAgriExportPrice+UnemploymentRate+dummy_legislation, data=TAB)
fn <- lm(HerbicideUseperHaCmr~wAgriExportPrice+UnemploymentRate+dummy_legislation, data=TAB)
waldtest(fs, fn, vcov = vcovHC(fs, type="HC0"))$F[2]


#### ---- Figures

pctChinaPestExp<-data_elasticity_estimation$`Share China PesticideExport`[lbeg:lend]
RealPestPriceCmr<-data_elasticity_estimation$`Pesticide price`[lbeg:lend]
df1<-data.frame(Year,pctChinaPestExp,RealPestPriceCmr )

filename <- "ShareChinaPestExp.png"
png(filename,width = 750, height = 600)
ggplot(data=df1, aes(x=Year, y=pctChinaPestExp)) +
  geom_line(color="red",size =1)+
  scale_y_continuous(labels = percent_format(accuracy = 1))+
  labs(x="Year",y=" ")+
  theme_minimal()+
  theme(axis.text.x = element_text(face="bold", size=14),
        axis.text.y = element_text(face="bold", size=14),
        axis.title=element_text(size=18,face="bold"))
dev.off()

filename <- "PesticidesUseperHa.png"
png(filename,width = 750, height = 600)
ggplot(data=df1, aes(x=Year, y=data_elasticity_estimation$`Pesticide Use Per Ha`[lbeg:lend])) +
  geom_line(color="red",size =1)+
  labs(x="Year",y="Kg AI per ha")+
  theme_minimal()+
  theme(axis.text.x = element_text(face="bold", size=18),
        axis.text.y = element_text(face="bold", size=18),
        axis.title=element_text(size=18,face="bold"))
dev.off()


filename <- "PesticidesUseTot.png"
png(filename,width = 750, height = 600)
ggplot(data=df1, aes(x=Year, y=data_elasticity_estimation$`Pesticides Utilisation (tons)`[lbeg:lend])) +
  geom_line(color="red",size =1)+
  labs(x="Year",y="Tons per year")+
  theme_minimal()+
  theme(axis.text.x = element_text(face="bold", size=18),
        axis.text.y = element_text(face="bold", size=18),
        axis.title=element_text(size=18,face="bold"))
dev.off()


filename <- "PestPriceCameroun.png"
png(filename,width = 750, height = 600)
ggplot(data=df1, aes(x=Year, y=RealPestPriceCmr)) +
  geom_line(color="blue", size =1)+
  labs(x="Year",y="$/kg of FP")+
  theme_minimal()+
  theme(axis.text.x = element_text(face="bold", size=18),
        axis.text.y = element_text(face="bold", size=18),
        axis.title=element_text(size=18,face="bold"))
dev.off()
