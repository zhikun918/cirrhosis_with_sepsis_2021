# =======Load required r packages
library(tidyverse)
library(haven)
library(gtsummary)
library(ggplot2)
library(rms)
library(glmnet)
library(MASS)
library(pROC)
library(plotROC)

1# ======Set working directory
setwd("C:/Users/Xu/Desktop/R工作站/DXY-CPM-Case1")

# =======Read the data
# ===Derivation cohort
gustoW<- read_sav("./01-Data/GustoW.sav") # library(haven)
# ===Validation cohort
gustoS<- read_sav("./01-Data/sample4.sav") # library(haven)


# =======1. Check the data
dim(gustoW)
dim(gustoS)

head(gustoW)
head(gustoS)

summary(gustoW)
summary(gustoS)

str(gustoW)
str(gustoS)

#========Convert the data type

fctlist <- colnames(gustoW [, -c(3, 12, 13, 19)]) # 筛选分类变量
gustoW[fctlist] <- lapply(gustoW[fctlist], factor) # 转为分类变量

fctlist <- colnames(gustoS [, -c(3, 12, 13, 19)])
gustoS[fctlist] <- lapply(gustoS[fctlist], factor)

# ======Check data type
summary(gustoW)
summary(gustoS)


# ========2. describe the data

describe(gustoW)
describe(gustoS)

# =========显示更好 library(gtsummary)无法安装此包 或类似与tableone
gustoW %>%  
  tbl_summary(by = DAY30)

gustoS %>%  
  tbl_summary(by = DAY30)

# ======连续变量，作图
ggplot(gustoW, aes(x = AGE, fill = DAY30)) +
  geom_histogram()

ggplot(gustoW, aes(x = HEI, fill = DAY30)) +
  geom_histogram()

ggplot(gustoW, aes(x = WEI, fill = DAY30)) +
  geom_histogram()


# ======3. Explore association

# Nonlinear association between weight and death? 

dddt <- datadist(gustoW) # 加载rms包
options(datadist="dddt")

# ====age
# lmfit <- lrm(DAY30 ~ AGE, data = gustoW)

rcsfit <- lrm(DAY30 ~ rcs(AGE, 4), data = gustoW)

rcsOR <- Predict(rcsfit, AGE, fun=exp,  ref.zero = TRUE)
ggplot(rcsOR)
anova(rcsfit) # age非线性；p小于0.05 age：65截断值

# ===plot OR
# ggplot(rcsOR)
# ggplot(rcsOR)+
#  geom_line(aes(x = WEI, y = yhat), color = "red") +
#  geom_ribbon(aes(x= WEI, ymin = lower, ymax = upper) , alpha = 0.1)

# === height 
lmfit <- lrm(DAY30 ~ HEI, data = gustoW)
rcsfit <- lrm(DAY30 ~ rcs(HEI, 4), data = gustoW)
rcsOR<-Predict(rcsfit, HEI, fun=exp,  ref.zero = TRUE)
ggplot(rcsOR)
anova(rcsfit) # hei线性；p大于0.05

# === weight
lmfit <- lrm(DAY30 ~ WEI, data = gustoW)
rcsfit <- lrm(DAY30 ~ rcs(WEI, 3), data = gustoW)
rcsOR<-Predict(rcsfit, WEI, fun=exp,  ref.zero = TRUE)
ggplot(rcsOR)
anova(rcsfit) # wei非线性；p小于0.05


gustoW <- mutate(gustoW, bmi = WEI/(HEI/100)**2) # tidyverse包生存新的列

dddt <- datadist(gustoW)
options(datadist="dddt")

lmfit <- lrm(DAY30 ~ bmi, data = gustoW)
rcsfit <- lrm(DAY30 ~ rcs(bmi, 3), data = gustoW)
rcsOR<-Predict(rcsfit, bmi, fun=exp,  ref.zero = TRUE)
ggplot(rcsOR)
anova(rcsfit) # bmi非线性；p小于0.05 bmi：30截断值

gustoW <- mutate(gustoW, bmi30up = ifelse(bmi >= 30, 1,0))

# ===== selected variables

str(gustoW)
gustoW <- gustoW[,-c(3, 5, 12, 13, 19, 22:27)]

# === new var
gustoS <- mutate(gustoS, bmi = WEI/(HEI/100)**2)
gustoS <- mutate(gustoS, bmi30up = ifelse(bmi >= 30, 1,0))

gustoS<- gustoS[,-c(3, 5, 12, 13, 19, 22:27)] # 删除不要的变量

str(gustoW)
str(gustoS)


# ====4. Predictor selection

# === define formula
outcome = "DAY30"
xlist = colnames(gustoW[,-1])

formula <- formula(paste(paste(outcome,"~", collapse=" "), 
                         paste(xlist, collapse=" + ")))
formula


# === fit full model
lrm_full <- lrm(formula, data = gustoW ) # 加载rms包
lrm_full

# 或者 方便展示 推荐glm法 gtsummary包
glm_full <- glm(formula, data = gustoW, family = binomial)
glm_full %>% 
  tbl_regression(exponentiate =  TRUE) 
glm_full


# ===selection
# =====backward   相对简单

glm_back <- stepAIC(glm_full, direction="backward")
glm_back %>% 
  tbl_regression(exponentiate =  TRUE) 

# ====LASSO 高级方法
tmp_y <- gustoW$DAY30
xlist = colnames(gustoW[,-1])
tmp_x <- model.matrix(~.,gustoW[xlist]) 

model_lasso <-  glmnet(tmp_x, tmp_y, family="binomial", nlambda=30, alpha=1, standardize=TRUE)
plot(model_lasso,xvar="lambda",label=TRUE)

model_lasso #显示lamada

# =====find the optimal model via cross-validation 筛选最优lamada
glm_cv <- cv.glmnet(tmp_x, tmp_y, family="binomial", nlambda=30, alpha=1, standardize=TRUE)
plot(glm_cv)
glm_cv$lambda.min
coef(glm_cv, s=glm_cv$lambda.min) 


# increase lambda for further shrinkage 7个变量
glm_cv$lambda.1se
coef(glm_cv, s=glm_cv$lambda.1se) 


# ====== 4. Fit the final model

outcome = "DAY30"
xlist = c("SEX", "A65", "SHO", "HYP","HRT", "ANT", "PMI") # lasso筛选出的7个变量

formula <- formula(paste(paste(outcome,"~", collapse=" "), 
                         paste(xlist, collapse=" + ")))

# ===lrm way 方法一
cpm_final_lrm <- lrm(formula, data = gustoW , x = TRUE, y = TRUE)
summary(cpm_final_lrm)

# ===glm way 方法二；方便显示
cpm_final_glm <- glm(formula, data = gustoW, family="binomial")
cpm_final_glm %>% 
  tbl_regression(exponentiate =  TRUE) 

# ====5. display the prediction model
pdf("./04-Out/nomogram1.pdf")
plot(nomogram(cpm_final_lrm, fun = plogis, lp = FALSE), xfrac=.15, cex.axis=.7)
dev.off()


# ==== 6. internal validation 内部验证

# C statistics C 统计量
cv <- validate(cpm_final_lrm, method="crossvalidation", B=10)

# C-stat = (Dxy+1)/2

str(cv)
# C-stat 
(cv[1, 1] +1)/2
(cv[1, 5] +1)/2 # 矫正后


# === Calibration plot  校准曲线
pdf("./04-Out/calibration plot1.pdf")
plot(calibrate(cpm_final_lrm, method="crossvalidation", B=10))
dev.off()

# ==== 7. external validation 外部验证

# ==== linear predictor 
gustoS$lp <- predict(cpm_final_glm, newdata = gustoS, type="link") # newdata设置为外部验证数据集

# ====predicted probability
gustoS$phat <- predict(cpm_final_glm, newdata = gustoS, type="response")  



# ==== 7.1 C statistics C 统计量
roc_ext <- roc(gustoS$DAY30, gustoS$phat)

auc_ext<- roc_ext$auc 

# ===AUC and 95% CI roc曲线
auc_ext
ci.auc(roc_ext) 

pdf("./04-Out/roc plot1.pdf")
plot(roc_ext)
dev.off()



# ==== 7.3 Calibration 校准

# === Brier score 
# gustoS$DAY30 <- as.factor(gustoS$DAY30) 转为数字变量（1/2）
gustoS$DAY30 <- as.integer(gustoS$DAY30) -1
brier_ext <- mean((gustoS$DAY30 - gustoS$phat)^2)
brier_ext # 越小越好，0.05


#  === calibration
model.calibration <- glm(DAY30 ~ lp, data = gustoS, family = binomial)
# calibration intercept 接近0越好
model.calibration$coefficients[1]  

# calibration slope 接近1越好
model.calibration$coefficients[2]  


# ====7.4 Calibration slope and calibration curve

pdf("./04-Out/calibration plot1.pdf")
val.prob(p = gustoS$phat, 
         y = gustoS$DAY30, 
         logistic.cal=F,  
         # legendloc = FALSE, 
         statloc=F)
dev.off()

# 推荐使用psfmi包多重插补