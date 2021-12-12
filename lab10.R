library(survival)
library(ggplot2)
leukemia # status = 1 (no censored), 0 (censored) wh. x = prognostic factor
attach(leukemia)

# Exploratory analysis
table(status) # ok
table(status)/length(status)
table(x)
Surv(time, status)
fit0 = survfit(Surv(time, status)~1)
summary(fit0)
quantile(fit0)
plot(fit0)
plot(fit0,fun="cloglog")
table(status,x) # number of deceased based on chemoterapic treatment type
chisq.test(status,x) # Warning: f.q.<5
fisher.test(status,x)
fit1 = survfit(Surv(time,status)~x)
summary(fit1)
plot(fit1,col=1:2,cex=2)
legend=("topright",col=1:2,c("Manteined","Nonmanteined"),lty=c(1,1))
plot(fit1,fun="cloglog") # H0: lambda_gr1 = lambda_gr2
survdiff(Surv(time,status)~x)

# Cox regression
cox_fit=coxph(Surv(time,status)~x)
summary(cox_fit)
scatter.smooth(residuals(fit))
qqnorm(residuals(fit))
qqline(residuals(fit))
shapiro.test(residuals(fit))
cox.zph(fit)
plot(cox.zph(fit))

# Second dataset
ovarian
attach(ovarian)
R.ds = as.factor(resid.ds) # dicotomization
Rx = as.factor(rx)
e.ps = as.factor(ecog.ps)

# Exploratory analysis
table(e.ps)
table(e.ps)/26
table(Rx)
table(R.ds)
table(R.ds)/26
boxplot(age)
summary(age)
sd(age)
table(fustat)
table(fustat)/26
 
# Kaplan-Meier estimations
fit0 = survfit(Surv(futime,fustat)~1)
summary(fit0)
quantile(fit0)
plot(fit0)
plot(fit0, fun="cloglog")

fit_e = survfit(Surv(futime,fustat)~e.ps)
plot(fit_e) 
plot(fit_e,fun="cloglog")

# Log-Rank test on prognostic factors
survdiff(Surv(futime,fustat)~e.ps)
fit_rx = survfit(Surv(futime,fustat)~Rx)
plot(fit_rx)
plot(fit_rx,fun="cloglog")
survdiff(Surv(futime,fustat)~Rx)
fit_rds = survfit(Surv(futime,fustat)~R.ds)
plot(fit_rds)
plot(fit_rds,fun="cloglog")
survdiff(Surv(futime,fustat)~R.ds)
summary(fit_rds)
age1=age<median(age)
fit_age = survfit(Surv(futime,fustat)~age1)
plot(fit_age)
plot(fit_age,fun="cloglog")
survdiff(Surv(futime,fustat)~age1)
summary(fit_age)
boxplot(age~e.ps)
wilcox.test(age~e.ps)

# Cox regression: backward modeling
fit1 = coxph(Surv(futime,fustat)~age+R.ds+Rx+e.ps)
summary(fit1)
fit2 = coxph(Surv(futime,fustat)~age+R.ds+Rx)
fit3 = coxph(Surv(futime,fustat)~age+Rx)
fit4 = coxph(Surv(futime,fustat)~age)
summary(fit4)
fit5 = coxph(Surv(futime,fustat)~age1)
summary(fit5)

# Validation after selecting model no.4
qqnorm(residuals(fit4))
qqline(residuals(fit4))
shapiro.test(residuals(fit4))
cox.zph(fit4)
plot(cox.zph(fit4))
