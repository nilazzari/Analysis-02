# LABORATORIO 10: SURVIVAL ANALYSIS (NON PAR. (LOGRANK) E SEMIPAR. (COX))


library(survival)
leukemia
# time = survival
# status = 1 (no censored), 0 (censored)
# x = prognostic factor
attach(leukemia)

# EXPLORATORY ANALYSIS
# univariate
plot(time) # NO!!
boxplot(time)# NO !! NEI DATASET CON CENSURA NON ANALIZZARE MARGINALMENTE
# NEANCHE PER VIA GRAFICA LA VARIABILE TEMPO
table(status) # ok
table(status)/length(status)
# i dati censurati sono il 21%
table(x) # lo studio è bilanciato nei due tipi di chemioterapia

# (Y,d) è la variabile sopravvivenza avendo censura
# Y=T se d=1(il pz. ha manifestato l'evento durante l'osservazione), Y=c se d=0 se no

Surv(time, status)
# il + a fianco indica la censura

# analisi grafica della sopravvivenza
fit0 = survfit(Surv(time, status)~1)
# è la stima non parametrica di kaplan meier
summary(fit0)
# time dice a che mese si è verificato l'evento, nrisk dice il numero di soggetti a 
# rischio all'evento, nevent quanti hanno presentato l'evento, survival è la probabilità
# di sopravvivere rispetto a quella quantità di mesi(time), std.err è la stima di kaplan meier
# e infine viene data la stima intervallare
# la sopravvivenza media NON VA CALCOLATA coi dati censurati
quantile(fit0)
plot(fit0) # la curva è la sopravvivenza stimata con km e le bande di confidenza(greenwood)
plot(fit0,fun="cloglog")

# bivariate
table(status,x) # ndecessi a seconda del tipo di trattamento chemioterapico
chisq.test(status,x) #freq.attese<5
fisher.test(status,x)
fit1 = survfit(Surv(time,status)~x)
summary(fit1)
# vengono restituite le due stime di entrambe le chemioterapie
plot(fit1,col=1:2,cex=2)
legend=("topright",col=1:2,c("Manteined","Nonmanteined"),lty=c(1)) # correggiiiiiiiiiii
#ggplot2 è meglio
# il fatto che le due curve non si intereschino ci indica che il 
# logrank test è adatto per valutare i rischi proporzionali
plot(fit1,fun="cloglog")
# se in questo grafico le curve non si intersecano, allora è adatto per l'assunzione
# di rischi proporzionali (logranl, regr. di cox e modello weibull)
#test logrank: confronta due funzioni di rischio (il rischio nel 1 gruppo è uguale al 2?)
# H0: lambda_gr1 = lambda_gr2
survdiff(Surv(time,status)~x)

##############################

# modello di regressione di Cox (semiparametrico)
# lambda(t;x)=l0(t)*exp(b*x)
# lambda(t;x=1)(lambda(t;x=0)) --> = exp(b)
# l'inferenza avviene per verosimiglianza parziale
cox_fit=coxph(Surv(time,status)~x)
summary(cox_fit)
# l'expcoef dà il rapporto dei rischi, 2 volte più grande se la chemioterapia è nonmanteined
# validazione
scatter.smooth(residuals(fit))
# nessun andamento sistematico apprezzabile
qqnorm(residuals(fit))
qqline(residuals(fit))
shapiro.test(residuals(fit))
cox.zph(fit) # valuta un pvalue rispetto al regressore considerato
plot(cox.zph(fit))
# sull'asse orizzontale il tempo, sull'asse verticale i residui associati al tempo

# conclusioni: al 5% non c'è differenza tra i due trattamenti, al 10% si,
# rispetto al numero di pazienti esaminati

# se x è un fattore usiamo il test logrank e modello di cox
# se x è quantitativo NO LOGRANK e usiamo cox marginale,
# a meno che non la dicotomizziamo

#####################################################

ovarian
# fustat evento 1 morto 0 vivo, age età, resid.ds (disturbo, fattore prognostico1
# rx (tratttamento 1 o 2,fattore prognostico2), #ecog.ps (esito esame ecografico, fattore prognostico3)
# obiettivo: valutare da cosa dipende la sopravvivenza
attach(ovarian)
# diciamo a r che le ultime 3 variabili sono qualitative:
R.ds = as.factor(resid.ds)
Rx = as.factor(rx)
e.ps = as.factor(ecog.ps)

# analisi univariate
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
 
#stimatore di kaplan meier
fit0 = survfit(Surv(futime,fustat)~1) #  stiamo stimando marginalmente
summary(fit0)
# lettura dell'output come prima
quantile(fit0)
plot(fit0) # curva di sopravvivenza, diventa piatta dalla censura in poi
plot(fit0, fun="cloglog")

#analisi bivariate
# primo fattore prognostico: e.ps
fit_e = survfit(Surv(futime,fustat)~e.ps)
plot(fit_e) # non vi è proprio intersezione perche seguono lo stesso andamento, sono simili
plot(fit_e,fun="cloglog")
#logranktest:
survdiff(Surv(futime,fustat)~e.ps)
# no significativo

# secondo fattore prognostico: Rx
fit_rx = survfit(Surv(futime,fustat)~Rx)
plot(fit_rx) # non vi è proprio intersezione perche seguono lo stesso andamento, sono simili
plot(fit_rx,fun="cloglog")
#logranktest:
survdiff(Surv(futime,fustat)~Rx)
# no significativo

# terzo fattore prognostico: R.ds
fit_rds = survfit(Surv(futime,fustat)~R.ds)
plot(fit_rds) # non vi è proprio intersezione perche seguono lo stesso andamento, sono simili
plot(fit_rds,fun="cloglog")
#logranktest:
survdiff(Surv(futime,fustat)~R.ds)
# no intersezione ma..
summary(fit_rds)

# dicotomizzazione di age
age1=age<median(age)
fit_age = survfit(Surv(futime,fustat)~age1)
plot(fit_age) # non vi è proprio intersezione perche seguono lo stesso andamento, sono simili
plot(fit_age,fun="cloglog")
#logranktest:
survdiff(Surv(futime,fustat)~age1)
# no intersezione ma..
summary(fit_age)

# bivariate tra fattori prognostici (di rischio)
boxplot(age~e.ps)
wilcox.test(age~e.ps)

# modellazione (regr. di cox)
# o la facciamo onward o backward
# io preferisco backward
# trv parziale (annidati)
# età quantitativa possiamo lasciarla cosi nel modello di cox
fit1 = coxph(Surv(futime,fustat)~age+R.ds+Rx+e.ps)
summary(fit1)
# leviamo quella col massimo pvalue --> e.ps
fit2 = coxph(Surv(futime,fustat)~age+R.ds+Rx)
# leviamo quella col massimo pvalue --> R.ds
fit3 = coxph(Surv(futime,fustat)~age+Rx)
# leviamo quella col massimo pvalue e che non è significativo --> Rx
fit4 = coxph(Surv(futime,fustat)~age)
summary(fit4)
# preferibile rispetto al modello con la sola intercetta

# proviamo con l'età dicotomizzata
fit5 = coxph(Surv(futime,fustat)~age1)
summary(fit5)
# GRANDE PERDITA DI INFORMAZIONE

# validiamo fit4
qqnorm(residuals(fit4))
qqline(residuals(fit4))
shapiro.test(residuals(fit4))
#test per rischi proporzionali
cox.zph(fit4)
# ipotesi accettata
plot(cox.zph(fit4))
#andamento quasi lineare, ok rischi proporzionali
