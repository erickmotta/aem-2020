# Load packages
library(ggfortify)
library(survival)

# Figure 3E, Serratia challenge at site 1, 2018
tab = read.table("tables/Serratia_challenge_site1_2018.txt",h=T)

start = tab$Start
stop = tab$Stop 
status = tab$Death

S = Surv(start,stop,status)

summary(S)

model = coxph(S ~ Treatment, data = tab)

cups = coxph(S ~ Cup , data = tab)

summary(model)

fit <- survfit(S ~ Treatment, data = tab)
autoplot(fit)


# Figure 4C, Serratia challenge at site 1, 2019
tab = read.table("tables/Serratia_challenge_site1_2019.txt",h=T)

start = tab$Start
stop = tab$Stop 
status = tab$Death

S = Surv(start,stop,status)

summary(S)

model = coxph(S ~ Treatment, data = tab)

cups = coxph(S ~ Cup , data = tab)

summary(model)

fit <- survfit(S ~ Treatment, data = tab)
autoplot(fit)


# Figure 5C, Serratia challenge at site 2, 2019
tab = read.table("tables/Serratia_challenge_site2_2019.txt",h=T)

start = tab$Start
stop = tab$Stop 
status = tab$Death

S = Surv(start,stop,status)

summary(S)

model = coxph(S ~ Treatment, data = tab)

cups = coxph(S ~ Cup , data = tab)

summary(model)

fit <- survfit(S ~ Treatment, data = tab)
autoplot(fit)

# Figure 6B, Survival curve for sprayed bees
tab = read.table("tables/topical_exp_survival_curve.txt",h=T)

start = tab$Start
stop = tab$Stop 
status = tab$Death

S = Surv(start,stop,status)

summary(S)

model = coxph(S ~ Treatment, data = tab)

cups = coxph(S ~ Cup , data = tab)

summary(model)

fit <- survfit(S ~ Treatment, data = tab)
autoplot(fit)
