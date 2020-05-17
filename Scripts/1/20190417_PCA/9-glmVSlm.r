
# simulate data
x1 = rbinom(1000,1,.5)
x2 = rpois(1000,10)
simul <- data.frame(
    u = rGA(1000,1+x1*5,1),
    x1 = x1,
    x2 = x2)
    
# the fits of a gamma (GA) and gaussian (NO) regression 
# should lead to the same estimates for x1 and x2 
# (when using the same link function) 
summary(fit.GA <- gamlss(u ~ x1+x2, data = simul, 
            family = GA(mu.link="identity")))
summary(fit.NO <- gamlss(log(u) ~ x1+x2, data = simul, 
            family = NO(mu.link="identity")))

# the fits of a gamma and gaussian regression 
# should lead to the same estimates (when using the same link function)             
library(gamlss)    
Rsq(gamlss(u ~ x1+x2, data = simul, 
            family = GA(mu.link="identity")))
Rsq(gamlss(u ~ x1+x2, data = simul, 
            family = GA(mu.link="log")))
Rsq(gamlss(log(u) ~ x1+x2, data = simul,
            family = NO))
