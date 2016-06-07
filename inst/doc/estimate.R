## ------------------------------------------------------------------------
library(IRTpp)
tst = simulateTest(model="3PL")

## Calibrate this test with a 2PL model.
est = irtpp(tst$test,"1PL")

## ------------------------------------------------------------------------
names(est)
est$z
est$LL

## ------------------------------------------------------------------------
test.plot(est$z)

## ------------------------------------------------------------------------
## Calibrating the same test under a 3PL model and displaying the AIC and BIC statistics
est = irtpp(tst$test,"3PL", loglikflag=T)

## ------------------------------------------------------------------------
test.plot(est$z)

## ------------------------------------------------------------------------
test.plot(est$z,2)

## ------------------------------------------------------------------------
zz = parameter.matrix(est$z,byrow = F)
th = individual.traits(model="3PL", itempars = zz,method = "EAP",dataset = tst$test, probability_matrix = est$prob_mat)

##The latent traits.
hist(th[,ncol(th)],breaks=40)

