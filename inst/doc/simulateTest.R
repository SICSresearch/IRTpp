## ------------------------------------------------------------------------
library(IRTpp)
test <- simulateTest(model="2PL",items=10,individuals=1000)

## ------------------------------------------------------------------------
test$itempars

## ------------------------------------------------------------------------
responses <-test$test
responses[1:10,]

## ------------------------------------------------------------------------
t = simulateTest(model="1PL",items=4,individuals=5)
length(t$test)

## ------------------------------------------------------------------------
t3 = simulateTest(model="3PL",items=500,individuals=10);
summary(t3$itempars$c)

## ------------------------------------------------------------------------
bd = list(c_lower=0.2)
t3 = simulateTest(model="3PL",items=500,individuals=10,boundaries=bd);
summary(t3$itempars$c)

## ------------------------------------------------------------------------
t3 = simulateTest(model="3PL",items=10,individuals=100,threshold=0.2);

## ------------------------------------------------------------------------
response <- t3$test
summary(rowSums(response))

