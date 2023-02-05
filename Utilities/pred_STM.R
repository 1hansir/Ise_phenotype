pred.STM = function(fit, test.Z, tgrid)
{
  expit(outer(drop(test.Z%*% fit$beta),
              fit$hfun(tgrid), '+'))
}