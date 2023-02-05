pred.STM = function(fit, test.Z, tgrid,
                    glink = expit)
{
  glink(outer(drop(test.Z%*% fit$beta),
              fit$hfun(tgrid), '+'))
}