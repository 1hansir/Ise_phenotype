pak.list = c("data.table","MAP",'dplyr', 'Rcpp',
             'flexmix',"Matrix","pROC",'bit','ff',
             'survival','reticulate','survivalmodels',
             'ggplot2','ggpubr'
             )

for (pak in pak.list)
{
  yo = require(pak, character.only = T)
  if(!yo)
  {
    install.packages(pak, type="binary",repos = "http://cran.r-project.org")
    require(pak, character.only = T)
  }
}