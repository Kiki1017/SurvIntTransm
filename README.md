# SurvIntTransm
Influence of intensity of surveillance for preventing sustained transmission of covid19 in new locations

## 1. Running the app
You will need to have shiny package installed:
```r
library(shiny)
runGitHub("SurvIntTransm", "rretkute")
```
Note that this app depends on different R packages, including:
* [ggplot2](https://ggplot2.tidyverse.org/)
* [grid](https://www.rdocumentation.org/packages/grid/versions/3.6.2)
* [gridExtra](https://cran.r-project.org/web/packages/gridExtra/index.html) 
* [rstan](https://github.com/stan-dev/rstan/wiki/RStan-Getting-Started)
* [rootSolve](https://cran.r-project.org/web/packages/rootSolve/index.html)
* [loo](https://cran.r-project.org/web/packages/loo/index.html)

## 2. Methods
###  2.1.Model 
The app fits to the observed days from symptom onset to hispitalisatio either exponental [1], or gamma distrubution [2].

###  2.2.Leave-one-out cross-validation  
[3]


### References
[1]  Thompson, R.N. (2020) J. Clin. Med.  9(2), 498. https://doi.org/10.3390/jcm9020498 

[2] Thompson, R.N., Jalava, K., Obolski, U. (2019) Sustained transmission of Ebola in new locations: more likely than previously thought. Lancet Infect Dis. 19: 1058–27859. https://doi.org/10.1016/S1473-3099(19)30483-9

[3] Vehtari, A., Gelman, A., and Gabry, J. (2017). Practical Bayesian model evaluation using leave-one-out cross-validation and WAIC. Statistics and Computing. 27(5), 1413–1432. https://doi.org/10.1007/s11222-016-9696-4.

