# PopGenet_HWE_Shiny
Hardy-Weinberg Equilibrium Test + Shiny


## Required libraries
- adegenet
- pegas
- HardyWeinberg
- shiny

To install/update the packages, run the following commands into the R console:
```{r }
install.packages("adegenet", dependencies = TRUE)
install.packages("pegas", dependencies = TRUE)
install.packages("HardyWeinberg", dependencies = TRUE)
install.packages("shiny", dependencies = TRUE)
```


## Run from R
```{r }
library("shiny")
runGitHub(repo = "PopGenet_HWE_Shiny", username = "bruno-toupance", ref = "main")
```
