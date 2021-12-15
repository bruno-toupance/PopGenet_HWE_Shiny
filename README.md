# PopGenet_HWE_Shiny
Hardy-Weinberg Equilibrium Test + Shiny


## Required libraries
- adegenet
- pegas
- HardyWeinberg
- shiny

To install/update the packages, run the following commands into the R console:
```{r }
install.packages("adegenet")
install.packages("pegas")
install.packages("HardyWeinberg")
install.packages("shiny")
```


## Run from R
```{r }
library(shiny)
runGitHub(repo = "PopGenet-HWE_Shiny", username = "bruno-toupance", ref = "main")
```
