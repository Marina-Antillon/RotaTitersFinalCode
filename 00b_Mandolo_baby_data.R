library(readxl)
library(dplyr)
library(Epi)

Mandolo <- read_excel("./Mandolo_data_new.xlsx")

Mandolo$IgGtitre = 10^Mandolo$IgG6wk_log10
Mandolo$Seropos = Mandolo$IgAStatus=="Seropositives"

OR_sc=glm(Seropos ~ log(IgGtitre, base=2), family=binomial(link = "logit"), data=Mandolo)
tmp = summary(OR_sc)

exp(tmp$coefficients[2,1])
exp(tmp$coefficients[2,1]-1.96*tmp$coefficients[2,2])
exp(tmp$coefficients[2,1]+1.96*tmp$coefficients[2,2])

summary(Mandolo$IgGtitre[Mandolo$Seropos]) # Paper: 5745 (IQR=2323-13350)
summary(Mandolo$IgGtitre[!Mandolo$Seropos]) # Paper: 9689 (IQR=5512-23616)
