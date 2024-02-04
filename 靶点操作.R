library(tidyverse)
`(8,11,14.Docosatrienoic.acid,.methyl.ester)`  -> x1
x1 %>% filter(x1$Probability.!=0) %>%select("Common.name","Probability.") -> x1
`[(2R).2.[[[(2R).2.(benzoylamino).3.phenylpropanoyl]amino]methyl].3.phenylpropyl].acetate` ->x2
x2 %>% filter(x2$Probability.!=0) %>%select("Common.name","Probability.") -> x2
epicampesterol -> x3
x3 %>% filter(x3$Probability.!=0) %>%select("Common.name","Probability.") -> x3
bsitosterol -> x4
x4 %>% filter(x4$Probability.!=0) %>%select("Common.name","Probability.") -> x4
sitosterol -> x5
x5 %>% filter(x5$Probability.!=0) %>%select("Common.name","Probability.") -> x5
Stigmasterol -> x6
x6 %>% filter(x6$Probability.!=0) %>%select("Common.name","Probability.") -> x6
CLR -> x7
x7 %>% filter(x7$Probability.!=0) %>%select("Common.name","Probability.") -> x7
final_targets <- list(x1,x2,x3,x4,x5,x6,x7)
do.call(rbind,final_targets) -> targets
unique(targets) -> f_targets
unique(genes) -> targets_qu
xlsx::write.xlsx(targets,"targets_p.xlsx")
xlsx::write.xlsx(genes,"targets.xlsx")
xlsx::write.xlsx(targets_qu,"targets_quchong.xlsx")
