dt = read.table("score.csv",sep = ",", header = T)
dt = dt[,2:3]
library(stringr)
dt_new = as.data.frame(str_split_fixed(dt$poches, ";", 2))
df = cbind(dt_new,dt[,2])
colnames(df)[3] = 'scores'

