Malawi.vacc = subset(vacc.rates, vacc.rates$Country == 'Malawi')
Malawi.cases = subset(Af_cases, Af_cases$Cname == 'Malawi')
Malawi.br  =  subset(Birth.rates, Birth.rates$Country == 'Malawi')



Malawi.vacc.from.1998 =  as.numeric(Malawi.vacc[20:32])
Malawi.BR.from.1998   =  as.numeric(Malawi.br[39:51])
Malawi.SIA.number     =  c(156154, 1583664, 1851176, 2120557)
SIA.year              =  c(2, 5, 8, 11)

cl <- makeCluster(4, outfile="")
registerDoParallel(cl)
num.repeats  =  4
a = numerous.sims.diff.br.vacc(num.repeats,  disease.state.3, 
                           initial.start.time = 280, multi.vacc = Malawi.vacc.from.1998, multi.br = Malawi.BR.from.1998, input.infections.all.times = 0,
                           init.pop.size = 10700180, sia.years = SIA.year, sia.numbers = Malawi.SIA.number, vacc.success = 0.9)

malawi.sim = matrix(0, num.repeats, length(Malawi.BR.from.1998))
for(i in 1 : num.repeats){
  malawi.sim[i, ] = unlist(a[i])
}


disease.state.input = disease.state.3
vacc.success = 0.9
for ( i in 1: length(Malawi.BR.from.1998)){
  vacc.prop   =   Malawi.vacc.from.1998[i]*vacc.success / 100
  birth.rate  =   Malawi.BR.from.1998[i]
  disease.state.input   =   age.population.with.demographics(disease.state.input, mixing.matrix, infectious.indices, 
                                                             time.step , infectious.period, 
                                                             demographic.ages, num.comps, maternal.indices, 
                                                             mat.immunity.loss, vacc.prop,  365)  
}
