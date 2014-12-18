library(RColorBrewer)
library(akima)
library(fields)
source("Deterministic_functions.R")
source("Deterministic_step_for_measles.R")
source("Initial_conditions.R")
tanzania.age.groups = read.csv("Tanzania_pop.csv")
tanzania.ages  =  matrix(0, 98, 2)
tanzania.ages[, 1]  =  0:97
for (i in 1 : length(tanzania.age.groups[,2])){
  tanzania.ages[((i-1) * 5) + 1, 2]  =  tanzania.age.groups[i, 2] * 1000
}

for ( i in 1 : length(tanzania.ages[, 2])){
  if ( tanzania.ages[i, 2]  ==  0){
    l = which(tanzania.ages[1 : i, 2]  >  0)
    a = tanzania.ages[l[length(l)], 2]
    k  = which(tanzania.ages[(i + 1) : length(tanzania.ages[, 2]) , 2]  >  0)[1]
    b  =  tanzania.ages[(i + k), 2]
    tanzania.ages[i, 2]  =  a + (b - a) * (5 - k) / 5
  }
}
tanzania.ages[97, 2] = 500
tanzania.ages[ , 2] = tanzania.ages[ , 2] * (1000 * sum(tanzania.age.groups[, 2])) / sum(tanzania.ages[ , 2] )

birth.rate = 42
death.rate = 13.5



vacs  =  c(0.89, 0.9, 0.94, 0.91, 0.89, 0.9, 0.88, 0.91, 0.92, 0.95, 0.97, 0.99)
sia   =  c(4, 7, 10)
sia2  =  c(2, 5, 8, 11)
vacc.succ  =  0.9
infs  =  matrix(0, length(vacs), 1)
infs2  =  matrix(0, length(vacs), 1)

v = 0.89*0.9
count =  1
weighted.mixing  =  0
for ( i in 1 : 97){
  disease.state[seq((1+ (i-1) * num.comps), i * num.comps)]    =  disease.state.i[seq((1+ (i-1) * num.comps), i * num.comps)] * tanzania.ages[i, 2] /  sum(disease.state.i[seq((1+ (i-1) * num.comps), i * num.comps)])
}
disease.state[all.exposed.and.infectious.2]  =  0

disease.state[sus]  =  disease.state[sus] * 0.042
disease.state[recovered.indices]  =  disease.state[recovered.indices] + (disease.state[sus] * (0.958))
disease.state[infectious.indices[13]]  =  1
disease.state1   =    Calculate.infected.by.age.given.sus.profile(disease.state, beta_0, beta_1, time.step, birth.rate/(1000*365) ,death.rate/(1000*365), migrant.indices, av.migrants.per.age.per.day, all.infectious,
                                                                  mixing.matrix, all.exposed.and.infectious, disease.time, all.exposed.and.infectious.aged, next.step.exposed.and.infectious,
                                                                  recovered.indices, maternal.indices, mat.immunity.loss,  susceptible.indices, sus, sus2, v, infected.by.age, 0, R_0, infectious.period, weighted.mixing,
                                                                  365, max.age.of.transmission, supp.vac.period, supp.vac, adjust.infecteds, min.age.supp.vacc,
                                                                  max.age.sup.vacc, list.of.ages, include.vac.for.infection, years.per.infection.vacc)

infs.by.age = unlist(disease.state1[6])
sum(infs.by.age)
infs[1] = sum(infs.by.age)
disease.state = unlist(disease.state1[1])
#disease.state[recovered.indices]   =   disease.state[recovered.indices]  -  (disease.state[sus] * .12)
#disease.state[sus]  =  disease.state[sus] - (disease.state[sus] * .12)

#disease.state[13]   =  number.initial.infecteds

for ( i in 2 : length(vacs) ){
  v  =  vacs[i] * vacc.succ
  b  =  41.8
  disease.state[all.exposed.and.infectious.2]  =  0
  disease.state[13]  =   number.initial.infecteds  +  disease.state[13]
  disease.state1     =   Calculate.infected.by.age.given.sus.profile(disease.state, beta_0, beta_1, time.step, birth.rate ,death.rate, migrant.indices, av.migrants.per.age.per.day, all.infectious,
                                                                     mixing.matrix, all.exposed.and.infectious, disease.time, all.exposed.and.infectious.aged, next.step.exposed.and.infectious,
                                                                     recovered.indices, maternal.indices, mat.immunity.loss,  susceptible.indices, sus, sus2, v, infected.by.age, 0, R_0, infectious.period, weighted.mixing,
                                                                     365 , max.age.of.transmission, supp.vac.period, supp.vac, adjust.infecteds, min.age.supp.vacc,
                                                                     max.age.sup.vacc, list.of.ages, include.vac.for.infection, years.per.infection.vacc)

  infs[i]   =  sum( unlist(disease.state1[4]) )
  disease.state  =  unlist(disease.state1[1])

  if (sia[count] == i){
    if (count == 2){
      disease.state  =  reduce.susceptibles(6/12, 10, disease.state, 0.7 , num.comps, susceptible.indices, recovered.indices, list.of.ages)
    }
    else{
      disease.state  =  reduce.susceptibles(9/12, 5, disease.state, 0.7 , num.comps, susceptible.indices, recovered.indices, list.of.ages)
    }
    count  =  count + 1
    if( count > length(sia)){
      count  =  length(sia)
    }
  }
}

infs
count = 1

for ( i in 1 : length(vacs) ){
  v  =  vacs[i] * vacc.succ
  b  =  41.8
  disease.state[13]  =   number.initial.infecteds  +  disease.state[13]
  disease.state1     =   Calculate.infected.by.age.given.sus.profile(disease.state, beta_0, beta_1, time.step, birth.rate ,death.rate, migrant.indices, av.migrants.per.age.per.day, all.infectious,
                                                                     mixing.matrix, all.exposed.and.infectious, disease.time, all.exposed.and.infectious.aged, next.step.exposed.and.infectious,
                                                                     recovered.indices, maternal.indices, mat.immunity.loss,  susceptible.indices, sus, sus2, v, infected.by.age, 0, R_0, infectious.period, weighted.mixing,
                                                                     365 , max.age.of.transmission, supp.vac.period, supp.vac, adjust.infecteds, min.age.supp.vacc,
                                                                     max.age.sup.vacc, list.of.ages, include.vac.for.infection, years.per.infection.vacc)

  infs2[i]  =  sum( unlist(disease.state1[4]) )
  disease.state  =  unlist(disease.state1[1])
  if (sia2[count] == i){
    #if (count == 2){
    #  disease.state= reduce.susceptibles(6/12, 10, disease.state, 0.7 , num.comps, susceptible.indices, recovered.indices, list.of.ages)
    #}
    #else{
      disease.state  =  reduce.susceptibles(9/12, 5, disease.state, 0.7 , num.comps, susceptible.indices, recovered.indices, list.of.ages)
    #}
    count  =  count + 1
    if( count > length(sia2)){
      count  =  length(sia2)
    }
  }
}



#dd2    =   disease.state
#dd     =   dd2
#count  =   1
#for ( i in 1 : length(vacs) ){
#  v  =  .99 * vacc.succ
#  dd  =  Age.population(dd, b/(1000*365), 15/(1000*365), all.exposed.and.infectious, disease.time,
#                        all.exposed.and.infectious.aged, time.step, next.step.exposed.and.infectious,
#                        recovered.indices, maternal.indices, mat.immunity.loss,  susceptible.indices,
#                        sus, sus2, v, 365, supp.vac.period, supp.vacc, supp.vac.prop, min.age.sup.vacc, max.age.sup.vacc)
#  if (sia[count] == i){
#    if (count == 2){
#      dd  =  reduce.susceptibles(6/12, 10, dd, 0.7 , num.comps, susceptible.indices, recovered.indices, list.of.ages)
#    }
#    else{
#      dd  =  reduce.susceptibles(9/12, 5, dd, 0.7 , num.comps, susceptible.indices, recovered.indices, list.of.ages)
#    }
#    count = count + 1
#    if( count > length(sia)){
#      count = length(sia)
#    }
#  }
#}
#dd[13]   =   number.initial.infecteds  +  disease.state[13]
#disease.state1   =   Calculate.infected.by.age.given.sus.profile(dd, beta_0, beta_1, time.step, birth.rate ,death.rate, migrant.indices, av.migrants.per.age.per.day, all.infectious,
#                                                                 mixing.matrix, all.exposed.and.infectious, disease.time, all.exposed.and.infectious.aged, next.step.exposed.and.infectious,
#                                                                 recovered.indices, maternal.indices, mat.immunity.loss,  susceptible.indices, sus, sus2, v, infected.by.age, 0, R_0, infectious.period, weighted.mixing,
#                                                                 365 , max.age.of.transmission, supp.vac.period, supp.vac, adjust.infecteds, min.age.supp.vacc,
#                                                                 max.age.sup.vacc, list.of.ages, include.vac.for.infection, years.per.infection.vacc)

#sum( unlist(disease.state1[4]) )
