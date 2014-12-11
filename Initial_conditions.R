demographic.ages              =       read.csv("Ages.csv")
demographic.ages              =       round(read.csv("Pop_Zimbabwe.csv"))
contacts                      <-      read.csv("Contacts.csv")
uk.contacts                   =       read.csv("uk_contacts.csv")
conts                         =       read.csv("contact_rate.csv")
initial.pop                   =       sum(demographic.ages[, 2])


mixing.matrix                 <-      full.mixing.matrix(contacts, demographic.ages)      # average number of people met of each age grooup, stratified by age
mixing.matrix               =       matrix(0, length(demographic.ages[,1]) + 11, length(demographic.ages[,1])+11)
pop.by.age  =  matrix(0, 109,1)
for ( i in 1:12){
  mixing.matrix  [i, ]  =  demographic.ages[1, 2] / 12       # this produces a mixing matrix which is equivalent to uniform mixing
  pop.by.age[i]         =  demographic.ages[1, 2] / 12
}
for ( i in 13:(length(demographic.ages[,1]) + 11)){
  mixing.matrix  [i, ]  =  demographic.ages[i - 11, 2]       # this produces a mixing matrix which is equivalent to uniform mixing
  pop.by.age[i]         =  demographic.ages[i - 11, 2]
}

list.of.ages                  =       matrix(0, length(demographic.ages[ , 1]) + 11, 1)
list.of.ages[1 : 12]          =       (seq(0,11) / 12)
list.of.ages[13 : length(list.of.ages)]  =  seq(1 , max(demographic.ages[, 1]))

num.comps = 5

maternal.indices              =       seq(1, (length(list.of.ages) - 2) * num.comps, num.comps)
susceptible.indices           =       seq(2  , length(list.of.ages) * num.comps ,  num.comps )
exposed.indices               =       seq(3, length(list.of.ages) * num.comps, num.comps)
infectious.indices            =       seq(4, length(list.of.ages) * num.comps, num.comps)
recovered.indices             =       seq(num.comps, length(list.of.ages) * num.comps, num.comps)

time.step = 1

R_0  = 15




incubation.period             =       10         # length of incubation period on average
mu                            =       min(1, time.step/incubation.period)   # probability of moving from exposed to infectious class during a timestep
infectious.period             =       8          # number of days spent in the infected class on average
rho                           =       min(1,time.step/infectious.period)    # probability of losing infectiousness. not necessarily recovered from the infection, but no longer infectious.

mat.immunity.loss             =      matrix(1, length(maternal.indices), 1)
mat.immunity.loss[1:9]             =      c(0, 0, 0.1, 0.2, 0.3, 0.5, 0.8, 0.9, 1)
mat.immunity.loss[1:9]             =      (1/(1+exp(-(seq(0,8,1) - 5))))

for (i in 1 : 8){
  mat.immunity.loss[i]             =    1  -  exp(-i * 0.03)
}

