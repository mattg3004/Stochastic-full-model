demographic.ages              =       read.csv("Ages.csv")
demographic.ages              =       round(read.csv("Pop_Zimbabwe.csv"))
contacts                      <-      read.csv("Contacts.csv")
uk.contacts                   =       read.csv("uk_contacts.csv")
conts                         =       read.csv("contact_rate.csv")
vacc.rates                    =       read.csv("African_vaccination_rates.csv")
cases.by.year                 =       read.csv("Measles_cases_by_year.csv")
Af_cases                      =       subset(cases.by.year, cases.by.year$WHO_REGION == "AFR")
Birth.rates                   =       read.csv("Birth_rates.csv")
initial.pop                   =       sum(demographic.ages[, 2])

All.African.Countries = Af_cases$Cname

########################
# Construct the mixing matrix which will be used for simulation. 
# Here we produce a mixing matrix equivalent to uniform mixing
########################
# mixing.matrix                 <-      full.mixing.matrix(contacts, demographic.ages)      # average number of people met of each age grooup, stratified by age
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

initial.prop.susceptible  =  0.05
beta_1   =   0.6
########################
# Set up the disease state matrix, with MSEIR compartments
########################
num.comps = 5

maternal.indices              =       seq(1, (length(list.of.ages) - 2) * num.comps, num.comps)
susceptible.indices           =       seq(2  , length(list.of.ages) * num.comps ,  num.comps )
exposed.indices               =       seq(3, length(list.of.ages) * num.comps, num.comps)
infectious.indices            =       seq(4, length(list.of.ages) * num.comps, num.comps)
recovered.indices             =       seq(num.comps, length(list.of.ages) * num.comps, num.comps)


oldest.migrant       =      15 
migrant.indices      =      seq(4,  num.comps * (oldest.migrant + 13), num.comps)


time.step = 1

R_0  = 15




incubation.period             =       10         # length of incubation period on average
mu                            =       min(1, time.step/incubation.period)   # probability of moving from exposed to infectious class during a timestep
infectious.period             =       8          # number of days spent in the infected class on average
rho                           =       min(1,time.step/infectious.period)    # probability of losing infectiousness. not necessarily recovered from the infection, but no longer infectious.

mat.immunity.loss             =      matrix(1, length(maternal.indices), 1)
#mat.immunity.loss[1:9]             =      c(0, 0, 0.1, 0.2, 0.3, 0.5, 0.8, 0.9, 1)
mat.immunity.loss[1:8]             =      (1/(1+exp(-(seq(0,7,1) - 5))))


min.age.sia  =  9/12
max.age.sia  =  5
sia.proportion = 0.6
do.plots = 0
sia.period  =  3
av.migrants.per.year = 1
#for (i in 1 : 8){
#  mat.immunity.loss[i]             =    1  -  exp(-i * 0.03)
#}


