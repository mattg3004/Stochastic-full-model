source("Initial_conditions.R")
source("Functions.R")
num.steps  =  25000
birth.rate =  40
death.rate = 15
vacc.prop = 1
l = 1
vacc.prop = 0.85 * 0.9
#vacc.prop = 0
new.infections  =  matrix(0, num.steps, 1)
all.infections  =  matrix(0, num.steps, 1)
prop.sus        =  matrix(0, num.steps, 1)
initial.prop.susceptible  =  0
disease.state  =  initial.disease.state(demographic.ages  ,   initial.prop.susceptible  ,  num.comps,  list.of.ages, mat.immunity.loss)
t  =  0
quartz()

list[d2,  new.infs, all.infss, prop.sus] = Run.simulations (num.steps, disease.state, mixing.matrix, infectious.indices,
                                                            birth.rate, death.rate,
                                                            time.step , infectious.period, beta_1,
                                                            demographic.ages, num.comps, maternal.indices, 
                                                            mat.immunity.loss, vacc.prop, l, do.plots = 1, initial.time = 280)

disease.state.2 = d2
disease.state.2[infectious.indices[10] ] = 1
list[d2,  new.infs, all.infss, prop.sus]  = Run.simulations (365, disease.state.2, mixing.matrix, infectious.indices,
                 birth.rate, death.rate,
                 time.step , infectious.period, beta_1,
                 demographic.ages, num.comps, maternal.indices, 
                 mat.immunity.loss, vacc.prop, l, do.plots = 1, initial.time = 280)
