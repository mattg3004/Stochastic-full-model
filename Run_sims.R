source("Initial_conditions.R")
num.steps  =  25000
birth.rate =  40
death.rate = 15
vacc.prop = 1
l = 1
vacc.prop = 0.85 * 0.9
new.infections  =  matrix(0, num.steps, 1)
all.infections  =  matrix(0, num.steps, 1)
prop.sus        =  matrix(0, num.steps, 1)
initial.prop.susceptible  =  0
disease.state  =  initial.disease.state(demographic.ages  ,   initial.prop.susceptible  ,  num.comps,  list.of.ages, mat.immunity.loss)
t  =  0
quartz()
for (i in 1 : num.steps){
  disease.state  =  demographics(disease.state, birth.rate, death.rate, time.step)
  disease.state  =  migrants.l.per.year(disease.state, migrant.indices, time.step, l)
  beta_0  =  calibrate.beta(mixing.matrix, disease.state, max.age, time.step, infectious.period, R_0, list.of.ages, num.comps)
  beta  =  beta_0 * (1 + beta_1 * cos(2 * pi * t /365))
  
  A  =  draw.next.step.all(disease.state, mixing.matrix, infectious.indices, 
                                      time.step , infectious.period, beta,
                                      demographic.ages, num.comps, maternal.indices, 
                                      mat.immunity.loss, vacc.prop)
  disease.state      =   unlist(A[1])
  new.infections[i]  =   unlist(A[2])
  all.infections[i]  =   sum(disease.state[infectious.indices])
  prop.sus[i]        =   sum(disease.state[susceptible.indices]) / sum(disease.state)
  if (i %% 500 == 0){
    par(mfrow = c(3,1))
    plot(1:i, prop.sus[1:i], type = "l", col = 'seagreen3')
    plot(1:i, new.infections[1:i], type = "l", col = 'dodgerblue3')
    plot(1:i, all.infections[1:i], type = "l", col = 'red')
  }
  t  =  t  +  time.step
}