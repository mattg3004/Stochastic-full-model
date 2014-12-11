
initial.prop.susceptible  =  0
disease.state  =  initial.disease.state(demographic.ages  ,   initial.prop.susceptible  ,  num.comps,  list.of.ages, mat.immunity.loss)
t  =  0
for (i in 1 : num.steps){
  disease.state  =  demographics(disease.state, birth.rate, death.rate, time.step)
  disease.state  =  migrants.l.per.year(disease.state, migrant.indices, time.step, l)
  beta_0  =  calibrate.beta(mixing.matrix, disease.state, max.age, time.step, infectious.period, R_0, list.of.ages, num.comps)
  beta  =  beta * (1 + beta_1 * cos(2 * pi * t /365))
  
  updated.state  =  draw.next.step.all(disease.state, mixing.matrix, infectious.indices, 
                                      time.step , infectious.period, beta,
                                      demographic.ages, num.comps)
  
  t  =  t  +  time.step
}