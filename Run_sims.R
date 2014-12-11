
initial.prop.susceptible  =  0
disease.state  =  initial.disease.state(demographic.ages  ,   initial.prop.susceptible  ,  num.comps,  list.of.ages, mat.immunity.loss)

for (i in 1 : num.steps){
  disease.state  =  demographics(disease.state, birth.rate, death.rate, time.step)
  beta  =  calibrate.beta(mixing.matrix, disease.state, max.age, time.step, infectious.period, R_0, list.of.ages, num.comps)
}