num.sims = 1000
total = 0
disease.state[infectious.indices[10]]  =  1
pop.by.age  =  population.by.age(disease.state, num.comps)
for ( i in 1:length(pop.by.age)){
  mixing.matrix  [i, ]  =  pop.by.age[i]       # this produces a mixing matrix which is equivalent to uniform mixing
}
beta_0  =  calibrate.beta(mixing.matrix, disease.state, 80, time.step, infectious.period, R_0, list.of.ages, num.comps)
for (i in 1 : (infectious.period * num.sims)){
  A  =  draw.next.step.all(disease.state, mixing.matrix, infectious.indices, 
                           time.step , infectious.period, beta_0,
                           demographic.ages, num.comps, maternal.indices, 
                           mat.immunity.loss, vacc.prop)
  
  total  =  total + unlist(A[2])
}
total / (num.sims * (sum(disease.state[susceptible.indices]) / sum(disease.state)*15))
