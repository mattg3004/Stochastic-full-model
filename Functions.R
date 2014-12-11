#############################
# This function creates the initial disease state that the population is in. 
# Can specify the proportion of people who are susceptible to begin with using initial.prop.susceptible variable
#############################
initial.disease.state <- function(demographic.ages  ,   initial.prop.susceptible  ,  num.comps,  list.of.ages, mat.immunity.loss){
  disease.state                 =      matrix(0,num.comps*length(list.of.ages),1)
  
  for (i in 1:12){
    disease.state[(((i-1)*num.comps)+1):((i)*num.comps)]      =    round(c( (demographic.ages[1,2]*(initial.prop.susceptible / 12) * (1 - mat.immunity.loss[i])) , 
                                                                            (demographic.ages[1,2]*(initial.prop.susceptible / 12) * ( mat.immunity.loss[i])), 0, 0, 
                                                                            (demographic.ages[1,2]*(1-initial.prop.susceptible)/12)))
  }
  for (i in 13:length(list.of.ages)){
    disease.state[(((i-1)*num.comps)+1):((i)*num.comps)]      =    round(c(0,  demographic.ages[(i - 11),2]*initial.prop.susceptible  , 0, 0, 
                                                                           (demographic.ages[(i - 11),2]*(1-initial.prop.susceptible))))
  }
  return(disease.state)
}


#############################
# This function reduces the proportion of people who are susceptible between min.age and max.age
# The proportion who are removed is given by the proportion.sus.to.remove variable
# Can be used to do supplementary vaccination
#############################
reduce.susceptibles <- function(min.age, max.age, disease.state, proportion.sus.to.remove, num.comps, susceptible.indices, recovered.indices, list.of.ages){
  k = which(list.of.ages == min.age)
  l = which(list.of.ages == max.age)
  susceptibles  =  susceptible.indices[k: l]
  disease.state [ recovered.indices[k:l]]  = round(disease.state [ recovered.indices[k:l]]   +    ( disease.state[ susceptibles  ] * ( proportion.sus.to.remove ) ))
  disease.state [ susceptibles ]  =   round((disease.state[ susceptibles ] * ( 1 - proportion.sus.to.remove ) ) )
  return(disease.state)
}

#############################
# Calculate the population by age group
#############################
population.by.age  <-  function(disease.state, num.comps){
  num.age.brackets  =  length(disease.state) / num.comps
  pop.by.age.all.brackets  =  matrix(0, num.age.brackets, 1)
  for (i in 1 : num.age.brackets){
    pop.by.age.all.brackets[i]  =  sum(disease.state[seq((i-1) * num.comps + 1 , i * num.comps)])
  }
  return(pop.by.age.all.brackets)
}


#############################
# This function scales the transmission rate so that the given value of R_0 would be achieved with the current mixing matrix.
# Can also specify the maximum age of transmission that we wish to consider (max.age parameter). This assumes that as measles is traditionally a childhood disease, that the 
# people who were spreading it were children, therefore any R_0's which were calculated previously were from children only. Debatable whether it really makes sense to do this.
#############################
calibrate.beta <- function (mixing.matrix, disease.state, max.age, time.step, infectious.period, R_0, list.of.ages, num.comps){
  k  =  which(list.of.ages == max.age)
  num.contacts  =  matrix(0, k, 1)
  for ( i in 1 : k ){
    num.contacts[i]  =  sum(mixing.matrix[, i]) * min(1, time.step / infectious.period)
  }
  
  pop.by.age = population.by.age(disease.state, num.comps)
  g = 0
  for (i in 1 : k){
    g = g + (pop.by.age[i] / sum(pop.by.age[1:k])) * num.contacts[i]
  }
  beta  = R_0 * min(1, time.step / infectious.period) / g
  return(beta)
}



###############################
# Calculate the force of infection by age given the current state, beta and the mixing matrix
###############################
foi.by.next.gen <- function ( mixing.matrix, disease.state, infectious.indices, time.step , infectious.period, beta, demographic.ages, num.comps){
  infectious.by.age  =  disease.state[infectious.indices]
  pop.by.age  =  population.by.age (disease.state, num.comps)
  foi.by.age  =  (colSums((infectious.by.age * t( mixing.matrix) * beta) ) / pop.by.age) * min( 1, time.step / infectious.period)
  return(foi.by.age)  
}


###############################
# Add births and remove deaths from population
###############################
demographics <- function(disease.state, birth.rate, death.rate, time.step){
  N = sum(disease.state)
  if( birth.rate > 1 ){
    birth.rate = birth.rate  / (1000 * 365)
  }
  births.average    =   birth.rate * N * time.step
  births.total      =   rpois(1, births.average)
  disease.state[1]  =   disease.state[1]  +  births.total
  if( death.rate > 1 ){
    death.rate = death.rate  / (1000 * 365)
  }
  prob.survival   =  1-(death.rate)*(time.step)
  disease.state   =  round(disease.state * prob.survival) 
  
  return(disease.state)
}



###############################
# Add infected migrants to the population to introduce infection. On average, there will be l introductions per year, where l is an input parameter
###############################
migrants.l.per.year <- function(disease.state, migrant.indices, time.step, l){
  k  =  length(migrant.indices)
  migrant.infecteds       =     rpois(k  , l * time.step/ (365 * k))
  disease.state[migrant.indices]           =     disease.state[migrant.indices] + migrant.infecteds
  return(disease.state)
}



###############################
# Functions which draw the number of people in each age and disease specific state at the next time step
###############################
draw.maternal.under1 <- function(x, time.step){
  u = time.step / 30
  mat.immune.loss  =  x[1]
  numbers = x[2]
  rmultinom(1, numbers , c(  (1-u) * (1 - mat.immune.loss * u),  (1-u) * mat.immune.loss * u , 0 ,0 , 0,
                             u * (1 - mat.immune.loss * u),  u * mat.immune.loss * u , 0 , 0, 0))
  
}


draw.sus.under1 <- function(x){
  u = time.step / 30
  foi = min(1, x[1])
  numbers = x[2]
  rmultinom(1, numbers , c( 0, (1-u)*(1-foi) , (1-u)*foi ,0 , 0, 0, u*(1-foi) , u*foi , 0, 0))
}


draw.exposed.under1 <- function(x){
  u = time.step / 30
  numbers = x[1]
  mu = min(1, time.step / incubation.period)
  rmultinom(1, numbers, c(0 , 0, (1-u) * (1- mu), (1-u ) * mu,  0  , 0 , 0 ,  u  * (1 - mu)  ,  u  *  mu, 0))
}


draw.infecteds.under1 <- function(x){
  u = time.step / 30
  rho = min(1, time.step / infectious.period)
  numbers = x[1]
  rmultinom(1, numbers , c(0 , 0 , 0,  (1-u) * (1- rho) , (1-u ) * rho , 0 , 0 , 0,  u  * (1 - rho) ,  u  *  rho))
}


draw.recovered.under1 <- function(x){
  u = time.step / 30
  numbers = x[1]
  rmultinom(1, numbers , c(0 , 0 ,  0, 0 ,  (1-u ) , 0 , 0 , 0 , 0 ,  u  ))
}


draw.sus <- function(x){
  u = time.step / 365
  foi = min(1, x[1])
  numbers = x[2]
  rmultinom(1, numbers , c( 0, (1-u)*(1-foi) , (1-u)*foi ,0 , 0, 0, u*(1-foi) , u*foi , 0, 0))
}

draw.exposed <- function(x){
  u = time.step / 365
  numbers = x[1]
  mu = min(1, time.step / incubation.period)
  rmultinom(1, numbers, c(0 , 0, (1-u) * (1- mu), (1-u ) * mu,  0  , 0 , 0 ,  u  * (1 - mu)  ,  u  *  mu, 0))
}


draw.infecteds <- function(x){
  u = time.step / 365
  rho = min(1, time.step / infectious.period)
  numbers = x[1]
  rmultinom(1, numbers , c(0 , 0 , 0,  (1-u) * (1- rho) , (1-u ) * rho , 0 , 0 , 0,  u  * (1 - rho) ,  u  *  rho))
}


draw.recovered <- function(x){
  u = time.step / 365
  numbers = x[1]
  rmultinom(1, numbers , c(0 , 0 ,  0, 0 ,  (1-u ) , 0 , 0 , 0 , 0 ,  u  ))
}



###############################
# Function which draws the populations for the next stage after a timestep for under 1's, as these are stratified monthly, whilst the rest of the population is stratified yearly.
###############################
draw.next.step.under1 <- function(disease.state, mixing.matrix, infectious.indices, 
                           time.step , infectious.period, beta,
                           demographic.ages, num.comps){
  
  updated.state =  matrix(0, length(disease.state), 1)
  
  foi.ages  <-   foi.by.next.gen( mixing.matrix, disease.state, 
                                  infectious.indices, time.step,
                                  infectious.period, beta,
                                  demographic.ages, num.comps)
  
  x              =      matrix(0, 12, 3)
  x[ , 1]        =      mat.immunity.loss[1 : 12]
  x[ , 2]        =      disease.state[maternal.indices[1:12]]
  mat.outs       =      apply(x, 1, draw.maternal.under1, time.step = time.step)
  
  x              =      matrix(0, 12, 2)
  x[ , 1]        =      foi.ages[1 : 12]
  x[ , 2]        =      disease.state[susceptible.indices[1:12]]
  sus.outs       =      apply(x, 1, draw.sus.under1)
  
  x1             =      matrix(0, 12, 2)
  x1[ , 1]       =      disease.state[exposed.indices[1:12]]
  exposed.out    =      apply(x1, 1, draw.exposed.under1)
  
  
  x2             =      matrix(0, 12, 2)
  x2[ , 1]       =      disease.state[infectious.indices[1:12]]
  inf.out        =      apply(x2, 1, draw.infecteds.under1)
  
  x3             =      matrix(0, 12, 2)
  x3[ , 1]       =      disease.state[recovered.indices[1:12]]
  recovered.out   =     apply(x3, 1, draw.recovered.under1)
  
  new.infected       =   sum(sus.outs[3, ])  +  sum(sus.outs[8, ])
  number.infectious  =   sum(exposed.out[4, ]) + sum(exposed.out[9, ])
  
  for (p in 1 : 12){
    updated.state[seq(((num.comps*( p-1 )) + 1) ,num.comps*(p+1))]   =   updated.state[seq(((num.comps*( p-1 )) + 1) , num.comps*(p+1))]  +
      sus.outs[ , p]  +  exposed.out[ , p]  +  inf.out[ , p] + recovered.out[ , p] + mat.outs[ , p]
  }
  return(list(updated.state, new.infected))
}



###############################
# Function which draws the populations for the next stage after a timestep for over 1's
###############################
draw.next.step.over1 <- function(disease.state, mixing.matrix, infectious.indices, 
                                  time.step , infectious.period, beta,
                                  demographic.ages, num.comps){
  
  updated.state = matrix(0, length(disease.state), 1)
  foi.ages  <-   foi.by.next.gen( mixing.matrix, disease.state,
                                  infectious.indices, time.step,
                                  infectious.period, beta,
                                  demographic.ages, num.comps)
  
  x              =      matrix(0, length(foi.ages) - 12, 2)
  x[ , 1]        =      foi.ages[-(1:12)]
  x[ , 2]        =      disease.state[susceptible.indices[-(1:12)]]
  sus.outs       =      apply(x, 1, draw.sus)
  
  x1             =      matrix(0, length(foi.ages) - 12, 2)
  x1[ , 1]       =      disease.state[exposed.indices[-(1:12)]]
  exposed.out    =      apply(x1, 1, draw.exposed)
  
  
  x2             =      matrix(0, length(foi.ages) - 12, 2)
  x2[ , 1]       =      disease.state[infectious.indices[-(1:12)]]
  inf.out        =      apply(x2, 1, draw.infecteds)
  
  x3             =      matrix(0, length(foi.ages) - 12, 2)
  x3[ , 1]       =      disease.state[recovered.indices[-(1:12)]]
  recovered.out   =     apply(x3, 1, draw.recovered)
  
  new.infected       =   sum(sus.outs[3, ])  +  sum(sus.outs[8, ])
  number.infectious  =   sum(exposed.out[4, ]) + sum(exposed.out[9, ])
  
  for (p in 13 : ((length(disease.state) / num.comps) - 1)){
    updated.state[seq(((num.comps*( p-1 )) + 1) ,num.comps*(p+1))]   =   updated.state[seq(((num.comps*( p-1 )) + 1) , num.comps*(p+1))]  +  
      sus.outs[ , p-12]  +  exposed.out[ , p-12]  +  inf.out[ , p-12] + recovered.out[ , p-12]
  }
  p = length(disease.state) / num.comps
  updated.state[seq(((num.comps*( p-1 )) + 1) ,num.comps*(p))]   =   updated.state[seq(((num.comps*( p-1 )) + 1) , num.comps*(p))]  +  
    sus.outs[1:num.comps , p-12]  +  exposed.out[1:num.comps , p-12]  +  inf.out[1:num.comps , p-12]  +  recovered.out[ 1:num.comps, p-12]
  
  return(list(updated.state, new.infected))
}



###############################
# Function which draws the populations for all age groups for the next timestep
###############################
draw.next.step.all <- function(disease.state, mixing.matrix, infectious.indices, 
                                 time.step , infectious.period, beta,
                                 demographic.ages, num.comps){
  
  output.under1  =  draw.next.step.under1(disease.state, mixing.matrix, infectious.indices, 
                                           time.step , infectious.period, beta,
                                           demographic.ages, num.comps)
  updated.state.under1  =  unlist(output.under1[1])
  new.infecteds.under1  =  unlist(output.under1[2])
  output.over1   =  draw.next.step.over1(disease.state, mixing.matrix, infectious.indices, 
                                          time.step , infectious.period, beta,
                                          demographic.ages, num.comps)
  updated.state.over1  =  unlist(output.over1[1])
  new.infecteds.over1  =  unlist(output.over1[2])
  
  updated.state  =  updated.state.under1   +   updated.state.over1
  new.infecteds  =  new.infecteds.under1   +   new.infecteds.over1
  return(list(updated.state, new.infecteds))
}




