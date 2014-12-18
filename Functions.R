#############################
# This function creates the initial disease state that the population is in. 
# Can specify the proportion of people who are susceptible to begin with using initial.prop.susceptible variable
#############################
initial.disease.state <- function(demographic.ages  ,   initial.prop.susceptible  ,  num.comps,  list.of.ages, mat.immunity.loss){
  disease.state                 =      matrix(0,num.comps*length(list.of.ages),1)
  
  for (i in 1:12){
    disease.state[(((i-1)*num.comps)+1):((i)*num.comps)]      =    round(c( (demographic.ages[1,2]*(initial.prop.susceptible / 12) * (1 - mat.immunity.loss[i])) , 
                                                                            (demographic.ages[1,2]*(initial.prop.susceptible / 12) * ( mat.immunity.loss[i])), 
                                                                              0, 
                                                                              0, 
                                                                            (demographic.ages[1,2]*(1-initial.prop.susceptible)/12)))
  }
  for (i in 13:length(list.of.ages)){
    disease.state[(((i-1)*num.comps)+1):((i)*num.comps)]      =    round(c(0,  
                                                                           demographic.ages[(i - 11),2]*initial.prop.susceptible,
                                                                           0, 
                                                                           0, 
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




#############################
# Calculate the proportion of the each age group who is susceptible
#############################
calc.sus.by.age <- function(disease.state, sus, num.comps){
  p = disease.state[sus]
  p1 = matrix(0, 98,1)
  p1[1] = sum(p[1:12])
  p1[2:98] = p[13:109]
  l = population.by.age(disease.state, num.comps)
  k2 = matrix(0, 98, 1)
  k2 = sum(l[1:12])
  k2[2:98] = l[13:109]
  return(p1/k2)
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
  v  =  x[3]
  rmultinom(1, numbers , c(  (1-u) * (1 - mat.immune.loss * u),  (1-u) * mat.immune.loss * u , 0 ,0 , 0,
                             u * (1 - mat.immune.loss * u) * (1 - v),  u * mat.immune.loss * u * (1 - v) , 0 , 0, u * v))
  
}


draw.sus.under1 <- function(x){
  u = time.step / 30
  foi = min(1, x[1])
  numbers = x[2]
  v = x[3]
  rmultinom(1, numbers , c( 0, (1 - u) * (1 - foi) , (1 - u) * foi ,0 , 0,
                            0, u * (1 - foi) * (1 - v) , u * foi * (1-v) , 0, u * v))
}


draw.exposed.under1 <- function(x){
  u = time.step / 30
  numbers = x[1]
  mu = min(1, time.step / incubation.period)
  rmultinom(1, numbers, c(0 , 0, (1-u) * (1- mu), (1-u ) * mu,  0  ,
                          0 , 0 ,  u  * (1 - mu)  ,  u  *  mu, 0))
}


draw.infecteds.under1 <- function(x){
  u = time.step / 30
  rho = min(1, time.step / infectious.period)
  numbers = x[1]
  rmultinom(1, numbers , c(0 , 0 , 0,  (1-u) * (1- rho) , (1-u ) * rho ,
                           0 , 0 , 0,  u  * (1 - rho) ,  u  *  rho))
}


draw.recovered.under1 <- function(x){
  u = time.step / 30
  numbers = x[1]
  rmultinom(1, numbers , c(0 , 0 ,  0, 0 ,  (1-u ) ,
                           0 , 0 , 0 , 0 ,  u  ))
}


draw.sus <- function(x){
  u = time.step / 365
  foi = min(1, x[1])
  numbers = x[2]
  rmultinom(1, numbers , c( 0, (1-u)*(1-foi) , (1-u)*foi ,0 , 0, 
                            0, u*(1-foi) , u*foi , 0, 0))
}

draw.exposed <- function(x){
  u = time.step / 365
  numbers = x[1]
  mu = min(1, time.step / incubation.period)
  rmultinom(1, numbers, c(0 , 0, (1-u) * (1- mu), (1-u ) * mu,  0  ,
                          0 , 0 ,  u  * (1 - mu)  ,  u  *  mu, 0))
}


draw.infecteds <- function(x){
  u = time.step / 365
  rho = min(1, time.step / infectious.period)
  numbers = x[1]
  rmultinom(1, numbers , c(0 , 0 , 0,  (1-u) * (1- rho) , (1-u ) * rho ,
                           0 , 0 , 0,  u  * (1 - rho) ,  u  *  rho))
}


draw.recovered <- function(x){
  u = time.step / 365
  numbers = x[1]
  rmultinom(1, numbers , c(0 , 0 ,  0, 0 ,  (1-u ) ,
                           0 , 0 , 0 , 0 ,  u  ))
}



###############################
# Function which draws the populations for the next stage after a timestep for under 1's, as these are stratified monthly, whilst the rest of the population is stratified yearly.
###############################
draw.next.step.under1 <- function(disease.state, mixing.matrix, infectious.indices, 
                           time.step , infectious.period, beta,
                           demographic.ages, num.comps, maternal.indices, 
                           mat.immunity.loss, vacc.prop){
  
  updated.state =  matrix(0, length(disease.state), 1)
  
  foi.ages  <-   foi.by.next.gen( mixing.matrix, disease.state, 
                                  infectious.indices, time.step,
                                  infectious.period, beta,
                                  demographic.ages, num.comps)
  
  x              =      matrix(0, 12, 3)
  x[ , 1]        =      mat.immunity.loss[1 : 12]
  x[ , 2]        =      disease.state[maternal.indices[1:12]]
  x[9, 3]        =      vacc.prop
  mat.outs       =      apply(x, 1, draw.maternal.under1, time.step = time.step)
  
  x              =      matrix(0, 12, 3)
  x[ , 1]        =      foi.ages[1 : 12]
  x[ , 2]        =      disease.state[susceptible.indices[1:12]]
  x[9, 3]        =      vacc.prop
  sus.outs       =      apply(x, 1, draw.sus.under1)
  
  x1             =      matrix(0, 12, 1)
  x1[ , 1]       =      disease.state[exposed.indices[1:12]]
  exposed.out    =      sapply(x1, draw.exposed.under1)
  
  
  x2             =      matrix(0, 12, 1)
  x2[ , 1]       =      disease.state[infectious.indices[1:12]]
  inf.out        =      sapply(x2, draw.infecteds.under1)
  
  x3             =      matrix(0, 12, 1)
  x3[ , 1]       =      disease.state[recovered.indices[1:12]]
  recovered.out   =     sapply(x3, draw.recovered.under1)
  
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
  
  x1             =      matrix(0, length(foi.ages) - 12, 1)
  x1[ , 1]       =      disease.state[exposed.indices[-(1:12)]]
  exposed.out    =      sapply(x1, draw.exposed)
  
  
  x2             =      matrix(0, length(foi.ages) - 12, 1)
  x2[ , 1]       =      disease.state[infectious.indices[-(1:12)]]
  inf.out        =      sapply(x2, draw.infecteds)
  
  x3             =      matrix(0, length(foi.ages) - 12, 1)
  x3[ , 1]       =      disease.state[recovered.indices[-(1:12)]]
  recovered.out   =     sapply(x3, draw.recovered)
  
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
                               demographic.ages, num.comps, maternal.indices, 
                               mat.immunity.loss, vacc.prop){
  
  output.under1  =  draw.next.step.under1(disease.state, mixing.matrix, infectious.indices, 
                                          time.step , infectious.period, beta,
                                          demographic.ages, num.comps, maternal.indices, 
                                          mat.immunity.loss, vacc.prop)
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





###############################
# Function which draws the populations for the next stage after a timestep for under 1's, with no infection taking place.
# This is separate as these are stratified monthly, whilst the rest of the population is stratified yearly. 
###############################
age.under1 <- function(disease.state, mixing.matrix, infectious.indices, 
                                  time.step , infectious.period, beta,
                                  demographic.ages, num.comps, maternal.indices, 
                                  mat.immunity.loss, vacc.prop){
  
  updated.state =  matrix(0, length(disease.state), 1)
  
  foi.ages  <-   matrix(0, length(disease.state[infectious.indices]), 1)
  
  
  x              =      matrix(0, 12, 3)
  x[ , 1]        =      mat.immunity.loss[1 : 12]
  x[ , 2]        =      disease.state[maternal.indices[1:12]]
  x[9, 3]        =      vacc.prop
  mat.outs       =      apply(x, 1, draw.maternal.under1, time.step = time.step)
  
  
  x              =      matrix(0, 12, 3)
  x[ , 1]        =      foi.ages[1 : 12]
  x[ , 2]        =      disease.state[susceptible.indices[1:12]]
  x[9, 3]        =      vacc.prop
  sus.outs       =      apply(x, 1, draw.sus.under1)
  
  
  x3             =      matrix(0, 12, 1)
  x3[ , 1]       =      disease.state[recovered.indices[1:12]]
  recovered.out   =     sapply(x3, draw.recovered.under1)
 
  
  for (p in 1 : 12){
    updated.state[seq(((num.comps*( p-1 )) + 1) ,num.comps*(p+1))]   =   updated.state[seq(((num.comps*( p-1 )) + 1) , num.comps*(p+1))]  +
      sus.outs[ , p]  +  exposed.out[ , p]  +  inf.out[ , p] + recovered.out[ , p] + mat.outs[ , p]
  }
  
  return(updated.state)
}



###############################
# Function which draws the populations for the next stage after a timestep for over 1's, with no infection taking place
###############################
age.over1 <- function(disease.state, mixing.matrix, infectious.indices, 
                                 time.step , infectious.period, beta,
                                 demographic.ages, num.comps){
  
  updated.state = matrix(0, length(disease.state), 1)
  foi.ages  <-   matrix(0, length(disease.state[infectious.indices]), 1)
  
  x              =      matrix(0, length(foi.ages) - 12, 2)
  x[ , 1]        =      foi.ages[-(1:12)]
  x[ , 2]        =      disease.state[susceptible.indices[-(1:12)]]
  sus.outs       =      apply(x, 1, draw.sus)
  
  x3             =      matrix(0, length(foi.ages) - 12, 2)
  x3[ , 1]       =      disease.state[recovered.indices[-(1:12)]]
  recovered.out   =     sapply(x3, draw.recovered)
  
  for (p in 13 : ((length(disease.state) / num.comps) - 1)){
    updated.state[seq(((num.comps*( p-1 )) + 1) ,num.comps*(p+1))]   =   updated.state[seq(((num.comps*( p-1 )) + 1) , num.comps*(p+1))]  +  
      sus.outs[ , p-12]  +  exposed.out[ , p-12]  +  inf.out[ , p-12] + recovered.out[ , p-12]
  }
  p = length(disease.state) / num.comps
  updated.state[seq(((num.comps*( p-1 )) + 1) ,num.comps*(p))]   =   updated.state[seq(((num.comps*( p-1 )) + 1) , num.comps*(p))]  +  
    sus.outs[1:num.comps , p-12]  +  exposed.out[1:num.comps , p-12]  +  inf.out[1:num.comps , p-12]  +  recovered.out[ 1:num.comps, p-12]
  
  return(updated.state)
}



###############################
# Function which ages the population by a given length of time, with no infection occurring in that period
###############################
age.population <- function(disease.state, mixing.matrix, infectious.indices, 
                               time.step , infectious.period, beta,
                               demographic.ages, num.comps, maternal.indices, 
                               mat.immunity.loss, vacc.prop){
  
  updated.state.under1  =  age.under1(disease.state, mixing.matrix, infectious.indices, 
                                          time.step , infectious.period,
                                          demographic.ages, num.comps, maternal.indices, 
                                          mat.immunity.loss, vacc.prop)
  
  updated.state.over1   =  age.over1(disease.state, mixing.matrix, infectious.indices, 
                                         time.step , infectious.period,
                                         demographic.ages, num.comps)
  
  updated.state  =  updated.state.under1   +   updated.state.over1
  
  return(updated.state)
}



###############################
# Function which will perform aging of the population along with taking care of demographics for a given length of time
###############################
age.population.with.demographics <- function(disease.state, mixing.matrix, infectious.indices, 
                                              time.step , infectious.period, 
                                              demographic.ages, num.comps, maternal.indices, 
                                              mat.immunity.loss, vacc.prop, num.steps){
  
  for (i in 1 : num.steps){
    disease.state  =  demographics(disease.state, birth.rate, death.rate, time.step)
    disease.state  =  age.population(disease.state, mixing.matrix, infectious.indices, 
                                     time.step , infectious.period, beta,
                                     demographic.ages, num.comps, maternal.indices, 
                                     mat.immunity.loss, vacc.prop)
    
    
  }
  
  return(disease.state)
}



###############################
# Produce plots of various state descriptors
###############################
produce.plots <-function(i, prop.sus, new.infections, all.infections, disease.state, susceptible.indices, num.comps){
  
  par(mfrow = c(2,2))
  plot(1:i, 1- prop.sus[1:i], type = "l", col = 'seagreen3', xlab = 'Time', ylab = 'Population Immunity')
  plot(0:97, 1 - calc.sus.by.age(disease.state, susceptible.indices, num.comps), type = "l", col = 'purple', ylab = 'Immunity by age', xlab = 'Age')
  plot(1:i, new.infections[1:i], type = "l", col = 'dodgerblue3', xlab = 'Time', ylab = 'Infected by day')
  plot(1:i, all.infections[1:i], type = "l", col = 'red', xlab = 'Time', ylab = 'All infectious')  
}





###############################
# Make it so that it is possible to receive the output of a function which is in a list, and automatically split it into components using list[..]
###############################
list <- structure(NA,class="result")
"[<-.result" <- function(x,...,value) {
  args <- as.list(match.call())
  args <- args[-c(1:2,length(args))]
  length(value) <- length(args)
  for(i in seq(along=args)) {
    a <- args[[i]]
    if(!missing(a)) eval.parent(substitute(a <- v,list(a=a,v=value[[i]])))
  }
  x
}





###############################
# Run simulations using previously defined functions
###############################
Run.simulations <- function(num.steps, disease.state, mixing.matrix, infectious.indices,
                            birth.rate, death.rate,
                            time.step , infectious.period, seasonal.scaling,
                            demographic.ages, num.comps, maternal.indices, 
                            mat.immunity.loss, vacc.prop, l, do.plots, initial.time){
  
  # Use initial.time to specify where in the seasonal process the disease is. 
  # initial.time = 0 or 280 gives peak initial transmission which decreases or average initial transmission which increases
  
  t  =  initial.time
  # seasonal.scaling is the increase and decrease seen in the transmission rates during peak and lowest transmission. call this beta_1 for ease of use.
  beta_1 = seasonal.scaling
  new.infections  =  matrix(0, num.steps, 1)
  all.infections  =  matrix(0, num.steps, 1)
  prop.sus        =  matrix(0, num.steps, 1)
  for (run in 1 : num.steps){
    pop.by.age  =  population.by.age(disease.state, num.comps)
    for ( i in 1:length(pop.by.age)){
      mixing.matrix  [i, ]  =  pop.by.age[i]       # this produces a mixing matrix which is equivalent to uniform mixing
    }
    disease.state  =  demographics(disease.state, birth.rate, death.rate, time.step)
    disease.state  =  migrants.l.per.year(disease.state, migrant.indices, time.step, l)
    
    beta  =  calibrate.beta(mixing.matrix, disease.state, max.age, time.step, infectious.period, R_0, list.of.ages, num.comps) * (1 + beta_1 * cos(2 * pi * t /365))
    
    list[disease.state, new.infections[run]]  =  draw.next.step.all(disease.state, mixing.matrix, infectious.indices, 
                                                                time.step , infectious.period, beta,
                                                                demographic.ages, num.comps, maternal.indices, 
                                                                mat.immunity.loss, vacc.prop)
    
    all.infections[run]  =   sum(disease.state[infectious.indices])
    prop.sus[run]        =   sum(disease.state[susceptible.indices]) / sum(disease.state)
    if ((run %% 100 == 0) & (do.plots == 1)){
      produce.plots(run, prop.sus, new.infections, all.infections, disease.state, susceptible.indices, num.comps)
    }
    t  =  t  +  time.step
  }
  return(list(disease.state, new.infections, all.infections, prop.sus))
}




###################################
# Run simulations using an initial disease.state for a given number of years, where one infected individual is introduced at the beginning of each year, and the population is then aged between.
# The population prior to the introductionInclude sia's in this 
###################################
Run.continuous.sims <- function(total.years, number.of.replicates, initial.disease.state, 
                                mixing.matrix, infectious.indices, time.step, 
                                infectious.period, beta_1, demographic.ages,
                                num.comps, maternal.indices, 
                                mat.immunity.loss, vacc.prop, 
                                av.migrants.per.year, do.plots, 
                                initial.time, min.age.sia, max.age.sia, 
                                sia.proportion, list.of.ages,
                                birth.rate, death.rate){
  #
  total.infections.by.year  =  matrix(0, number.of.replicates, total.years)
  for (i in 1 : number.of.replicates){
  #  print(paste('Started', i, 'of', number.of.replicates))
    disease.state = initial.disease.state
    disease.state[infectious.indices]  =  0
   # print(sum(disease.state[susceptible.indices]) / sum(disease.state))
    for ( j in 1 : total.years){
      disease.state[infectious.indices[20]]  =  disease.state[infectious.indices[20]] + 1 
      
      list[disease.state, new.infections, all.infections, prop.sus] = Run.simulations(365, disease.state, mixing.matrix, 
                                                                                      infectious.indices, birth.rate, death.rate,
                                                                                      time.step , infectious.period,
                                                                                      beta_1, demographic.ages, num.comps, 
                                                                                      maternal.indices,  mat.immunity.loss, vacc.prop, 
                                                                                      av.migrants.per.year, do.plots, initial.time)
      
      total.infections.by.year[i, j]  =  sum(new.infections)
      
      if ( j %% sia.period == 0){
       # print(sum(disease.state[susceptible.indices]) / sum(disease.state))
        disease.state  =  reduce.susceptibles(min.age.sia, max.age.sia, disease.state, sia.proportion, num.comps, susceptible.indices, recovered.indices, list.of.ages)
      #  print(sum(disease.state[susceptible.indices]) / sum(disease.state))
      }
    }
   # print(paste('Finished', i, 'of', number.of.replicates))
  }
  return(total.infections.by.year)
}





###################################
# Run simulations in parallel using an initial disease.state for a given number of years, where one infected individual is introduced at the beginning of each year, and the population is then aged between.
# The population prior to the introductionInclude sia's in this 
###################################
numerous.sims <- function(num.replicates, num.years, initial.disease.state1, initial.start.time){
  infections.by.year = matrix(0, num.replicates, num.years)
  state  =  matrix(0, num.replicates, length(initial.disease.state))
  infections.by.year <- foreach (i = 1 : num.replicates, .export=ls(envir=globalenv())) %dopar%{
    print(paste('Started', i, 'of', num.replicates))
    infections.by.year[i,  ] <- Run.continuous.sims(num.years, 1, initial.disease.state1, mixing.matrix,
                                                    infectious.indices, time.step , infectious.period, beta_1,
                                                    demographic.ages, num.comps, maternal.indices, 
                                                    mat.immunity.loss, vacc.prop, av.migrants.per.year, 
                                                    do.plots, initial.start.time , min.age.sia, max.age.sia,
                                                    sia.proportion, list.of.ages, birth.rate, death.rate)  
  }
  return(infections.by.year)
}




###################################
# Run simulations where we can change the birth rate and vaccination rates multiple times
###################################
diff.br.vacc.sims <- function(initial.disease.state, 
                              initial.start.time, multi.vacc, multi.br, 
                              sia.years, sia.numbers, vacc.success,
                              mixing.matrix, infectious.indices, input.infections.all.times, 
                              death.rate, time.step , infectious.period, seasonal.scaling = beta_1,
                              demographic.ages, num.comps, maternal.indices, init.pop.size,
                              mat.immunity.loss, l  =  0, do.plots , initial.time = initial.start.time){
  
  sim.infs.by.year  =  matrix(0, length(multi.br), 1)
  count = 1
  input.disease.state = round(initial.disease.state * init.pop.size / sum(initial.disease.state))
  for (i in 1: length(multi.br)){
    #print(paste('On', i, 'of', length(multi.br)))
    birth.rate  =  multi.br[i]
    vacc  =  vacc.success*(multi.vacc[i]/100)
    input.infections = 0
    if( input.infections.all.times == 1){
      input.infections = 1
    }
    if( i == length(multi.br)){
      input.infections = 1
    }
    input.disease.state[infectious.indices[12]]  =   input.disease.state[infectious.indices[12]]   +  input.infections
    list[disease.state, new.infections, all.infections, prop.sus]  =  Run.simulations(365, input.disease.state, mixing.matrix, infectious.indices,
                                                                                      birth.rate, death.rate,
                                                                                      time.step , infectious.period, seasonal.scaling = beta_1,
                                                                                      demographic.ages, num.comps, maternal.indices, 
                                                                                      mat.immunity.loss, vacc, l  =  0, do.plots , initial.time = initial.start.time)
    if( count < (length(sia.years) + 1)){
      if ( i == sia.years[count] ){
        num  =  sia.numbers[count]
        l = which(list.of.ages == 9/12)
        k = which(list.of.ages == 5)
        target.pop  =  sum(disease.state[((l * num.comps) + 1):(k * num.comps)])
        total.prop.SIA  =  num / target.pop
        print(total.prop.SIA)
        disease.state  =  reduce.susceptibles(9/12, 5, disease.state, total.prop.SIA, num.comps, susceptible.indices, recovered.indices, list.of.ages)
        count  =  count + 1
      }  
    }
    
    input.disease.state   =   disease.state
    sim.infs.by.year[i]  =  sum(new.infections)
  }
  return(sim.infs.by.year)
}




###################################
# In this function we can change the vaccination and birth rates multiple times
# Run simulations in parallel using an initial disease.state for a given number of years, where one infected individual is introduced at the beginning of each year, and the population is then aged between.
# The population prior to the introductionInclude sia's in this. 
###################################
numerous.sims.diff.br.vacc <- function(num.replicates, initial.disease.state1, 
                                       initial.start.time, multi.vacc, multi.br, input.infections.all.times,
                                       init.pop.size, sia.years, sia.numbers, vacc.success){
  
  sim.infs.by.year   =  matrix(0 , num.replicates, length(multi.br))
  state  =  matrix(0, num.replicates, length(initial.disease.state))
  sim.infs.by.year <- foreach (run = 1 : num.replicates, .export=ls(envir=globalenv()))    %dopar%  {

    print(paste('Started', run, 'of', num.replicates))
    sim.infs.by.year[run, ] = diff.br.vacc.sims (initial.disease.state1, 
                                  initial.start.time, multi.vacc, multi.br, 
                                  sia.years, sia.numbers, vacc.success, 
                                  mixing.matrix, infectious.indices, input.infections.all.times,
                                  death.rate, time.step , infectious.period, seasonal.scaling = beta_1,
                                  demographic.ages, num.comps, maternal.indices, init.pop.size,
                                  mat.immunity.loss, l  =  0, do.plots , initial.time = initial.start.time)
  }
  return(sim.infs.by.year)
}




###################################
# Output coefficient of variation for 10 year periods, with 1 year moving window. 
# Also output mean birth and vaccination rates along with cumulative incidence in those 10 years
###################################
output.coeff.of.var.and.other.variables <- function(window.length){

  cases.by.country.by.year = read.csv("Measles_cases_by_year.csv")
  Birth.rates = read.csv("Birth_rates.csv")
  pop.by.year = read.csv("Population_by-country_by_year.csv")
  
  
  
  African.pop.by.year = subset(pop.by.year, pop.by.year$WHO_REGION == "AFR")
  African_vaccination = read.csv("African_vaccination_rates.csv", stringsAsFactors = FALSE)
  African_birth_rates = subset(Birth.rates, Birth.rates$WHO_REGION == "AFR")
  African_data = subset(cases.by.country.by.year, cases.by.country.by.year$WHO_REGION == "AFR")
  
  
  num.windows = 25 + (10 - window.length)
  
  coeff.var = matrix(0, length(African_data[ , 1]), num.windows)
  for ( backwards in 0 : (num.windows - 1)){
    for ( i in 1 : length(African_data[ , 1])){
      coeff.var[i, num.windows - backwards]  =  sd(African_data[i, (5 + backwards):(5 + window.length - 1 + backwards )], na.rm = TRUE) /  mean(as.numeric(African_data[i, (5 + backwards):(5 + window.length - 1 + backwards)]), na.rm = TRUE)
    }
  }
  
  
  
  incidence.per.1000 = matrix(0, length(African_data[ , 1]), num.windows)
  for ( backwards in 0 : (num.windows - 1)){
    for ( i in 1 : length(African_data[ , 1])){
      incidence.per.1000[i, num.windows - backwards]  =  sum(as.numeric(African_data[i,  (5 + backwards):(5 + window.length - 1 + backwards)]) / as.numeric(African.pop.by.year[i, (54 - backwards) : (54 - (window.length - 1) - backwards)]), na.rm = TRUE) * 1000
    }
  }
  
  
  
  
  mean.vac  =  matrix(0, length(African_vaccination[, 1]), num.windows) 
  for (j in 0 :  (num.windows - 1)){
    for ( i in 1 : length(African_vaccination[, 1])){
      mean.vac[i, j + 1]  =  mean(as.numeric(African_vaccination[i, (2 + j):(2 + window.length - 1  + j )]), na.rm = TRUE)
    }
  }
  
  African_birth_rates$X2013 = African_birth_rates$X2012
  African_birth_rates2 = African_birth_rates[, -54]
  African_birth_rates2 =African_birth_rates2 [ ,-(1:20)]
  mean.br  =  matrix(0, length(African_birth_rates[, 1]), num.windows) 
  for (j in 0 :  (num.windows - 1)){
    for ( i in 1 : length(African_birth_rates[, 1])){
      mean.br[i, j+1]  =  mean(as.numeric(African_birth_rates2[i, (1 + j) : (1 + window.length - 1  + j)]), na.rm = TRUE)
    }
  }
  
  return(list(coeff.var, incidence.per.1000, mean.vac, mean.br))
}







###################################
# Plot coeff of variance against incidence for a given period of time
###################################
plot.coeff.var <- function(African.countries, countries.to.plot, 
                           start.year, window.length, 
                           scaling, text.size){
  
  list[coeff.var, incidence.per.1000, mean.vac, mean.br] = output.coeff.of.var.and.other.variables(window.length)
  
  k = NULL
  if (length(countries.to.plot)  <  length(African.countries)){
    for ( i in 1 : length(countries.to.plot)){
      p = which(African.countries == countries.to.plot[i])
      k = c(k, p)      
    }
  } else {k = seq(1, 47)}
  j = start.year - 1979
  African.cov = data.frame(labels1 = countries.to.plot, Coeff.of.var = coeff.var[k, j], Incidence = incidence.per.1000[k, j], BR = mean.br[k, j], mean.vacc = mean.vac[k, j], Inverse.BR = 1/mean.br[k, j])
  #quartz()
  to.plot = ggplot(African.cov, aes(x = Coeff.of.var, y = Incidence, label = labels1)) 
  if(scaling == 'none')
  {
    qqq1 = to.plot + geom_point(aes(size = BR, colour = mean.vacc)) + scale_color_gradient(high = "blue", low = "red") +
      scale_y_continuous(limits = c(0, 100) ) + scale_x_continuous(limits = c(0, 3.5) ) +
      labs(x = "Coefficient of variation", y = paste('Incidence per 1000:',toString(start.year),'-',toString(start.year + window.length - 1)))  +
      theme(axis.text.x = element_text(colour="black"), axis.text.y = element_text(colour="black"))  + 
      geom_text(size = text.size, angle = 45) 
  }
  
  
  if(scaling == 'log')
  {
    qqq1 = to.plot + geom_point(aes(size = BR, colour = mean.vacc)) + scale_color_gradient(high = "blue", low = "red") +
      scale_y_continuous(limits = c(0.001, 100), trans= 'log10' ) + scale_x_continuous(limits = c(0, 3.5) ) +
      labs(x = "Coefficient of variation", y = paste('Incidence per 1000:',toString(start.year),'-',toString(start.year + window.length - 1)))  +
      theme(axis.text.x = element_text(colour="black"), axis.text.y = element_text(colour="black"))  + 
      geom_text(size = 3, angle = 45) 
  }
  
  if(scaling == 'sqrt')
  {
    qqq1 = to.plot + geom_point(aes(size = BR, colour = mean.vacc)) + scale_color_gradient(high = "blue", low = "red") +
      scale_y_continuous(limits = c(0, 100), trans= 'sqrt' ) + scale_x_continuous(limits = c(0, 3.5) ) +
      labs(x = "Coefficient of variation", y = paste('Incidence per 1000:',toString(start.year),'-',toString(start.year + window.length - 1)))  +
      theme(axis.text.x = element_text(colour="black"), axis.text.y = element_text(colour="black"))  + 
      geom_text(size = 3, angle = 45) 
  }
  
  
  
  print(qqq1  )
}
