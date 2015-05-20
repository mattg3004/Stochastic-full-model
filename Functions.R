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
reduce.susceptibles <- function(min.age, max.age, disease.state, proportion.sus.to.remove, num.comps, 
                                susceptible.indices, recovered.indices, list.of.ages, vacc.success){
  k = which(list.of.ages == min.age)
  l = which(list.of.ages == max.age)
  susceptibles  =  susceptible.indices[k: l]
  disease.state [ recovered.indices[k:l]]  = round(disease.state [ recovered.indices[k:l]]   +    vacc.success * ( disease.state[ susceptibles  ] * ( proportion.sus.to.remove ) )) 
  disease.state [ susceptibles ]  =   round(vacc.success * (disease.state[ susceptibles ] * ( 1 - proportion.sus.to.remove ) ) )
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
                                          time.step , infectious.period, beta,
                                          demographic.ages, num.comps, maternal.indices, 
                                          mat.immunity.loss, vacc.prop)
  
  updated.state.over1   =  age.over1(disease.state, mixing.matrix, infectious.indices, 
                                         time.step , infectious.period, beta,
                                         demographic.ages, num.comps)
  
  updated.state  =  updated.state.under1   +   updated.state.over1
  
  return(updated.state)
}



###############################
# Function which will perform aging of the population along with taking care of demographics for a given length of time
###############################
age.population.with.demographics <- function(disease.state, mixing.matrix, infectious.indices, 
                                              time.step , infectious.period, beta,
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
    #  print(sum(disease.state[susceptible.indices]) / sum(disease.state))
    # print(vacc.prop)
    if(run %% 50 == 0){
      #   plot(calc.sus.by.age(disease.state, susceptible.indices, num.comps))
    }
    pop.by.age  =  population.by.age(disease.state, num.comps)
    for ( i in 1:length(pop.by.age)){
      mixing.matrix  [i, ]  =  pop.by.age[i]       # this produces a mixing matrix which is equivalent to uniform mixing
    }
    disease.state  =  demographics(disease.state, birth.rate, death.rate, time.step)
    disease.state  =  migrants.l.per.year(disease.state, migrant.indices, time.step, l)
    
    beta  =  calibrate.beta(mixing.matrix, disease.state, max.age, time.step, infectious.period, R_0, list.of.ages, num.comps) * (1 + beta_1 * cos(2 * pi * t /365))
    if (vacc.prop > 1){
      vacc.prop  =  vacc.prop / 100
    }
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
# The population prior to the introduction Include sia's in this 
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
# Run simulations in parallel using an initial disease.state for a given number of years, 
# where one infected individual is introduced at the beginning of each year, and the population is then aged between.
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
                              mat.immunity.loss, l  , do.plots , 
                              initial.time = initial.start.time, unreached.pop,
                              vacc.efficiency, additional.campaigns.in.large.outbreak,
                              large.outbreak.threshold, additional.campaign.size,
                              constant.sia.prop, sia.prop){
  
  
  #print(l)
  sim.infs.by.year  =  matrix(0, length(multi.br), 1)
  incidence.by.year  =  matrix(0, length(multi.br), 1)
  count = 1
  input.disease.state = round(initial.disease.state * init.pop.size / sum(initial.disease.state))
  for (i in 1: length(multi.br)){
    
    #print(paste('On', i, 'of', length(multi.br)))
    birth.rate  =  multi.br[i]
    vacc  =  max((vacc.success*(multi.vacc[i]/100)) - unreached.pop, 0)
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
                                                                                      mat.immunity.loss, vacc, l , do.plots , initial.time = initial.start.time)
    if( count < (length(sia.years) + 1)){
      if ( i == sia.years[count] ){
        num  =  sia.numbers[count]
        p = which(list.of.ages == 9/12)
        k = which(list.of.ages == 5)
        target.pop  =  sum(disease.state[((p * num.comps) + 1):(k * num.comps)])
        total.prop.SIA  =  max(0, (vacc.success * (num / target.pop) * vacc.efficiency) - unreached.pop)
        total.prop.SIA  =  min(1, total.prop.SIA)
        if(constant.sia.prop == 1){
          total.prop.SIA  =  sia.prop
        }
        # print(paste("population =",sum(disease.state),"target population =", target.pop))
        #  print(total.prop.SIA)
        disease.state  =  reduce.susceptibles(9/12, 5, disease.state, total.prop.SIA, num.comps, 
                                              susceptible.indices, recovered.indices, list.of.ages, vacc.success)
        count  =  count + 1
      }  
    }
    
    input.disease.state   =  disease.state
    sim.infs.by.year[i]   =  sum(new.infections)
    incidence.by.year[i]  =  1000*( sum(new.infections) / sum(disease.state) )
    if ( (sum(new.infections) > large.outbreak.threshold) & (additional.campaigns.in.large.outbreak == 1)){
      disease.state  =  reduce.susceptibles(9/12, 5, disease.state, additional.campaign.size, num.comps, 
                                            susceptible.indices, recovered.indices, list.of.ages, vacc.success)
    }
  }
  return(list(sim.infs.by.year, incidence.by.year))
}




###################################
# In this function we can change the vaccination and birth rates multiple times
# Run simulations in parallel using an initial disease.state for a given number of years, 
# where one infected individual is introduced at the beginning of each year, and the population is then aged between.
###################################
numerous.sims.diff.br.vacc <- function(num.replicates, initial.disease.state1, 
                                       initial.start.time, multi.vacc, multi.br, 
                                       input.infections.all.times, init.pop.size, 
                                       sia.years, sia.numbers, vacc.success,
                                       num.migrants.per.year, unreached.pop,
                                       vacc.efficiency, additional.campaigns.in.large.outbreak,
                                       large.outbreak.threshold, additional.campaign.size,
                                       constant.sia.prop, sia.prop){
  
  sim.infs.by.year   = list()
  state  =  matrix(0, num.replicates, length(initial.disease.state))
  sim.infs.by.year <- foreach (run = 1 : num.replicates, .export=ls(envir=globalenv()))    %dopar%  {
    
    print(paste('Started', run, 'of', num.replicates))
    sim.infs.by.year[[run ]] = diff.br.vacc.sims (initial.disease.state1, 
                                                  initial.start.time, multi.vacc, multi.br, 
                                                  sia.years, sia.numbers, vacc.success, 
                                                  mixing.matrix, infectious.indices, input.infections.all.times,
                                                  death.rate, time.step , infectious.period, seasonal.scaling = beta_1,
                                                  demographic.ages, num.comps, maternal.indices, init.pop.size,
                                                  mat.immunity.loss, l = num.migrants.per.year, do.plots , 
                                                  initial.time = initial.start.time, unreached.pop,
                                                  vacc.efficiency, additional.campaigns.in.large.outbreak,
                                                  large.outbreak.threshold, additional.campaign.size,
                                                  constant.sia.prop, sia.prop)
    
    
  }
  return(sim.infs.by.year)
}




###################################
# Output coefficient of variation for 10 year periods, with 1 year moving window. 
# Also output mean birth and vaccination rates along with cumulative incidence in those 10 years
###################################
output.coeff.of.var.and.other.variables <- function(window.length){

 # cases.by.country.by.year = read.csv("Measles_cases_by_year.csv")
 # Birth.rates = read.csv("Birth_rates.csv")
 # pop.by.year = read.csv("Population_by-country_by_year.csv")
  
  
  
  African.pop.by.year = subset(pop.by.year, pop.by.year$WHO_REGION == "AFR")
  #African_vaccination = read.csv("African_vaccination_rates.csv", stringsAsFactors = FALSE)
  African_birth_rates = subset(Birth.rates, Birth.rates$WHO_REGION == "AFR")
  African_data = subset(cases.by.country.by.year, cases.by.country.by.year$WHO_REGION == "AFR")
  
  
  num.windows = 25 + (10 - window.length)
  
  coeff.var = matrix(0, length(African_data[ , 1]), num.windows)
  for ( j in 0 : (num.windows - 1)){
    for ( i in 1 : length(African_data[ , 1])){
      coeff.var[i, num.windows - j]  =  sd(African_data[i, (5 + j):(5 + window.length - 1 + j )], na.rm = TRUE) /  
        mean(as.numeric(African_data[i, (5 + j):(5 + window.length - 1 + j)]), na.rm = TRUE)
      if(is.na( mean(as.numeric(African_data[i, (5 + j):(5 + window.length - 1 + j)]), na.rm = TRUE)) == FALSE ){
        if( mean(as.numeric(African_data[i, (5 + j):(5 + window.length - 1 + j)]), na.rm = TRUE) == 0) {
          coeff.var[i, num.windows - j]  =  0
        } 
      }
    }
  }
  
  
  
  incidence.per.1000 = matrix(0, length(African_data[ , 1]), num.windows)
  for ( j in 0 : (num.windows - 1)){
    for ( i in 1 : length(African_data[ , 1])){
      incidence.per.1000[i, num.windows - j]  =  sum(as.numeric(African_data[i,  (5 + j):(5 + window.length - 1 + j)]) / 
                                                               as.numeric(African.pop.by.year[i, (54 - j) : (54 - (window.length - 1) - j)]), na.rm = TRUE) * 1000
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
                           start.year, coeff.var, incidence.per.1000, 
                           mean.vac, mean.br, window.length, 
                           scaling, text.size, upper.limit, column.number ){
  
  
  
  k = NULL
  if (length(countries.to.plot)  <  length(African.countries)){
    for ( i in 1 : length(countries.to.plot)){
      p = which(African.countries == countries.to.plot[i])
      k = c(k, p)      
    }
  } else {k = seq(1, 47)}
  j = start.year - 1979
  j = column.number
  African.cov = data.frame(labels1 = countries.to.plot, Coeff.of.var = coeff.var[k, j], Incidence = incidence.per.1000[k, j], BR = mean.br[k, j], mean.vacc = mean.vac[k, j], Inverse.BR = 1/mean.br[k, j])
  #quartz()
  to.plot = ggplot(African.cov, aes(x = Coeff.of.var, y = Incidence, label = labels1)) 
  if(scaling == 'none')
  {
    qqq1 = to.plot + geom_point(aes(size = BR, colour = mean.vacc)) + scale_color_gradient(high = "blue", low = "red") +
      scale_y_continuous(limits = c(0, upper.limit) ) + scale_x_continuous(limits = c(0, 3.5) ) +
      labs(x = "Coefficient of variation", y = paste('Incidence per 1000:',toString(start.year),'-',toString(start.year + window.length - 1)))  +
      theme(axis.text.x = element_text(colour="black"), axis.text.y = element_text(colour="black"))  + 
      geom_text(size = text.size, angle = 45) 
  }
  
  
  if(scaling == 'log')
  {
    qqq1 = to.plot + geom_point(aes(size = BR, colour = mean.vacc)) + scale_color_gradient(high = "blue", low = "red") +
      scale_y_continuous(limits = c(0.001, upper.limit), trans= 'log10' ) + scale_x_continuous(limits = c(0, 3.5) ) +
      labs(x = "Coefficient of variation", y = paste('Incidence per 1000:',toString(start.year),'-',toString(start.year + window.length - 1)))  +
      theme(axis.text.x = element_text(colour="black"), axis.text.y = element_text(colour="black"))  + 
      geom_text(size = 3, angle = 45) 
  }
  
  if(scaling == 'sqrt')
  {
    qqq1 = to.plot + geom_point(aes(size = BR, colour = mean.vacc)) + scale_color_gradient(high = "blue", low = "red") +
      scale_y_continuous(limits = c(0, upper.limit), trans = 'sqrt' ) + scale_x_continuous(limits = c(0, 3.5) ) +
      labs(x = "Coefficient of variation", y = paste('Incidence per 1000:',toString(start.year),'-',toString(start.year + window.length - 1)))  +
      theme(axis.text.x = element_text(colour="black"), axis.text.y = element_text(colour="black"))  + 
      geom_text(size = 3, angle = 45) 
  }
  
  
  
  print(qqq1)
}




###################################
# Plot coeff of variance against incidence for a given period of time
###################################
plot.coeff.var.with.inputs <- function(African.countries, countries.to.plot, 
                                       start.year, window.length, 
                                       scaling, text.size){
  
 
  
  k = NULL
  if (length(countries.to.plot)  <  length(African.countries)){
    for ( i in 1 : length(countries.to.plot)){
      p = which(African.countries == countries.to.plot[i])
      k = c(k, p)      
    }
  } else {k = seq(1, 47)}
  j = start.year - 1979
  African.cov = data.frame(labels1 = countries.to.plot[1:length(k)], Coeff.of.var = coeff.var[k, j], Incidence = incidence.per.1000[k, j], BR = mean.br[k, j], mean.vacc = mean.vac[k, j], Inverse.BR = 1/mean.br[k, j])
  #quartz()
  to.plot = ggplot(African.cov, aes(x = Coeff.of.var, y = Incidence, label = labels1)) 
  if(scaling == 'none')
  {
    qqq1 = to.plot + geom_point(aes(size = BR, colour = mean.vacc, alpha = .7)) + scale_color_gradient2(low="red", mid="yellow", high="blue", limits=c(0, 100),  midpoint=50, guide = FALSE) +
      scale_y_continuous(limits = c(0, 100) ) + scale_x_continuous(limits = c(0, 3.5) ) + scale_size_continuous(range=c(10, 30), guide = FALSE) +
      labs(x = "Coefficient of variation", y = paste('Incidence per 1000:',toString(start.year),'-',toString(start.year + window.length - 1)))  +
      theme(axis.text.x = element_text(colour="black"), axis.text.y = element_text(colour="black"))  + 
      geom_text(size = text.size)  + theme_bw() + scale_alpha(guide = "none")
  }
  
  
  if(scaling == 'log')
  {
    qqq1 = to.plot + geom_point(aes(size = BR, colour = mean.vacc, alpha = .7)) + scale_color_gradient2(low="red", mid="yellow", high="blue", limits=c(0, 100),  midpoint=50, guide = FALSE) +
      scale_y_continuous(limits = c(0.001, 100), trans= 'log10' ) + scale_x_continuous(limits = c(0, 3.5) ) +  scale_size_continuous(range=c(10, 30), guide = FALSE) +
      labs(x = "Coefficient of variation", y = paste('Incidence per 1000:',toString(start.year),'-',toString(start.year + window.length - 1)))  +
      theme(axis.text.x = element_text(colour="black"), axis.text.y = element_text(colour="black"))  + 
      geom_text(size = text.size) + theme_bw() + scale_alpha(guide = "none")
  }
  
  if(scaling == 'sqrt')
  {
    qqq1 = to.plot + geom_point(aes(size = BR, colour = mean.vacc, alpha = .8)) + scale_color_gradientn(colours = colorRamps::matlab.like2(100)) +
      scale_y_continuous(limits = c(0, 160), trans= 'sqrt' ) + scale_x_continuous(limits = c(0, 3.5), breaks = c(seq(from = 0, to = 3.5, by = 0.5) )) +  scale_size_continuous(range=c(10, 30), guide = FALSE) +
      labs(x = "Coefficient of variation", y = paste('Incidence per 1000:',toString(start.year),'-',toString(start.year + window.length - 1)))  +
      theme(axis.text.x = element_text(colour="black"), axis.text.y = element_text(colour="black"))  + 
      geom_text(size = text.size) + theme_bw() + scale_alpha(guide = "none") + theme(panel.grid.major = element_line(size = 1))  
  }
  
  
  
  print(qqq1  )
}






###################################
# Plots animation
###################################
coeff.vs.incidence.plot <- function(countries.to.plot, scaling, text.size, ani.h, ani.w ){
  library(igraph)
  library(animation)
  library(ggplot2)
  
  ani.record(reset=TRUE)
  start.year  =  1980
  if (countries.to.plot == 'All'){
    countries.to.plot = African.countries  
  }
  for ( i in 1 : (length(coeff.var[1, ])-1)){
    plot.coeff.var.with.inputs(African.countries, countries.to.plot, start.year+i, 
                                window.length = 10, scaling, text.size )
    ani.record()
  }
  oopts = ani.options(interval = 0.5)
  ani.replay()
  
  saveHTML(ani.replay(), ani.height = ani.h,  ani.width = ani.w)

}




###################################
# Output data for animation
###################################
output.coeff.vs.incidence.data <- function(window.length, interp.resolution ){
  list[coeff.var2, incidence.per.10002, mean.vac2, mean.br2] = output.coeff.of.var.and.other.variables(window.length)
  coeff.var = matrix(0, length(coeff.var2[, 1]), interp.resolution * ( length(coeff.var2[1, ])-1) + 1)
  incidence.per.1000 = matrix(0, length(coeff.var2[, 1]), interp.resolution * (length(coeff.var2[1, ])-1) + 1)
  mean.vac = matrix(0, length(coeff.var2[, 1]), interp.resolution * (length(coeff.var2[1, ])-1) + 1)
  mean.br = matrix(0, length(coeff.var2[, 1]), interp.resolution * (length(coeff.var2[1, ])-1) + 1)
  
  for( i in 1 : length(coeff.var2[1, ])){
    coeff.var[, ((i-1)* interp.resolution) + 1]  =   coeff.var2[, i]
    incidence.per.1000[, ((i-1)* interp.resolution) + 1]  =   incidence.per.10002[, i]
    mean.vac[, ((i-1)* interp.resolution) + 1]  =   mean.vac2[, i]
    mean.br[, ((i-1)* interp.resolution) + 1]  =   mean.br2[, i]
  }
  
  for(i in 1 : (length(coeff.var2[1, ]) - 1 )){
    for(j in 1 : (interp.resolution - 1)){
      coeff.var[, ((i-1)* interp.resolution) + 1 + j]  =  coeff.var2[, i]  + (coeff.var2[, i + 1] - coeff.var2[, i]) * j / interp.resolution
      incidence.per.1000[, ((i-1)* interp.resolution) + 1 + j]  =  incidence.per.10002[, i]  + (incidence.per.10002[, i + 1] - incidence.per.10002[, i]) * j / interp.resolution
      mean.vac[, ((i-1)* interp.resolution) + 1 + j]  =  mean.vac2[, i]  + (mean.vac2[, i + 1] - mean.vac2[, i]) * j / interp.resolution
      mean.br[, ((i-1)* interp.resolution) + 1 + j]  =  mean.br2[, i]  + (mean.br2[, i + 1] - mean.br2[, i]) * j / interp.resolution
    }
  }
  
  return( list(coeff.var, incidence.per.1000, mean.vac, mean.br))
  
}






###################################
# Output data for google bubble chart
###################################
output.data.for.animation <- function(window.length, interp.resolution){
  
  list[coeff.var, incidence.per.1000, mean.vac, mean.br]  <- output.coeff.vs.incidence.data(window.length, interp.resolution )
  
  
  count = 1
  anim.data = matrix(0, (length(coeff.var[, 1])) * length(coeff.var[1, ]), 6)
  anim.data  =  data.frame(anim.data)
  colnames(anim.data) = c("Country", "Coefficient.of.Variation", "Incidence", "Mean.vaccination", "Mean.birth.rate", "Year")
  
  for(i in 1:length(coeff.var[1, ])){
    anim.data$Country[((count - 1) * length(coeff.var[, 1]) + 1): ((count ) * length(coeff.var[, 1]))]  =  African_vaccination$Country
    anim.data$Coefficient.of.Variation[((count - 1) * length(coeff.var[, 1]) + 1): ((count ) * length(coeff.var[, 1]))]  =  round(as.numeric(coeff.var[, i]),2)
    anim.data$Incidence[((count - 1) * length(coeff.var[, 1]) + 1): ((count ) * length(coeff.var[, 1]))]  =  round(as.numeric(incidence.per.1000[, i]),2)
    anim.data$Mean.vaccination[((count - 1) * length(coeff.var[, 1]) + 1): ((count ) * length(coeff.var[, 1]))]  =  round(as.numeric(mean.vac[, i]),  2)
    anim.data$Mean.birth.rate[((count - 1) * length(coeff.var[, 1]) + 1): ((count ) * length(coeff.var[, 1]))]  =  round(as.numeric(mean.br[, i]), 2)
    anim.data$Year[((count - 1) * length(coeff.var[, 1]) + 1): ((count ) * length(coeff.var[, 1]))]  =  (1980 + (i-1) / interp.resolution)
    
    count  =  count + 1
  }
  anim.data2 = matrix(0, (length(coeff.var[, 1])) * length(coeff.var[1, ]) + 2 * length(unique(anim.data$Year)), 6)
  anim.data2  =  data.frame(anim.data2)
  colnames(anim.data2) = c("Country", "Coefficient.of.Variation", "Incidence", "Mean.vaccination", "Mean.birth.rate", "Year")
  anim.data2[1:length(anim.data[, 1]), ]  =  anim.data[]
  anim.data2[(length(anim.data[, 1]) + 1) : length(anim.data2[, 1]), 6]  =  c(unique(anim.data$Year),unique(anim.data$Year))
  anim.data2[(length(anim.data[, 1]) + 1) : length(anim.data2[, 1]), 1]  =  ""
  for(i in 1 : length(unique(anim.data$Year))){
    t  =  subset(anim.data, anim.data$Year ==  unique(anim.data$Year)[i])
    anim.data2[(length(anim.data[, 1]) + i) , 4]  =  0
    anim.data2[(length(anim.data[, 1]) + i)  , 5]  =  min(t$Mean.birth.rate)
  }
  
  for(i in 1 : length(unique(anim.data$Year))){
    t  =  subset(anim.data, anim.data$Year ==  unique(anim.data$Year)[i])
    anim.data2[(length(anim.data[, 1]) + length(unique(anim.data$Year)) +  i) , 4]  =  100
    anim.data2[(length(anim.data[, 1]) + length(unique(anim.data$Year)) +  i) , 5]  =  min(t$Mean.birth.rate)
  }
  return(anim.data2)
}





###################################
# Output coefficient of variation for a variable year periods, with 1 year moving window. 
# Can specify which regions we are interested in using this function
# Also output mean birth and vaccination rates along with cumulative incidence in those 10 years
###################################
output.coeff.of.var.and.other.variables.specify.regions <- function(window.length, regions, official.vacc){
  
  cases.by.country.by.year = read.csv("Measles_cases_by_year.csv", stringsAsFactors = FALSE)
  Birth.rates = read.csv("Birth_rates.csv", stringsAsFactors = FALSE)
  pop.by.year = read.csv("All_populations.csv", stringsAsFactors = FALSE)
  if(official.vacc == 1){
    vacc.rates = read.csv("Vacc_rates_official.csv", stringsAsFactors = FALSE)
    vacc.rates = cbind(vacc.rates[, 2], vacc.rates[, seq(36, 3, -1)], vacc.rates$WHO_REGION)
    colnames(vacc.rates) = c("Country",paste(seq(1980, 2013)), "WHO_REGION")
  } else{
    vacc.rates = read.csv("Measles_vac_all.csv", stringsAsFactors = FALSE)
  }
  
  
  subset.pop.by.year = subset(pop.by.year, pop.by.year$WHO_REGION %in% regions)
  subset.vaccination = subset(vacc.rates, vacc.rates$WHO_REGION %in% regions)
  subset.birth.rates = subset(Birth.rates, Birth.rates$WHO_REGION %in% regions)
  subset.data = subset(cases.by.country.by.year, cases.by.country.by.year$WHO_REGION %in% regions)
  
  missing1 = setdiff(subset.vaccination$Country, subset.pop.by.year$Country.Name)
  if(length(missing1) > 0){
    j = which(subset.vaccination$Country %in% missing1)
    subset.vaccination = subset.vaccination[-j, ]
  }
  missing2 = setdiff(subset.birth.rates$Country, subset.vaccination$Country)
  if(length(missing2) > 0){
    j = which(subset.birth.rates$Country %in% missing2)
    subset.birth.rates = subset.birth.rates[-j, ]
  }
  missing3 = setdiff(subset.data$Cname, subset.vaccination$Country)
  if(length(missing3) > 0){
    j = which(subset.data$Cname %in% missing3)
    subset.data = subset.data[-j, ]
  }
  
  
  p1  =  subset.pop.by.year
  p2  =  subset.vaccination
  p3  =  subset.birth.rates
  p4  =  subset.data
  
  num.windows = 25 + (10 - window.length)
  for ( i in 1 : length(subset.vaccination[, 1])){
    C  =  subset.vaccination$Country[i]
    p2[i, ]  =  subset.vaccination[i, ]
    j = which(subset.pop.by.year$Country.Name == C)
    p1[i, ]  =  subset.pop.by.year[j, ]
    j = which(subset.birth.rates$Country == C)
    p3[i, ]  =  subset.birth.rates[j, ]
    j = which(subset.data$Cname == C)
    p4[i, ]  =  subset.data[j, ]
  }
  subset.pop.by.year = p1
  subset.vaccination = p2
  subset.birth.rates = p3
  subset.data = p4
  
  coeff.var = matrix(0, length(subset.data[ , 1]), num.windows)
  for ( j in 0 : (num.windows - 1)){
    for ( i in 1 : length(subset.data[ , 1])){
      coeff.var[i, num.windows - j]  =  sd(subset.data[i, (5 + j):(5 + window.length - 1 + j )], na.rm = TRUE) /  
        mean(as.numeric(subset.data[i, (5 + j):(5 + window.length - 1 + j)]), na.rm = TRUE)
      if(is.na( mean(as.numeric(subset.data[i, (5 + j):(5 + window.length - 1 + j)]), na.rm = TRUE)) == FALSE ){
        if( mean(as.numeric(subset.data[i, (5 + j):(5 + window.length - 1 + j)]), na.rm = TRUE) == 0) {
          coeff.var[i, num.windows - j]  =  0
        } 
      }
    }
  }
  
  
  
  incidence.per.1000 = matrix(0, length(subset.data[ , 1]), num.windows)
  for ( j in 0 : (num.windows - 1)){
    for ( i in 1 : length(subset.data[ , 1])){
      incidence.per.1000[i, num.windows - j]  =  sum(as.numeric(subset.data[i,  (5 + j):(5 + window.length - 1 + j)]) / 
                                                               as.numeric(subset.pop.by.year[i, (54 - j) : (54 - (window.length - 1) - j)]), na.rm = TRUE) * 1000
    }
  }
  
  
  
  
  
  mean.vac  =  matrix(0, length(subset.vaccination[, 1]), num.windows) 
  for (j in 0 :  (num.windows - 1)){
    for ( i in 1 : length(subset.vaccination[, 1])){
      mean.vac[i, j + 1]  =  mean(as.numeric(subset.vaccination[i, (2 + j):(2 + window.length - 1  + j )]), na.rm = TRUE)
    }
  }
  
  subset.birth.rates$X2013 = subset.birth.rates$X2012
  subset.birth.rates2 = subset.birth.rates[, - (length(subset.birth.rates[1, ]) - 1)]
  subset.birth.rates2 =subset.birth.rates2 [ ,-(1:20)]
  mean.br  =  matrix(0, length(subset.birth.rates[, 1]), num.windows) 
  for (j in 0 :  (num.windows - 1)){
    for ( i in 1 : length(subset.birth.rates[, 1])){
      mean.br[i, j+1]  =  mean(as.numeric(subset.birth.rates2[i, (1 + j) : (1 + window.length - 1  + j)]), na.rm = TRUE)
    }
  }
  
  return(list(coeff.var, incidence.per.1000, mean.vac, mean.br, p1$WHO_REGION, p1$Country.Name))
}






###################################
# Output coefficient of variation for a variable year periods, with 1 year moving window. 
# Can specify which regions we are interested in using this function
# Also output mean birth and vaccination rates along with cumulative incidence in those 10 years
###################################
output.coeff.vs.incidence.data.specify.region <- function(window.length, interp.resolution, regions, official.vacc ){
  # output all data in matrix form
  list[coeff.var2, incidence.per.10002, mean.vac2, mean.br2, region.code, countries] = output.coeff.of.var.and.other.variables.specify.regions(window.length,regions, official.vacc)
  
  # the matrix form of the data is interpolated so that the animation will be smooth
  coeff.var = matrix(0, length(coeff.var2[, 1]), interp.resolution * ( length(coeff.var2[1, ])-1) + 1)
  incidence.per.1000 = matrix(0, length(coeff.var2[, 1]), interp.resolution * (length(coeff.var2[1, ])-1) + 1)
  mean.vac = matrix(0, length(coeff.var2[, 1]), interp.resolution * (length(coeff.var2[1, ])-1) + 1)
  mean.br = matrix(0, length(coeff.var2[, 1]), interp.resolution * (length(coeff.var2[1, ])-1) + 1)
  
  for( i in 1 : length(coeff.var2[1, ])){
    coeff.var[, ((i-1)* interp.resolution) + 1]  =   coeff.var2[, i]
    incidence.per.1000[, ((i-1)* interp.resolution) + 1]  =   incidence.per.10002[, i]
    mean.vac[, ((i-1)* interp.resolution) + 1]  =   mean.vac2[, i]
    mean.br[, ((i-1)* interp.resolution) + 1]  =   mean.br2[, i]
  }
  
  for(i in 1 : (length(coeff.var2[1, ]) - 1 )){
    for(j in 1 : (interp.resolution - 1)){
      coeff.var[, ((i-1)* interp.resolution) + 1 + j]  =  coeff.var2[, i]  + (coeff.var2[, i + 1] - coeff.var2[, i]) * j / interp.resolution
      incidence.per.1000[, ((i-1)* interp.resolution) + 1 + j]  =  incidence.per.10002[, i]  + (incidence.per.10002[, i + 1] - incidence.per.10002[, i]) * j / interp.resolution
      mean.vac[, ((i-1)* interp.resolution) + 1 + j]  =  mean.vac2[, i]  + (mean.vac2[, i + 1] - mean.vac2[, i]) * j / interp.resolution
      mean.br[, ((i-1)* interp.resolution) + 1 + j]  =  mean.br2[, i]  + (mean.br2[, i + 1] - mean.br2[, i]) * j / interp.resolution
    }
  }
  
  return( list(coeff.var, incidence.per.1000, mean.vac, mean.br, region.code, countries))
  
}







###################################
# Output data for google bubble chart where we can specify which regions we are interested in seeing
# Setting official.vac to 1 uses official WHO vaccination rates, whilst any other value uses WHO estimates
###################################
output.data.for.animation.specify.region <- function(window.length, interp.resolution, regions, official.vacc){
  
  list[coeff.var, incidence.per.1000, mean.vac, mean.br, region.code, countries]  <- output.coeff.vs.incidence.data.specify.region(window.length, interp.resolution, regions, official.vacc)
  
  
  count = 1
  anim.data = matrix(0, (length(coeff.var[, 1])) * length(coeff.var[1, ]), 7)
  anim.data  =  data.frame(anim.data)
  colnames(anim.data) = c("Country", "Coefficient.of.Variation", "Incidence", "Mean.vaccination", "Mean.birth.rate", "Year", "WHO_REGION")
  
  for(i in 1:length(coeff.var[1, ])){
    anim.data$Country[((count - 1) * length(coeff.var[, 1]) + 1): ((count ) * length(coeff.var[, 1]))]  =  countries
    anim.data$Coefficient.of.Variation[((count - 1) * length(coeff.var[, 1]) + 1): ((count ) * length(coeff.var[, 1]))]  =  round(as.numeric(coeff.var[, i]),2)
    anim.data$Incidence[((count - 1) * length(coeff.var[, 1]) + 1): ((count ) * length(coeff.var[, 1]))]  =  round(as.numeric(incidence.per.1000[, i]),2)
    anim.data$Mean.vaccination[((count - 1) * length(coeff.var[, 1]) + 1): ((count ) * length(coeff.var[, 1]))]  =  round(as.numeric(mean.vac[, i]),  2)
    anim.data$Mean.birth.rate[((count - 1) * length(coeff.var[, 1]) + 1): ((count ) * length(coeff.var[, 1]))]  =  round(as.numeric(mean.br[, i]), 2)
    anim.data$Year[((count - 1) * length(coeff.var[, 1]) + 1): ((count ) * length(coeff.var[, 1]))]  =  (1980 + (i-1) / interp.resolution)
    anim.data$WHO_REGION[((count - 1) * length(coeff.var[, 1]) + 1): ((count ) * length(coeff.var[, 1]))]  =  region.code
    count  =  count + 1
  }
  
  anim.data2 = matrix(0, (length(coeff.var[, 1])) * length(coeff.var[1, ]) + (2 * length(regions) * length(unique(anim.data$Year))), 7)
  anim.data2  =  data.frame(anim.data2)
  colnames(anim.data2) = c("Country", "Coefficient.of.Variation", "Incidence", "Mean.vaccination", "Mean.birth.rate", "Year", "WHO_REGION")
  anim.data2[1:length(anim.data[, 1]), ]  =  anim.data
  
  year.mins = matrix(0, length(unique(anim.data$Year)), 2)
  
  for(i in 1 : length(unique(anim.data$Year))){
    t  =  subset(anim.data, anim.data$Year ==  unique(anim.data$Year)[i])
    year.mins[i, 1] = unique(anim.data$Year)[i]
    year.mins[i, 2] = min(t$Mean.birth.rate,na.rm = T)  
  }
  
  l =  expand.grid("", -1, 0, c(0,100), 0,unique(anim.data$Year), regions)
  for(i in 1 : (2 * length(regions) * length(unique(anim.data$Year)))){
    y = l[i, 6]
    j = which(year.mins[, 1] == y)
    l[i, 5]  =  year.mins[j, 2]
  }
  
  anim.data2[(length(anim.data[, 1]) + 1 ): length(anim.data2[, 1]), ]  =  l
  for( i in 1 : (2 * length(regions) * length(unique(anim.data$Year)))){
    anim.data2[length(anim.data[, 1]) + i, 7] = regions[as.numeric(anim.data2[length(anim.data[, 1]) + i,7])]
    anim.data2[length(anim.data[, 1]) + i, 1] = ""
  }
  return(anim.data2)
}





###################################
# Make overall susceptibility of the population a given value
###################################
new.susceptibility.level <- function(disease.state, sus.level, susceptible.indices, recovered.indices){
  a = sum(disease.state[susceptible.indices])/ sum (disease.state)  
  disease.state[susceptible.indices]  =  disease.state[susceptible.indices] - (1- sus.level / a) * disease.state[susceptible.indices]
  disease.state[recovered.indices]  =  disease.state[recovered.indices] + (1-sus.level / a) * disease.state[susceptible.indices]
  return(disease.state)
}



###################################
# For specified susceptible proportions by age, create a disease state that matches these proportions
###################################
construct.specified.disease.state <- function( num.comps, sus.props.by.age){
  disease.state = initial.disease.state(demographic.ages  ,   1  ,  num.comps,  list.of.ages, mat.immunity.loss)
  if (length(sus.props.by.age) == 98){
    j  =  which(list.of.ages == 1)
    for ( k in 1 : (j-1)){
      num.people = sum(disease.state[c(susceptible.indices[k], recovered.indices[k])])
      disease.state[susceptible.indices[k]] = sus.props.by.age[1] * num.people
      disease.state[recovered.indices[k]] = ( 1 - sus.props.by.age[1]) * num.people
    }
    for (k in j : length(list.of.ages)){
      num.people = sum(disease.state[c(susceptible.indices[k], recovered.indices[k])])
      disease.state[susceptible.indices[k]] = sus.props.by.age[k - (j-2)] * num.people
      disease.state[recovered.indices[k]] = ( 1 - sus.props.by.age[k - (j-2)]) * num.people
    }
  } else{
    for ( k in 1 : length(sus.props.by.age)){
      num.people = sum(disease.state[c(susceptible.indices[k], recovered.indices[k])])
      disease.state[susceptible.indices[k]] = sus.props.by.age[k] * num.people
      disease.state[recovered.indices[k]] = ( 1 - sus.props.by.age[k]) * num.people
    }
  }
  return(disease.state)
}





###################################
# For initial disease state, calculate the expected coefficient of variation and 
# incidence over a given time period
###################################
model.of.coeff.var.incidence <- function(num.repeats, disease.state.input, initial.sus.level, 
                                         input.sus.props.by.age, sus.props.by.age,
                                         vacc.by.year, br.by.year, initial.pop.size, 
                                         vacc.success, av.migrants.per.year, sia.years,
                                         sia.numbers, start.time.for.spread,
                                         unreached.population, vacc.efficieny,
                                         new.figure, colour){
  
  if (input.sus.props.by.age == 1){
    disease.state.input = construct.specified.disease.state( num.comps, sus.props.by.age)
  }
  # change the overall susceptibility level to a specified proportion of the population
  count = 1
  while(((sum(disease.state.input[susceptible.indices]) / sum(disease.state.input) > (initial.sus.level + 0.002)) | 
          (sum(disease.state.input[susceptible.indices]) / sum(disease.state.input) < (initial.sus.level - 0.002))) &
          count < 10){
    disease.state.input  =  new.susceptibility.level(disease.state.input, sus.level = initial.sus.level, susceptible.indices, recovered.indices )
    count = count + 1
  }
  
  
  a = numerous.sims.diff.br.vacc(num.repeats,  disease.state.input, initial.start.time = start.time.for.spread, multi.vacc = vacc.by.year, 
                                 multi.br = br.by.year, input.infections.all.times = 1, init.pop.size = initial.pop.size, 
                                 sia.years = sia.years, sia.numbers = sia.numbers, vacc.success = vacc.success, 
                                 num.migrants.per.year = av.migrants.per.year,
                                 unreached.pop = unreached.population, vacc.efficieny)
  
  c.var = matrix(0, num.repeats, 1)
  cumulative.inc = matrix(0, num.repeats, 1)
  total.cases = matrix(0, num.repeats, 1)
  for (i in 1 : num.repeats){
    c.var[i] =  sd(unlist(a[i]))/mean(unlist(a[i]))
    cumulative.inc[i]  =  sum(unlist(a[[i]][[2]]))
    total.cases[i]   =  sum(unlist(a[[i]][[1]]))
  }
  if(new.figure == 1){
    plot(mean(c.var), mean(cumulative.inc), col = colour, pch = 16, cex = 3)
  }else{
    points(mean(c.var), mean(cumulative.inc), col = colour, pch = 16, cex = 3)
  }
  
  return(list(a, c.var, cumulative.inc, total.cases))
}




###################################
# Output coefficient of variation for a variable year periods, with 1 year moving window. 
# Can specify which regions we are interested in using this function
# Also output mean birth and vaccination rates along with cumulative incidence in those 10 years
###################################
output.coeff.of.var.and.other.variables.specify.regions.gaussian.weighting <- function(window.length, regions, gaussian.st.dev, cutoff = 10){
   
  cases.by.country.by.year = read.csv("Measles_cases_by_year.csv", stringsAsFactors = FALSE)
  cases.by.country.by.year = read.csv("Measles_cases_by_year2.csv", stringsAsFactors = FALSE)
  Birth.rates = read.csv("Birth_rates.csv", , stringsAsFactors = FALSE)
  pop.by.year = read.csv("All_populations.csv", stringsAsFactors = FALSE)
  vacc.rates = read.csv("Measles_vac_all.csv", stringsAsFactors = FALSE)
  
  subset.pop.by.year = subset(pop.by.year, pop.by.year$WHO_REGION %in% regions)
  subset.vaccination = subset(vacc.rates, vacc.rates$WHO_REGION %in% regions)
  subset.birth.rates = subset(Birth.rates, Birth.rates$WHO_REGION %in% regions)
  subset.data = subset(cases.by.country.by.year, cases.by.country.by.year$WHO_REGION %in% regions)
  
  missing1 = setdiff(subset.vaccination$Country, subset.pop.by.year$Country.Name)
  if(length(missing1) > 0){
    j = which(subset.vaccination$Country %in% missing1)
    subset.vaccination = subset.vaccination[-j, ]
  }
  missing2 = setdiff(subset.birth.rates$Country, subset.vaccination$Country)
  if(length(missing2) > 0){
    j = which(subset.birth.rates$Country %in% missing2)
    subset.birth.rates = subset.birth.rates[-j, ]
  }
  missing3 = setdiff(subset.data$Cname, subset.vaccination$Country)
  if(length(missing3) > 0){
    j = which(subset.data$Cname %in% missing3)
    subset.data = subset.data[-j, ]
  }
  
  
  p1  =  subset.pop.by.year
  p2  =  subset.vaccination
  p3  =  subset.birth.rates
  p4  =  subset.data
  
  num.windows = 25 + (10 - window.length)
  for ( i in 1 : length(subset.vaccination[, 1])){
    C  =  subset.vaccination$Country[i]
    p2[i, ]  =  subset.vaccination[i, ]
    j = which(subset.pop.by.year$Country.Name == C)
    p1[i, ]  =  subset.pop.by.year[j, ]
    j = which(subset.birth.rates$Country == C)
    p3[i, ]  =  subset.birth.rates[j, ]
    j = which(subset.data$Cname == C)
    p4[i, ]  =  subset.data[j, ]
  }
  subset.pop.by.year = p1
  subset.vaccination = p2
  subset.birth.rates = p3
  subset.data = p4
  
  require(stats)
  x = seq(1980, 2013)
  xout = seq(1980, 2013, 1)
  list[interp.subset.data, interp.subset.vacc, interp.subset.br, interp.subset.pop] = interp.datasets(subset.data, 
                                                                                                      subset.vaccination, 
                                                                                                      subset.birth.rates, 
                                                                                                      subset.pop.by.year,
                                                                                                      x,
                                                                                                      xout)
  mean.cases = matrix(0, length(interp.subset.data[, 1]), length(interp.subset.data[1, ]))
  mean.cases[, 1] = subset.data$Cname
  mean.cases[, 2] = subset.data$WHO_REGION
  
  coeff.var.cases = matrix(0, length(interp.subset.data[, 1]), length(interp.subset.data[1, ]))
  coeff.var.cases[, 1] = subset.data$Cname
  coeff.var.cases[, 2] = subset.data$WHO_REGION
  
  incidence.per.1000 = matrix(0, length(interp.subset.data[, 1]), length(interp.subset.data[1, ]))
  incidence.per.1000[, 1] = subset.data$Cname
  incidence.per.1000[, 2] = subset.data$WHO_REGION
  
  mean.br = matrix(0, length(interp.subset.data[, 1]), length(interp.subset.data[1, ]))
  mean.br[, 1] = subset.data$Cname
  mean.br[, 2] = subset.data$WHO_REGION
  
  mean.vac = matrix(0, length(interp.subset.data[, 1]), length(interp.subset.data[1, ]))
  mean.vac[, 1] = subset.data$Cname
  mean.vac[, 2] = subset.data$WHO_REGION
  
  for(i in 3: length(interp.subset.data[1, ])){
    
    w.input = xout - xout[i-2]
    w = output.weights.gaussian.with.cutoff(w.input, gaussian.st.dev, cutoff)
    # w is the weights that we use to calculate incidence and coefficient of variation over time
    #w = dnorm(w.input, mean = 0, sd = gaussian.st.dev)
    for (j in 1 : length(interp.subset.data[, 1])){  
      
      mean.cases[j, i] = weighted.mean(as.numeric(interp.subset.data[j ,-(1:2)]), w)
      if(mean.cases[j, i] == 0){
        coeff.var.cases[j, i] = 0
      } else{
      coeff.var.cases[j, i] = sqrt(weighted.var(as.numeric(interp.subset.data[j ,-(1:2)]), w)) / as.numeric(mean.cases[j,i])
      }
      incidence.per.1000[j, i]  =  1000 * as.numeric(mean.cases[j, i])  /  (as.numeric(sum(as.numeric(interp.subset.pop[j ,-(1:2)], na.rm = T) * w) / sum(w)))
      mean.br[j, i] = as.numeric(sum(as.numeric(interp.subset.br[j ,-(1:2)]) * w , na.rm = T) / sum(w))
      mean.vac[j, i] = as.numeric(sum(as.numeric(interp.subset.vacc[j ,-(1:2)]) * w , na.rm = T) / sum(w))   
      
    }
  }
  xout = seq(1980, 2013, 1/interp.resolution)
  mean.cases = cbind(mean.cases[,(1:2)], interpolate.give.dataset(mean.cases[, -(1:2)], x, xout))
  coeff.var.cases = cbind(coeff.var.cases[,(1:2)], interpolate.give.dataset(coeff.var.cases[, -(1:2)], x, xout))
  incidence.per.1000 = cbind(incidence.per.1000[,(1:2)], interpolate.give.dataset(incidence.per.1000[, -(1:2)], x, xout))
  mean.br = cbind(mean.br[,(1:2)], interpolate.give.dataset(mean.br[, -(1:2)], x, xout))
  mean.vac = cbind(mean.vac[,(1:2)], interpolate.give.dataset(mean.vac[, -(1:2)], x, xout))
  
  mean.vac[, -(1:2)] = round(as.numeric(mean.vac[, -(1:2)]), 2)
  mean.cases[, -(1:2)] = round(as.numeric(mean.cases[, -(1:2)]), 2)
  mean.br[, -(1:2)] = round(as.numeric(mean.br[, -(1:2)]), 2)
  incidence.per.1000[, -(1:2)] = round(as.numeric(incidence.per.1000[, -(1:2)]), 2)
  coeff.var.cases[, -(1:2)] = round(as.numeric(coeff.var.cases[, -(1:2)]), 2)
  

  output.data = matrix(0, (length(coeff.var.cases[, 1])) * length(coeff.var.cases[1, -(1:2)]) + (2 * length(regions) * length(xout)), 7)
  output.data  =  data.frame(output.data)
  colnames(output.data) = c("Country", "Coefficient.of.Variation", "Incidence", "Mean.vaccination", "Mean.birth.rate", "Year", "WHO_REGION")
  
  output.data$Country[seq(1, (length(coeff.var.cases[, 1])) * length(coeff.var.cases[1, -(1:2)]))] = rep(mean.cases[, 1], length(coeff.var.cases[1, -(1:2)]))
  output.data$WHO_REGION[seq(1, (length(coeff.var.cases[, 1])) * length(coeff.var.cases[1, -(1:2)]))] = rep(mean.cases[, 2], length(coeff.var.cases[1, -(1:2)]))
  count = 1
  for(i in 3 : length(coeff.var.cases[1, ])){
    output.data$Coefficient.of.Variation[seq(((count - 1) * (length(coeff.var.cases[, 1]))) + 1, count * (length(coeff.var.cases[, 1])))] = coeff.var.cases[, i]
    output.data$Incidence[seq(((count - 1) * (length(coeff.var.cases[, 1]))) + 1, count * (length(coeff.var.cases[, 1])))] = incidence.per.1000[, i]
    output.data$Mean.vaccination[seq(((count - 1) * (length(coeff.var.cases[, 1]))) + 1, count * (length(coeff.var.cases[, 1])))] = mean.vac[, i]
    output.data$Mean.birth.rate[seq(((count - 1) * (length(coeff.var.cases[, 1]))) + 1, count * (length(coeff.var.cases[, 1])))] = mean.br[, i]
    output.data$Year[seq(((count - 1) * (length(coeff.var.cases[, 1]))) + 1, count * (length(coeff.var.cases[, 1])))] = xout[i - 2]  
    count = count + 1
  }


  year.mins = matrix(0, length(xout), 2)

  for(i in 1 : length(xout)){
    t  =  subset(output.data, output.data$Year ==  unique(xout)[i])
    year.mins[i, 1] = xout[i]
    year.mins[i, 2] = as.numeric(min(t$Mean.birth.rate,na.rm = T)  )
  }

  l =  expand.grid("", -1, 0, c(0,100), 0, xout, regions)
  
  for(i in 1 : (2 * length(regions) * length(xout))){
    y = l[i, 6]
    j = which(year.mins[, 1] == y)
    l[i, 5]  =  year.mins[j, 2]
  }

  output.data[((length(coeff.var.cases[, 1])) * length(coeff.var.cases[1, -(1:2)]) + 1 ): length(output.data[, 1]), ]  =  l
  
  for( i in 1 : (2 * length(regions) * length(xout))){
    output.data[ (length(coeff.var.cases[, 1])) * length(coeff.var.cases[1, -(1:2)]) + i, 7] = regions[as.numeric(output.data[ (length(coeff.var.cases[, 1])) * length(coeff.var.cases[1, -(1:2)]) + i, 7] )]
    output.data[ (length(coeff.var.cases[, 1])) * length(coeff.var.cases[1, -(1:2)])  + i, 1] = ""
  }
  output.data$Coefficient.of.Variation = as.numeric(output.data$Coefficient.of.Variation)
  output.data$Incidence  =  as.numeric(output.data$Incidence)
  output.data$Mean.vaccination  =  as.numeric(output.data$Mean.vaccination)
  output.data$Mean.birth.rate  =  as.numeric(output.data$Mean.birth.rate)
  output.data$Year   =  as.numeric(output.data$Year)
  output.data$Coefficient.of.Variation[which(output.data$Coefficient.of.Variation == "Inf")] = 0
  return(output.data)
}





#####################

interpolate.give.dataset <- function(data,  
                                     x,
                                     xout){
  interp.data = matrix(0, length(data[, 1]), length(xout))
  for ( i in 1 : length(data[, 1])){
    y = data[i, ]
    if(length(which(!is.na(y)) == F) < 2){
      interp.data[i, ]  =  0} else{
        list[qq,ww] =  approx (as.numeric(x), as.numeric(y),  method = "linear", xout )
        interp.data[i, ]  =  round(ww, 2)
      }
  }
  
  return(interp.data)
  
}





###############
# Calculate weighted mean of data, with input weights
###############
weighted.mean <- function(x, w, na.rm = T) {
  if (na.rm) {
    w <- w[i <- !is.na(x)]
    x <- x[i]
  }
  sum.w <- sum(w, na.rm = na.rm)
  mean.w <- sum(x * w, na.rm = na.rm) / sum(w)
  return(mean.w)
}




###############
# Calculate weighted variance of data, with input weights
###############
weighted.var <- function(x, w, na.rm = T) {
  if (na.rm) {
    w <- w[i <- !is.na(x)]
    x <- x[i]
  }
  #sum.w <- sum(w, na.rm = na.rm)
  #sum.w2 <- sum(w^2, na.rm = na.rm)
  mean.w <- sum(x * w, na.rm = na.rm) / sum(w)
  #(sum.w / (sum.w^2 - sum.w2)) * sum(w * (x - mean.w)^2, na.rm =
   #                                    na.rm)
  sum(w * (x - mean.w) ^2 , na.rm = na.rm) / sum(w, na.rm = na.rm)
}





#################
# output the gaussian weight for the averaging process
#################
output.weights.gaussian.with.cutoff <- function(x, st.dev, cutoff, neg.only = F){
  
  weights = dnorm(x, mean = 0, sd = st.dev)
  if(neg.only == T){
    weights[which(x > cutoff)] = 0
  } else{
    weights[which(abs(x) > cutoff)] = 0
  }
  return(weights)
}






###########################
# separate output of simulation informatively
###########################
decompose.simulation.output <- function(x, y, z, window.length){
  
  num.windows = y - window.length + 1
  cases.by.year = matrix(0, x, y)
  incidence.by.year = matrix(0, x, y)
  c.var = matrix(0, x, num.windows)
  cumulative.inc = matrix(0, x, num.windows)
  total.cases = matrix(0, x, 1)
  for (i in 1 : x){ 
    cases.by.year[i, ] = unlist(z[[i]][[1]])
    total.cases[i]   =  sum(unlist(z[[i]][[1]]))
    incidence.by.year[i, ] = unlist(z[[i]][[2]])
  }
  
  for(i in 1 : x){
    for(j in 1 : num.windows){
      c.var[i,j] =  sd(unlist(z[[i]][[2]])[seq(j, j + window.length - 1) ])/mean(unlist(z[[i]][[2]])[seq(j, j + window.length - 1)])
      cumulative.inc[i,j]  =  sum(unlist(z[[i]][[2]])[seq(j, j + window.length - 1)])
    }
  }
  
  
  return (list(cases.by.year, c.var, cumulative.inc, total.cases, incidence.by.year))
}  





###########################
# output plots for paper, showing the change in the average incidence and coefficient of variation over time for a chosen region and time period
###########################
plots.for.paper <- function(data = anim.data, year = 1990, region = 'AFR', previously.plotted.year = 1990, scaling = 'sqrt', arrow.color = 'red', text.size = 5, region.text){
  library(grid)
  max.incidence = max(data$Incidence[which(data$Year == 1980 & data$WHO_REGION == region)], na.rm = T)
  a = subset(data, data$WHO_REGION == region & anim.data$Year == year & anim.data$Country != "")
  b = data.frame(labels1 = a$Country, Coeff.of.var =  a$Coefficient.of.Variation, Incidence =  a$Incidence, BR =  a$Mean.birth.rate, Mean_vaccination =  a$Mean.vaccination, Inverse.BR = 1/a$Mean.birth.rate)
  to.plot = ggplot(b, aes(x = Coeff.of.var, y = Incidence, label = labels1)) 
  
  if(scaling == 'sqrt')
  {
    qqq1 = to.plot + geom_point(aes(size = BR, colour = 100 - Mean_vaccination, alpha = .8)) + scale_color_gradientn(limits = c(0, 100),colours = colorRamps::matlab.like2(256)) + 
      scale_y_continuous(limits = c(0, round(max.incidence + 2)), trans= 'sqrt' ) + scale_x_continuous(limits = c(0, 5), breaks = c(seq(from = 0, to = 5, by = 0.5) )) +
      scale_size_continuous(range=c(10, 30), guide = FALSE) +
      labs(x = "Coefficient of variation", y = paste("Incidence in", region.text, toString(year)), color = "Non vaccinated" )  +
      theme(axis.text.x = element_text(colour="black"), axis.text.y = element_text(colour="black"))  + 
      geom_text(size = text.size) + theme_bw() + scale_alpha(guide = "none") +
      theme(panel.grid.major = element_line(size = 1), axis.text=element_text(size=24), axis.title=element_text(size=28)  )
  }
  if(scaling == 'none')
  {
    qqq1 = to.plot + geom_point(aes(size = BR, colour = 100 - Mean_vaccination, alpha = .8)) + scale_color_gradientn(limits = c(0, 100),colours = colorRamps::matlab.like2(256)) + 
      scale_y_continuous(limits = c(0, round(max.incidence + 2)) ) + scale_x_continuous(limits = c(0, 5), breaks = c(seq(from = 0, to = 5, by = 0.5) )) +
      scale_size_continuous(range=c(10, 30), guide = FALSE) +
      labs(x = "Coefficient of variation", y = paste("Incidence in", region.text, toString(year)), color = "Non vaccinated" )  +
      theme(axis.text.x = element_text(colour="black"), axis.text.y = element_text(colour="black"))  + 
      geom_text(size = text.size) + theme_bw() + scale_alpha(guide = "none") +
      theme(panel.grid.major = element_line(size = 1), axis.text=element_text(size=14), axis.title=element_text(size=18)  )
  }
  
  
  
  Mean_X = matrix(0, year - previously.plotted.year + 1, 1)
  Mean_Y = matrix(0, year - previously.plotted.year + 1, 1)
  for ( i in previously.plotted.year : year){
    Mean_X[i - previously.plotted.year + 1] = mean(anim.data$Coefficient.of.Variation[which(anim.data$Year == i & anim.data$WHO_REGION == region & anim.data$Country != "")], na.rm = T)
    Mean_Y[i - previously.plotted.year + 1] = mean(anim.data$Incidence[which(anim.data$Year == i & anim.data$WHO_REGION == region & anim.data$Country != "")], na.rm = T)
  }
  df <- data.frame(x = Mean_X, y  = Mean_Y)
  if(year != 1980){
    qqq1 + geom_path(data = df, aes(x, y, label = NULL ), size = 3, arrow = arrow(), color = arrow.color)
  }else{
    qqq1
  }
}








###########################
# apply gaussian weights to case data to average over a given window
###########################
gaussian.averaging <- function(y, xout, st.dev, cutoff){
  
  coeff.var.cases = y
  incidence.per.1000 = y
  for(i in 1: length(y[1, ])){
    
    w.input = xout - xout[i]
    w = output.weights.gaussian.with.cutoff(w.input, st.dev, cutoff)
    # w is the weights that we use to calculate incidence and coefficient of variation over time
    #w = dnorm(w.input, mean = 0, sd = gaussian.st.dev)
    for (j in 1 : length(y[, 1])){  
      
      mean.cases[j, i] = weighted.mean(as.numeric(y[j ,]), w)
      if(mean.cases[j, i] == 0){
        coeff.var.cases[j, i] = 0
      } else{
        coeff.var.cases[j, i] = sqrt(weighted.var(as.numeric(y[j ,]), w)) / as.numeric(mean.cases[j,i])
      }
      incidence.per.1000[j, i]  =  sum(y[j, ] * w, na.rm = T) 
    }
  }
  return(list(coeff.var.cases, incidence.per.1000))
}




###################################
# Output coefficient of variation for a variable year periods, with 1 year moving window. 
# Can specify which regions we are interested in using this function
# Also output mean birth and vaccination rates along with cumulative incidence in those 10 years. 
# Do by selecting case numbers according to a given distribution and then taking sample mean and variance
###################################
output.coeff.of.var.and.other.variables.specify.regions.gaussian.weighting.random.selection <- function(window.length, regions, 
                                                                                                        gaussian.st.dev, 
                                                                                                        cutoff = 10,
                                                                                                        rand.length = 1000,
                                                                                                        neg.averaging.only= F){
  
  cases.by.country.by.year = read.csv("Measles_cases_by_year.csv", stringsAsFactors = FALSE)
  cases.by.country.by.year = read.csv("Measles_cases_by_year2.csv", stringsAsFactors = FALSE)
  Birth.rates = read.csv("Birth_rates.csv", stringsAsFactors = FALSE)
  pop.by.year = read.csv("All_populations.csv", stringsAsFactors = FALSE)
  vacc.rates = read.csv("Measles_vac_all.csv", stringsAsFactors = FALSE)
  
  subset.pop.by.year = subset(pop.by.year, pop.by.year$WHO_REGION %in% regions)
  subset.vaccination = subset(vacc.rates, vacc.rates$WHO_REGION %in% regions)
  subset.birth.rates = subset(Birth.rates, Birth.rates$WHO_REGION %in% regions)
  subset.data = subset(cases.by.country.by.year, cases.by.country.by.year$WHO_REGION %in% regions)
  
  missing1 = setdiff(subset.vaccination$Country, subset.pop.by.year$Country.Name)
  if(length(missing1) > 0){
    j = which(subset.vaccination$Country %in% missing1)
    subset.vaccination = subset.vaccination[-j, ]
  }
  missing2 = setdiff(subset.birth.rates$Country, subset.vaccination$Country)
  if(length(missing2) > 0){
    j = which(subset.birth.rates$Country %in% missing2)
    subset.birth.rates = subset.birth.rates[-j, ]
  }
  missing3 = setdiff(subset.data$Cname, subset.vaccination$Country)
  if(length(missing3) > 0){
    j = which(subset.data$Cname %in% missing3)
    subset.data = subset.data[-j, ]
  }
  
  
  p1  =  subset.pop.by.year
  p2  =  subset.vaccination
  p3  =  subset.birth.rates
  p4  =  subset.data
  
  num.windows = 25 + (10 - window.length)
  for ( i in 1 : length(subset.vaccination[, 1])){
    C  =  subset.vaccination$Country[i]
    p2[i, ]  =  subset.vaccination[i, ]
    j = which(subset.pop.by.year$Country.Name == C)
    p1[i, ]  =  subset.pop.by.year[j, ]
    j = which(subset.birth.rates$Country == C)
    p3[i, ]  =  subset.birth.rates[j, ]
    j = which(subset.data$Cname == C)
    p4[i, ]  =  subset.data[j, ]
  }
  subset.pop.by.year = p1
  subset.vaccination = p2
  subset.birth.rates = p3
  subset.data = p4
  
  require(stats)
  x = seq(1980, 2013)
  xout = seq(1980, 2013, 1)
  list[interp.subset.data, interp.subset.vacc, interp.subset.br, interp.subset.pop] = interp.datasets(subset.data, 
                                                                                                      subset.vaccination, 
                                                                                                      subset.birth.rates, 
                                                                                                      subset.pop.by.year,
                                                                                                      x,
                                                                                                      xout)
  mean.cases = matrix(0, length(interp.subset.data[, 1]), length(interp.subset.data[1, ]))
  mean.cases[, 1] = subset.data$Cname
  mean.cases[, 2] = subset.data$WHO_REGION
  
  coeff.var.cases = matrix(0, length(interp.subset.data[, 1]), length(interp.subset.data[1, ]))
  coeff.var.cases[, 1] = subset.data$Cname
  coeff.var.cases[, 2] = subset.data$WHO_REGION
  
  incidence.per.1000 = matrix(0, length(interp.subset.data[, 1]), length(interp.subset.data[1, ]))
  incidence.per.1000[, 1] = subset.data$Cname
  incidence.per.1000[, 2] = subset.data$WHO_REGION
  
  mean.br = matrix(0, length(interp.subset.data[, 1]), length(interp.subset.data[1, ]))
  mean.br[, 1] = subset.data$Cname
  mean.br[, 2] = subset.data$WHO_REGION
  
  mean.vac = matrix(0, length(interp.subset.data[, 1]), length(interp.subset.data[1, ]))
  mean.vac[, 1] = subset.data$Cname
  mean.vac[, 2] = subset.data$WHO_REGION
  
  for(i in 3: length(interp.subset.data[1, ])){
    
    w.input = xout - xout[i-2]
    w = output.weights.gaussian.with.cutoff(w.input, gaussian.st.dev, cutoff, neg.averaging.only)
    # w is the weights that we use to calculate incidence and coefficient of variation over time
    #w = dnorm(w.input, mean = 0, sd = gaussian.st.dev)
    for (j in 1 : length(interp.subset.data[, 1])){  
      k = runif(rand.length)
      rand.cases = matrix(0, rand.length, 1)
      b = cumsum(w) / sum(w)
      for(q in 1 : rand.length){
        rand.cases[q] = as.numeric(interp.subset.data[j, which(b > k[q])[1] + 2])
      }
      mean.cases[j, i] = mean(rand.cases, na.rm = T)
      if(mean.cases[j, i] == 0){
        coeff.var.cases[j, i] = 0
      } else{
        coeff.var.cases[j, i] = as.numeric(sqrt(var(rand.cases, na.rm = T)) / as.numeric(mean.cases[j,i]))
      }
      incidence.per.1000[j, i]  =  1000 * as.numeric(mean.cases[j, i])  /  (as.numeric(sum(as.numeric(interp.subset.pop[j ,-(1:2)], na.rm = T) * w) / sum(w)))
      mean.br[j, i] = as.numeric(sum(as.numeric(interp.subset.br[j ,-(1:2)]) * w , na.rm = T) / sum(w))
      mean.vac[j, i] = as.numeric(sum(as.numeric(interp.subset.vacc[j ,-(1:2)]) * w , na.rm = T) / sum(w))   
      
    }
  }
  xout = seq(1980, 2013, 1/interp.resolution)
  mean.cases = cbind(mean.cases[,(1:2)], interpolate.give.dataset(mean.cases[, -(1:2)], x, xout))
  coeff.var.cases = cbind(coeff.var.cases[,(1:2)], interpolate.give.dataset(coeff.var.cases[, -(1:2)], x, xout))
  incidence.per.1000 = cbind(incidence.per.1000[,(1:2)], interpolate.give.dataset(incidence.per.1000[, -(1:2)], x, xout))
  mean.br = cbind(mean.br[,(1:2)], interpolate.give.dataset(mean.br[, -(1:2)], x, xout))
  mean.vac = cbind(mean.vac[,(1:2)], interpolate.give.dataset(mean.vac[, -(1:2)], x, xout))
  
  mean.vac[, -(1:2)] = round(as.numeric(mean.vac[, -(1:2)]), 2)
  mean.cases[, -(1:2)] = round(as.numeric(mean.cases[, -(1:2)]), 2)
  mean.br[, -(1:2)] = round(as.numeric(mean.br[, -(1:2)]), 2)
  incidence.per.1000[, -(1:2)] = round(as.numeric(incidence.per.1000[, -(1:2)]), 2)
  coeff.var.cases[, -(1:2)] = round(as.numeric(coeff.var.cases[, -(1:2)]), 2)
  
  
  output.data = matrix(0, (length(coeff.var.cases[, 1])) * length(coeff.var.cases[1, -(1:2)]) + (2 * length(regions) * length(xout)), 7)
  output.data  =  data.frame(output.data)
  colnames(output.data) = c("Country", "Coefficient.of.Variation", "Incidence", "Mean.vaccination", "Mean.birth.rate", "Year", "WHO_REGION")
  
  output.data$Country[seq(1, (length(coeff.var.cases[, 1])) * length(coeff.var.cases[1, -(1:2)]))] = rep(mean.cases[, 1], length(coeff.var.cases[1, -(1:2)]))
  output.data$WHO_REGION[seq(1, (length(coeff.var.cases[, 1])) * length(coeff.var.cases[1, -(1:2)]))] = rep(mean.cases[, 2], length(coeff.var.cases[1, -(1:2)]))
  count = 1
  for(i in 3 : length(coeff.var.cases[1, ])){
    output.data$Coefficient.of.Variation[seq(((count - 1) * (length(coeff.var.cases[, 1]))) + 1, count * (length(coeff.var.cases[, 1])))] = coeff.var.cases[, i]
    output.data$Incidence[seq(((count - 1) * (length(coeff.var.cases[, 1]))) + 1, count * (length(coeff.var.cases[, 1])))] = incidence.per.1000[, i]
    output.data$Mean.vaccination[seq(((count - 1) * (length(coeff.var.cases[, 1]))) + 1, count * (length(coeff.var.cases[, 1])))] = mean.vac[, i]
    output.data$Mean.birth.rate[seq(((count - 1) * (length(coeff.var.cases[, 1]))) + 1, count * (length(coeff.var.cases[, 1])))] = mean.br[, i]
    output.data$Year[seq(((count - 1) * (length(coeff.var.cases[, 1]))) + 1, count * (length(coeff.var.cases[, 1])))] = xout[i - 2]  
    count = count + 1
  }
  
  
  year.mins = matrix(0, length(xout), 2)
  
  for(i in 1 : length(xout)){
    t  =  subset(output.data, output.data$Year ==  unique(xout)[i])
    year.mins[i, 1] = xout[i]
    year.mins[i, 2] = as.numeric(min(t$Mean.birth.rate,na.rm = T)  )
  }
  
  l =  expand.grid("", -1, 0, c(0,100), 0, xout, regions)
  
  for(i in 1 : (2 * length(regions) * length(xout))){
    y = l[i, 6]
    j = which(year.mins[, 1] == y)
    l[i, 5]  =  year.mins[j, 2]
  }
  
  output.data[((length(coeff.var.cases[, 1])) * length(coeff.var.cases[1, -(1:2)]) + 1 ): length(output.data[, 1]), ]  =  l
  
  for( i in 1 : (2 * length(regions) * length(xout))){
    output.data[ (length(coeff.var.cases[, 1])) * length(coeff.var.cases[1, -(1:2)]) + i, 7] = regions[as.numeric(output.data[ (length(coeff.var.cases[, 1])) * length(coeff.var.cases[1, -(1:2)]) + i, 7] )]
    output.data[ (length(coeff.var.cases[, 1])) * length(coeff.var.cases[1, -(1:2)])  + i, 1] = ""
  }
  output.data$Coefficient.of.Variation = as.numeric(output.data$Coefficient.of.Variation)
  output.data$Incidence  =  as.numeric(output.data$Incidence)
  output.data$Mean.vaccination  =  as.numeric(output.data$Mean.vaccination)
  output.data$Mean.birth.rate  =  as.numeric(output.data$Mean.birth.rate)
  output.data$Year   =  as.numeric(output.data$Year)
  output.data$Coefficient.of.Variation[which(output.data$Coefficient.of.Variation == "Inf")] = 0
  return(output.data)
}




###################################
# Output coefficient of variation for a variable year periods, with 1 year moving window. 
# Can specify which regions we are interested in using this function
# Also output mean birth and vaccination rates along with cumulative incidence in those 10 years
###################################
output.coeff.of.var.and.other.variables.specify.regions.col.forward.5.years <- function(back.length, forward.length, regions, official.vacc){
  
  ##' load required data
  
  cases.by.country.by.year = read.csv("Measles_cases_by_year.csv", stringsAsFactors = FALSE)
  Birth.rates = read.csv("Birth_rates.csv", stringsAsFactors = FALSE)
  pop.by.year = read.csv("All_populations.csv", stringsAsFactors = FALSE)
  
  ##' Specify if we want to consider the official vaccination proportions or estimates. Estimates are probably more accurate than official data.
  
  if(official.vacc == 1){
    vacc.rates = read.csv("Vacc_rates_official.csv", stringsAsFactors = FALSE)
    vacc.rates = cbind(vacc.rates[, 2], vacc.rates[, seq(36, 3, -1)], vacc.rates$WHO_REGION)
    colnames(vacc.rates) = c("Country",paste(seq(1980, 2013)), "WHO_REGION")
  } else{
    vacc.rates = read.csv("Measles_vac_all.csv", stringsAsFactors = FALSE)
  }
  
  ##' Only consider data for the WHO regions that we specify, i.e. only include data from African countries, if the only region is "AFR"
  
  subset.pop.by.year = subset(pop.by.year, pop.by.year$WHO_REGION %in% regions)
  subset.vaccination = subset(vacc.rates, vacc.rates$WHO_REGION %in% regions)
  subset.birth.rates = subset(Birth.rates, Birth.rates$WHO_REGION %in% regions)
  subset.data = subset(cases.by.country.by.year, cases.by.country.by.year$WHO_REGION %in% regions)
  
  ##' Make sure that all the data sets only consider the same set of countries.
  
  missing1 = setdiff(subset.vaccination$Country, subset.pop.by.year$Country.Name)
  if(length(missing1) > 0){
    j = which(subset.vaccination$Country %in% missing1)
    subset.vaccination = subset.vaccination[-j, ]
  }
  missing2 = setdiff(subset.birth.rates$Country, subset.vaccination$Country)
  if(length(missing2) > 0){
    j = which(subset.birth.rates$Country %in% missing2)
    subset.birth.rates = subset.birth.rates[-j, ]
  }
  missing3 = setdiff(subset.data$Cname, subset.vaccination$Country)
  if(length(missing3) > 0){
    j = which(subset.data$Cname %in% missing3)
    subset.data = subset.data[-j, ]
  }
  ##' Order datasets in alphabetical order according to country.
  
  subset.pop.by.year = subset.pop.by.year[order(subset.pop.by.year$Country.Name), ]
  subset.vaccination = subset.vaccination[order(subset.vaccination$Country), ]
  subset.birth.rates = subset.birth.rates[order(subset.birth.rates$Country), ]
  subset.data = subset.data[order(subset.data$Cname), ]
  
  ##' Re-order and rename columns
  
  population  =  cbind(subset.pop.by.year$Country.Name, subset.pop.by.year$WHO_REGION, subset.pop.by.year[paste("X", seq(1980,2013), sep = "")])
  colnames(population) = c("Country", "WHO_REGION", seq(1980, 2013))
  vacc.rates  =  cbind(subset.vaccination$Country, subset.vaccination$WHO_REGION, subset.vaccination[paste("X", seq(1980,2013), sep = "")])
  colnames(vacc.rates) = c("Country", "WHO_REGION", seq(1980, 2013))
  birth.rates  =  cbind(subset.birth.rates$Country, subset.birth.rates$WHO_REGION, subset.birth.rates[paste("X", c(seq(1980,2012), 2012), sep = "")])
  colnames(birth.rates) = c("Country", "WHO_REGION", seq(1980, 2013))
  cases =  cbind(subset.data$Cname, subset.data$WHO_REGION, subset.data[paste("X", seq(1980,2013), sep = "")])
  colnames(cases) = c("Country", "WHO_REGION", seq(1980, 2013))
  
  ##' the number of data points that we will end up with is the number of years in the dataset minus the number of years that we use for retrospective analysis, given by back.length and also minus the number of years we look forward by, given by forward.length
  num.windows = length(seq(1980, 2013)) - back.length - forward.length
  initial.year = 1980
  year = initial.year
  
  
  ##' coeff.var is the coefficient variation over previous 10 years.
  ##' past.incidence will hold the incidence per 1000 over the previous 10 years. Also work out the mean birth rate in these years.
  
  coeff.var = matrix(0, length(cases[ , 1]) , num.windows)
  past.incidence = matrix(0, length(cases[ , 1]), num.windows)
  mean.br = matrix(0, length(cases[ , 1]), num.windows)
  
  for ( j in 1 : num.windows){
    for ( i in 1 : length(cases[ , 1])){
      
      ##' set the i,j entry to be the coefficient of variation of country i over the given 10 year period in the past. Entry which ultimately corresponds to 1990 in this matrix will be the coefficient of variation of cases from 1980-1989
      
      coeff.var[i, j]  =  sd(cases[i, paste(seq(year, year + back.length - 1))], na.rm = TRUE) /  
        mean(as.numeric(cases[i, paste(seq(year, year + back.length - 1))]), na.rm = TRUE)
      
      ##' it is possible that the denominator of this could be 0 or NA, so prevent this from giving strange results.
      if(is.na( mean(as.numeric(cases[i, paste(seq(year, year + back.length - 1))]), na.rm = TRUE)) == FALSE ){
        if( mean(as.numeric(cases[i, paste(seq(year, year + back.length - 1))]), na.rm = TRUE) == 0) {
          coeff.var[i, j]  =  0
        } 
      }
      past.incidence[i, j]  =  sum(as.numeric(cases[i,  paste(seq(year, year + back.length - 1))]) / 
                                     as.numeric(population[i, paste(seq(year, year + back.length - 1))]), na.rm = TRUE) * 1000
      if(length(which(is.na(as.numeric(birth.rates[i, paste(seq(year, year + back.length - 1))])))) < back.length){
        mean.br[i, j]  =  mean(as.numeric(birth.rates[i, paste(seq(year, year + back.length - 1))]), na.rm = TRUE)
      }
    }
    year = year + 1
  }

  colnames(coeff.var) = paste("X",seq(1990, 2013-forward.length), sep = "")
  colnames(past.incidence) =  paste("X",seq(1990, 2013-forward.length), sep = "")
  colnames(mean.br) = paste("X",seq(1990, 2013-forward.length), sep = "")
  
  
  
  ##' Calculate the incidence over the next number of years given by the input forward.length
  
  year = 1990
  future.incidence  =  matrix(0, length(cases[, 1]), num.windows) 
  for (j in 1 : num.windows){
    for ( i in 1 : length(cases[, 1])){
      future.incidence[i, j]  =  sum(as.numeric(cases[i,  paste(seq(year,year + forward.length))]) / 
                                       as.numeric(population[i, paste(seq(year,year + forward.length))]), na.rm = TRUE) * 1000
    }
    year = year + 1
  }
  colnames(future.incidence) = paste("X",seq(1990, 2013-forward.length), sep = "")
  
  
  return(list(coeff.var, past.incidence, future.incidence, mean.br, cases$Country, cases$WHO_REGION))
}




###################################
##' Interpolate a matrix
###################################
interpolate.data <- function(x, y, xout){
  Z = matrix(0, length(y[, 1]), length(xout))
  for(i in 1 : length(y[, 1])){
    if(length(which(!is.na(((y[i, ]))))) > 1){
      list[ A, B] = approx(x, y = as.numeric((y[i, ])), xout)
      Z[i,] = round(B, 2)
    }
  }
  return(Z)
}




###################################
# Output coefficient of variation for a variable year periods, with 1 year moving window. 
# Can specify which regions we are interested in using this function
# Also output mean birth and vaccination rates along with cumulative incidence in those 10 years
###################################
output.coeff.vs.incidence.data.specify.region.col.forward.5.years <- function(back.length, forward.length, interp.resolution, regions, official.vacc ){
  # output all data in matrix form
  
  list[coeff.var2, past.incidence2, future.incidence2, mean.br2, countries, region.code] = output.coeff.of.var.and.other.variables.specify.regions.col.forward.5.years(back.length, forward.length, regions, official.vacc)
  num.years = length(coeff.var2[1, ]) 
  num.countries = length(countries)
  
  # the matrix form of the data is interpolated so that the animation will be smooth. specify the amount of smoothing that we want via interp.resolution. If this is 10, then it means we will output 10 data points for each year.
  x = 1990:2008
  xout = seq(x[1], tail(x, 1), 1/interp.resolution)
  
  ##' interpolate the coefficient of variation matrix, add country and region data and name the columns appropriately
  coeff.var3 = interpolate.data(x, coeff.var2, xout)
  coeff.var =  data.frame((coeff.var3))
  colnames(coeff.var) = paste("X", xout, sep = "")
  
  ##' Do the same for past incidence, mean birth rate and future incidence
  past.inc3 = interpolate.data(x, past.incidence2, xout)
  past.incidence = data.frame(past.inc3)
  colnames(past.incidence) = paste("X", xout, sep = "")
  
  mean.br3 = interpolate.data(x, mean.br2, xout)
  mean.br = data.frame(mean.br3)
  colnames(mean.br) = paste("X", xout, sep = "")
  
  future.inc3 = interpolate.data(x, future.incidence2, xout)
  future.incidence = data.frame(future.inc3)
  colnames(future.incidence) = paste("X", xout, sep = "")
  

  return( list(coeff.var, past.incidence, future.incidence, mean.br, countries, region.code))
  
}







###################################
# Output data for google bubble chart where we can specify which regions we are interested in seeing
# Setting official.vac to 1 uses official WHO vaccination rates, whilst any other value uses WHO estimates
###################################
output.data.for.animation.specify.region.col.forward.5.years <- function(back.length, forward.length, interp.resolution, regions, official.vacc){
  
  list[coeff.var, past.incidence, future.incidence, mean.br, countries, region.code]  <- output.coeff.vs.incidence.data.specify.region.col.forward.5.years(back.length, forward.length, interp.resolution, regions, official.vacc)
  
  ##' for the animation, we output line list data, rather than a matrix. Therefore each country will have a line in the data set for every year the data covers.
  count = 1
  num.years = length(coeff.var[1, ])
  num.countries = length(countries)
  anim.data = matrix(0, num.years * num.countries, 7)
  anim.data  =  data.frame(anim.data)
  colnames(anim.data) = c("Country", "Coefficient.of.Variation", "Past.Incidence", "Future.Incidence", "Mean.birth.rate", "Year", "WHO_REGION")
  
  anim.data$Country = rep(countries, num.years)
  anim.data$WHO_REGION = rep(region.code, num.years)
  for(i in 1:length(coeff.var[1, ])){
    ##' Specify the set of indices that we want to fill in, these will be used for each named column
    indices = seq(((count - 1) * num.countries + 1), (count  * num.countries))
    
    anim.data$Coefficient.of.Variation[indices]  =  round(as.numeric(coeff.var[, i]), 2)
    anim.data$Past.Incidence[indices]  =  round(as.numeric(past.incidence[, i]),2)
    anim.data$Future.Incidence[indices]  =  round(as.numeric(future.incidence[, i]),  2)
    anim.data$Mean.birth.rate[indices]  =  round(as.numeric(mean.br[, i]), 2)
    anim.data$Year[indices]  =  (1990 + (i-1) / interp.resolution)
    count  =  count + 1
  }
  
 # anim.data2 = matrix(0, (length(coeff.var[, 1])) * length(coeff.var[1, ]) + (2 * length(regions) * length(unique(anim.data$Year))), 7)
#  anim.data2  =  data.frame(anim.data2)
#  colnames(anim.data2) = c("Country", "Coefficient.of.Variation", "Past.Incidence", "Future.Incidence", "Mean.birth.rate", "Year", "WHO_REGION")
#  anim.data2[1:length(anim.data[, 1]), 2:6]  =  anim.data[, 2:6]
#  anim.data2[1:length(anim.data[, 1]), 1]  =  as.character(anim.data[, 1])
#  anim.data2[1:length(anim.data[, 1]), 7]  =  as.character(anim.data[, 7])
  
#  year.mins = matrix(0, length(unique(anim.data$Year)), 2)
  
#  for(i in 1 : length(unique(anim.data$Year))){
#    t  =  subset(anim.data, anim.data$Year ==  unique(anim.data$Year)[i])
#    year.mins[i, 1] = unique(anim.data$Year)[i]
#    year.mins[i, 2] = min(t$Mean.birth.rate,na.rm = T)  
#  }
#  
#  l =  expand.grid("", -1, 0, c(0,max(anim.data$Future.Incidence)), 0,unique(anim.data$Year), regions)
#  for(i in 1 : (2 * length(regions) * length(unique(anim.data$Year)))){
#    y = l[i, 6]
#    j = which(year.mins[, 1] == y)
#    l[i, 5]  =  year.mins[j, 2]
#  }
  
#  anim.data2[(length(anim.data[, 1]) + 1 ): length(anim.data2[, 1]), ]  =  l
#  for( i in 1 : (2 * length(regions) * length(unique(anim.data$Year)))){
#    anim.data2[length(anim.data[, 1]) + i, 7] = regions[as.numeric(anim.data2[length(anim.data[, 1]) + i,7])]
#    anim.data2[length(anim.data[, 1]) + i, 1] = ""
#  }
#  return(anim.data2)
  return(anim.data)
}




################################
### Sort data for the animation process
################################
get.data.for.animation <- function(regions){
  
  #cases.by.country.by.year = read.csv("Measles_cases_by_year.csv", stringsAsFactors = FALSE)
  cases.by.country.by.year = read.csv("Measles_cases_by_year2.csv", stringsAsFactors = FALSE)
  Birth.rates = read.csv("Birth_rates.csv", , stringsAsFactors = FALSE)
  pop.by.year = read.csv("All_populations.csv", stringsAsFactors = FALSE)
  vacc.rates = read.csv("Measles_vac_all.csv", stringsAsFactors = FALSE)
  
  subset.pop.by.year = subset(pop.by.year, pop.by.year$WHO_REGION %in% regions)
  subset.vaccination = subset(vacc.rates, vacc.rates$WHO_REGION %in% regions)
  subset.birth.rates = subset(Birth.rates, Birth.rates$WHO_REGION %in% regions)
  subset.data = subset(cases.by.country.by.year, cases.by.country.by.year$WHO_REGION %in% regions)
  
  missing1 = setdiff(subset.vaccination$Country, subset.pop.by.year$Country.Name)
  if(length(missing1) > 0){
    j = which(subset.vaccination$Country %in% missing1)
    subset.vaccination = subset.vaccination[-j, ]
  }
  missing2 = setdiff(subset.birth.rates$Country, subset.vaccination$Country)
  if(length(missing2) > 0){
    j = which(subset.birth.rates$Country %in% missing2)
    subset.birth.rates = subset.birth.rates[-j, ]
  }
  missing3 = setdiff(subset.data$Cname, subset.vaccination$Country)
  if(length(missing3) > 0){
    j = which(subset.data$Cname %in% missing3)
    subset.data = subset.data[-j, ]
  }
  
  
  p1  =  subset.pop.by.year
  p2  =  subset.vaccination
  p3  =  subset.birth.rates
  p4  =  subset.data
  
  for ( i in 1 : length(subset.vaccination[, 1])){
    C  =  subset.vaccination$Country[i]
    p2[i, ]  =  subset.vaccination[i, ]
    j = which(subset.pop.by.year$Country.Name == C)
    p1[i, ]  =  subset.pop.by.year[j, ]
    j = which(subset.birth.rates$Country == C)
    p3[i, ]  =  subset.birth.rates[j, ]
    j = which(subset.data$Cname == C)
    p4[i, ]  =  subset.data[j, ]
  }
  subset.pop.by.year = p1
  subset.vaccination = p2
  subset.birth.rates = p3
  subset.data = p4
  
  return(list(subset.pop.by.year, subset.vaccination, subset.birth.rates, subset.data))
}



prepare.matrices.for.animation <- function(interp.subset.data, subset.data){
  
  mean.cases = matrix(0, length(interp.subset.data[, 1]), length(interp.subset.data[1, ]))
  mean.cases[, 1] = subset.data$Cname
  mean.cases[, 2] = subset.data$WHO_REGION
  
  coeff.var.cases = matrix(0, length(interp.subset.data[, 1]), length(interp.subset.data[1, ]))
  coeff.var.cases[, 1] = subset.data$Cname
  coeff.var.cases[, 2] = subset.data$WHO_REGION
  
  incidence.per.1000 = matrix(0, length(interp.subset.data[, 1]), length(interp.subset.data[1, ]))
  incidence.per.1000[, 1] = subset.data$Cname
  incidence.per.1000[, 2] = subset.data$WHO_REGION
  
  mean.br = matrix(0, length(interp.subset.data[, 1]), length(interp.subset.data[1, ]))
  mean.br[, 1] = subset.data$Cname
  mean.br[, 2] = subset.data$WHO_REGION
  
  mean.vac = matrix(0, length(interp.subset.data[, 1]), length(interp.subset.data[1, ]))
  mean.vac[, 1] = subset.data$Cname
  mean.vac[, 2] = subset.data$WHO_REGION
  
  return(list(mean.cases, coeff.var.cases, incidence.per.1000, mean.br, mean.vac))
}



############################
###
############################

calculate.all.weighted.means <- function(x, st.dev, cutoff, xout){
  N = length(x[, 1])
  M = length(x[1, ])
  means = matrix(0, N, M)
  for (j in 1 : M){
    w.input = xout - xout[j]
    w = output.weights.gaussian.with.cutoff(w.input, gaussian.st.dev, cutoff)
    for (i in 1 : N){
      means[i, j] = weighted.mean(as.numeric(x[i ,]), w)
    }
  }
  return(means)
}

###################################
# Output coefficient of variation for a variable year periods, with 1 year moving window. 
# Can specify which regions we are interested in using this function
# Also output mean birth and vaccination rates along with cumulative incidence in those 10 years
###################################
animation.gaussian.weighting.method <- function(window.length, regions, gaussian.st.dev, cutoff = 10){
  
  require(stats)
  list[subset.pop.by.year, subset.vaccination, subset.birth.rates, subset.data] = get.data.for.animation(regions)

  x = seq(1980, 2013)
  xout = seq(1980, 2013, 1)
  
  ### interpolate the datasets to have entries for all points in time once the interpolation is done.
  list[interp.subset.data, interp.subset.vacc, interp.subset.br, interp.subset.pop] = interp.datasets(subset.data, 
                                                                                                      subset.vaccination, 
                                                                                                      subset.birth.rates, 
                                                                                                      subset.pop.by.year,
                                                                                                      x,
                                                                                                      xout)
  
  ### output matrices the correct size for our animation
  
  list[mean.cases, coeff.var.cases, incidence.per.1000, mean.br, mean.vac] = prepare.matrices.for.animation(interp.subset.data, subset.data)
  
  mean.cases = calculate.all.weighted.means(interp.subset.data[, -(1:2)], st.dev = 7, cutoff = 50, xout)
  w = matrix(0, length(xout), length(xout))
  for(i in 1: length(xout)){
    w.input = xout - xout[i]
    w[i, ] = output.weights.gaussian.with.cutoff(w.input, gaussian.st.dev, cutoff)
  }
  
#   num.windows = 25 + (10 - window.length)
#   year = 1980
#   coeff.var = matrix(0, length(subset.data[ , 1]), num.windows)
#   incidence.per.1000 = matrix(0, length(subset.data[ , 1]), num.windows)
#   mean.br = matrix(0, length(subset.data[ , 1]), num.windows)
#   mean.vacc = matrix(0, length(subset.data[ , 1]), num.windows)
#   for ( j in 1 : num.windows){
#     for ( i in 1 : length(subset.data[ , 1])){
#       coeff.var[i, j]  =  sd(interp.subset.data[i, paste(seq(year, year + window.length - 1))], na.rm = TRUE) /  
#         mean(as.numeric(interp.subset.data[i, paste(seq(year, year + window.length - 1))]), na.rm = TRUE)
#       if(is.na( mean(as.numeric(interp.subset.data[i, paste(seq(year, year + window.length - 1))]), na.rm = TRUE)) == FALSE ){
#         if( mean(as.numeric(interp.subset.data[i, paste(seq(year, year + window.length - 1))]), na.rm = TRUE) == 0) {
#           coeff.var[i, j]  =  0
#         } 
#       }
#       incidence.per.1000[i, j]  =  sum(as.numeric(interp.subset.data[i,  paste(seq(year, year + window.length - 1))]) / 
#                                          as.numeric(interp.subset.pop[i, paste(seq(year, year + window.length - 1))]), na.rm = TRUE) * 1000
#       if(length(which(is.na(as.numeric(interp.subset.br[i, paste(seq(year, year + window.length - 1))])))) < window.length){
#         mean.br[i, j]  =  mean(as.numeric(interp.subset.br[i, paste(seq(year, year + window.length - 1))]), na.rm = TRUE)
#       }
#       if(length(which(is.na(as.numeric(interp.subset.vacc[i, paste(seq(year, year + window.length - 1))])))) < window.length){
#         mean.vacc[i, j]  =  mean(as.numeric(interp.subset.vacc[i, paste(seq(year, year + window.length - 1))]), na.rm = TRUE)
#       }
#     }
#     year = year + 1
#   }
#   
#   x1 = seq(1, length(coeff.var[1, ]))
#   w1 = matrix(0, length(x1), length(x1))
#   for (i in 1 : length(x1)){
#     w.input = x1 - x1[i]
#     w1[i, ] = output.weights.gaussian.with.cutoff(w.input, gaussian.st.dev, cutoff)
#   }
#   
#   coeff.2 = coeff.var
#   incidence.2 = incidence.per.1000
#   mbr2 = mean.br
#   mvacc2 = mean.vacc
#   for(i in 1 : length(coeff.var[1, ])){
#     for(j in 1 : length(coeff.var[, 1])){
#       coeff.2[j, i] = sum(coeff.var[j, ] * w1[i, ])
#       incidence.2[j, i] = sum(incidence.per.1000[j, ] * w1[i, ])
#       mbr2[j, i] = sum(mean.br[j, ] * w1[i, ])
#       mvacc2[j, i] = sum(mean.vacc[j, ] * w1[i, ])
#     }
#   }
#   coeff.var.cases = coeff.2
#   
  for(i in 1: length(xout)){
    for (j in 1 : length(interp.subset.data[, 1])){  
      
       if(mean.cases[j, i] == 0){
         coeff.var.cases[j, i + 2] = 0
       } else{
         coeff.var.cases[j, i + 2] = coeff.var.calc(as.numeric(interp.subset.data[j ,-(1:2)]), w[i, ], mean.cases[j, ]) 
       }
      incidence.per.1000[j, i + 2]  =  1000 * as.numeric(mean.cases[j, i])  /  (as.numeric(sum(as.numeric(interp.subset.pop[j ,-(1:2)], na.rm = T) * w) / sum(w)))
      mean.br[j, i + 2] = as.numeric(sum(as.numeric(interp.subset.br[j ,-(1:2)]) * w , na.rm = T) / sum(w))
      mean.vac[j, i + 2] = as.numeric(sum(as.numeric(interp.subset.vacc[j ,-(1:2)]) * w , na.rm = T) / sum(w))
    }
  }
  xout = seq(1980, 2013, 1/interp.resolution)
  mean.cases = cbind(mean.vac[,(1:2)], interpolate.give.dataset(mean.cases, x, xout))
  coeff.var.cases = cbind(coeff.var.cases[,(1:2)], interpolate.give.dataset(coeff.var.cases[, -(1:2)], x, xout))
  incidence.per.1000 = cbind(incidence.per.1000[,(1:2)], interpolate.give.dataset(incidence.per.1000[, -(1:2)], x, xout))
  mean.br = cbind(mean.br[,(1:2)], interpolate.give.dataset(mean.br[, -(1:2)], x, xout))
  mean.vac = cbind(mean.vac[,(1:2)], interpolate.give.dataset(mean.vac[, -(1:2)], x, xout))
  
  mean.vac[, -(1:2)] = round(as.numeric(mean.vac[, -(1:2)]), 2)
  mean.cases[, -(1:2)] = round(as.numeric(mean.cases[, -(1:2)]), 2)
  mean.br[, -(1:2)] = round(as.numeric(mean.br[, -(1:2)]), 2)
  incidence.per.1000[, -(1:2)] = round(as.numeric(incidence.per.1000[, -(1:2)]), 2)
  coeff.var.cases[, -(1:2)] = round(as.numeric(coeff.var.cases[, -(1:2)]), 2)
  
  
  output.data = matrix(0, (length(coeff.var.cases[, 1])) * length(coeff.var.cases[1, -(1:2)]) + (2 * length(regions) * length(xout)), 7)
  output.data  =  data.frame(output.data)
  colnames(output.data) = c("Country", "Coefficient.of.Variation", "Incidence", "Mean.vaccination", "Mean.birth.rate", "Year", "WHO_REGION")
  
  output.data$Country[seq(1, (length(coeff.var.cases[, 1])) * length(coeff.var.cases[1, -(1:2)]))] = rep(mean.cases[, 1], length(coeff.var.cases[1, -(1:2)]))
  output.data$WHO_REGION[seq(1, (length(coeff.var.cases[, 1])) * length(coeff.var.cases[1, -(1:2)]))] = rep(mean.cases[, 2], length(coeff.var.cases[1, -(1:2)]))
  count = 1
  for(i in 3 : length(coeff.var.cases[1, ])){
    output.data$Coefficient.of.Variation[seq(((count - 1) * (length(coeff.var.cases[, 1]))) + 1, count * (length(coeff.var.cases[, 1])))] = coeff.var.cases[, i]
    output.data$Incidence[seq(((count - 1) * (length(coeff.var.cases[, 1]))) + 1, count * (length(coeff.var.cases[, 1])))] = incidence.per.1000[, i]
    output.data$Mean.vaccination[seq(((count - 1) * (length(coeff.var.cases[, 1]))) + 1, count * (length(coeff.var.cases[, 1])))] = mean.vac[, i]
    output.data$Mean.birth.rate[seq(((count - 1) * (length(coeff.var.cases[, 1]))) + 1, count * (length(coeff.var.cases[, 1])))] = mean.br[, i]
    output.data$Year[seq(((count - 1) * (length(coeff.var.cases[, 1]))) + 1, count * (length(coeff.var.cases[, 1])))] = xout[i - 2]  
    count = count + 1
  }
  
  
  year.mins = matrix(0, length(xout), 2)
  
  for(i in 1 : length(xout)){
    t  =  subset(output.data, output.data$Year ==  unique(xout)[i])
    year.mins[i, 1] = xout[i]
    year.mins[i, 2] = as.numeric(min(t$Mean.birth.rate,na.rm = T)  )
  }
  
  l =  expand.grid("", -1, 0, c(0,100), 0, xout, regions)
  
  for(i in 1 : (2 * length(regions) * length(xout))){
    y = l[i, 6]
    j = which(year.mins[, 1] == y)
    l[i, 5]  =  year.mins[j, 2]
  }
  
  output.data[((length(coeff.var.cases[, 1])) * length(coeff.var.cases[1, -(1:2)]) + 1 ): length(output.data[, 1]), ]  =  l
  
  for( i in 1 : (2 * length(regions) * length(xout))){
    output.data[ (length(coeff.var.cases[, 1])) * length(coeff.var.cases[1, -(1:2)]) + i, 7] = regions[as.numeric(output.data[ (length(coeff.var.cases[, 1])) * length(coeff.var.cases[1, -(1:2)]) + i, 7] )]
    output.data[ (length(coeff.var.cases[, 1])) * length(coeff.var.cases[1, -(1:2)])  + i, 1] = ""
  }
  output.data$Coefficient.of.Variation = as.numeric(output.data$Coefficient.of.Variation)
  output.data$Incidence  =  as.numeric(output.data$Incidence)
  output.data$Mean.vaccination  =  as.numeric(output.data$Mean.vaccination)
  output.data$Mean.birth.rate  =  as.numeric(output.data$Mean.birth.rate)
  output.data$Year   =  as.numeric(output.data$Year)
  output.data$Coefficient.of.Variation[which(output.data$Coefficient.of.Variation == "Inf")] = 0
  return(output.data)
}





#####################

interpolate.give.dataset <- function(data,  
                                     x,
                                     xout){
  interp.data = matrix(0, length(data[, 1]), length(xout))
  for ( i in 1 : length(data[, 1])){
    y = data[i, ]
    if(length(which(!is.na(y)) == F) < 2){
      interp.data[i, ]  =  0} else{
        list[qq,ww] =  approx (as.numeric(x), as.numeric(y),  method = "linear", xout )
        interp.data[i, ]  =  round(ww, 2)
      }
  }
  
  return(interp.data)
  
}



####################
# interpolate datasets to a given grid
####################

####################
# interpolate datasets to a given grid
####################
interp.datasets <- function(subset.data, 
                            subset.vaccination, 
                            subset.birth.rates, 
                            subset.pop.by.year, 
                            x,
                            xout){
  
  interp.subset.data = matrix(0, length(subset.data[, 1]), length(xout) + 2)
  interp.subset.data[, 1] = subset.data$Cname
  interp.subset.data[, 2] = subset.data$WHO_REGION
  
  interp.subset.vacc = matrix(0, length(subset.vaccination[, 1]), length(xout) + 2)
  interp.subset.vacc[, 1] = subset.data$Cname
  interp.subset.vacc[, 2] = subset.data$WHO_REGION
  
  interp.subset.br = matrix(0, length(subset.birth.rates[, 1]), length(xout) + 2)
  interp.subset.br[, 1] = subset.data$Cname
  interp.subset.br[, 2] = subset.data$WHO_REGION
  
  interp.subset.pop = matrix(0, length(subset.pop.by.year[, 1]), length(xout) + 2)
  interp.subset.pop[, 1] = subset.data$Cname
  interp.subset.pop[, 2] = subset.data$WHO_REGION
  
  for ( i in 1 : length(subset.pop.by.year[, 1])){
    y = subset.data[i, paste("X", seq(1980, 2013), sep = "")]
    list[qq,ww] =  approx (as.numeric(x), as.numeric(y),  method = "linear", xout )
    interp.subset.data[i, 3: length(interp.subset.data[1, ])]  =  round(ww, 2)
    colnames(interp.subset.data) = c("Country", "WHO_REGION", seq(1980, 2013))
    
    y1 = subset.vaccination[i, paste("X", seq(1980, 2013), sep = "")]
    list[qq,ww] =  approx (as.numeric(x), as.numeric(y1),  method = "linear", xout )
    interp.subset.vacc[i, 3: length(interp.subset.vacc[1, ])]  =  round(ww, 2)
    colnames(interp.subset.vacc) = c("Country", "WHO_REGION", seq(1980, 2013))
    
    y2 =as.numeric( c(subset.birth.rates[i, paste("X", seq(1980, 2012), sep = "")], subset.birth.rates[i, paste("X", 2012, sep = "")]))
    if(length(which(is.na(y2) == F)) < 2) 
    {interp.subset.br[i, 3: length(interp.subset.data[1, ])]  = 0} else{
      list[qq,ww] =  approx (as.numeric(x), as.numeric(y2),  method = "linear", xout )
      interp.subset.br[i, 3: length(interp.subset.data[1, ])]  =  round(ww, 2)
    }
    colnames(interp.subset.br) = c("Country", "WHO_REGION", seq(1980, 2013))
    
    
    y3 = subset.pop.by.year[i, paste("X", seq(1980, 2013), sep = "")]
    list[qq,ww] =  approx (as.numeric(x), as.numeric(y3),  method = "linear", xout )
    interp.subset.pop[i, 3: length(interp.subset.data[1, ])]  =  round(ww, 2)
    colnames(interp.subset.pop) = c("Country", "WHO_REGION", seq(1980, 2013))
  }
  
  return(list(interp.subset.data, interp.subset.vacc, interp.subset.br, interp.subset.pop))
  
}





###############
# Calculate weighted mean of data, with input weights
###############
weighted.mean <- function(x, w, na.rm = T) {
  if (na.rm) {
    w <- w[i <- !is.na(x)]
    x <- x[i]
  }
  sum.w <- sum(w, na.rm = na.rm)
  mean.w <- sum(x * w, na.rm = na.rm) / sum(w)
  return(mean.w)
}




###############
# Calculate weighted variance of data, with input weights
###############
coeff.var.calc <- function(x, w, means, na.rm = T) {
 # if (na.rm) {
#    w <- w[i <- !is.na(x)]
 #   x <- x[i]
#  }
 # print(paste("w =", length(w),"x =", length(x), "means =", length(means)))
  sqrt(sum((w * (x - means)^2) / sum(w, na.rm = na.rm), na.rm = na.rm)) / sum(w * means, na.rm = na.rm)
  
  #mean.w <- sum(x * w, na.rm = na.rm) / sum(w)
  #(sum.w / (sum.w^2 - sum.w2)) * sum(w * (x - mean.w)^2, na.rm =
  #                                    na.rm)
  #sum(w * (x - mean.w) ^2 , na.rm = na.rm) / sum(w, na.rm = na.rm)
}

