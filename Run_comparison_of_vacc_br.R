require("foreach")
source("Initial_conditions.R")
source("Functions.R")

require("doParallel")
cl <- makeCluster(4, outfile="")
registerDoParallel(cl)
list[out4, d.out] = numerous.sims(num.replicates = 4, num.years =  2, initial.disease.state =  disease.state.2, initial.start.time = 280)

sim.results = matrix(0, replicates, replicates)
for (i in 1 : replicates){
  sim.results[i, ]  =  unlist(out1[i])
}

small.outbreak = 1000
mid.outbreak = 5000
large.outbreak = 25000
br = seq(50, 20, -7.5)
vaccs = seq(0.6, 0.9, 0.05)
count = 1
pr.small.outbreak2 = matrix(0, length(br) * length(vaccs), 24)
pr.large.outbreak2 = matrix(0, length(br) * length(vaccs), 24)
pr.mid.outbreak2 = matrix(0, length(br) * length(vaccs), 24)
large.outbreak.in.3.year.periods = matrix(0, length(br) * length(vaccs), 8)
mid.outbreak.in.3.year.periods = matrix(0, length(br) * length(vaccs), 8)
small.outbreak.in.3.year.periods = matrix(0, length(br) * length(vaccs), 8)
inputs = matrix(0, length(br) * length(vaccs), 2)
cl <- makeCluster(4, outfile="")
replicates = 24
years = 24
sim.results = matrix(0, replicates, years)
registerDoParallel(cl)
count = 1
for(i in 1:length(vaccs)){
  for (j in 1 : length(br)){
    birth.rate =  br[j]
    vacc.prop  =  vaccs[i]
    output1   =   numerous.sims(num.replicates = replicates, num.years = years, initial.disease.state =  disease.state.2, initial.start.time = 280)
    for (k in 1 : years){
      sim.results[k, ]  =  unlist(output1[k])
    }
    for ( k in 1 : replicates){
      pr.small.outbreak2[count, k]   =   length(which(sim.results[, k ] > small.outbreak)) / replicates
      pr.large.outbreak2[count, k]   =   length(which(sim.results[, k ] > large.outbreak)) / replicates  
      pr.mid.outbreak2[count, k]     =   length(which(sim.results[, k ] > mid.outbreak)) / replicates
    }
    for(k in 1 : 8){
      large.outbreak.in.3.year.periods[count, k]   =   length(which( (sim.results[, (k-1)*3 + 1 ] > large.outbreak) |  (sim.results[, (k-1)*3 + 2 ] > large.outbreak) | (sim.results[, (k*3) ] > large.outbreak)) ) / replicates  
      mid.outbreak.in.3.year.periods[count, k]     =   length(which( (sim.results[, (k-1)*3 + 1 ] > mid.outbreak) |  (sim.results[, (k-1)*3 + 2 ] > mid.outbreak) | (sim.results[, (k*3) ] > mid.outbreak)) ) / replicates
      small.outbreak.in.3.year.periods[count, k]   =   length(which( (sim.results[, (k-1)*3 + 1 ] > small.outbreak) |  (sim.results[, (k-1)*3 + 2 ] > small.outbreak) | (sim.results[, (k*3) ] > small.outbreak)) ) / replicates
    }
    
    inputs[count, 1]  =  birth.rate
    inputs[count, 2]  =  vacc.prop
    count  =  count+1
    print(paste('Done', count, 'of',  length(br)*length(vaccs)))
  }  
}