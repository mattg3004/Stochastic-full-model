require(ggplot2)
frames = 64


for(i in 1 :frames){
  setwd("/Users/mgraham/Dropbox/Measles/Stochastic-full-model (1)/Animations_Lower_Gaussian_st_dev")
  if (i < 10) {name = paste('000',i,'plot.png',sep='')}
  
  if (i < 100 && i >= 10) {name = paste('00',i,'plot.png', sep='')}
  png(name, height = 1000, width = 1000)
  
  if(i < 33){
    plots.for.gif(data = anim.data, year = min(anim.data$Year)+i-1, region = 'AFR',
                  previously.plotted.year = 1990, scaling = 'sqrt', 
                  arrow.color = 'red', text.size = 5, region.text = "Africa", 
                  window.length=10, empty = 0)
  }
  
#   if(i > 30 & i < 36){
#     plots.for.gif(data = anim.data, year = 2013, region = 'AFR',
#                   previously.plotted.year = 1990, scaling = 'sqrt', 
#                   arrow.color = 'red', text.size = 5, region.text = "Africa", 
#                   window.length=10, empty = 0)
#   }
#   
  if( i > 32 & i < 65){
    plots.for.gif(data = anim.data, year = min(anim.data$Year)+i-33, region = 'AMR',
                  previously.plotted.year = 1990, scaling = 'sqrt', 
                  arrow.color = 'red', text.size = 5, region.text = "Americas", 
                  window.length=10)
  }
#   if(i > 64){
#     plots.for.gif(data = anim.data, year = 2013, region = 'AMR',
#                   previously.plotted.year = 1990, scaling = 'sqrt', 
#                   arrow.color = 'red', text.size = 5, region.text = "Americas", 
#                   window.length=10)
#   }
  dev.off()
  setwd("/Users/mgraham/Dropbox/Measles/Stochastic-full-model (1)")
}



frames = 30

for(i in 1 :frames){
  setwd("/Users/mgraham/Dropbox/Measles/Stochastic-full-model (1)/Animation_pngs_Africa")
#  if (i < 10) {name = paste('000',i),'plot.png',sep='')}
  
  #if (i < 100 && i >= 10) 
  name = paste('00',i + 40,'plot.png', sep='')
  png(name, height = 1000, width = 1000)
  plots.for.gif(data = anim.data, year = 1984+i-1, region = 'AMR',
                previously.plotted.year = 1990, scaling = 'sqrt', 
                arrow.color = 'red', text.size = 5, region.text = "Americas", 
                window.length=10)
  
  dev.off()
  setwd("/Users/mgraham/Dropbox/Measles/Stochastic-full-model (1)")
}
