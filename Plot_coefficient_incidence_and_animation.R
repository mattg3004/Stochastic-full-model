#load libraries
library(igraph)
library(animation)
library(ggplot2)


ani.record(reset=TRUE)
quartz()
window.length = 10
start.year  =  1980
text.size   =  2
countries.to.plot = c("Malawi", "Niger", "Tanzania", "Zambia","Algeria", "Namibia", "Senegal")
#countries.to.plot = African.countries
upper.limit = 25
scaling = 'none'

interp.resolution = 10

list[coeff.var, incidence.per.1000, mean.vac, mean.br]  <- output.coeff.vs.incidence.data(window.length, interp.resolution )

                        
count = 1
anim.data = matrix(0, length(coeff.var[, 1]) * length(coeff.var[1, ]), 6)
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



for ( i in 1 : (length(coeff.var[1, ]))){
  plot.coeff.var(African.countries, countries.to.plot, start.year + floor((i-1) / interp.resolution), 
                 coeff.var, incidence.per.1000, mean.vac, mean.br,
                 window.length = window.length, scaling,
                 text.size, upper.limit, column.number = i )
  ani.record()
}
#oopts = ani.options(interval = 0.5)
#ani.replay()

saveHTML(ani.replay(), ani.height = 500, ani.width = 500)

specific.countries = c("Malawi", "Niger", "Tanzania", "Zambia","Algeria", "Namibia", "Senegal")
coeff.vs.incidence.plot(countries.to.plot = 'All', scaling = 'log', text.size = 2, ani.h = 500,  ani.w = 500 )

#saveVideo(ani.replay(), video.name = "graphanimation.mp4", ffmpeg = "/Users/Matthew/Downloads/SnowLeopard_Lion_Mountain_Lion_Mavericks_Yosemite_17.12.2014/ffmpeg")
