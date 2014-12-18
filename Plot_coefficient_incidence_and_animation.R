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
countries.to.plot = African.countries
for ( i in 1 : (length(coeff.var[1, ])-1)){
  plot.coeff.var(African.countries, countries.to.plot, start.year+i, 
                 window.length = window.length,'none', text.size )
  ani.record()
}
oopts = ani.options(interval = 0.5)
ani.replay()

saveHTML(ani.replay())


#saveVideo(ani.replay(), video.name = "graphanimation.mp4", ffmpeg = "/Users/Matthew/Downloads/SnowLeopard_Lion_Mountain_Lion_Mavericks_Yosemite_17.12.2014/ffmpeg")
