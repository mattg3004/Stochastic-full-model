library(RColorBrewer)
library(rgl)

plot3d.space.given.country.trajectory <- function(Countries, new.plot, manual.scale, xscale, yscale, zscale){
  
  
  years = unique(round(anim.data$Year))
  sub.data  =  subset(anim.data, anim.data$Country %in% Countries & anim.data$Year %in% years)
  cols <- colorRampPalette(brewer.pal(11,"RdYlGn"))(length(sub.data[, 1]))
  if(new.plot == 1){
    open3d()
  }
  if( manual.scale == 1){
    plot3d(sub.data$Coefficient.of.Variation, sub.data$Mean.vaccination, sub.data$Incidence,
           xlab = "Coefficient of variation", ylab = "Immunity", zlab = "Incidence", type = "l", col =  cols,
           xlim = xscale, ylim = yscale, zlim = zscale, main = Countries, forceClipregion = T, lwd =3)
    print(zscale)
  } else{
    plot3d(sub.data$Coefficient.of.Variation, sub.data$Mean.vaccination, sub.data$Incidence,
           xlab = "Coefficient of variation", ylab = "Immunity", zlab = "Incidence", type = "l", col =  cols, main = Countries,lwd =3)
  }
  
}


plot3d.space.given.time <- function(Countries, year, new.plot, manual.scale, xscale, yscale, zscale){
  
  sub.data  =  subset(anim.data, anim.data$Year == year & anim.data$Country %in% Countries)
  cols <- colorRampPalette(brewer.pal(11,"RdYlGn"))(length(Countries))
  if(new.plot == 1){
    open3d()
  }
  if( manual.scale == 1){
    plot3d(sub.data$Coefficient.of.Variation, sub.data$Mean.vaccination, sub.data$Incidence,
           xlab = "Coefficient of variation", ylab = "Immunity", zlab = "Incidence", type = "s", col =  cols,
           xlim = xscale, ylim = yscale, zlim = zscale, main = year)
  } else{
    plot3d(sub.data$Coefficient.of.Variation, sub.data$Mean.vaccination, sub.data$Incidence,
           xlab = "Coefficient of variation", ylab = "Immunity", zlab = "Incidence", type = "s", col =  cols, main = year)
  }
  
}








plot3d.space.given.time(Countries = All.African.Countries, year = 1980, new.plot = 1, manual.scale = 1, xscale = c(0, 3), yscale = c(0, 100), zscale = c(0, 150))
plot3d.space.given.time(Countries = All.African.Countries, year = 1990, new.plot = 1, manual.scale = 1, xscale = c(0, 3), yscale = c(0, 100), zscale = c(0, 150))
plot3d.space.given.time(Countries = All.African.Countries, year = 2000, new.plot = 1, manual.scale = 1, xscale = c(0, 3), yscale = c(0, 100), zscale = c(0, 150))




plot3d.space.given.country.trajectory(Countries = "Malawi" , new.plot = 1, manual.scale = 0, xscale = c(0, 3), yscale = c(0, 100), zscale = c(0, 150))
plot3d.space.given.country.trajectory(Countries = "Cabo Verde" , new.plot = 1, manual.scale = 0, xscale = c(0, 3), yscale = c(0, 100), zscale = c(0, 150))
plot3d.space.given.country.trajectory(Countries = "Tanzania" , new.plot = 1, manual.scale = 0, xscale = c(0, 1), yscale = c(0, 100), zscale = c(0, 1))



