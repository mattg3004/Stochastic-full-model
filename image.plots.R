library(RColorBrewer)
library(akima)
library(fields)


quartz()
x1 = 1 / inputs[, 1]
y1 = inputs[, 2]
prob.outbreak = small.outbreak.in.3.year.periods[, 2]
s <- interp(x1,y1,prob.outbreak, xo=seq(min(x1), max(x1), length = 400),
            yo=seq(min(y1), max(y1), length = 400))

image.plot(s, cex.lab = 2, cex.axis = 2, cex.main = 2, cex.sub = 2)
contour(s, add = TRUE)
points(1/42,0.9 * 0.95, cex = 5, col = rgb(1,1,1), lwd = 3)
points(1/43, .88*.9, cex = 5, pch = 0, col = rgb(1,1,1), lwd = 3)
points(1/27,0.8 * 0.9, cex = 5, pch = 2, col = rgb(1,1,1),lwd = 3)
points(1/31,0.72 * 0.9, cex = 5, pch = 5, col = rgb(1,1,1), lwd = 3)


quartz()
x1 = 1 / inputs[, 1]
y1 = inputs[, 2]
prob.outbreak = mid.outbreak.in.3.year.periods[, 2]
s <- interp(x1,y1,prob.outbreak, xo=seq(min(x1), max(x1), length = 400),
            yo=seq(min(y1), max(y1), length = 400))

image.plot(s, cex.lab = 2, cex.axis = 2, cex.main = 2, cex.sub = 2)
contour(s, add = TRUE)
points(1/42,0.9 * 0.95, cex = 5, col = rgb(1,1,1), lwd = 3)
points(1/43, .88*.9, cex = 5, pch = 0, col = rgb(1,1,1), lwd = 3)
points(1/27,0.8 * 0.9, cex = 5, pch = 2, col = rgb(1,1,1),lwd = 3)
points(1/31,0.72 * 0.9, cex = 5, pch = 5, col = rgb(1,1,1), lwd = 3)




quartz()
x1 = 1 / inputs[, 1]
y1 = inputs[, 2]
prob.outbreak = large.outbreak.in.3.year.periods[, 2]
s <- interp(x1,y1,prob.outbreak, xo=seq(min(x1), max(x1), length = 400),
            yo=seq(min(y1), max(y1), length = 400))

image.plot(s, cex.lab = 2, cex.axis = 2, cex.main = 2, cex.sub = 2)
contour(s, add = TRUE)
points(1/42,0.9 * 0.95, cex = 5, col = rgb(1,1,1), lwd = 3)
points(1/43, .88*.9, cex = 5, pch = 0, col = rgb(1,1,1), lwd = 3)
points(1/27,0.8 * 0.9, cex = 5, pch = 2, col = rgb(1,1,1),lwd = 3)
points(1/31,0.72 * 0.9, cex = 5, pch = 5, col = rgb(1,1,1), lwd = 3)


quartz()
x1 = 1 / inputs[, 1]
y1 = inputs[, 2]
prob.outbreak = small.outbreak.in.3.year.periods[, 4]
s <- interp(x1,y1,prob.outbreak, xo=seq(min(x1), max(x1), length = 400),
            yo=seq(min(y1), max(y1), length = 400))

image.plot(s, cex.lab = 2, cex.axis = 2, cex.main = 2, cex.sub = 2)
contour(s, add = TRUE)
points(1/42,0.9 * 0.95, cex = 5, col = rgb(1,1,1), lwd = 3)
points(1/43, .88*.9, cex = 5, pch = 0, col = rgb(1,1,1), lwd = 3)
points(1/27,0.8 * 0.9, cex = 5, pch = 2, col = rgb(1,1,1),lwd = 3)
points(1/31,0.72 * 0.9, cex = 5, pch = 5, col = rgb(1,1,1), lwd = 3)


quartz()
x1 = 1 / inputs[, 1]
y1 = inputs[, 2]
prob.outbreak = mid.outbreak.in.3.year.periods[, 4]
s <- interp(x1,y1,prob.outbreak, xo=seq(min(x1), max(x1), length = 400),
            yo=seq(min(y1), max(y1), length = 400))

image.plot(s, cex.lab = 2, cex.axis = 2, cex.main = 2, cex.sub = 2)
contour(s, add = TRUE)
points(1/42,0.9 * 0.95, cex = 5, col = rgb(1,1,1), lwd = 3)
points(1/43, .88*.9, cex = 5, pch = 0, col = rgb(1,1,1), lwd = 3)
points(1/27,0.8 * 0.9, cex = 5, pch = 2, col = rgb(1,1,1),lwd = 3)
points(1/31,0.72 * 0.9, cex = 5, pch = 5, col = rgb(1,1,1), lwd = 3)




quartz()
x1 = 1 / inputs[, 1]
y1 = inputs[, 2]
prob.outbreak = large.outbreak.in.3.year.periods[, 4]
s <- interp(x1,y1,prob.outbreak, xo=seq(min(x1), max(x1), length = 400),
            yo=seq(min(y1), max(y1), length = 400))

image.plot(s, cex.lab = 2, cex.axis = 2, cex.main = 2, cex.sub = 2)
contour(s, add = TRUE)
points(1/42,0.9 * 0.95, cex = 5, col = rgb(1,1,1), lwd = 3)
points(1/43, .88*.9, cex = 5, pch = 0, col = rgb(1,1,1), lwd = 3)
points(1/27,0.8 * 0.9, cex = 5, pch = 2, col = rgb(1,1,1),lwd = 3)
points(1/31,0.72 * 0.9, cex = 5, pch = 5, col = rgb(1,1,1), lwd = 3)

quartz()
x1 = 1 / inputs[, 1]
y1 = inputs[, 2]
prob.outbreak = small.outbreak.in.3.year.periods[, 8]
s <- interp(x1,y1,prob.outbreak, xo=seq(min(x1), max(x1), length = 400),
            yo=seq(min(y1), max(y1), length = 400))

image.plot(s, cex.lab = 2, cex.axis = 2, cex.main = 2, cex.sub = 2)
contour(s, add = TRUE)
points(1/42,0.9 * 0.95, cex = 5, col = rgb(1,1,1), lwd = 3)
points(1/43, .88*.9, cex = 5, pch = 0, col = rgb(1,1,1), lwd = 3)
points(1/27,0.8 * 0.9, cex = 5, pch = 2, col = rgb(1,1,1),lwd = 3)
points(1/31,0.72 * 0.9, cex = 5, pch = 5, col = rgb(1,1,1), lwd = 3)



quartz()
x1 = 1 / inputs[, 1]
y1 = inputs[, 2]
prob.outbreak = mid.outbreak.in.3.year.periods[, 8]
s <- interp(x1,y1,prob.outbreak, xo=seq(min(x1), max(x1), length = 400),
            yo=seq(min(y1), max(y1), length = 400))

image.plot(s, cex.lab = 2, cex.axis = 2, cex.main = 2, cex.sub = 2)
contour(s, add = TRUE)
points(1/42,0.9 * 0.95, cex = 5, col = rgb(1,1,1), lwd = 3)
points(1/43, .88*.9, cex = 5, pch = 0, col = rgb(1,1,1), lwd = 3)
points(1/27,0.8 * 0.9, cex = 5, pch = 2, col = rgb(1,1,1),lwd = 3)
points(1/31,0.72 * 0.9, cex = 5, pch = 5, col = rgb(1,1,1), lwd = 3)



quartz()
x1 = 1 / inputs[, 1]
y1 = inputs[, 2]
prob.outbreak = large.outbreak.in.3.year.periods[, 8]
s <- interp(x1,y1,prob.outbreak, xo=seq(min(x1), max(x1), length = 400),
            yo=seq(min(y1), max(y1), length = 400))

image.plot(s, cex.lab = 2, cex.axis = 2, cex.main = 2, cex.sub = 2)
contour(s, add = TRUE)
points(1/42,0.9 * 0.95, cex = 5, col = rgb(1,1,1), lwd = 3)
points(1/43, .88*.9, cex = 5, pch = 0, col = rgb(1,1,1), lwd = 3)
points(1/27,0.8 * 0.9, cex = 5, pch = 2, col = rgb(1,1,1),lwd = 3)
points(1/31,0.72 * 0.9, cex = 5, pch = 5, col = rgb(1,1,1), lwd = 3)






quartz()
x1 = 1 / inputs[, 1]
y1 = inputs[, 2]
prob.outbreak = pr.large.outbreak2[, 24]
s <- interp(x1,y1,prob.outbreak, xo=seq(min(x1), max(x1), length = 400),
            yo=seq(min(y1), max(y1), length = 400))

image.plot(s, cex.lab = 2, cex.axis = 2, cex.main = 2, cex.sub = 2)
contour(s, add = TRUE)
points(1/42,0.9 * 0.95, cex = 5, col = rgb(1,1,1), lwd = 3)
points(1/43, .88*.9, cex = 5, pch = 0, col = rgb(1,1,1), lwd = 3)
points(1/27,0.8 * 0.9, cex = 5, pch = 2, col = rgb(1,1,1),lwd = 3)
points(1/31,0.72 * 0.9, cex = 5, pch = 5, col = rgb(1,1,1), lwd = 3)




