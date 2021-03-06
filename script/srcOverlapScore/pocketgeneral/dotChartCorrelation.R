#!/usr/bin/env Rscript

require(plotrix)


radial.plot = function (lengths, radial.pos = NULL, labels = NA, label.pos = NULL, 
    radlab = FALSE, start = 0, clockwise = FALSE, rp.type = "r", 
    label.prop = 1.15, main = "", xlab = "", ylab = "", line.col = par("fg"), 
    lty = par("lty"), lwd = 4, mar = c(6, 6, 6, 6), 
    show.grid = TRUE, show.grid.labels = 3, show.radial.grid = TRUE, 
    grid.col = "#333333", grid.bg = "transparent", grid.left = FALSE, 
    grid.unit = NULL, point.symbols = NULL, point.col = NULL, 
    show.centroid = FALSE, radial.lim = NULL, radial.labels = NULL, 
    boxed.radial = TRUE, poly.col = NULL, add = FALSE, ...) 
{

    #modif environement
    par(cex.axis = 2.5)
    par(cex.lab = 2.5)	
    #print (par())

    if (is.null(radial.lim)) 
        radial.lim <- range(lengths)
    length.dim <- dim(lengths)
    if (is.null(length.dim)) {
        npoints <- length(lengths)
        nsets <- 1
        lengths <- matrix(lengths, nrow = 1)
    }
    else {
        npoints <- length.dim[2]
        nsets <- length.dim[1]
        lengths <- as.matrix(lengths)
    }
    lengths <- lengths - radial.lim[1]
    lengths[lengths < 0] <- NA
    if (is.null(radial.pos[1])) 
        radial.pos <- seq(0, pi * (2 - 2/npoints), length.out = npoints)
    radial.pos.dim <- dim(radial.pos)
    if (is.null(radial.pos.dim)) 
        radial.pos <- matrix(rep(radial.pos, nsets), nrow = nsets, 
            byrow = TRUE)
    else radial.pos <- as.matrix(radial.pos)
    if (clockwise) 
        radial.pos <- -radial.pos
    if (start) 
        radial.pos <- radial.pos + start
    if (show.grid) {
        if (length(radial.lim) < 3) 
            grid.pos <- pretty(radial.lim)
        else grid.pos <- radial.lim
        if (grid.pos[1] < radial.lim[1]) 
            grid.pos <- grid.pos[-1]
        maxlength <- max(grid.pos - radial.lim[1])
        angles <- seq(0, 1.96 * pi, by = 0.04 * pi)
    }
    else {
        grid.pos <- NA
        maxlength <- diff(radial.lim)
    }
    oldpar <- par("xpd", "mar", "pty")
    if (!add) {
        par(mar = mar, pty = "s")
        plot(c(-maxlength, maxlength), c(-maxlength, maxlength), 
            type = "n", axes = FALSE, main = main, xlab = xlab, 
            ylab = ylab)
        if (show.grid) {
            for (i in seq(length(grid.pos), 1, by = -1)) {
                xpos <- cos(angles) * (grid.pos[i] - radial.lim[1])
                ypos <- sin(angles) * (grid.pos[i] - radial.lim[1])
                polygon(xpos, ypos, border = grid.col, col = grid.bg)
            }
        }
    }
    par(xpd = TRUE)
    if (length(line.col) < nsets) 
        line.col <- 1:nsets
    if (length(rp.type) < nsets) 
        rp.type <- rep(rp.type, length.out = nsets)
    if (length(point.symbols) < nsets) 
        point.symbols <- rep(point.symbols, length.out = nsets)
    if (length(point.col) < nsets) 
        point.col <- rep(point.col, length.out = nsets)
    if (length(poly.col) < nsets) 
        poly.col <- rep(poly.col, length.out = nsets)
    if (length(lty) < nsets) 
        lty <- rep(lty, length.out = nsets)
    if (length(lwd) < nsets) 
        lwd <- rep(lwd, length.out = nsets)
    for (i in 1:nsets) {
        if (nsets > 1) {
            linecol <- line.col[i]
            polycol <- poly.col[i]
            pointcol <- point.col[i]
            pointsymbols <- point.symbols[i]
            ltype <- lty[i]
            lwidth <- lwd[i]
        }
        else {
            linecol <- line.col
            polycol <- poly.col
            pointcol <- point.col
            pointsymbols <- point.symbols
            ltype <- lty
            lwidth <- lwd
        }
        rptype <- unlist(strsplit(rp.type[i], ""))
        if (match("s", rptype, 0)) {
            if (is.null(pointsymbols)) 
                pointsymbols <- i
            if (is.null(pointcol)) 
                pointcol <- i
        }
        xpos <- cos(radial.pos[i, ]) * lengths[i, ]
        ypos <- sin(radial.pos[i, ]) * lengths[i, ]
        if (match("r", rptype, 0)) 
            segments(0, 0, xpos, ypos, col = linecol, lty = ltype, 
                lwd = lwidth, ...)
        if (match("p", rptype, 0)){ 
            polygon(xpos, ypos, border = linecol, col = polycol, 
                lty = ltype, lwd = lwidth, ...)
		print(ypos[1])
	}
        if (match("s", rptype, 0)) 
            points(xpos, ypos, pch = pointsymbols, col = pointcol, 
                ...)
        if (show.centroid) 
            if (match("p", rptype, 0)) {
                nvertices <- length(xpos)
                polygonarea <- xpos[nvertices] * ypos[1] - xpos[1] * 
                  ypos[nvertices]
                for (vertex in 1:(nvertices - 1)) polygonarea <- polygonarea + 
                  xpos[vertex] * ypos[vertex + 1] - xpos[vertex + 
                  1] * ypos[vertex]
                polygonarea <- polygonarea/2
                centroidx <- (xpos[nvertices] + xpos[1]) * (xpos[nvertices] * 
                  ypos[1] - xpos[1] * ypos[nvertices])
                centroidy <- (ypos[nvertices] + ypos[1]) * (xpos[nvertices] * 
                  ypos[1] - xpos[1] * ypos[nvertices])
                for (vertex in 1:(nvertices - 1)) {
                  centroidx <- centroidx + (xpos[vertex] + xpos[vertex + 
                    1]) * (xpos[vertex] * ypos[vertex + 1] - 
                    xpos[vertex + 1] * ypos[vertex])
                  centroidy <- centroidy + (ypos[vertex] + ypos[vertex + 
                    1]) * (xpos[vertex] * ypos[vertex + 1] - 
                    xpos[vertex + 1] * ypos[vertex])
                }
                points(centroidx/(6 * polygonarea), centroidy/(6 * 
                  polygonarea), col = point.col[i], pch = point.symbols[i], 
                  cex = 2, ...)
            }
            else points(mean(xpos), mean(ypos), col = pointcol, 
                pch = pointsymbols, cex = 2, ...)
    }
    if (!add) {
        if (is.na(labels[1])) {
            label.pos <- seq(0, 1.8 * pi, length = 9)
            labels <- as.character(round(label.pos, 2))
        }
        if (is.null(label.pos[1])) {
            lablen <- length(labels)
            label.pos <- seq(0, pi * (2 - 2/lablen), length.out = lablen)
        }
        if (clockwise) 
            label.pos <- -label.pos
        if (start) 
            label.pos <- label.pos + start
        xpos <- cos(label.pos) * maxlength
        ypos <- sin(label.pos) * maxlength
        if (show.radial.grid) 
            segments(0, 0, xpos, ypos, col = grid.col)
        xpos <- cos(label.pos) * maxlength * label.prop
        ypos <- sin(label.pos) * maxlength * label.prop
        if (radlab) {
            for (label in 1:length(labels)) {
                labelsrt <- (180 * label.pos[label]/pi) + 180 * 
                  (label.pos[label] > pi/2 && label.pos[label] < 
                    3 * pi/2)
		
                text(xpos[label], ypos[label], labels[label], 
                  cex = par("cex.axis"), srt = labelsrt)
            }
        }
        else{
	ypos[as.integer(length (ypos)/4) + 2]=ypos[as.integer(length (ypos)/4) +2] + 0.01
	ypos[as.integer(length (ypos)/4)+1 ]=ypos[as.integer(length (ypos)/4)+1 ] + 0.025
	ypos[as.integer(length (ypos)/4)+3 ]=ypos[as.integer(length (ypos)/4)+3 ] + 0.01


	ypos[as.integer(length (ypos)/1.3) ]=ypos[as.integer(length (ypos)/1.3) ] - 0.02
	ypos[as.integer(length (ypos)/1.3) +1 ]=ypos[as.integer(length (ypos)/1.3) ] + 0.015
	ypos[as.integer(length (ypos)/1.3) -1 ]=ypos[as.integer(length (ypos)/1.3) ] + 0.015

	xpos[as.integer(length (xpos)/1.3) +1 ]=xpos[as.integer(length (xpos)/1.3) ] + 0.2
	xpos[as.integer(length (xpos)/1.3) -1 ]=xpos[as.integer(length (xpos)/1.3) ]- 0.2


	boxed.labels(xpos, ypos, labels, ypad = 0.7, border = FALSE, 
            cex = par("cex.axis"))
		}
        if (show.grid.labels) {
            if (show.grid.labels%%2) {
                ypos <- grid.pos - radial.lim[1]
                xpos <- rep(0, length(grid.pos))
                if (show.grid.labels == 1) 
                  ypos <- -ypos
            }
            else {
                xpos <- grid.pos - radial.lim[1]
                ypos <- rep(0, length(grid.pos))
                if (show.grid.labels == 2) 
                  xpos <- -xpos
            }
            if (is.null(radial.labels)) 
                radial.labels = as.character(grid.pos)
            if (!is.null(grid.unit)) 
                radial.labels[length(grid.pos)] <- paste(radial.labels[length(grid.pos)], 
                  grid.unit)
            if (boxed.radial) 
                boxed.labels(xpos, ypos, radial.labels, border = FALSE, 
                  cex = par("cex.lab"))
            else text(xpos, ypos, radial.labels, cex = par("cex.lab"))
        }
    }
    invisible(oldpar)
}








###############################
###############################

args <- commandArgs(TRUE)
path_dataset = args[1]


d = read.table (path_dataset, sep = "\t", header = TRUE)


for (i in 1:dim(d)[2]){
	if (d[2,i] < 0.001){
		colnames (d)[i] = paste (colnames (d)[i], "***", sep = "")
	}
	else if (d[2,i] < 0.01){
		colnames (d)[i] = paste (colnames (d)[i], "**", sep = "")
	}
	else if (d[2,i] < 0.05){
		colnames (d)[i] = paste (colnames (d)[i], "*", sep = "")
	}

}


d = d[1,]
d = abs (d)
d = d[order(d)]


png (paste (path_dataset, "_radar.png", sep = ""), 2500, 2000)
par(mar=c(25,25,0,0))
radial.plot(d, labels=colnames (d),rp.type="p",main="", radial.lim=c(0,0.8),line.col="black", mar=c (10,10,10,10) )

dev.off()


svg (paste (path_dataset, "_radar.svg", sep = ""), 2500, 2000)
#par(mar=c(25,25,0,0))
radial.plot(d, labels=colnames (d),rp.type="p",main="", radial.lim=c(0,0.8),line.col="black", mar=c (10,10,10,10) )

dev.off()




print (mean(as.numeric(d[1,])))
print (sd(as.numeric(d[1,])))
print (min (as.numeric(d[1,])))
print (max (as.numeric(d[1,])))

print (mean(as.numeric(d[1,1:17])))
print (sd(as.numeric(d[1,1:17])))
print (max(as.numeric(d[1,1:17])))
print (min(as.numeric(d[1,1:17])))

print (mean(as.numeric(d[1,17:length(d[1,])])))
print (sd(as.numeric(d[1,17:length(d[1,])])))
print (max(as.numeric(d[1,17:length(d[1,])])))
print (min(as.numeric(d[1,17:length(d[1,])])))

png(paste (path_dataset, "_hist_global.png", sep = ""))
plot (density(as.numeric(d[1,])))
dev.off()

png(paste (path_dataset, "_hist_geo.png", sep = ""))
plot (density(as.numeric(d[1,1:17])))
dev.off()


png(paste (path_dataset, "_hist_physico.png", sep = ""))
plot (density(as.numeric(d[1,17:length(d[1,])])))
dev.off()

png(paste (path_dataset, ".png", sep = ""), 2000, 3508)
dotchart(as.matrix(d), cex = 1.6, pch = 4, main = "Coef correlation")
dev.off()	


png (paste (path_dataset, "_radar.png", sep = ""), 2500, 2000)
par(mar=c(25,25,0,0))
radial.plot(d, labels=colnames (d),rp.type="p",main="", radial.lim=c(0,0.8),line.col="black", mar=c (10,10,10,10) )

dev.off()


postscript (file = paste (path_dataset, "_radar.ps", sep = ""), width = 1000, height = 1000)
#par(mar=c(25,25,0,0))
radial.plot(d, labels=colnames (d),rp.type="p",main="", radial.lim=c(0,0.8),line.col="black")

dev.off()




