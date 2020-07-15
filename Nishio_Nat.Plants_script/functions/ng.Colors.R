
#####################
##### ng.Colors #####
#####################

#original color functions
#(originally from tim.colors)
#
#ng.colors.helper()
#   make intermediate color names by cubic spline interpolation
#ng.tim.colors()
#   black & white added to tim.colors
#ng.po.colors()
#   black -> purple -> red -> yellow -> white
#
#<Example>
#x<- 1:40; y<- 1:60; z<- outer( x,y,"+") 
#image.plot(x, y, z, col=ng.po.colors(64)) 
#

library(fields) #for splint()

ng.colors.helper <- function(orig, n){
  rgb.tim <- t(col2rgb(orig))
  temp <- matrix(NA, ncol = 3, nrow = n)
  x <- seq(0, 1, , length(orig))
  xg <- seq(0, 1, , n)
  for (k in 1:3) {
    hold <- splint(x, rgb.tim[, k], xg)
    hold[hold < 0] <- 0
    hold[hold > 255] <- 255
    temp[, k] <- round(hold)
  }
  return(temp)
}

ng.tim.colors <- function (n = 64)
{
  #   orig <- c("#000000", "#880088", "#FF0000", "#FFFF00", "#FFFFFF")
  orig <- c("#000000", "#00001D", "#000039", "#000055", "#000070",
            "#00008A", "#0000A3", "#0000BB", "#0000D1", "#0000E5", "#0000F7",
            "#0109FF", "#021EFF", "#0437FF", "#0552FF", "#066EFF", "#078BFF",
            "#07A7FF", "#07C2FF", "#06DAFF", "#03EFFF", "#00FFFF", "#00FFEA",
            "#00FFD2", "#00FFB8", "#00FF9C", "#00FF80", "#00FF64", "#00FF49",
            "#00FF30", "#00FF1A", "#00FF07", "#08FF00", "#1CFF00", "#34FF00",
            "#4EFF00", "#6AFF00", "#87FF00", "#A3FF00", "#BEFF00", "#D7FF00",
            "#EDFF00", "#FFFF00", "#FFEB00", "#FFD200", "#FFB600", "#FF9700",
            "#FF7800", "#FF5A00", "#FF3D00", "#FF2500", "#FF1100", "#FF0400",
            "#FD0006", "#FA0114", "#F80B25", "#F61B3A", "#F63151", "#F64C6A",
            "#F76B86", "#F98DA3", "#FBB2C1", "#FDD8E0", "#FFFFFF")
  if (n == 64)
    return(orig)
  temp <- ng.colors.helper(orig, n)
  rgb(temp[, 1], temp[, 2], temp[, 3], maxColorValue = 255)
}
cat("ng.tim.colors() loaded\n")

ng.po.colors <- function (n = 64)
{
  orig <- c("#000000", "#880088", "#FF0000", "#FFFF00", "#FFFFFF")
  temp <- ng.colors.helper(orig, n)
  rgb(temp[, 1], temp[, 2], temp[, 3], maxColorValue = 255)
}
cat("ng.po.colors() loaded\n")

#memo to develop another color function
if(0){
  n <- 128
  orig <- c("#000000", "#880088", "#FF0000", "#FFFF00", "#FFFFFF")
  rgb.tim <- t(col2rgb(orig))
  temp <- matrix(NA, ncol = 3, nrow = n)
  x <- seq(0, 1, , length(orig))
  xg <- seq(0, 1, , n)
  for (k in 1:3) {
    hold <- splint(x, rgb.tim[, k], xg)
    hold[hold < 0] <- 0
    hold[hold > 255] <- 255
    temp[, k] <- round(hold)
  }
  plot(rep(1, n), col=rgb(temp[, 1], temp[, 2], temp[, 3], maxColorValue = 255), pch=16)
  #    rgb(temp[, 1], temp[, 2], temp[, 3], maxColorValue = 255)
}
