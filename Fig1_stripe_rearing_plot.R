## Code for plot Figure 1 
## 2020.10.21
rm(list = ls())
setwd(choose.dir()) # set working directory to the folder of the project (~/EmergeProbV1)
source("utils.R")
library(circular)

dt = 1 # ms
ne = 80

par(mfrow = c(2,4), mar = c(4,4,4,3), lwd  =2)
## A schematic plot of the network (left blank)
plot.new()

## Figure 1B: STDP learning rule
tm_plst = 34 # tau - (ms)
tp_plst = 14 # tau + (ms)
mp_ratio = 1.05 # ratio between LTD and LTP
Bm_plst =  exp((-1./tm_plst) *dt)
Bp_plst =  exp((-1./tp_plst) *dt)
A_ltp = 5.e-3 # amplitude for LTP
A_ltd = mp_ratio*tp_plst*A_ltp/tm_plst # amplitude for LTD
x_STDP_minus = seq(-100,0,.1)
x_STDP_plus = seq(0,100,.1)
y_STDP_minus = -A_ltd*exp(x_STDP_minus/tm_plst)
y_STDP_plus = A_ltp*exp(-x_STDP_plus/tp_plst)
plot(x_STDP_minus, y_STDP_minus, ylim = c(-0.003, 0.005), xlim = c(-100,100), type = "l", 
     xlab = "Time", ylab = "Weight change", main = "B STDP learning rule", bty = "n")
abline(h = 0, col = "grey50")
abline(v = 0, col = "grey50")
lines(x_STDP_plus, y_STDP_plus)

## Figure 1A: Continuous Inputt Pattern
# plot rate for 90 degree input
pre_n=100
ff.bg = 10 # back ground firing rate (gaussian function)
ff.A = 40 # peak firing rate (gaussian function)
ff.sigma = 10 # sigma (gaussian function)
g.x=1:pre_n
ff_rate <- function(ff.mu){
  f=ff.bg + ff.A * exp( -pmin((g.x-ff.mu)^2, (g.x-pre_n-ff.mu)^2, (g.x+pre_n-ff.mu)^2)/(2*ff.sigma^2) )
}
plot(g.x, ff_rate(pre_n/2), type = "l", xlab = "Input neuron", ylab = "Spike rates (Hz)",ylim = c(0,60),
     main = "F Continuous orientation input pattern", bty = "n")
lines(g.x, ff_rate(pre_n/4), lty = 3)
lines(g.x,ff_rate(pre_n/4*3), lty = 3)

## Figure 1C: prior density distribution (von mises distribution, center at -90, -45, 0, 45 degree)
lens_set = c(45,90,135, 0)
lens_names = c("-45", "0", "45","90")
lens_col = c("chartreuse3", "blue3","firebrick3","gold1")
# prior density distribution (von mises distribution, center at 90 degree)
for( i in 1:length(lens_names)){
  lens_name = lens_set[i]# in the code 90 is the horizontal orientation
  lens_ori = pi * lens_name/180 *2 #pi*2/4 *2 # 90
  lens_kappa = 0.5
  if(i == 1){
    test = dvonmises(seq(0,2*pi,by = 0.01),  circular(lens_ori), lens_kappa)*2
    plot(seq(0,pi,by = 0.005)/pi*180-90, test*pi/180, main ="G Prior of orientations under stripe rearing",
         type = "l", col = lens_col[i], xaxt="n", bty = "n",
         ylim = c(0,0.01), xlab = "orientation", ylab = "density")
    axis(side = 1, at = c(-90,-45,0,45), labels = c(-90,-45,0,45))
  }else{
    test = dvonmises(seq(0,2*pi,by = 0.01),  circular(lens_ori), lens_kappa)*2
    lines(seq(0,pi,by = 0.005)/pi*180-90, test*pi/180, col = lens_col[i])
  }
  #legend("topleft",legend = lens_names, col = lens_col, lty = 1)
  
}

## Figure 1F: Distribution of preferred orientation before and after stripe rearing
repeat_n = 25
hist_n = 9
hist_bin_width = 20
hist_record_init = matrix(0, nrow = hist_n, ncol = repeat_n)
hist_record_final = matrix(0, nrow = hist_n, ncol = repeat_n)
## -90 degree
load("data/lens_results_uniforminit_0.RData")
for(repeat_i in 1:repeat_n){
  results = results_list[[repeat_i]]
  # init
  hist_record_init[,repeat_i] = hist((results$analyze_results.init$tuninge[1,]+10) %% 180 , breaks = seq(0,180,length.out = 10), plot = F)$counts/ne
  # final
  hist_record_final[,repeat_i] = hist((results$analyze_results$tuninge[1,]+10) %% 180 , breaks = seq(0,180,length.out = 10), plot = F)$counts/ne
}
hist_sd_init = apply(hist_record_init, MARGIN = 1, FUN = sd)
hist_sd_final = apply(hist_record_final, MARGIN = 1, FUN = sd)
barCenters = barplot(rowMeans(hist_record_final), ylim = c(0,0.2), col = "gold1", xaxt = "n",width = 1, space  = 0,
                     xlab = "Preferred Orientation", ylab = "Percentage", main = "H. Distrbution of PO after stripe rearing")
axis(side = 1, at = c(0.5,2,3.5, 5, 6.5, 8), labels = c(-90,-60,-30,0,30,60))
barplot(rowMeans(hist_record_init), ylim = c(0,0.3), border = "grey50", col = rgb(0,0,0,0),width = 1, space  = 0, add = T)
text(1, 0.18, labels = "-90",cex = 1.5)

## -45 degree
load("data/lens_results_uniforminit_45.RData")
for(repeat_i in 1:repeat_n){
  results = results_list[[repeat_i]]
  # init
  hist_record_init[,repeat_i] = hist((results$analyze_results.init$tuninge[1,]+10) %% 180 , breaks = seq(0,180,length.out = 10), plot = F)$counts/ne
  # final
  hist_record_final[,repeat_i] = hist((results$analyze_results$tuninge[1,]+10) %% 180 , breaks = seq(0,180,length.out = 10), plot = F)$counts/ne
}

hist_sd_init = apply(hist_record_init, MARGIN = 1, FUN = sd)
hist_sd_final = apply(hist_record_final, MARGIN = 1, FUN = sd)
barCenters = barplot(rowMeans(hist_record_final), ylim = c(0,0.2), col = "chartreuse3", xaxt = "n",width = 1, space  = 0,
                     xlab = "Preferred Orientation")
axis(side = 1, at = c(0.5,2,3.5, 5, 6.5, 8), labels = c(-90,-60,-30,0,30,60))
barplot(rowMeans(hist_record_init), ylim = c(0,0.3), border = "grey50", col = rgb(0,0,0,0),width = 1, space  = 0, add = T)
text(1, 0.18, labels = "-45",cex = 1.5)

## 0 degree Horizontal 
load("data/lens_results_uniforminit_90.RData")
for(repeat_i in 1:repeat_n){
  results = results_list[[repeat_i]]
  # init
  hist_record_init[,repeat_i] = hist((results$analyze_results.init$tuninge[1,]+10) %% 180 , breaks = seq(0,180,length.out = 10), plot = F)$counts/ne
  # final
  hist_record_final[,repeat_i] = hist((results$analyze_results$tuninge[1,]+10) %% 180 , breaks = seq(0,180,length.out = 10), plot = F)$counts/ne
}

hist_sd_init = apply(hist_record_init, MARGIN = 1, FUN = sd)
hist_sd_final = apply(hist_record_final, MARGIN = 1, FUN = sd)
barCenters = barplot(rowMeans(hist_record_final), ylim = c(0,0.2), col = "blue3", xaxt = "n",width = 1, space  = 0,
                     xlab = "Preferred Orientation")
axis(side = 1, at = c(0.5,2,3.5, 5, 6.5, 8), labels = c(-90,-60,-30,0,30,60))
barplot(rowMeans(hist_record_init), ylim = c(0,0.3), border = "grey50", col = rgb(0,0,0,0),width = 1, space  = 0, add = T)
text(1, 0.18, labels = "0",cex = 1.5)

## 45
load("data/lens_results_uniforminit_135.RData")
for(repeat_i in 1:repeat_n){
  results = results_list[[repeat_i]]
  # init
  hist_record_init[,repeat_i] = hist((results$analyze_results.init$tuninge[1,]+10) %% 180 , breaks = seq(0,180,length.out = 10), plot = F)$counts/ne
  # final
  hist_record_final[,repeat_i] = hist((results$analyze_results$tuninge[1,]+10) %% 180 , breaks = seq(0,180,length.out = 10), plot = F)$counts/ne
}

hist_sd_init = apply(hist_record_init, MARGIN = 1, FUN = sd)
hist_sd_final = apply(hist_record_final, MARGIN = 1, FUN = sd)
barCenters = barplot(rowMeans(hist_record_final), ylim = c(0,0.2), col = "firebrick3", xaxt = "n",width = 1, space  = 0,
                     xlab = "Preferred Orientation")
axis(side = 1, at = c(0.5,2,3.5, 5, 6.5, 8), labels = c(-90,-60,-30,0,30,60))
barplot(rowMeans(hist_record_init), ylim = c(0,0.3), border = "grey50", col = rgb(0,0,0,0),width = 1, space  = 0, add = T)
text(1, 0.18, labels = "45",cex = 1.5)
