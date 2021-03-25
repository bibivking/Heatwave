# MMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM
# Author  : Anna Ukkola
# Version : 1.0 (22.03.2021)"
# WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

library(raster)

#Year to read in
years <- list(2019)


#File paths
path_gw <- "/g/data/w35/mm3972/model/cable/runs/AWAP_pumping/CTL-30x3+30yr/outputs-daily"

path_nogw <- "/g/data/w35/mm3972/model/cable/runs/AWAP_pumping/FREE_DRAIN/outputs-spinup30x3yr"


#Get iveg
iveg <- raster(paste0(path_gw, "/cable_out_2001_SE_Aus.nc"), varname="iveg")
iveg[iveg >4] <- NA


#Get HW index
hw_ind <- brick("/g/data/w35/mm3972/scripts/ehfheatwaves/HW_Event_Indicator_1970-2019.nc")
start <- which(getZ(hw_ind) == "2019-01-01")
hw_ind <- hw_ind[[start:(start+30)]]
hw_ind[hw_ind < 1] <- NA


#function to calculate VPD -----------

qair_to_vpd <- function(Qair_ts, Tair_ts) {

  # unit convert
  Tair_ts  = Tair_ts -  273.15

  #Calculate VPD
  Purf     = 1000.*100.
  Es = 100.0 * 6.112 * exp((17.67 * Tair_ts) / (243.5 + Tair_ts))

  #vapor pressure
  Ea = (Qair_ts * Purf) / (0.622 + (1.0 - 0.622) * Qair_ts)
  VPD = (Es - Ea) * 0.001 # in kPa

  return(VPD)
}

#-------------------



##########
### GW ###
##########

deltat <- vector()
fwsoil <- vector()
vpd    <- vector()
tair   <- vector()

#Get data for Janury 2019

for (l in 1:length(years)) {

  for(y in 1){#length(years[[l]])) {

    file <- paste0(path_gw, "/cable_out_", years[[l]][y], "_SE_Aus.nc")

    leafT_data  <- mask(mask(brick(file, varname="VegT"), iveg)[[1:31]], hw_ind)
    fwsoil_data <- mask(mask(brick(file, varname="Fwsoil"), iveg)[[1:31]], hw_ind)
    qair_data   <- mask(mask(brick(file, varname="Qair"), iveg)[[1:31]], hw_ind)
    tair_data   <- mask(mask(brick(file, varname="Tair"), iveg)[[1:31]], hw_ind)


    #Get values
    deltat  <- append(deltat, as.vector(values(leafT_data))[which(!is.na(as.vector(values(leafT_data))))])
    fwsoil  <- append(fwsoil, as.vector(values(fwsoil_data))[which(!is.na(as.vector(values(fwsoil_data))))])
    vpd     <- append(vpd, as.vector(values(qair_data))[which(!is.na(as.vector(values(qair_data))))])
    tair    <- append(tair, as.vector(values(tair_data))[which(!is.na(as.vector(values(tair_data))))])

  }

}


#Convert Qair to VPD
vpd <- qair_to_vpd(vpd, tair)



#############
### no GW ###
#############

deltat_nogw <- vector()
fwsoil_nogw <- vector()

for (l in 1:length(years)) {

  for(y in 1){#length(years[[l]])) {

    file <- paste0(path_nogw, "/cable_out_", years[[l]][y], "_SE_Aus.nc")

    leafT_data  <- mask(mask(brick(file, varname="VegT"), iveg)[[1:31]], hw_ind)
    fwsoil_data <- mask(mask(brick(file, varname="Fwsoil"), iveg)[[1:31]], hw_ind)

    #Get values
    deltat_nogw  <- append(deltat_nogw, as.vector(values(leafT_data))[which(!is.na(as.vector(values(leafT_data))))])
    fwsoil_nogw  <- append(fwsoil_nogw, as.vector(values(fwsoil_data))[which(!is.na(as.vector(values(fwsoil_data))))])

  }

}



############
### Plot ###
############

#Difference GW - no GW

deltaT_difference <- deltat - deltat_nogw

fwsoil_difference <- fwsoil - fwsoil_nogw

#Convert Tair from K to C
tair <- tair - 273.15


#Initialise plot
png("./plots/Fig7_DeltaT_vs_fwsoil_vpd_tair.png",
    height=2.70, width=8, units="in", res=400)

par(mai=c(0.55, 0.2, 0.2, 0.2))
par(omi=c(0, 0.3, 0.1, 0.1))

par(mfcol=c(1,3))


### Fwsoil ###

#Set density colours
dens_cols  <- densCols(fwsoil_difference, deltaT_difference,
                       colramp=colorRampPalette(c("#deebf7", "#3182bd")),
                       nbin=1000)

#Plot
plot(fwsoil_difference, deltaT_difference, pch=20, cex=0.3, col=dens_cols,
     ylab="", xlab="")

#Calculate correlation
cor <- cor.test(fwsoil_difference, deltaT_difference)$estimate
mtext(side=1, adj=0.05, line=-2, paste0("r = ", round(cor, digits=2)), cex=0.8)

#Labels
mtext(side=1, "Δβ (-)", line=2.5, cex=0.8)
mtext(side=2, "ΔT difference (°C)", line=2.5, cex=0.8)
mtext(side=3, "a)", adj=0.05, line=-1.5, cex=0.8)


### VPD ###

#Set density colours
dens_cols  <- densCols(vpd, deltaT_difference,
                       colramp=colorRampPalette(c("#deebf7", "#3182bd")),
                       nbin=1000)

#Plot
plot(vpd, deltaT_difference, pch=20, cex=0.3, col=dens_cols,
     ylab="", xlab="")

#Calculate correlation
cor <- cor.test(vpd, deltaT_difference)$estimate
mtext(side=1, adj=0.05, line=-2, paste0("r = ", round(cor, digits=2)), cex=0.8)

#Labels
mtext(side=1, "D (kPa)", line=2.5, cex=0.8)
mtext(side=3, "b)", adj=0.05, line=-1.5, cex=0.8)



### Tair ###

#Set density colours
dens_cols  <- densCols(tair, deltaT_difference,
                       colramp=colorRampPalette(c("#deebf7", "#3182bd")),
                       nbin=1000)

#Plot
plot(tair, deltaT_difference, pch=20, cex=0.3, col=dens_cols,
     ylab="", xlab="")

#Calculate correlation
cor <- cor.test(tair, deltaT_difference)$estimate
mtext(side=1, adj=0.05, line=-2, paste0("r = ", sprintf("%.2f", cor)), cex=0.8)

#Labels
mtext(side=1, "Tair (°C)", line=2.5, cex=0.8)
mtext(side=3, "c)", adj=0.05, line=-1.5, cex=0.8)


dev.off()
