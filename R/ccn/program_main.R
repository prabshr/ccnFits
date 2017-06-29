Copyright 2017 Prabhakar Shrestha 
# program_main.R is part of ccnFits.

#    ccnFits is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.

#    ccnFits is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.

#    You should have received a copy of the GNU General Public License
#    along with ccnFits.  If not, see <http://www.gnu.org/licenses/>.

# SGN--program_main.R
# Prabhakar Shrestha , pshrestha@uni-bonn.de
# Description:
# Generate polynomial fits for Nccn as a function of (w,T) for fixed pressures (Pfix);
# Functions Used:  Psimet, Gdtp, cohard_smax, cohard                      %
# Output:              mat files, pfits0X, X=0, 1, ... dimsize(fixedpressure)           
# All units are in cgs
# Piece-wise linear model fitting due to non-linear dependency on wind speed.

source("$HOME/ccnFits/R/functions_met/Psimet.R")
#
source("$HOME/ccnFits/R/functions_met/Gdtp.R")
source("$HOME/ccnFits/R/functions_met/supsat.R")
#
source("$HOME/ccnFits/R/functions_met/cohard.R")
source("$HOME/ccnFits/R/functions_met/cohard_smax.R")
source("$HOME/ccnFits/R/functions_met/hyprgeo2.R")
source("$HOME/ccnFits/R/functions_met/hyprsr2.R")
source("$HOME/ccnFits/R/functions_met/cgamaf.R")

# CCN spectra, C=x(1), [cm-3]
#x    = c(3270, 1.56, 0.7, 136) # CCNType3 (Default)
x    = c(22489, 2.26, 1.41, 43.76) #May 15 2009
x = c(54106, 2.84, 1.83, 35.50)

# Discrete Fits at:
fixedpressure     = c(1000,980,960,940,920,900)                       # [mb]
fixedtemperature  = c(-20, -15, -10, -5, 0, 5, 10, 15, 20, 25, 30)    # [C]
fixedwvel         = c(seq(0.01,1,0.01), seq(1.01,10,0.1), seq(10.1,100,1), seq(101,1000,10),seq(1010,10000,100))

#fixedwvel         =  seq(0.01,1,0.01)                                  # [cm/s]
#fixedwvel         =  seq(1,10,0.1)                                  # [cm/s]
#fixedwvel         = seq(10,100,1)                                     # [cm/s]
#fixedwvel         = seq(100,1000,10)                                  # [cm/s]
#fixedwvel         = seq(1000,10000,100)                               # [cm/s]

#

#Create Matrix to Store Data for Fits
nwvel = length(fixedwvel)
nwtem = length(fixedtemperature)
nwprs = length(fixedpressure)

nrows = nwvel*nwtem*nwprs
#ncols = 4
#Na    = matrix(rep(1,ncols*nrows), nrow = nrows, ncol = ncols)
Nv  = vector(mode="numeric", length=nrows)
Sv  = vector(mode="numeric", length=nrows)
#
Wv  = vector(mode="numeric", length=nrows)
Pv  = vector(mode="numeric", length=nrows)
Tv  = vector(mode="numeric", length=nrows)
#
Wv2 = vector(mode="numeric", length=nrows)
Pv2  = vector(mode="numeric", length=nrows)
Tv2  = vector(mode="numeric", length=nrows)
#
Wv3  = vector(mode="numeric", length=nrows)
Pv3  = vector(mode="numeric", length=nrows)
Tv3  = vector(mode="numeric", length=nrows)
#
Wv4  = vector(mode="numeric", length=nrows)
Pv4  = vector(mode="numeric", length=nrows)
Tv4  = vector(mode="numeric", length=nrows)

count = 1

for (iprs in 1:length(fixedpressure)) {      # PRESSURE LOOP [mb]

# Constant
Rd   = 287.0e4                               # Gas constant of dry air, [erg/g/K]
rhow = 1.                                    # Density of water, [g/cm3]
P    = fixedpressure[iprs]                   # mb
B    = beta(x[2]/2,3/2)                      # Beta function
w    = fixedwvel                             # VECTOR of vertical velocity for fits.

for (T in fixedtemperature) {                # TEMPERATURE LOOP [C]
    T=T+273.15
# Derived
rhoa = 1000*P/(Rd*T)

# Estimated parameters
Psi  = Psimet(T,P)                           # [cm-1] [ ]
Psi1 = Psi[1]
Psi2 = Psi[2]
G    = Gdtp(T,P)                             # [cm2/s]

# Max super-saturation iteration
for (iter in 1:5) {                          # ITERATIVE LOOP
	
    if (iter == 1) {
    	hyprgeo = rep(1,length(w))
    } else {
        hyprgeo = cohard_smax(x,smax)
        
    }
    numer = rhoa*(Psi1*w)^1.5
    denom = hyprgeo*(2*x[2]*(x[1])*pi*rhow*Psi2*(G^1.5)*B)
    smax  = (1e4*numer/denom)^(1/(x[2]+2))                   #1e4 factor converts from [fraction] to [%]
    Nccn = cohard(x,smax)
}                                            # ITERATIVE LOOP

istart = count
iend   = count+nwvel-1
Nv[istart:iend] = log10(Nccn*1E6)
Sv[istart:iend] = log10(smax)
#
Wv[istart:iend] = log10(w)
Tv[istart:iend] = T
Pv[istart:iend] = P
#
Wv2[istart:iend] = log10(w)*log10(w)
Tv2[istart:iend] = T*T
Pv2[istart:iend] = P*P
#
Wv3[istart:iend] = log10(w)*log10(w)*log10(w)
Tv3[istart:iend] = T*T*T
Pv3[istart:iend] = P*P*P
#
Wv4[istart:iend] = log10(w)*log10(w)*log10(w)*log10(w)
Tv4[istart:iend] = T*T*T*T
Pv4[istart:iend] = P*P*P*P
count = count + nwvel
}                                            # TEMPERATURE LOOP [C]

}                                            # PRESSURE LOOP [mb]

# Save vectors in data frame
#df = data.frame(Nv,Sv,Wv,Tv,Pv)
# Use multiple linear regression for Nccn = f(w,T,P)
#fitNv = lm(Nv ~ (Wv + Tv + Pv)^3, data=df)
# Use multiple linear regression for Smax = f(w,T,P)
#fitSv = lm(Sv ~ (Wv + Tv + Pv)^3, data=df)

#
# Save vectors in data frame
df = data.frame(Nv,Sv,Wv,Tv,Pv,Wv2,Tv2,Pv2,Wv3,Tv3,Pv3,Wv4,Tv4,Pv4)
#
# Use multiple linear regression for Nccn = f(w,T,P)
fitNv = lm(Nv ~ (Tv4 + Tv3 + Tv2 + Tv)*Wv4 + (Tv4 + Tv3 + Tv2 + Tv)*Wv3 + (Tv4 + Tv3 + Tv2 + Tv)*Wv2 + (Tv4 + Tv3 + Tv2 + Tv)*Wv + (Tv4 + Tv3 + Tv2 + Tv) + (Pv2 + Pv)*(Wv4 + Wv3 + Wv2 + Wv)*(Tv3 + Tv2 + Tv) - Tv:Wv3 - Tv:Pv - Tv:Wv4:Pv - Tv:Wv3:Pv -Tv:Wv2:Pv -Tv:Wv:Pv , data=df)
# Use multiple linear regression for Smax = f(w,T,P)
fitSv = lm(Sv ~ (Tv4 + Tv3 + Tv2 + Tv)*Wv4 + (Tv4 + Tv3 + Tv2 + Tv)*Wv3 + (Tv4 + Tv3 + Tv2 + Tv)*Wv2 + (Tv4 + Tv3 + Tv2 + Tv)*Wv + (Tv4 + Tv3 + Tv2 + Tv) + (Pv2 + Pv)*(Wv4 + Wv3 + Wv2 + Wv)*(Tv3 + Tv2 + Tv) - Tv:Wv3 - Tv:Pv - Tv:Wv4:Pv - Tv:Wv3:Pv -Tv:Wv2:Pv -Tv:Wv:Pv , data=df)
#

#diagnostic plots
#layout(matrix(c(1,2,3,4),2,2))
#plot(fit)

#ALGORITHM FOR FORTRAN EQUATION
#
writeLines(paste("CCN SPECTRA", x[1]," ",x[2]," ",x[3]," ",x[4],"\n"))
writeLines(paste("%s/:/*/g","\n"))
writeLines(paste("%s/Tv/T/g","\n"))
writeLines(paste("%s/Pv/P/g","\n"))
writeLines(paste("%s/Wv/W/g","\n"))

writeLines("fitNv\n")
attrs = attributes(coef(fitNv))
coefs = as.matrix(coef(fitNv),dimnames="")
for (i in seq(1,length(coefs),by=5))
{
 if (i == 1) {
   lhs="1E-6*10**("
 } else {
   lhs=""
 }
 ji = i+5-1
 if (ji > length(coefs)) {
    ji = length(coefs)
 }
 for (j in i:ji)
{
 if (coefs[j] < 0 ) {
   lhs <- cat(lhs,coefs[j],"*",attrs$names[j])
 } else {
   lhs <- cat(lhs,"+",coefs[j],"* ", attrs$names[j])
 }
}
 if (ji < length(coefs)) {
   writeLines(paste(lhs,"&\n"))
 } else {
   writeLines(paste(lhs,")\n"))
 }
}
#
#
writeLines("fitSv\n")
attrs = attributes(coef(fitSv))
coefs = as.matrix(coef(fitSv),dimnames="")
for (i in seq(1,length(coefs),by=5))
{
 if (i == 1) {
   lhs="10**("
 } else {
   lhs=""
 }
 ji = i+5-1
 if (ji > length(coefs)) {
    ji = length(coefs)
 }
 for (j in i:ji)
{
 if (coefs[j] < 0 ) {
   lhs <- cat(lhs,coefs[j],"*",attrs$names[j])
 } else {
   lhs <- cat(lhs,"+",coefs[j],"* ", attrs$names[j])
 }
}
 if (ji < length(coefs)) {
   writeLines(paste(lhs,"&\n"))
 } else {
   writeLines(paste(lhs,")\n"))
 }
}
#
#
# PLOTTING/VALIDATION
Tval = Tv
Wval = Wv
Pval = Pv

Na =1E-6*(10^(predict(fitNv,list(Tv=Tval,Wv=Wval,Pv=Pval,Tv2=Tval*Tval,Wv2=Wval*Wval,Pv2=Pval*Pval,Tv3=Tval*Tval*Tval,Wv3=Wval*Wval*Wval,Pv3=Pval*Pval*Pval,Tv4=Tval*Tval*Tval*Tval,Wv4=Wval*Wval*Wval*Wval,Pv4=Pval*Pval*Pval*Pval))))
Sa =10^(predict(fitSv,list(Tv=Tval,Wv=Wval,Pv=Pval,Tv2=Tval*Tval,Wv2=Wval*Wval,Pv2=Pval*Pval,Tv3=Tval*Tval*Tval,Wv3=Wval*Wval*Wval,Pv3=Pval*Pval*Pval,Tv4=Tval*Tval*Tval*Tval,Wv4=Wval*Wval*Wval*Wval,Pv4=Pval*Pval*Pval*Pval)))

#coefs = coef(fitNv)
#Na = 1E-6*(10^(coefs['(Intercept)'] + coefs['Wv']*Wval + coefs['Tv']*Tval + coefs['Pv']*Pval + coefs['Wv:Tv']*Wval*Tval + coefs['Wv:Pv']*Wval*Pval + coefs['Tv:Pv']*Tval*Pval + coefs['Wv:Tv:Pv']*Wval*Tval*Pval))
#

#coefs = coef(fitSv)
#Sa = (10^(coefs['(Intercept)'] + coefs['Wv']*Wval + coefs['Tv']*Tval + coefs['Pv']*Pval + coefs['Wv:Tv']*Wval*Tval + coefs['Wv:Pv']*Wval*Pval + coefs['Tv:Pv']*Tval*Pval + coefs['Wv:Tv:Pv']*Wval*Tval*Pval))

#Debug plot
Nc = 1E-6*10^Nv
Sc = 10^Sv
plot(Sc, Nc, log="xy",col="red",ylim=range(c(Nc,Na)), xlim=range(c(Sc,Sa)))
par(new = TRUE)
plot(Sa,Na,log="xy",col="green", axes= FALSE, ylim=range(c(Nc,Na)),xlim=range(c(Sc,Sa)))

#Clear Memory
rm(list=ls())

