Copyright 2017 Prabhakar Shrestha
# wtpFit.R is part of ccnFits.

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

source("$HOME/ccnFits/R/functions_met/Gdtp.R")
source("$HOME/ccnFits/R/functions_met/supsat.R")

xdata = c(273 ,284, 294, 300)
ydata = Gdtp(xdata,900)



# look at it
plot(xdata,ydata)
# some starting values
p1 = 1
p2 = 0.2

# do the fit
fit = nls(ydata ~ p1*cos(p2*xdata) + p2*sin(p1*xdata), start=list(p1=p1,p2=p2))

# summarise
summary(fit)
new = data.frame(xdata = seq(min(xdata),max(xdata),len=200))
lines(new$xdata,predict(fit,newdata=new))
