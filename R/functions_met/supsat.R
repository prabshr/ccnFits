Copyright 2017 Prabhakar Shrestha
# supsat.R is part of ccnFits.

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

# Saturation partial pressure based on Clausius Clapeyron equation
# Temperature in [K] , esat is in mb
# Also saturation specific humidity is estimate given pressure
# P is in [mb]

supsat = function(T,P) {
eps  = 0.622            
T    = T-273.14                            # [C]

esat = 6.112*exp(17.67*T/(T+243.5))     # [mb]

qvsat= 1000*eps*esat/P                  # [g/kg]

xsat  = c(esat, qvsat)
return(xsat)
}
