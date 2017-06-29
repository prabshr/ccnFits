Copyright 2017 Prabhakar Shrestha
# Psimet.R is part of ccnFits.

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

Psimet = function(T,P) { 

# T is air temperature, [K], and P is pressure, [hPa]
# Units: Psi1 [cm-1], Psi2 [ ]
# All units are in cgs

# CONSTANTS
Rd  = 287.0e4           # gas constant of dry air, [erg/g/K]
Rv  = 461.5e4           # gas constatn of water vapor, [erg/g/K]
rhow= 1.                # density of water, [g/cm3]
eps = 0.622             # mol. wt. ratio, Rd/Rv;
cp  = 1005.0e4          # sp. heat of air, [erg/g/K]
cpv = 1870.0e4          # sp. heat capacity at constant pressure, [erg/g/K]
c   = 4187.0e4          # sp. heat capacity of liquid water, [erg/g/K]
Lvo = 2.5e010           # latent heat of vaporization at 0 C, [erg/g]
To  = 273.15            # reference temperature, [K]
g   = 981.0             # gravitational accelearation , [cm/s2]

# DERIVED
xsat = supsat(T,P)     # esat is saturated vapor pressure, [hPa]
esat = xsat[1]
qvsat = xsat[2]
Lvt =Lvo -(c-cpv)*(T-To)       # latent heat of vaporization, [erg/g]

# ESTIMATE
Psi1 = (g/(T*Rd))*(eps*Lvt/(cp*T)-1)

Psi2 = P/(eps*esat) + (eps*Lvt^2)/(Rd*cp*T^2)

Psi  = c(Psi1, Psi2)

return(Psi)
}
