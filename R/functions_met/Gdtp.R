# Gdtp.R is part of ccnFits.

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

Gdtp = function(T,P) {
# T is air temperature, [K], and P is pressure, [hPa]
# G has the final units of [cm2/s]
# All units are in cgs


# CONSTANTS
Rd  = 287.0e4     #gas constant of dry air, [erg/g/K]
Rv  = 461.5e4     #gas constatn of water vapor, [erg/g/K]
rhow= 1.          #density of water, [g/cm3]
cpv = 1870.0e4    #sp. heat capacity at constant pressure, [erg/g/K]
c   = 4187.0e4    #sp. heat capacity of liquid water, [erg/g/K]
Lvo = 2.5e010     #latent heat of vaporization at 0 C, [erg/g]
To  = 273.15      #reference temperature, [K]
Po  = 1000.       #reference pressure [hPa]

# DERIVED
xsat = supsat(T,P)               #esat is saturated vapor pressure, [hPa]
esat  = xsat[1]
qvsat = xsat[2]
Lvt =Lvo -(c-cpv)*(T-To)                 #latent heat of vaporization, [erg/g]
Dv  = 0.211e4*((T/To)^1.94)*(Po/P)/1e4   #cm2/s Hall and Prauppacher (1976). Pg,.503 book
ka = 1.4132E7*(1.49628E-5*T^1.5)/(T+120) #erg/cm/s/K from CLCM

# ESTIMATE
esat=esat*1000.                              #Need to convert pressure into Pascal,into g/cm/s2

Fd= rhow*Rv*T/(esat*Dv)                      #associated with diffusion
Fk= rhow*(Lvt/(ka*T))*((Lvt/(Rv*T))-1)       #associated with conduction
G = 1/(Fd+Fk)                                #[cm2/s]
return(G)
}
