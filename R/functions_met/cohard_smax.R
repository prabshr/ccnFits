Copyright 2017 Prabhakar Shrestha
# cohard_smax.R is part of ccnFits.

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

cohard_smax = function(x,s) {
# ESTIMATE Nccn given supersaturation and parameters
# Cohard et al. 1998
# All units are in cgs, C [cm-3] s-supersaturation [%]
C    = x[1]
k    = x[2]
mu   = x[3]
beta = x[4]

# THIS USES THE HYPERGEOMETRIC FUNCTION EQ.10 Cohard et al. 1998
hyprgeo = vector(mode="numeric", length=length(s))
for (i in 1:length(s)) {

    hyprgeo[i]=Re(hyprgeo2(mu,k/2,k/2+3/2,-beta*s[i]*s[i]))
    
}
return(hyprgeo)
}
