Copyright 2017 Prabhakar Shrestha
# hyprsr2.R is part of ccnFits.

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

hyprsr2 = function(A,B,C,Z) {
#Evaluates hypergeometric function via series-real arguments
#       REAL A,B,C,Z
#       INTEGER N,NSV
#       REAL AA,BB,CC,ZZ,FAC,TEMP,SRES,EPS
EPS=1.E-6
FAC=1.0
TEMP=FAC
AA=A
BB=B
CC=C
for (N in 1:100) {
    NSV=N
    FAC=((AA*BB)/CC)*FAC*Z/N
    SRES=TEMP+FAC
    if (abs((SRES-TEMP)/TEMP)<EPS) {
    	HYPRSR2=SRES
    	return(HYPRSR2)
    	break
    	 #GOTO 20
    } 
    TEMP=SRES
    AA=AA+1.
    BB=BB+1.
    CC=CC+1.
}
#disp( 'HYPRSR2 FAILS TO CONVERGE  A,B,C,Z='), disp([A,B,C,Z])
#       STOP
#    20 HYPRSR2=SRES  
#      PRINT *, 'HYPRSR2  A,B,C,Z,N,SRES=', A,B,C,Z,NSV,SRES
#       RETURN
#       END
return(HYPRSR2)
}
