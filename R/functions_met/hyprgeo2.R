Copyright 2017 Prabhakar Shrestha
# hyprgeo2.R is part of ccnFits.

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

# Generalized hypergeometric function evaluation for real arguments
# Applies transformed series when applicable.  Uses Numerical
# Recipes in FORTRAN as last resort (complex arguments)
#       REAL A,B,C,Z
#       COMPLEX AC,BC,CC,ZC,ANSC,HYPRGEO
#       REAL ZZ,BR1,BR2,HP1,HP2
#       REAL GA,GB,GC

#      PRINT *, 'HYPRGEO2 INPUT  A,B,C,Z=', A,B,C,Z
hyprgeo2 = function(A,B,C,Z) {
	BR1=0.75
	BR2=2.0
	
	if (abs(Z) <= BR1) {
		HYPRGEO2 =hyprsr2(A,B,C,Z)
		} else if ((Z > BR1) & (Z < 1.0)) {
    		ZZ=1.-Z
    		GA=cgamaf(A)
    		GB=cgamaf(B)
    		GC=cgamaf(C)
    		HP1=GC*cgamaf(C-A-B)/(cgamaf(C-A)*cgamaf(C-B))*hyprsr2(A,B,A+B-C+1,ZZ)
    		HP2=GC*cgamaf(A+B-C)/(GA*GB)*(ZZ^(C-A-B))*hyprsr2(C-A,C-B,C-A-B+1,ZZ)
    		HYPRGEO2=HP1+HP2
    		#    PRINT *, 'HYPRGEO2 1-Z TRANSFORM  HP1,HP2,RES=', HP1,HP2,HYPRGEO2
    		} else if ((Z < -BR1) & (Z > -BR2)) {
    		ZZ=1./(1.-Z)
    		GA=cgamaf(A)
    		GB=cgamaf(B)
    		GC=cgamaf(C)
    		HP1=GC*cgamaf(B-A)/(GB*cgamaf(C-A))*(ZZ^A)*hyprsr2(A,C-B,A-B+1,ZZ)
    		HP2=GC*cgamaf(A-B)/(GA*cgamaf(C-B))*(ZZ^B)*hyprsr2(B,C-A,B-A+1,ZZ)
    		HYPRGEO2=HP1+HP2
    		#     PRINT *, 'HYPRGEO2 1/(1-Z) TRANSFORM  HP1,HP2,RES=', HP1,HP2,
    		#     1        HYPRGEO2
    		} else if (Z <= -BR2) {
    		ZP=-Z
    		ZR=1./Z
    		ZPR=1./ZP
    		GA=cgamaf(A)
    		GB=cgamaf(B)
    		GC=cgamaf(C)
    		HP1=GC*cgamaf(B-A)/(GB*cgamaf(C-A))*(ZPR^A)*hyprsr2(A,1.-C+A,1.-B+A,ZR)
    		HP2=GC*cgamaf(A-B)/(GA*cgamaf(C-B))*(ZPR^B)*hyprsr2(B,1.-C+B,1.-A+B,ZR)
    		HYPRGEO2=HP1+HP2
    		#      PRINT *, 'HYPRGEO2 1/Z ALT. TRANSFORM  HP1,HP2,RES=', HP1,HP2,
    		#     1        HYPRGEO2
    		} else if (Z >= 1.0) {
    		AC=complex(A,0.0)
    		BC=complex(B,0.0)
    		CC=complex(C,0.0)
    		ZC=complex(Z,0.0)
    		# CPS    [ANSC]=hyprgeo(AC,BC,CC,ZC);
    		# CPS    HYPRGEO2=real(ANSC);
    		#    disp('Code Needs to be updated:::'
    		} else {
		print("hyprgeo2, no code written")
	}
	 return(HYPRGEO2)
}
