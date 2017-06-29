# cgamaf.R is part of ccnFits.

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

# Generalized gamma function evaluation, not restricted to X.GE.1
#       REAL X,XP1,PI,GMALN
cgamaf = function(X) {
	
	PI=3.1415926535898
	
	if (X >= 1.0) {
		GMALN=log(gamma(X))
		CGAMAF=exp(GMALN)
		} else if (X > 0.0) {
		XP1=X+1.0
		GMALN=log(gamma(XP1))
		CGAMAF=exp(GMALN)/X
		} else if (X < 0.0) {
		XP1=1.0-X
		GMALN=log(gamma(XP1))
		CGAMAF=PI/(exp(GMALN)*sin(XP1*PI))
		} else {
		print("cgamaf, No code written")
		}
	return(CGAMAF)
}
