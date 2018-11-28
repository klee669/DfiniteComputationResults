load("../D-finite-core.sage")



# Define polynomial ring, ore algebra, and variables for the system.
Pols.<x> = PolynomialRing(QQ); Dops.<Dx> = OreAlgebra(Pols)
(x1, x2, x3, x4, x5, y1, y2) = var('x1 x2 x3 x4 x5 y1 y2')
# x1, x2, x3, x4, x5 : variables for polynomial system
# y1, y2, y3, y4 : variables for D-finite functions
A.<x1, x2, x3, x4, x5, y1, y2> = PolynomialRing(QQ)



exp1 = Dx^2 + 2*x*Dx
exp2 = x^2 * Dx^2 + x * Dx + (x^2 - 81)  # expression for an error function.
explist = [exp1, exp2]



# regular singular points of bessel functions
regularSingularPoints = (expForDeri(expForDeri(exp2))).leading_coefficient().roots(AA,multiplicities = False)


sys = [x1^2 + x2^2 - 61, 2 * x4 * x5 - 11 , x3 - (1/2) * x5 - x1 , x4 - y1, x5 - y2] 
varlist = [x1, x2, x3, [x4, x3], [x5, x2]]



ini1 = [0, 2/sqrt(pi)]
ini2 = [-20643840/pi, 1/185794560 - 7129/(468202291200 * pi) + euler_gamma/(92897280 * pi) - log(2)/(92897280 * pi)]
iniAt = 0
inilist = [ini1, ini2]


point = (627898967/100000000,464481310/100000000,-38382379/100000000, -41273856/100000000, -133256269/10000000) # given point (values for x1, x2, x3, x4 and x5)



b = beta(sys, varlist, explist, inilist, point, iniAt)


print ("This example runs alpha theory-based test on the system with error function and Bessel function over different radii")
print ("The radius of convergence of the Bessel function is " + RR(regularSingularPoints[4]).__repr__())
print table([(i*0.1,gammaValue(sys,varlist,explist,inilist,i*0.1*regularSingularPoints[4], point,iniAt,1e-40).n(),(b * gammaValue(sys,varlist,explist,inilist,i*0.1*regularSingularPoints[4], point,iniAt,1e-40)).n()) for i in [0.00001, 0.0001,0.001,0.01,0.1,0.3]], header_row=["radius", r"gamma",r"alpha"], frame=True)


