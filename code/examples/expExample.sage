load("../D-finite-core.sage")


# Define polynomial ring, ore algebra, and variables for the system.
Pols.<x> = PolynomialRing(QQ); Dops.<Dx> = OreAlgebra(Pols)
(x1, x2, y1) = var('x1 x2 y1')
# x1, x2 : variables for polynomial system
# y1 : variables for D-finite functions
A.<x1, x2, y1> = PolynomialRing(QQ)



expression = Dx - 4
explist = [expression]


sys = [x2 - 0.0183, x2 - y1]
varlist = [x1, [x2, x1]]



# approximate solution
point = (-1.000, 0.018316)


# define initial values at t=0
ini1 = [1]
iniAt = 0
inilist = [ini1]



print("This example shows the comparison with alphaCertified. The first table shows the gamma value depends on different radii, and the second graph shows the list plots from the data in the first table, and the actual value of 'alphaCertified'")
print table([(point,i*0.05,(gammaValue(sys,varlist,explist,inilist,i*0.05, point,iniAt,1e-40)).n()) for i in range(1,20)], header_row=["point",r"radius", r"$gamma$"], frame=True)

listOfPts = [(i*0.025, (gammaValue(sys,varlist,explist,inilist,i*0.025, point,iniAt,1e-40)).n()) for i in range(1,40)]
P = line([(0,84),(1,84)], color='red', legend_label='alphaCertified')
Q = line(listOfPts, color='blue')
show(P+Q)
