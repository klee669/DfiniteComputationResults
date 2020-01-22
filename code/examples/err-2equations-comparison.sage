# an example of system with D-finite function
# error function.

load("../D-finite-core.sage")



# Define polynomial ring, ore algebra, and variables for the system.
Pols.<x> = PolynomialRing(QQ); Dops.<Dx> = OreAlgebra(Pols)
(x1, x2, x3, x4, y1, y2) = var('x1 x2 x3 x4 y1 y2')
# x1, x2, x3, x4 : variables for polynomial system
# y1, y2 : variables for D-finite functions
A.<x1, x2, x3, x4, y1, y2> = PolynomialRing(QQ)




exp = Dx^2 + 2*x*Dx  # expression for an error function.
explist = [exp, exp]



sys = [x3*x4 - 1/2, x1^2 + x2^2 - 4, x3 - y1, x4 - y2] # system
varlist = [x1, x2, [x3, x1], [x4, x2]]




point1 = (480322/1000000,1941466/1000000,503038/1000000,993961/1000000) # given point (values for x1, x2, x3, and x4)
point2 = (48032/100000,194147/100000,50304/100000,99396/100000) 
point3 = (4803/10000,19415/10000,5030/10000,9940/10000) 
point4 = (480/1000,1941/1000,503/1000,994/1000)
point5 = (48/100,194/100,50/100,99/100)
point6 = (5/10, 19/10, 5/10, 1)
point7 = (1, 2, 1, 1)
point = [point1, point2, point3, point4, point5, point6, point7]




ini = [0, 2/sqrt(pi)]
inilist = [ini, ini]
iniAt = 0

rad = 0.4





interval1 = [RBF(RealInterval(0.480323,0.480321)),RBF(RealInterval(1.941467,1.941465)),RBF(RealInterval(0.503039,0.503037)),RBF(RealInterval(0.993962,0.993960))]
interval2 = [RBF(RealInterval(0.48033,0.48031)),RBF(RealInterval(1.94148,1.94146)),RBF(RealInterval(0.50305,0.50303)),RBF(RealInterval(0.99397,0.99395))]
interval3 = [RBF(RealInterval(0.4804,0.4802)),RBF(RealInterval(1.9416,1.9414)),RBF(RealInterval(0.5031,0.5029)),RBF(RealInterval(0.9941,0.9939))]
interval4 = [RBF(RealInterval(0.481,0.479)),RBF(RealInterval(1.942,1.940)),RBF(RealInterval(0.504,0.502)),RBF(RealInterval(0.995,0.993))]
interval5 = [RBF(RealInterval(0.49,0.47)),RBF(RealInterval(1.95,1.93)),RBF(RealInterval(0.51,0.49)),RBF(RealInterval(1,0.98))]
interval6 = [RBF(RealInterval(0.4,0.6)),RBF(RealInterval(1.8,2)),RBF(RealInterval(0.4,0.6)),RBF(RealInterval(0.9,1.1))]
interval7 = [RBF(RealInterval(0,2)),RBF(RealInterval(1,3)),RBF(RealInterval(0,2)),RBF(RealInterval(0,2))]
interval = [interval1, interval2, interval3, interval4, interval5, interval6, interval7]

k1 = krawczykoperator(sys, varlist, explist, inilist, interval[0],iniAt, 1e-40)
k2 = krawczykoperator(sys, varlist, explist, inilist, interval[1],iniAt, 1e-40)
k3 = krawczykoperator(sys, varlist, explist, inilist, interval[2],iniAt, 1e-40)
k4 = krawczykoperator(sys, varlist, explist, inilist, interval[3],iniAt, 1e-40)
k5 = krawczykoperator(sys, varlist, explist, inilist, interval[4],iniAt, 1e-40)
k6 = krawczykoperator(sys, varlist, explist, inilist, interval[5],iniAt, 1e-40)
k7 = krawczykoperator(sys, varlist, explist, inilist, interval[6],iniAt, 1e-40)
k = [k1,k2,k3,k4,k5,k6,k7]


print ("This example compares Krawczyk method and alpha theory based tests on different accuracies")
print ("alpha theory-based test runs over the radius r =" + rad.__repr__())
print table([(7-i, k[i][0][0].endpoints()[0]>interval[i][0].endpoints()[0] and k[i][0][0].endpoints()[1]<interval[i][0].endpoints()[1]
              and k[i][1][0].endpoints()[0]>interval[i][1].endpoints()[0] and k[i][1][0].endpoints()[1]<interval[i][1].endpoints()[1]
              and k[i][2][0].endpoints()[0]>interval[i][2].endpoints()[0] and k[i][2][0].endpoints()[1]<interval[i][2].endpoints()[1]
              and k[i][3][0].endpoints()[0]>interval[i][3].endpoints()[0] and k[i][3][0].endpoints()[1]<interval[i][3].endpoints()[1],(beta(sys, varlist, explist, inilist, point[i],iniAt) * gammaValue(sys,varlist,explist,inilist,rad, point[i],iniAt,1e-40)).n() < ((13-3*sqrt(17))/4).n())  for i in [0..6]], header_row=["decimal places", r"krawczyk",r"alpha"], frame=True)
