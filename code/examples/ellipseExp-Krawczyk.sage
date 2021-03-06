load("../D-finite-core.sage")



# Define polynomial ring, ore algebra, and variables for the system.
Pols.<x> = PolynomialRing(QQ); Dops.<Dx> = OreAlgebra(Pols)
(x1, x2, x3, x4, x5, x6, x7, a, b, c, d,  y1, y2, y3, y4) = var('x1 x2 x3 x4 x5 x6 x7 a b c d y1 y2 y3 y4')
# x1 ~ d : variables for polynomial system
# y1, y2, y3, y4 : variables for D-finite functions
A.<x1, x2, x3, x4, x5, x6, x7, a, b, c, d, y1, y2, y3, y4> = PolynomialRing(QQ)




expression1 = (x-x^3)*Dx^2 + (1-x^2)*Dx + x
expression2 = (-x^4+x^2)*Dx^2 + (-3*x^3 + x)*Dx - 1
explist = [expression1,expression1, expression2, expression2]



sys = [x3+x5*(2*x1), 2*x4+x6*(2*x2), x1+ 2*x3*x5+4*x7*c, 2*x2+8*x6*x4+8*x7*d, x3^2-1+x1^2, 4*x4^2-4+x2^2, a+2*b-17,a-y1,b-y2,c-y3,d-y4]
varlist = [x1,x2,x3,x4,x5,x6, x7, [a, x3], [b, x4], [c, x3], [d, x4]]



# approximate solution
point = (RBF((0.8337852,0.8337854)), RBF((1.5601132,1.5601134)), RBF((0.5520887,0.5520889)), RBF((0.6257088,0.6257090)), RBF((-0.3310738,-0.3310736)), RBF((-0.4010664,-0.4010662)), RBF((0.0590726,0.0590728)), RBF((5.7728450,5.7728452)), RBF((5.6135773,5.6135775)), RBF((-1.9815465,-1.9815463)), RBF((-2.3543458,-2.3543456)))
varlist = [x1,x2,x3,x4,x5,x6, x7, [a, x3], [b, x4], [c, x3], [d, x4]]



# define initial values at t=1/2
ini1 = [4*elliptic_ec(1/4), 8*(elliptic_ec(1/4)-elliptic_kc(1/4))]
ini2 = [8*(elliptic_ec(1/4)-elliptic_kc(1/4)),4*(-4*(elliptic_ec(1/4)-elliptic_kc(1/4))+2*(2*(elliptic_ec(1/4)-elliptic_kc(1/4))-(8/3)*(elliptic_ec(1/4)-(3/4)*elliptic_kc(1/4))))]
iniAt = 1/2
inilist = [ini1,ini1,ini2,ini2]





K = krawczykoperator(sys, varlist, explist, inilist, point, iniAt,1e-40)
print("This example shows the result of Krawczyk method-based test certification on the maximization problem about two ellipses")
print("box = "+ point.__repr__())
print("Krawczyk operator = "+ (K).__repr__())
for i in range(0, len(sys)):
    if (K[i][0].endpoints()[0] < point[i].endpoints()[0] or K[i][0].endpoints()[1] > point[i].endpoints()[1]):
        print "false"
print("The Krawczyk operator is contained in the interval. The unique solution is contained in the box.")
