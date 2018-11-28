
########## Description ###########
# SAGE functions for certifying square system with D-finite functions.
# This file consists of 3 parts, which are
# 1. functions for D-finite functions.
# 2. functions for Krawczyk method.
# 3. functions for alpha theory.


########## Setting ##########
# Consider the square system 'F' with n+m variables. The first n equations should be a polynomial of n+m varibles, and the last m equations should be a x_{n+i}- g(x_s(i)) where g is a D-finite function and s(i) is a number between 1 to n+i-1.



# We need ore_algebra package.
import math
from ore_algebra import OreAlgebra



########## functions for D-finite functions ##########



###### D-finite functions ######
# D-finite function is a solution of ordinary linear differential equation with polynomial coefficients.
# 'OreAlgebra' in 'ore_algebra' does computations for D-finite functions using analytic continuation.





#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#





def iniForDeri(expression, ini, iniAt):
    r"""
        Computes the initial values for the derivative of D-finite functions.

    INPUT:

        -  ``expression`` - An element in univariate Ore algebra in 'Dx' over univariate polynomial ring in 'x' over rational field
        -  ``ini`` - A list of numbers

    OUTPUT: A list of numbers

    EXAMPLES::

        sage: import math
        sage: from ore_algebra import OreAlgebra
        sage: Pols.<x> = PolynomialRing(QQ) ; Dops.<Dx> =  OreAlgebra(Pols)
        sage: iniForDeri(Dx^2 + 2*x*Dx, [0,2/sqrt(pi)])
        [2/sqrt(pi), 0]  ### an initial condition for 'Dx^2 + 2*x*Dx + 2'.

    AUTHORS:
    
        - Kisun Lee klee669@gatech.edu (2018-08-21)

    SEEALSO::

        :func:evaluateJacobianInterval
        :func:invmat
        :func:gammaVaue
    

    """

    # 1. extract the information about expression. (degree, leading term)
    degOfExp = order(expression)
    A.<x> = PolynomialRing(QQ)
    leadingTerm = A(expression.coefficients()[-1])


    # 2. make the list of initial values for the derivative.
    # set 'derivativeIni = ini'
    derivativeIni = []
    for i in range(0, len(ini)):
        derivativeIni.append(ini[i])

    # 3. If there is no root on the 'leadingTerm' (D-finite function has no singularity)  compute the new initial value using previous initial values.
    if leadingTerm.subs(x=iniAt) != 0:
        derive = leadingTerm * (Dx ** degOfExp) - expression
        deriveCoefflist = derive.coefficients(sparse=false)
        newInitialValue = 0
        for i in range(0, len(deriveCoefflist)):
            newInitialValue = newInitialValue + (1/leadingTerm.subs(x=iniAt)) * derivativeIni[i] * deriveCoefflist[i].subs(x=iniAt)
        derivativeIni.append(newInitialValue)
        derivativeIni = derivativeIni[1:]


    # 4. If the input  D-finite function has 'regular singular points', then use local basis to compute the asymptotic behavior of the function.
    else:
        localBasis = expression.local_basis_monomials(0)
        for i in range(0, len(localBasis)):
            derivativeIni[i] = derivativeIni[i] * diff(localBasis[i]).subs(x=1)
    return derivativeIni




#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#


def expForDeri(expression):
    r"""
        Computes the recursive expression for the derivative of D-finite functions

    INPUT:

        -  ``expression`` - An element in univariate Ore algebra in 'Dx' over univariate polynomial ring in 'x' over rational field

    OUTPUT: An element in univariate Ore algebra in 'Dx' over univariate polynomial ring in 'x' over rational field

    EXAMPLES::

        sage: import math
        sage: from ore_algebra import OreAlgebra
        sage: Pols.<x> = PolynomialRing(QQ) ; Dops.<Dx> =  OreAlgebra(Pols)
        sage: expForDeri(Dx^2 + 2*x*Dx)
        Dx^2 + 2*x*Dx + 2

    AUTHORS:
    
        - Kisun Lee klee669@gatech.edu (2018-08-21)

    SEEALSO::

        :func:iniForDeri
        :func:evaluateJacobianInterval
        :func:invmat
        :func:gammaValue
    

    """

    
    # 1. When the given expression has no constant term, then it is simple.
    # It is enough to divide by Dx from the left.
    coefflist = expression.coefficients(sparse=False)
    constTerm = coefflist[0]
    A.<x> = PolynomialRing(QQ)
    Aops.<Dx> = OreAlgebra(A)
    if constTerm == 0:
        return Aops(Dx * expression/Dx)

    # 2. When the expression has nonzero constant term, then 'Dx*exp*constTerm - (constTerm)'*exp' gives an expression for the derivative.
    # i.e., Let 'p_r * f^(r) + ... + p_1 * f' + p_0 * f = 0'. Then,
    # Dx * expression = p_r * f^(r+1) + (p'_r + p_{r-1}) * f^(r) + ... + (p'_1 + p_0) * f' + p'_0 * f.
    # Therefore, in order to eliminate f-term, we should consider the following.
    else:
        expressionForDeri = constTerm * (Dx*expression) - diff(constTerm) * expression
        leadingTerm = A(expressionForDeri.coefficients()[-1])
        if leadingTerm.degree() == 0:
            return Aops(A((1/leadingTerm)) * expressionForDeri/Dx)
        else:
            return Aops(expressionForDeri/Dx)






#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#



def position(testlist, entry):
    r"""
        Returns the position(entry number) of given element in the list.

    INPUT:

        -  ``testlist`` - A list
        -  ``entry`` - An element in 'testlist'

    OUTPUT: Integer

    EXAMPLES::

        sage: position([1,2,3],2)
        1

        sage: position([1,2,3,4,5],4)
        3
        
        sage: position([1,2,3,4,5],7)
        no such entry

    AUTHORS:
    
        - Kisun Lee klee669@gatech.edu (2018-08-21)

    SEEALSO::

        :func:evaluateSystemPoint
        :func:evaluateSystemAtMid
        :func:evaluateJacobianPoint
        :func:evaluateJacobianInterval
        :func:invmat
        :func:gammaValue
    

    """



    # 1. For loops check i-th entry of the list and count the number of iterations until i-th entry of the list is equal to the element that we want to find a position.
    # When the loop finds the same entry, then it breaks.
    length = len(testlist)
    k = 0
    for i in range(0, length):
        if testlist[i] != entry:
            k = k+1
        else:
            break

    # 2. If the for loop above finds the given element, then it returns the number of iteration(i.e., the position of the element in the list).
    # Otherwise, it returns error message.
    if k == length:
        print "no such entry"
    else:
        return k



########## functions for Krawczyk method ##########


###### Krawcyzk method ######
# We define the Krawczyk operator K(X) = y - Y * #F(y) + (I_n - Y * #F'(X)) * (X - y)
# X : an interval vector
# y : a point in X
# Y : an invertible matrix
# #F(y) : interval box containing F(y)
# I_n : n by n identity matrix
# #F'(X) : interval matrix containing F'(X)
# X - y : an interval containing the origin obtained by translating X by y.
# If K(X) is contained in X, then we conclude that X has a unique root of 'F'.



###### RBF function in SAGE ######
# In order to apply the Krawczyk method, we should use intervals.
# SAGE function 'RBF' can be used in order to make an interval. It converts a list (or sequence) of length 2 into a sequence.
# For example, the command "RBF([1,2])" constructs a closed interval [1,2].
# In order to define interval vector or matrix use a function matrix on the list of intervals.
# For example, matrix([[RBF([1,2])],[RBF([2,4])]]) (2 by 1 interval vector)




#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#



    
def makemidvec(X):
    r"""
        Produces a vector consists of midpoints of an input interval vector.

    INPUT:

        -  ``X`` - An interval vector

    OUTPUT: A vector

    EXAMPLES::

        sage: makemidvec([RBF([1,2]),RBF([3,4])])
        [3/2]
        [7/2]

    AUTHORS:
    
        - Kisun Lee klee669@gatech.edu (2018-08-21)

    SEEALSO::

        :func:zeroball
        :func:krawczykoperator
    

    """
    
    # 1. Construct a list of entries for a vector.
    # Convert a midpoint into a rational number.
    listofy = []
    for i in range(0,len(X)):
        listofy.append([Rational(X[i].mid())])

    # 2. Return a matrix of midpoints.
    return matrix(listofy)





#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#



def zeroball(X):
    r"""
        Compute a traslated interval vector which is centered at the origin

    INPUT:

        -  ``X`` - An interval vector

    OUTPUT: A vector

    EXAMPLES::

        sage: zeroball([RBF([1,2]),RBF([3,4])])
        [[+/- 0.501]]
        [[+/- 0.501]]

    AUTHORS:
    
        - Kisun Lee klee669@gatech.edu (2018-08-21)

    SEEALSO::

        :func:krawczykoperator
    

    """

    # 1. make a vector consists of midpoints of each interval in 'X'.
    matent = []
    midvec = makemidvec(X)


    # 2. translate 'X' using a vector of midpoints.
    for i in range(0, len(X)):
        matent.append(X[i] - RBF(midvec[i][0]))
    return transpose(matrix(matent))




#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#




def evaluateSystemPoint(sys, varlist, explist, inilist, point, iniAt):
    r"""
        Evaluates ``sys`` at ``point``.

    INPUT:

        -  ``sys`` : A list of elements of multivariate polynomial ring
        -  ``varlist`` : A list of lists (variables used in ``sys``)
        -  ``explist`` : A list of elements in univariate Ore algebra in 'Dx' over univariate polynomial ring in 'x' over rational field
        -  ``inilist`` : A list of lists (initial values)
        -  ``point`` : A vector (or list) of numbers

    OUTPUT: A list of numbers

    EXAMPLES::

        sage: Pols.<x> = PolynomialRing(QQ) ; Dops.<Dx> =  OreAlgebra(Pols)
        sage: (x1, x2, x3, x4, y1, y2) = var('x1 x2 x3 x4 y1 y2')
        sage: A.<x1, x2, x3, x4, y1, y2> = PolynomialRing(QQ)
        sage: exp = Dx^2 + 2*x*Dx
        sage: explist = [exp, exp]
        sage: sys = [x3 * x4 - 1/2, x1^2 + x2^2 - 4, x3 - y1, x4 - y2]
        sage: ini = [0, 2/sqrt(pi)]
        sage: inilist = [ini, ini]
        sage: varlist = [x1, x2, [x3, x1], [x4, x2]]
        sage: point = (4/10, 2/10, 4/10, 10/10)
        sage: evaluateSystemPoint(sys, varlist, explist, inilist, point)
        [ -0.10000000000   -3.8000000000 -0.028392355047   0.77729741079]

    AUTHORS:
    
        - Kisun Lee klee669@gatech.edu (2018-08-21)

    SEEALSO::

        :func:beta
    

    """
    

    # 1. If 'point' is not a list (for example, a sequence), convert that into a list type.
    # Also, 'valueList' will take a value of D-finite functions evaluated at 'point'.
    valueList = list(point)


    # 2. Separate 'D-finite function part' from 'varlist'
    variablesForDfinitePart = varlist[-len(explist):]


    # 3. Construct a list of variables for system. (Note that 'varlist' is a "partially" nested list.)
    # For example, this part makes '[x1, x2, [x3, x1], [x4, x2]]' into [x1, x2, x3, x4].
    variablesForSystem = varlist[0:(len(varlist) - len(explist))] 
    for i in range(0, len(variablesForDfinitePart)):
        variablesForSystem.append(variablesForDfinitePart[i][0])


    # 4. Evaluate D-finite functions at points in 'point' (= dFinitevaluesAtPoint), and put them into 'valueList'.
    # use 'numerical_solution' function and take a midpoint as 'dFinitevaluesAtPoint'.
    # Make 'valueList' as a tuple type element.
    for i in range(0, len(explist)):
        evaluateAt = RBF(point[position(variablesForSystem, variablesForDfinitePart[i][1])]).mid()
#        dFinitevaluesAtPoint = ((explist[i]).numerical_solution(inilist[i], [0, evaluateAt],1e-10)).mid()
        dFinitevaluesAtPoint = ((explist[i]).numerical_solution(inilist[i], [iniAt, evaluateAt],1e-10)).mid()
        valueList.append(dFinitevaluesAtPoint)
    valueList = tuple(valueList)


    # 5. Since in 'values' only has values of D-finite functions over 'point', we should put these variables ('y1' and 'y2' in the above example) in 'variablesForSystem'.
    # This part makes '[x1, x2, x3, x4]' into '[x1, x2, x3, x4, y1, y2]'
    listofvari = variablesForSystem
    for i in range(0, len(explist)):
        variablesForSystem.append(-(sys[(len(varlist) - len(explist)) + i] - listofvari[(len(varlist) - len(explist)) + i]))


    # 6. This part actually evaluates 'sys' over 'values'.
    # Convert 'sys' into a matrix.
    # for-loop substitutes 'values' into the converted matrix.
    mat = matrix(sys)
    for i in range(0, len(valueList)):
        mat = mat.subs({variablesForSystem[i] : valueList[i]})
    return mat



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#




def evaluateSystemAtMid(sys, varlist, explist, inilist, X, iniAt):
    r"""
        Evaluates ``sys`` at midpoint of ``X``.

    INPUT:

        -  ``sys`` : A list of elements of multivariate polynomial ring
        -  ``varlist`` : A list of lists (variables used in ``sys``)
        -  ``explist`` : A list of elements in univariate Ore algebra in 'Dx' over univariate polynomial ring in 'x' over rational field
        -  ``inilist`` : A list of lists (initial values)
        -  ``X`` : A vector (or list) of 'RBF's.

    OUTPUT: A list of numbers

    EXAMPLES::

        sage: Pols.<x> = PolynomialRing(QQ) ; Dops.<Dx> =  OreAlgebra(Pols)
        sage: (x1, x2, x3, x4, y1, y2) = var('x1 x2 x3 x4 y1 y2')
        sage: A.<x1, x2, x3, x4, y1, y2> = PolynomialRing(QQ)
        sage: exp = Dx^2 + 2*x*Dx
        sage: explist = [exp, exp]
        sage: sys = [x3 * x4 - 1/2, x1^2 + x2^2 - 4, x3 - y1, x4 - y2]
        sage: ini = [0, 2/sqrt(pi)]
        sage: inilist = [ini, ini]
        sage: varlist = [x1, x2, [x3, x1], [x4, x2]]
        sage: X = [RBF([3/10,5/10]), RBF([1/10,3/10]), RBF([3/10,5/10]), RBF([9/10,11/10])]
        sage: evaluateSystemAtMid(sys, varlist, explist, inilist, point)
        [ -0.10000000000   -3.8000000000 -0.028392355047   0.77729741079]

    AUTHORS:
    
        - Kisun Lee klee669@gatech.edu (2018-08-21)

    SEEALSO::

        :func:krawczykoperator
    

    """


    
    # 1. 'valueList' has midpoints of each interval in 'X'. 
    # Also, 'valueList' will take a value of D-finite functions evaluated at these midpoints.
    valueList = []
    for i in range(0,len(X)):
        valueList.append(X[i].mid())

    
    # 2. Separate 'D-finite function part' from 'varlist'
    variablesForDfinitePart = varlist[-len(explist):] 


    # 3. Construct a list of variables for system. (Note that 'varlist' is a "partially" nested list.)
    # For example, this part makes '[x1, x2, [x3, x1], [x4, x2]]' into [x1, x2, x3, x4].
    variablesForSystem = varlist[0:(len(varlist) - len(explist))] 
    for i in range(0, len(variablesForDfinitePart)):
        variablesForSystem.append(variablesForDfinitePart[i][0])


    # 4. Evaluate D-finite functions at points in 'point' (= dFinitevaluesAtPoint), and put them into 'valueList'.
    # use 'numerical_solution' function and take a midpoint as 'dFinitevaluesAtPoint'.
    # Make 'valueList' as a tuple type element.
    for i in range(0, len(explist)):
        evaluateAt = (X[position(variablesForSystem, variablesForDfinitePart[i][1])]).mid()
#        dFinitevaluesAtPoint = (explist[i]).numerical_solution(inilist[i],[0,evaluateAt],1e-10).mid()
        dFinitevaluesAtPoint = (explist[i]).numerical_solution(inilist[i],[iniAt,evaluateAt],1e-10).mid()
        valueList.append(dFinitevaluesAtPoint)
    valueList = tuple(valueList)


    # 5. Since in 'values' has only values of D-finite functions over midpoints of 'X', we should put these variables ('y1' and 'y2' in the above example) in 'variablesForSystem'.
    # This part makes '[x1, x2, x3, x4]' into '[x1, x2, x3, x4, y1, y2]'
    listofvari = variablesForSystem
    for i in range(0, len(explist)):
        variablesForSystem.append(-(sys[(len(varlist) - len(explist)) + i] - listofvari[(len(varlist) - len(explist)) + i]))

        
    # 6. This part actually evaluates 'sys' over 'values'.
    # Convert 'sys' into a matrix.
    # for-loop substitutes 'values' into the converted matrix.
    mat = matrix(sys)
    for i in range(0, len(valueList)):
        mat = mat.subs({variablesForSystem[i] : valueList[i]})


    # 7. Lastly, convert list of intervals into an interval vector (so, this should take 'transpose' at last).
    # Use 'CBF' to convert entries into Complex balls.
    matent = []
    for i in range(0, len(mat[0])):
        matent.append(CBF(mat[0][i]))
    return transpose(matrix(matent))





#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#




def evaluateJacobianPoint(sys, varlist, explist, inilist, point, iniAt):
    r"""
        Evaluates the derivative of ``sys`` at ``point``.

    INPUT:

        -  ``sys`` : A list of elements of multivariate polynomial ring
        -  ``varlist`` : A list of lists (variables used in ``sys``)
        -  ``explist`` : A list of elements in univariate Ore algebra in 'Dx' over univariate polynomial ring in 'x' over rational field
        -  ``inilist`` : A list of lists (initial values)
        -  ``point`` : A vector (or list) of numbers

    OUTPUT: A list of numbers

    EXAMPLES::

        sage: Pols.<x> = PolynomialRing(QQ) ; Dops.<Dx> =  OreAlgebra(Pols)
        sage: (x1, x2, x3, x4, y1, y2) = var('x1 x2 x3 x4 y1 y2')
        sage: A.<x1, x2, x3, x4, y1, y2> = PolynomialRing(QQ)
        sage: exp = Dx^2 + 2*x*Dx
        sage: explist = [exp, exp]
        sage: sys = [x3 * x4 - 1/2, x1^2 + x2^2 - 4, x3 - y1, x4 - y2]
        sage: ini = [0, 2/sqrt(pi)]
        sage: inilist = [ini, ini]
        sage: varlist = [x1, x2, [x3, x1], [x4, x2]]
        sage: point = (4/10, 2/10, 4/10, 10/10)
        sage: evaluateJacobianPoint(sys, varlist, explist, inilist, point)
        [                    0                     0                     1                   2/5]
        [                  4/5                   2/5                     0                     0]
        [-0.428392355046668455                     0                     1                     0]
        [                    0  -0.22270258921047845                     0                     1]

    AUTHORS:
    
        - Kisun Lee klee669@gatech.edu (2018-08-21)

    SEEALSO::

        :func:beta
        :func:mu
    

    """

    
    # 1. All equations in 'sys' should be converted into a SAGE Symbolic Expression
    # If 'point' is not a list (for example, a sequence), convert that into a list type.
    # Also, 'valueList' will take a value of D-finite functions evaluated at 'point'.
    m = matrix(SR,len(sys))
    valueList = list(point)

    
    # 2. Separate 'D-finite function part' from 'varlist'
    variablesForDfinite = varlist[-len(explist):] # extract the variables for D-finite functions.


    # 3. Construct a list of variables for system. (Note that 'varlist' is a "partially" nested list.)
    # For example, this part makes '[x1, x2, [x3, x1], [x4, x2]]' into [x1, x2, x3, x4].
    variablesForSystem = varlist[0:(len(varlist) - len(explist))] # 
    for i in range(0, len(variablesForDfinite)):
        variablesForSystem.append(variablesForDfinite[i][0])
    jacOfVariables = jacobian(sys,variablesForSystem)

        
    # 4. Since in 'values' only has values of D-finite functions over 'point', we should put these variables ('y1' and 'y2' in the above example) in 'variablesForSystem'.
    listofvari = variablesForSystem
    for i in range(0, len(explist)):
        variablesForSystem.append(-(sys[(len(varlist) - len(explist)) + i] - listofvari[(len(varlist) - len(explist)) + i]))

        
    # 5. Since 'sys' has x-variables 'x1 x2 x3 ..' and y-variables 'y1 y2 y3... ', we should compute Jacobian for variables ('jacOfVariables') and the matrix obtained by differentiating 'D-finite function part ('m') respectively.
    # Eventually, we add two matrices.
    for i in range(0, len(explist)):
        m[(len(varlist) - len(explist)) + i,position(variablesForSystem, variablesForDfinite[i][1])] = - variablesForSystem[len(varlist) + i]
    jac = jacOfVariables + m


    # 6. Evaluate D-finite functions at points in 'point', and put them into 'valueList'.
    for i in range(0, len(explist)):
        evaluateAt = point[position(variablesForSystem, variablesForDfinite[i][1])]
#        dFinitevaluesAtPoint = (expForDeri(explist[i])).numerical_solution(iniForDeri(explist[i],inilist[i]),[0,evaluateAt])
        dFinitevaluesAtPoint = (expForDeri(explist[i])).numerical_solution(iniForDeri(explist[i],inilist[i],iniAt),[iniAt,evaluateAt])
        valueList.append(CBF(dFinitevaluesAtPoint).real().mid()) # should convert values into 'CBF'. (Operations don't support between CBF and numbers)

    
    # 7. This part actually evaluates 'sys' over 'valueList'.
    # for-loop substitutes 'valueList' into the converted matrix.
    for i in range(0, len(valueList)):
        jac = jac.subs({variablesForSystem[i] : valueList[i]})
    return jac




#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#






def evaluateJacobianInterval(sys, varlist, explist, inilist, X, iniAt, targetAccuracy):
    r"""
        Evaluates the derivative of ``sys`` at ``X``.

    INPUT:

        -  ``sys`` : A list of elements of multivariate polynomial ring
        -  ``varlist`` : A list of lists (variables used in ``sys``)
        -  ``explist`` : A list of elements in univariate Ore algebra in 'Dx' over univariate polynomial ring in 'x' over rational field
        -  ``inilist`` : A list of lists (initial values)
        -  ``X`` : A vector (or list) of 'RBF's
        -  ``targetAccuracy`` : A number

    OUTPUT: A matrix of intervals

    EXAMPLES::

        sage: Pols.<x> = PolynomialRing(QQ) ; Dops.<Dx> =  OreAlgebra(Pols)
        sage: (x1, x2, x3, x4, y1, y2) = var('x1 x2 x3 x4 y1 y2')
        sage: A.<x1, x2, x3, x4, y1, y2> = PolynomialRing(QQ)
        sage: exp = Dx^2 + 2*x*Dx
        sage: explist = [exp, exp]
        sage: sys = [x3 * x4 - 1/2, x1^2 + x2^2 - 4, x3 - y1, x4 - y2]
        sage: ini = [0, 2/sqrt(pi)]
        sage: inilist = [ini, ini]
        sage: varlist = [x1, x2, [x3, x1], [x4, x2]]
        sage: X = [RBF([3/10,5/10]), RBF([1/10,3/10]), RBF([3/10,5/10]), RBF([9/10,11/10])]
        sage: targetAccuracy = 1e-40
        sage: evaluateJacobianInterval(sys, varlist, explist, inilist, X, targetAccuracy)
        [                0                 0  [1e+0 +/- 0.101]       [+/- 0.501]]
        [       [+/- 1.01]       [+/- 0.601]                 0                 0]
        [       [+/- 1.06]                 0                 1                 0]
        [                0 [-1.1 +/- 0.0893]                 0                 1]

    AUTHORS:
    
        - Kisun Lee klee669@gatech.edu (2018-08-21)

    SEEALSO::

        :func:krawczykoperator
    

    """


    
    # 1. All equations in 'sys' should be converted into a SAGE Symbolic Expression
    # Also, 'valueList' will take a value of D-finite functions evaluated at 'point'.    
    m = matrix(SR,len(sys))
    valueList = list(X)

    
    # 2. Separate 'D-finite function part' from 'varlist'
    variablesForDfinite = varlist[-len(explist):]


    # 3. Construct a list of variables for system. (Note that 'varlist' is a "partially" nested list.)
    # For example, this part makes '[x1, x2, [x3, x1], [x4, x2]]' into [x1, x2, x3, x4].
    variablesForSystem = varlist[0:(len(varlist) - len(explist))]  
    for i in range(0, len(variablesForDfinite)):
        variablesForSystem.append(variablesForDfinite[i][0])


    # 4. Compute the Jacobian of 'sys' with respect to variables 'x1, x2, ...'
    jac = jacobian(sys,variablesForSystem)


    # 5. Since in 'values' only has values of D-finite functions over 'point', we should put these variables ('y1' and 'y2' in the above example) in 'variablesForSystem'.
    listofvari = variablesForSystem
    for i in range(0, len(explist)):
        variablesForSystem.append(-(sys[(len(varlist) - len(explist)) + i] - listofvari[(len(varlist) - len(explist)) + i]))


    # 6. Since 'sys' has x-variables 'x1 x2 x3 ..' and y-variables 'y1 y2 y3... ', we should compute Jacobian for variables ('jacOfVariables') and the matrix obtained by differentiating 'D-finite function part ('m') respectively.
    # Eventually, we add two matrices.
    for i in range(0, len(explist)):
        m[(len(varlist) - len(explist)) + i,position(variablesForSystem, variablesForDfinite[i][1])] = - variablesForSystem[len(varlist) + i]
    jac = jac + m

    

    # 7. Evaluate the derivative of D-finite functions over 'X', and put them into 'valueList'.
    # 'numerical_solution' evaluates the derivative of D-finite functions with accuracy 'targetAccuracy'.
    # for-loop makes computes the bound for D-finite functions in 'sys'.
    for i in range(0, len(explist)):
        derivativeExpression = expForDeri(explist[i])
        derivativeInitialValues = iniForDeri(explist[i], inilist[i],iniAt)
        evaluateAt = X[position(variablesForSystem, variablesForDfinite[i][1])]
#        valueList.append((derivativeExpression).numerical_solution(derivativeInitialValues, [0, evaluateAt], targetAccuracy))
        valueList.append((derivativeExpression).numerical_solution(derivativeInitialValues, [iniAt, evaluateAt], targetAccuracy))


    # 8. This part evaluates 'sys' over 'valueList'.
    # for-loop substitutes 'valueList' into 'jac'.
    for i in range(0, len(valueList)):
        jac = jac.subs({variablesForSystem[i] : valueList[i]})
    return jac




#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#




def invmat(sys, varlist, explist, inilist, X, iniAt):
    r"""
        Approximates the inverse of ``sys`` at the midpoint of ``X``.

    INPUT:

        -  ``sys`` : A list of elements of multivariate polynomial ring
        -  ``varlist`` : A list of lists (variables used in ``sys``)
        -  ``explist`` : A list of elements in univariate Ore algebra in 'Dx' over univariate polynomial ring in 'x' over rational field
        -  ``inilist`` : A list of lists (initial values)
        -  ``X`` : A vector (or list) of 'RBF's

    OUTPUT: A matrix of intervals

    EXAMPLES::

        sage: Pols.<x> = PolynomialRing(QQ) ; Dops.<Dx> =  OreAlgebra(Pols)
        sage: (x1, x2, x3, x4, y1, y2) = var('x1 x2 x3 x4 y1 y2')
        sage: A.<x1, x2, x3, x4, y1, y2> = PolynomialRing(QQ)
        sage: exp = Dx^2 + 2*x*Dx
        sage: explist = [exp, exp]
        sage: sys = [x3 * x4 - 1/2, x1^2 + x2^2 - 4, x3 - y1, x4 - y2]
        sage: ini = [0, 2/sqrt(pi)]
        sage: inilist = [ini, ini]
        sage: varlist = [x1, x2, [x3, x1], [x4, x2]]
        sage: X = [RBF([3/10,5/10]), RBF([1/10,3/10]), RBF([3/10,5/10]), RBF([9/10,11/10])]
        sage: invmat(sys, varlist, explist, inilist, X)
        [ [10.61194084185438 +/- 5.80e-16] [-11.50477422535320 +/- 3.18e-15] [-10.61194084185438 +/- 5.80e-16] [-4.244776336741752 +/- 3.02e-16]]
        [[-21.22388168370876 +/- 1.16e-15]  [25.50954845070639 +/- 3.66e-15]  [21.22388168370876 +/- 1.16e-15]  [8.489552673483503 +/- 3.98e-16]]
        [ [10.20381938028256 +/- 1.48e-15] [-11.06231555149910 +/- 2.29e-15] [-9.203819380282559 +/- 4.74e-16] [-4.081527752113023 +/- 1.23e-16]]
        [[-23.00954845070640 +/- 2.80e-15]  [27.65578887874776 +/- 1.52e-16]  [23.00954845070640 +/- 2.80e-15]  [10.20381938028256 +/- 1.48e-15]]

    AUTHORS:
    
        - Kisun Lee klee669@gatech.edu (2018-08-21)

    SEEALSO::

        :func:krawczykoperator
    

    """

    
    # 1. All equations in 'sys' should be converted into a SAGE Symbolic Expression
    # Put midpoints of all intervals in 'X' into 'valueList'
    # Also, 'valueList' will take a value of D-finite functions evaluated at 'point'.
    m = matrix(SR,len(sys))
    valueList = []
    for i in range(0,len(X)):
        valueList.append(X[i].mid())

    
    # 2. Separate 'D-finite function part' from 'varlist'
    variablesForDfinite = varlist[-len(explist):] # extract the variables for D-finite functions.


    # 3. Construct a list of variables for system. (Note that 'varlist' is a "partially" nested list.)
    # For example, this part makes '[x1, x2, [x3, x1], [x4, x2]]' into [x1, x2, x3, x4].
    variablesForSystem = varlist[0:(len(varlist) - len(explist))] # 
    for i in range(0, len(variablesForDfinite)):
        variablesForSystem.append(variablesForDfinite[i][0])


    # 4. Compute the Jacobian of 'sys' with respect to variables 'x1, x2, ...'
    jac = jacobian(sys,variablesForSystem)


    # 5. Since in 'values' only has values of D-finite functions over 'point', we should put these variables ('y1' and 'y2' in the above example) in 'variablesForSystem'.
    listofvari = variablesForSystem
    for i in range(0, len(explist)):
        variablesForSystem.append(-(sys[(len(varlist) - len(explist)) + i] - listofvari[(len(varlist) - len(explist)) + i]))


    # 6. Since 'sys' has x-variables 'x1 x2 x3 ..' and y-variables 'y1 y2 y3... ', we should compute Jacobian for variables ('jacOfVariables') and the matrix obtained by differentiating 'D-finite function part ('m') respectively.
    # Eventually, we add two matrices.
    for i in range(0, len(explist)):
        m[(len(varlist) - len(explist)) + i,position(variablesForSystem, variablesForDfinite[i][1])] = - variablesForSystem[len(varlist) + i]
    jac = jac + m


    # 7. Evaluate the derivative of D-finite functions at 'm(X)', and put them into 'valueList'.
    for i in range(0, len(explist)):
        evaluateAt = (X[position(variablesForSystem, variablesForDfinite[i][1])]).mid()
        derivativeExpression = expForDeri(explist[i])
        derivativeInitialValues = iniForDeri(explist[i],inilist[i], iniAt)
#        valueList.append(CBF(derivativeExpression.numerical_solution(derivativeInitialValues,[0,evaluateAt])).real().mid())
        valueList.append(CBF(derivativeExpression.numerical_solution(derivativeInitialValues,[iniAt,evaluateAt])).real().mid())
    for i in range(0, len(valueList)):
        jac = jac.subs({variablesForSystem[i] : valueList[i]})
        

    # 8. compute the inverse of 'jac'.
    jacinv = jac.inverse()


    # 9. Convert entries in 'jac' into 'RBF'.
    # This process is needed because operations between "numbers" and 'RBF' are not supported.
    mats = []
    for i in range(0, jacinv.nrows()):
        row = []
        for j in range(0, jacinv.ncols()):
            row.append(CBF(jacinv[i][j]))
        mats.append(row)
    return matrix(mats)




#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#


def intervalMatrixNorm(mat):
    rowNormList = []
    for i in range(0, mat.nrows()):
        rowNorm = 0
        for j in range(0, mat.ncols()):
            rowNorm = rowNorm + max(abs(RBF(mat[i][j]).lower()),abs(RBF(mat[i][j]).upper()))
        rowNormList.append(rowNorm)
    return max(rowNormList)




def krawczykoperator(sys,varlist,explist,inilist,X, iniAt, targetAccuracy):
    r"""
        Computes the Krawczyk operator at ``X``.

    INPUT:

        -  ``sys`` : A list of elements of multivariate polynomial ring
        -  ``varlist`` : A list of lists (variables used in ``sys``)
        -  ``explist`` : A list of elements in univariate Ore algebra in 'Dx' over univariate polynomial ring in 'x' over rational field
        -  ``inilist`` : A list of lists (initial values)
        -  ``X`` : A vector (or list) of 'RBF's
        -  ``targetAccuracy`` : A number

    OUTPUT: A vector of intervals

    EXAMPLES::

        sage: Pols.<x> = PolynomialRing(QQ) ; Dops.<Dx> =  OreAlgebra(Pols)
        sage: (x1, x2, x3, x4, y1, y2) = var('x1 x2 x3 x4 y1 y2')
        sage: A.<x1, x2, x3, x4, y1, y2> = PolynomialRing(QQ)
        sage: exp = Dx^2 + 2*x*Dx
        sage: explist = [exp, exp]
        sage: sys = [x3 * x4 - 1/2, x1^2 + x2^2 - 4, x3 - y1, x4 - y2]
        sage: ini = [0, 2/sqrt(pi)]
        sage: inilist = [ini, ini]
        sage: varlist = [x1, x2, [x3, x1], [x4, x2]]
        sage: X = [RBF([3/10,5/10]), RBF([1/10,3/10]), RBF([3/10,5/10]), RBF([9/10,11/10])]
        sage: targetAccuracy = 1e-40
        sage: krawczykoperator(sys, varlist, explist, inilist, X, targetAccuracy)
        [[-4e+1 +/- 1.54]]
        [ [9e+1 +/- 2.68]]
        [[-4e+1 +/- 3.06]]
        [ [1e+2 +/- 5.34]]

    AUTHORS:
    
        - Kisun Lee klee669@gatech.edu (2018-08-21)

    

    """
    

    # 1. make a vector of midpoints of 'X'.
    y = makemidvec(X)

    
    # 2. approximate an invertible matrix.
    Y = invmat(sys,varlist,explist,inilist,X,iniAt)

    
    # 3. compute the Krawczyk operator.
    roughvalue =  y - Y * evaluateSystemAtMid(sys,varlist,explist,inilist,X, iniAt) + (matrix.identity(len(sys)) - Y * evaluateJacobianInterval(sys,varlist,explist,inilist,X,iniAt,targetAccuracy)) * zeroball(X)
#    print ("|I_n-Y#F'(I)| = " + str(intervalMatrixNorm((matrix.identity(len(sys)) - Y * evaluateJacobianInterval(sys,varlist,explist,inilist,X,targetAccuracy)))))
#    print ("Uniqueness = " + str(intervalMatrixNorm((matrix.identity(len(sys)) - Y * evaluateJacobianInterval(sys,varlist,explist,inilist,X,targetAccuracy))) < 1))
    

    # 4. convert entries in 'roughvalue' to intervals. (In order to make an interval matrix 'K(X)')
    oper = []
    for x in range(0, len(sys)):
        oper.append([RBF(roughvalue[x][0])])
    return matrix(oper)




########## functions for alpha theory ##########


###### alpha theory ######
# We define alpha(F,x) = beta(F,x) * gamma(F,x)
# beta(F,x) = \|(F'(X))^-1 * F(x)\|
# gamma(F,x) = sup_{k>1} \|(F'(x))^-1 * F^{(k)}(x) / k!\|^{1/k-1}
# If alpha(F,x) < (13 - 3* sqrt(17)) / 4, then we conclude that x in an apprximate solution to F.
# In this case, we may have a bound for gamma(F,x) <= mu(F,x) * (d^(3/2)/(2*\|(1,x)\|) + C) where C is a sum of bounds for sup_{k>1} \|(g'(x))^-1 * g^{(k)}(x) / k!\|^{1/k-1}.





#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#


def beta(sys, varlist, explist, inilist, point, iniAt):
    r"""
        Computes the value of beta for ``sys`` at ``point``.

    INPUT:

        -  ``sys`` : A list of elements of multivariate polynomial ring
        -  ``varlist`` : A list of lists (variables used in ``sys``)
        -  ``explist`` : A list of elements in univariate Ore algebra in 'Dx' over univariate polynomial ring in 'x' over rational field
        -  ``inilist`` : A list of lists (initial values)
        -  ``point`` : A vector (or list) of numbers

    OUTPUT: A number

    EXAMPLES::

        sage: Pols.<x> = PolynomialRing(QQ) ; Dops.<Dx> =  OreAlgebra(Pols)
        sage: (x1, x2, x3, x4, y1, y2) = var('x1 x2 x3 x4 y1 y2')
        sage: A.<x1, x2, x3, x4, y1, y2> = PolynomialRing(QQ)
        sage: exp = Dx^2 + 2*x*Dx
        sage: explist = [exp, exp]
        sage: sys = [x3 * x4 - 1/2, x1^2 + x2^2 - 4, x3 - y1, x4 - y2]
        sage: ini = [0, 2/sqrt(pi)]
        sage: inilist = [ini, ini]
        sage: varlist = [x1, x2, [x3, x1], [x4, x2]]
        sage: point = (4/10, 2/10, 4/10, 10/10)
        sage: beta(sys, varlist, explist, inilist, point)
        13.532009384573664

    AUTHORS:
    
        - Kisun Lee klee669@gatech.edu (2018-08-21)



    """

    return (matrix(QQ,(evaluateJacobianPoint(sys, varlist, explist, inilist, point,iniAt))).inverse() * transpose(evaluateSystemPoint(sys, varlist, explist, inilist, point, iniAt))).norm()






#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#





def pointnorm(point):
    r"""
        Computes the value of 'projectivized' norm of ``point``.

    INPUT:

        -  ``point`` : A vector (or list) of numbers

    OUTPUT: A number

    EXAMPLES::
 
        sage: pointnorm([1,2,4]) 
        sqrt(22)
        sage: pointnorm([1,2+I,4])
        sqrt(23)

    AUTHORS:
    
        - Kisun Lee klee669@gatech.edu (2018-08-21)


    SEEALSO::
        
        :func:deltaMatrix
        :func:gammaValue


    """

    
    entries = []
    for i in range(0, len(point)):
        entries.append((abs(point[i]))^2)
    norm = 1 + sum(entries)
    return sqrt(norm)






#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#



def polynorm(poly):
    r"""
        Computes the value of (Bombieri-Weyl) norm of ``poly``.

    INPUT:

        -  ``poly`` : An element of multivariate polynomial ring

    OUTPUT: A number

    EXAMPLES::
 
        sage: BB.<x1,x2> = PolynomialRing(QQ)
        sage: polynorm(2*x1^2 + 3*x1*x2)
        sqrt(17/2)

    AUTHORS:
    
        - Kisun Lee klee669@gatech.edu (2018-08-21)


    SEEALSO::
        
        :func:sysnorm


    """

    
    # 1. Note that the Bombieri-Weyl norm is defined from the factorial of degree of a polynomial and its monomials. Thus, we should collect information of degree and monomials of 'poly' first.
    d = poly.degree()
    listOfTerms = poly.monomials()


    # 2. 'main for-loop' This computes the value of norm for each term in 'poly', and finally sum all of them.
    summands = []
    for i in range(0, len(listOfTerms)):


        # 2.1. This part collects all variables used in each term of 'poly'.
        vari  = listOfTerms[i].variables()

        # 2.2. Computes degree of each term with respect to its variables used.
        deg = [0]
        for j in range(0, len(vari)):
            deg.append(listOfTerms[i].degree(vari[j]))

        # 2.3. 'reversedDegreesum' computes '(degree) - (sum of degrees for each term)'. 
        reversedDegreesum = d - sum(deg)

        # 2.4. This part computes factorials of degrees computed above.
        k = 0
        for j in range(0, len(deg)):
            deg[j] = factorial(deg[j])
            k = prod(deg)

        # 2.5. Put all computed value into the list 'summands'.
        summands.append(abs(poly.coefficients()[i])^2 * (factorial(reversedDegreesum) * k)/factorial(d))


    return sqrt(sum(summands))






#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#






def sysnorm(sys):
    r"""
        Computes the norm for ``sys``

        -  For system 'F', we define the norm of 'F' by '||F|| = Sum_i ||p_i||' where 'F = [ p_1, p_2, ..., p_n ]'.

    INPUT:

        -  ``sys`` : A list of elements of multivariate polynomial ring

    OUTPUT: A number

    EXAMPLES::
 
        sage: A.<x1, x2, x3, x4, y1, y2> = PolynomialRing(QQ)
        sage: sys = [x3*x4 - 1/2, x1^2 + x2^2 - 4, x3 - y1, x4 - y2] 
        sage: sysnorm(sys)
        1/2*sqrt(91)

    AUTHORS:
    
        - Kisun Lee klee669@gatech.edu (2018-08-21)


    SEEALSO::
        
        :func:mu


    """


    
    k = []
    for i in range(0, len(sys)):
        k.append((polynorm(sys[i]))^2)
    return sqrt(sum(k))






#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#






def deltaMatrix(sys,point):
    r"""
        Computes the Delta matrix of the polynomial system ``sys`` at ``point``.

        -  For the polynomial system 'P = [ p_1, p_2, ..., p_n ]', Delta matrix is a diagonal matrix with its diagonal entry defined by 'sqrt(deg(p_i)) * ||(1,x)||^{deg(p_i) - 1}.

    INPUT:

        -  ``sys`` : A list of elements of multivariate polynomial ring
        -  ``point`` : A vector (or list) of numbers


    OUTPUT: A number

    EXAMPLES::
 
        sage: A.<x1, x2> = PolynomialRing(QQ)
        sage: deltaMatrix([x1*x2, x1^2 + x2^3 - 3], [3,4])
        [sqrt(26)*sqrt(2)                0]
        [               0       26*sqrt(3)]

    AUTHORS:
    
        - Kisun Lee klee669@gatech.edu (2018-08-21)


    SEEALSO::
        
        :func:mu


    """

    
    k = []
    for i in range(0, len(sys)):
        k.append(sqrt(sys[i].degree()) * pointnorm(point)**(sys[i].degree() - 1))
    return diagonal_matrix(k)






#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#





def mu(sys, varlist, explist, inilist, point):
    r"""
        Computes the value of mu for ``sys`` at ``point``.

    INPUT:

        -  ``sys`` : A list of elements of multivariate polynomial ring
        -  ``varlist`` : A list of lists (variables used in ``sys``)
        -  ``explist`` : A list of elements in univariate Ore algebra in 'Dx' over univariate polynomial ring in 'x' over rational field
        -  ``inilist`` : A list of lists (initial values)
        -  ``point`` : A vector (or list) of numbers

    OUTPUT: A number

    EXAMPLES::

        sage: Pols.<x> = PolynomialRing(QQ) ; Dops.<Dx> =  OreAlgebra(Pols)
        sage: (x1, x2, x3, x4, y1, y2) = var('x1 x2 x3 x4 y1 y2')
        sage: A.<x1, x2, x3, x4, y1, y2> = PolynomialRing(QQ)
        sage: exp = Dx^2 + 2*x*Dx
        sage: explist = [exp, exp]
        sage: sys = [x3 * x4 - 1/2, x1^2 + x2^2 - 4, x3 - y1, x4 - y2]
        sage: ini = [0, 2/sqrt(pi)]
        sage: inilist = [ini, ini]
        sage: varlist = [x1, x2, [x3, x1], [x4, x2]]
        sage: point = (4/10, 2/10, 4/10, 10/10)
        sage: mu(sys, varlist, explist, inilist, point)
        96.87420782217045

    AUTHORS:
    
        - Kisun Lee klee669@gatech.edu (2018-08-21)

    SEEALSO::

        :func:gammaValue
    

    """
    
    polysys = sys[0:(len(sys) - len(explist))]
    polysysnorm = sysnorm(polysys)
    firstblock = deltaMatrix(polysys,point) * polysysnorm
    secondblock = identity_matrix(len(explist))
    delta = block_diagonal_matrix(firstblock,secondblock)
    jac = evaluateJacobianPoint(sys, varlist, explist, inilist, point, iniAt)
    return max(1, norm(matrix(QQ,jac).inverse() * delta))





#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#


def gammaValue(sys, varlist, explist, inilist, rad, point, iniAt, targetAccuracy):
    r"""
        Computes the value of gamma for ``sys`` at ``point``.

    INPUT:

        -  ``sys`` : A list of elements of multivariate polynomial ring
        -  ``varlist`` : A list of lists (variables used in ``sys``)
        -  ``explist`` : A list of elements in univariate Ore algebra in 'Dx' over univariate polynomial ring in 'x' over rational field
        -  ``inilist`` : A list of lists (initial values)
        -  ``rad`` : A number
        -  ``point`` : A vector (or list) of numbers
        -  ``targetAccuracy`` : A number

    OUTPUT: A number

    EXAMPLES::

        sage: Pols.<x> = PolynomialRing(QQ) ; Dops.<Dx> =  OreAlgebra(Pols)
        sage: (x1, x2, x3, x4, y1, y2) = var('x1 x2 x3 x4 y1 y2')
        sage: A.<x1, x2, x3, x4, y1, y2> = PolynomialRing(QQ)
        sage: exp = Dx^2 + 2*x*Dx
        sage: explist = [exp, exp]
        sage: sys = [x3 * x4 - 1/2, x1^2 + x2^2 - 4, x3 - y1, x4 - y2]
        sage: ini = [0, 2/sqrt(pi)]
        sage: inilist = [ini, ini]
        sage: varlist = [x1, x2, [x3, x1], [x4, x2]]
        sage: rad = 1
        sage: point = (4/10, 2/10, 4/10, 10/10)
        sage: targetAccuracy = 1e-40
        sage: g = gammaValue(sys, varlist, explist, inilist, rad, point, targetAccuracy)
        sage: g.n()
        17226.1733144475

    AUTHORS:
    
        - Kisun Lee klee669@gatech.edu (2018-08-21)


    """


    

    # 1. Split the system into polynomial part and D-finite function part.
    polysys = sys[0:(len(sys) - len(explist))]


    # 2. Compute the maximum degree of polynomials in 'sys'.
    d = []
    for i in range(0, len(polysys)):
        d.append(polysys[i].degree())
    d = max(d)


    # 3. Separate 'D-finite function part' from 'varlist'    
    variablesForDfinite = varlist[-len(explist):] # extract the variables for D-finite functions.


    # 4. Construct a list of variables for system. (Note that 'varlist' is a "partially" nested list.)
    # For example, this part makes '[x1, x2, [x3, x1], [x4, x2]]' into [x1, x2, x3, x4].
    variablesForSystem = varlist[0:(len(varlist) - len(explist))] # 
    for i in range(0, len(variablesForDfinite)):
        variablesForSystem.append(variablesForDfinite[i][0])

    
    # 5. Since in 'values' only has values of D-finite functions over 'point', we should put these variables ('y1' and 'y2' in the above example) in 'variablesForSystem'.
    listofvari = variablesForSystem
    for i in range(0, len(explist)):
        variablesForSystem.append(-(sys[(len(varlist) - len(explist)) + i] - listofvari[(len(varlist) - len(explist)) + i]))

        
    # 6. Compute the maximum values of the derivatives of each D-finite function, and put them into 'mlist'.
    mlist = []
    for i in range(0, len(explist)):

        
        # construct a cover for the disk 
        evaluateAt = point[position(variablesForSystem, variablesForDfinite[i][1])]
        circleCover = CBF([evaluateAt - (rad + rad*I), evaluateAt + (rad + rad*I)])


        # bound for D-finite function (no differentiation)
#        boundInterval0 = (explist[i].numerical_solution(inilist[i], [0, CBF(circleCover)],targetAccuracy))
        boundInterval0 = (explist[i].numerical_solution(inilist[i], [iniAt, CBF(circleCover)],targetAccuracy))
        upperbound0 = sqrt(max(abs(boundInterval0.real().upper()), abs(boundInterval0.real().lower()))^2 + max(abs(boundInterval0.imag().upper()), abs(boundInterval0.imag().lower()))^2)


        # bound for the derivative
        derivativeExpression = expForDeri(explist[i])
        derivativeInitialValues = iniForDeri(explist[i],inilist[i],iniAt)
#        boundInterval1 = ((derivativeExpression).numerical_solution(derivativeInitialValues, [0, CBF(circleCover)]))
        boundInterval1 = ((derivativeExpression).numerical_solution(derivativeInitialValues, [iniAt, CBF(circleCover)]))
        upperbound1 = sqrt(max(abs(boundInterval1.real().upper()), abs(boundInterval1.real().lower()))^2 + max(abs(boundInterval1.imag().upper()), abs(boundInterval1.imag().lower()))^2)


        # bound for the second derivative
        secondDerivativeExpression = expForDeri(derivativeExpression)
        secondDerivativeInitialValues = iniForDeri(derivativeExpression, derivativeInitialValues, iniAt)
#        boundInterval2 = ((secondDerivativeExpression).numerical_solution(secondDerivativeInitialValues, [0, CBF(circleCover)]))
        boundInterval2 = ((secondDerivativeExpression).numerical_solution(secondDerivativeInitialValues, [iniAt, CBF(circleCover)]))
        upperbound2 = sqrt(max(abs(boundInterval2.real().upper()), abs(boundInterval2.real().lower()))^2 + max(abs(boundInterval2.imag().upper()), abs(boundInterval2.imag().lower()))^2)

        
        mlist.append(min(max([1,upperbound0/rad]),max([1,upperbound1/2]),max([1,upperbound2 * rad/2]))/(rad))
        

    # 7. using 'mu' and information obtained from above, compute 'gamma' value and return it.
    return mu(sys,varlist,explist,inilist,point)*(sqrt(d^3)/(2*pointnorm(point))+sum(mlist))





