import galois
GF2 = galois.GF(2)
GF4 = galois.GF(4)
GF16 = galois.GF(16)
p=2
m= 2
bigfield = GF4
cero = bigfield(0)
one = bigfield(1) 
dos = bigfield(2)
poly_x =galois.Poly.Identity(bigfield)
maxdegree= 40
x_to=[]
x_to.append( galois.Poly([1], field=bigfield) )
for deg in range(maxdegree + 1):
    x_to.append( poly_x * x_to[deg] )
def sigma(a, power=1):
    if power == 0:
        return a
    else:
        return sigma(a, power-1)**p
def decompose_poly( polyA ):
    output = []
    for i in range(m):
        output.append(galois.Poly.Zero(bigfield))
    coeficientes = polyA.coeffs[::-1]
    for i in range(len(coeficientes)):
        modi = i % m
        output[modi] = output[modi] + (coeficientes[i] *\
         * x_to[ i - modi ] ) 
    return output
def createBi(poly):
    Bi = []
    aux1= poly
    for i in range(m):
        if i==0:
            Bi.append(poly)
        else:
            coeficientes = aux1.coeffs
            for j in range(len(coeficientes)):
                coeficientes[j]= sigma(coeficientes[j])
            aux2 =  galois.Poly(coeficientes)
            Bi.append(aux2)
            aux1=aux2
    return Bi
def products_aibi(Adecomp, Bi):
    output = []
    for index in range(m):
        output.append( Adecomp[index] * Bi[index] *\ 
        * x_to[index] )
    return output 
def add_aibi (vector_of_aibi):
    output = galois.Poly.Zero(bigfield)
    for index in range(m):
        output = output + vector_of_aibi[index]
    return output 
def skewmultiplication(poly1, poly2):
    aidecomp = decompose_poly(poly1)
    bi = createBi(poly2)
    aibi = products_aibi(aidecomp, bi)
    output = add_aibi(aibi)
    return output
def createBi_inverse(poly):
    Bi = []
    aux1= poly
    for i in range(m):
        if i==0:
            Bi.append(poly)
        else:
            coeficientes = aux1.coeffs
            for j in range(len(coeficientes)):
                coeficientes[j]= sigma(coeficientes[j],m-1)
            aux2 =  galois.Poly(coeficientes)
            Bi.append(aux2)
            aux1=aux2
    return Bi
def skewmultiplication_minusone(poly1, poly2):
    aidecomp = decompose_poly(poly1)
    bi = createBi_inverse(poly2)
    aibi = products_aibi(aidecomp, bi)
    output = add_aibi(aibi)
    return output
def B_to_n(poly, n):
    Bi = []
    aux1= poly
    for i in range(n+1):
        if i==0:
            Bi.append(poly)
        else:
            coeficientes = aux1.coeffs
            for j in range(len(coeficientes)):
                coeficientes[j]= sigma(coeficientes[j])
            aux2 =  galois.Poly(coeficientes)
            Bi.append(aux2)
            aux1=aux2
    return Bi[n]
def reciprocal(poly, delta): 
    coefficients = poly.coeffs
    coefficients = coefficients[::-1]
    output = galois.Poly(coefficients)
    if ( poly.degree == delta ):
        pass
    elif (poly.degree < delta):
        output = output * x_to[ delta - poly.degree ]
    else:
        print("Polynomial Degree is higher than expected")
        return None    
    return output
def polymod(poly, power): 
    res = poly % x_to[power + 1]
    return res
def skew_Reuclideandivision( polyA, polyB ):
    d1 = polyA.degree 
    d2 = polyB.degree 
    if (d1 < d2):
        return [galois.Poly.Zero(bigfield), polyA]
    Bn = B_to_n(polyB, d1-d2)
    Btilde = reciprocal(Bn , d2)
    Qtildecoef = one / (Btilde(cero))
    Qtilde = galois.Poly([Qtildecoef], field=bigfield)
    i=1
    while i < (d1 - d2 + 1):
        step1 = polymod( Btilde, 2*i -1 ) 
        step2 = Qtilde + Qtilde - ( *\
         skewmultiplication_minusone( *\
         skewmultiplication_minusone(Qtilde,step1),Qtilde))
        Qtilde = polymod( step2, 2*i )
        i = 2 * i
    Qtilde=polymod(skewmultiplication_minusone(polymod(*\
     reciprocal(polyA, d1) ,d1-d2) , Qtilde ), d1-d2 )
    Q = reciprocal( Qtilde, d1-d2 ) 
    R = polyA - skewmultiplication(Q, polyB)
    output = [Q,R]
    return output
def B_to_minusn(poly, n): 
    Bi = []
    aux1= poly
    for i in range(n+1):
        if i==0:
            Bi.append(poly)
        else:
            coeficientes = aux1.coeffs
            for j in range(len(coeficientes)):
                coeficientes[j]= sigma(coeficientes[j],m-1)
            aux2 =  galois.Poly(coeficientes)
            Bi.append(aux2)
            aux1=aux2
    return Bi[n]
def skew_Leuclideandivision( polyA, polyB ):
    d1 = polyA.degree 
    d2 = polyB.degree 
    if (d1 < d2):
        return [galois.Poly.Zero(bigfield), polyA]
    Btilde = reciprocal(polyB , d2)
    Qtildecoef = one / (Btilde(cero))
    Qtilde = galois.Poly([Qtildecoef], field=bigfield)
    i=1
    while i < (d1 - d2 + 1):
        step1 = polymod( Btilde, 2*i ) 
        step2 = Qtilde + Qtilde - ( *\
         skewmultiplication_minusone(*\
         skewmultiplication_minusone(*\
         Qtilde , step1) , Qtilde ) )
        Qtilde = polymod( step2, 2*i )
        i = 2 * i
    Qtilde = polymod(  skewmultiplication_minusone(Qtilde,*\
     polymod( reciprocal(polyA, d1) ,d1-d2) ), d1-d2 )
    Q = reciprocal( Qtilde, d1-d2 )
    Q = B_to_minusn(Q, d2)
    R = polyA - skewmultiplication( polyB, Q )
    output = [Q,R]
    return output
zeropoly= galois.Poly.Zero(bigfield)
onepoly= galois.Poly.One(bigfield)
def extended_grcd(a, b): 
        x0 = onepoly
        x1 = zeropoly
        y0 = zeropoly
        y1 = onepoly
        while b != zeropoly:
            tempa = a
            tempb = b
            q, b = skew_Reuclideandivision(tempa,tempb)
            a = tempb
            tempx0 = x0
            x0 = x1
            x1 = tempx0 - skewmultiplication(q,x0)
            tempy0 = y0
            y0 = y1
            y1 = tempy0 - skewmultiplication(q,y0)
        return [a, x0, y0]
def extended_glcd(a, b):
        x0 = onepoly
        x1 = zeropoly
        y0 = zeropoly
        y1 = onepoly
        while b != zeropoly:
            tempa = a
            tempb = b
            q, b = skew_Leuclideandivision(tempa,tempb)
            a = tempb
            tempx0 = x0
            x0 = x1
            x1 = tempx0 - skewmultiplication(x0,q)
            tempy0 = y0
            y0 = y1
            y1 = tempy0 - skewmultiplication(y0,q)
        return [a, x0, y0]
def isFtotaldivisorofG(f,g):
    n = f.degree
    xig = g
    for i in range(n):
        gxi = g * x_to[i]
        coeficientes = xig.coeffs
        for j in range(len(coeficientes)):
            coeficientes[j]= sigma(coeficientes[j])
        xig =  galois.Poly(coeficientes)
        q1, r1 = skew_Leuclideandivision( xig , f )
        q2, r2 = skew_Reuclideandivision( gxi , f )
        if r1 != zeropoly:
            return False
        if r2 != zeropoly:
            return False
    return True