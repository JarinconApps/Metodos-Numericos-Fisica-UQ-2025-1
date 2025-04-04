# Importar las librerias necesarias
import numpy as np
import matplotlib.pyplot as plt
import sympy

# Definimos la función: MultiplicacionSintetica
def MultiplicacionSintetica(raices, console=False):

    """
    ## ***Función***: MultiplicacionSintetica
    - **Descripción:** Dada una lista de puntos, las raíces de un polinomio, calcula los coeficientes del polinomio
    - **Parámetros:**
        - *raices:* Lista de raices del polinomio
        - *console:* Valor booleano por defecto en False. Si esta en True permite ver mensajes de los procesos en la consola
    - **Valor de Retorno:** Una lista con los coeficientes del polinomio
    """

    # Crear la lista de salida
    coeficientes = np.array([1, 0])

    # Transponer las raíces
    raices_ = np.array(raices)
    raices_ = -1*raices_

    # Crear un ciclo para recorrer las raíces del polinomio
    for r in raices_:
        
        # Multiplicar la raíz por la lista de coeficientes
        producto = r*coeficientes
        producto = np.insert(producto, 0, 0)
        producto = producto[:-1]

        # Sumar el producto con los coeficientes
        coeficientes = producto +  coeficientes

        # Agregar un nuevo cero a la lista de coeficientes
        coeficientes = np.append(coeficientes, 0)
        if console:
            print("Coeficientes: ", coeficientes)
            print("Producto: ", producto)

    # Eliminar el último elemento
    coeficientes = coeficientes[:-1]
    return coeficientes

# Definir la función que calcula el polinomio de Lagrange
def polyLagrange(nodos, k):
    
    # Crear un vector de nodos excepto el k
    headNodos = nodos[:k]
    tailNodos = nodos[k+1:]
    nodos_ = np.concatenate((headNodos, tailNodos))
    
    # Calcular la cantidad de nodos
    n = len(nodos_)
    
    # Calcular el denominador
    k_ = nodos[k]
    denominador = 1
    for i in range(n):
        denominador = denominador * (k_ - nodos_[i])            
    
    # Usar la función de Multiplicación Sintética para calcular el polinomio
    return MultiplicacionSintetica(nodos_) / denominador

# Definir la función que calcula el polinomio de Lagrange al cuadrado
def polyLagrangeCuadrado(nodos, k):
    
    # Crear un vector de nodos excepto el k
    headNodos = nodos[:k]
    tailNodos = nodos[k+1:]
    nodos_ = np.concatenate((headNodos, tailNodos))
    nodos_ = np.concatenate((nodos_, nodos_))
    
    # Calcular la cantidad de nodos
    n = len(nodos_)
    
    # Calcular el denominador
    k_ = nodos[k]
    denominador = 1
    for i in range(n):
        denominador = denominador * (k_ - nodos_[i])            
    
    # Usar la función de Multiplicación Sintética para calcular el polinomio
    return MultiplicacionSintetica(nodos_) / denominador

# Definir la funcion
def polyInterpoLagrange(nodos):
    
    # Obtener los vectores de x e y
    x, y = nodos
    
    # Obtener la cantidad de nodos
    nX = len(x)
    nY = len(y)
    
    if nX != nY:
        print("x no tiene la misma cantidad de elementos que y")
        return
    
    # Crear un ciclo para calcular cada polinomio de Lagrange
    polyInterpola = np.zeros(nX)
    
    for k in range(nX):
        polyInterpola += y[k] * polyLagrange(x, k)        
        
    # Devolver el resultado
    return polyInterpola

# Definir el método que calcula los números primos
def Primos(nMin, nMax):
    # Crear un ciclo para determinar los números primos
    Primos_ = np.array([])

    for n in range(nMin, nMax + 1):
        
        # Crear un ciclo interno para contar los divisisores del número n
        divisores = 0
        for m in range(1, n+1):
            
            # Preguntar si m es divisor de n
            if (n % m == 0):
                divisores += 1
            
            # Si divisores es mayor que 2, n no es primo
            if divisores > 2:
                break
        
        # Si la cantidad de divisores de n es dos es un número primo
        if divisores==2:
            Primos_ = np.append(Primos_, n)
    
    # Devolver los números primos
    return Primos_

# Definir la función evalPoly
def evalPoly(coeficientes, dominio, console=False):

    """
    # ***Función:*** evalPoly
    - **Descripcion:** Evalua un polinomio dado una lista de puntos y sus coeficientes
    - **Parámetros:**
        - *coeficientes:* Lista de coeficientes del polinomio de la forma $[a_n, a_{n-1}, \\ldots, a_1, a_0]$
        - *dominio:* Lista de valores a evaluar en el polinomio
        - *console:* Valor booleando. Permite mostrar los mensajes de procesos
    - **Valor de Retorno:** Una lista de imagenes del dominio evaluado en el polinomio
    """

    # Obtener el grado del polinomio
    grado = np.shape(coeficientes)[0] - 1

    # Crear la imagen    
    imagen = np.empty((0, grado))
    dominio_ = np.array(dominio)
    
    # Crear un ciclo para evaluar el polinomio
    imagen = np.empty((0, grado))

    for x in dominio_:
        suma = 0
        for n,c in enumerate(coeficientes):
            suma = suma + c*x**(grado-n)

        imagen = np.append(imagen, suma)

    # Devolver la imagen
    return imagen

# Definir la derivada de un polinomio
def DerivadaPolinomio(coeficientes, console= False):
    
    # Determinar el grado del polinomio
    grado = len(coeficientes) - 1
    
    #Polinomio Derivado
    polyDer = []
    
    # Crear un ciclo para derivar el polinomio
    for i,c in enumerate(coeficientes):
        polyDer.append(c*(grado - i))
        
    # Devolver el polinomio derivado
    return np.array(polyDer[:-1])

# Definir el método de Neville
def Neville(nodos, console=False):
    
    # Obtener los valores de X e Y de los nodos
    nodos_x, nodos_y = nodos
    
    # Definir el x simbólico
    x = sympy.symbols('x')
    
    # Definir la matriz de elementos y los polinomios de grado cero (0)
    Neville_ = np.array([nodos_y])
    
    # Calcular los polinomios de grado uno (1)
    lenX = len(nodos_x)
    for n in range(lenX-1):
        
        P_ = np.array([])
        P = Neville_[n]
        lenP = len(P) - n
        
        if console:
            print("--- Polinomios de Grado (",n + 1,") --------------------------------------------------------------------------")
        
        for k in range(lenP - 1):
            i = k
            j = i + n + 1
            denominador = nodos_x[j] - nodos_x[i]
            xi = nodos_x[i]
            xj = nodos_x[j]
            poly_ = ((x - xi)*P[j] - (x - xj)*P[j-1]) / denominador
            poly = sympy.expand(poly_)
            
            if console:
                # print("P(",i,",",j,") = (1 /",denominador,") ( x -",xi,")(",P[j],") - ( x -",xj,")(",P[j-1],") = ", poly)
                print("P(",i,",",j,") = ", poly)
            
            P_ = np.append(P_, poly)
            
        # Asegurar la dimensión de P_anterior
        dimP_ = len(P_)
        Q_ = np.pad(P_, (lenX - dimP_, 0), 'constant')
        Neville_ = np.vstack((Neville_, Q_)) 
        
    expMath = 'P_'+str(n) + '(x) = '+ sympy.latex(Neville_[-1][-1])
    
    return Neville_, expMath 


# Definir el método de las diferencias divididas
def DiferenciasDivididas(nodos, console=False):

    """
    #### ***Función:*** DiferenciasDivididas
    - **Descripción:** Esta función calcula un polinomio de grado $n$ dados un conjunto de nodos
    - **Parámetros:**
        - *nodos:* Es una dupla de la forma (x[...], y[...])
        - *console:* Valor de tipo booleano para mostrar mensajes por consola
    - **Valor de Retorno:** Devuelve un diccionario con todos los polinomios $P_0,\ldots,P_n$
    """
    
    # Obtener los nodos de X e Y
    x, y = nodos
    
    # Definir un diccionario para los polinomios y crear el primer polinomio
    Poly = {"P0": [y[0]]}
    
    # Crear un ciclo para recorrer todos los nodos
    cantNodos = len(x)
    for n in range(1, cantNodos):
        
        # Obtener el polinomio anterior
        if console:
            print("Iteración ",n," ----------------------------------------------------------------")
        P_ = Poly["P"+str(n-1)]
        x_n = np.array([x[n]])

        if console:
            print("P"+str(n), P_)
            print("x_n: ", x_n)
        
        # Evaluar el polinomio P_{n-1} en x_n
        y_0 = evalPoly(P_, x_n)
        
        # Calcular la productoria
        Prod = 1
        raices = np.array([])
        for k in range(n):
            Prod = Prod * (x[n] - x[k])
            raices = np.append(raices, x[k])

        # Definir el polinomio n-1
        ms_ = np.insert(P_, 0, 0)
        a_n = (y[n] - y_0) / Prod
        ms = a_n * MultiplicacionSintetica(raices)
        
        if console:
            print("a_n: ", a_n)

        # Sumar los dos polinomios
        P_n = ms_ + ms
        
        if console:
            print("P_n: ", P_n)
        
        # Agregar el polinomio encontrado al diccionario
        Poly["P"+str(n)] = P_n
    
    # Devolver el polinomio
    return Poly

# Definir el método de las diferencias divididas
def NodosEquiespaciados(nodos, x_0, h, console=False):

    """
    #### ***Función:*** NodosEquiespaciados
    - **Descripción:** Esta función calcula un polinomio de grado $n$ dados un conjunto de nodos
    - **Parámetros:**
        - *nodos:* Es una dupla de la forma (x[...], y[...])
        - *x_0:* Valor inicial en x
        - *h:* Espacio o distancia entre los nodos x_i y x_{i-1}
        - *console:* Valor de tipo booleano para mostrar mensajes por consola
    - **Valor de Retorno:** Devuelve un diccionario con todos los polinomios $P_0,\ldots,P_n$
    """
    
    # Obtener los nodos  Y
    y = nodos
    
    # Definir un diccionario para los polinomios y crear el primer polinomio
    Poly = {"P0": [y[0]]}
    
    # Crear un ciclo para recorrer todos los nodos
    cantNodos = len(y)
    for n in range(1, cantNodos):
        
        # Obtener el polinomio anterior
        if console:
            print("Iteración ",n," ----------------------------------------------------------------")
        P_ = Poly["P"+str(n-1)]
        x_n = np.array([x_0 + n*h])

        if console:
            print("P"+str(n-1), P_)
            print("x_n: ", x_n)
        
        # Evaluar el polinomio P_{n-1} en x_n
        y_0 = evalPoly(P_, x_n)
        
        # Calcular la productoria
        Factorial = 1
        raices = np.array([])
        for k in range(n):
            Factorial = Factorial * (k+1)
            raices = np.append(raices, x_0 + k*h)
            
        if console:
            print("Raíces: ", raices)
            print("Factorial: ",Factorial)

        # Definir el polinomio n-1
        ms_ = np.insert(P_, 0, 0)
        a_n = (y[n] - y_0) / (Factorial * h**(k+1))
        ms = a_n * MultiplicacionSintetica(raices)
        
        if console:
            print("a_n: ", a_n)

        # Sumar los dos polinomios
        P_n = ms_ + ms
        
        if console:
            print("P_n: ", P_n)
        
        # Agregar el polinomio encontrado al diccionario
        Poly["P"+str(n)] = P_n
    
    # Devolver el polinomio
    return Poly

# Definir el método de interpolación de Hermite
def Hermite(nodos, console=False):
    
    """
    #### ***Función:*** Hermite
    - **Descripción:** Permite calcular el polinomio interpolador de un conjunto de nodos (x[...], y[...], y'[...]) a través del método de los polinomios de Hermite
    - **Parámetros:** 
        - *nodos:* Tupla de los nodos que se quieren interpolar (x[...], y[...], y'[...]) 
        - *console:* Valor boolean que permite mostrar el proceso de cálculo si esta en True
    - **Valor de Retorno:** Un vector con los coeficientes del polinomio interpolador
    """
    
    # Obtener los nodos x,y y dy
    nodos_x, nodos_y, nodos_dy = nodos
    
    # Calcular la cantidad de nodos a trabajar y el grado del polinomio base
    cantNodos = len(nodos_x)
    grado = cantNodos - 1
    
    # Definir el resultado del polinomio interpolador 2*grado + 1
    PolyInterp = np.array([0*i for i in range(2*grado + 2)])
    
    # Crear un ciclo para calcular cada uno de los polinomios de Lagrange, H y \hat(H)
    for k in range(cantNodos):
        
        # Calcular el polinomio de lagrange en k
        if console:
            print("")
            print("--- Para el nodo k =",k,"--------------------------------------------")
        
        L = polyLagrange(nodos_x, k)
        L2 = polyLagrangeCuadrado(nodos_x, k)  
        
        if console:
            print("L(x) =",L)
            print("L2(x) =", L2)
        
        # Calcular la derivada del polinomio de Lagrange y su evaluación en el nodo x_k
        LDer_K = DerivadaPolinomio(L)
        if console:
            print("Dx(L) =", LDer_K)
            
        xk = np.array([nodos_x[k]])
        vLDer_K = evalPoly(LDer_K, xk)[0]
        
        if console:
            print("x_"+str(k)+" =",nodos_x[k],"; Dx(x_k) =", vLDer_K)
        
        # Usar multiplicación sintética para multiplicar L^2 por [1-2(x-xk)L'(x)] y Hallar H(x)
        r = -1/(2*vLDer_K) - nodos_x[k]
        H = r*L2
            
        H = np.insert(H, 0, 0)
        H = np.append(L2, 0) + H
        
        # Se debe multiplicar el resultado por -2*vLDer_K por la factorización realizada en la Multiplicación Sintética
        H = -2*vLDer_K*H
        
        if console:
            print("H: ", H)
        
        # Usar multiplicación sintética para multiplicar L^2 por (x-x_k) y Hallar \hat(H)(x)
        r = -nodos_x[k]
        H_ = r*L2
        H_ = np.insert(H_, 0, 0)
        H_ = np.append(L2, 0) + H_
        
        if console:
            print("H_: ", H_)
        
        # Crear la suma del polinomio interpolador y_k*H(x) + dy_k*\hat(H)(x)
        PolyInterp = PolyInterp + nodos_y[k]*H + nodos_dy[k]*H_
        
    # Devolver el resultado del polinomio interpolador
    return PolyInterp