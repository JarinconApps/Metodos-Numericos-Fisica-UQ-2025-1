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