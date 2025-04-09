# Definir la función evalPoly
import numpy as np

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