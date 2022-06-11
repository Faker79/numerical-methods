
import numpy as np


def RK4_ORDEN_2(f_1,f_2,ci_1,ci_2,a,b,n_p):
    """

    Parameters
    ----------
    f_1 : funcion
        Cambio de variable y´= x.
    f_2 : funcion
        Despeje del diferencial de edo de 2 orden.
    ci_1 : Int
        Condicion inicial.
    ci_2 : Int
        Condicion Inicial.
    a : Int
        Limite Inferior.
    b : Int
        Limite superior.
    n_p : Int
        Numero de iteraciones.

    Returns nodo_x,nodo_y,nodo_t
    -------
    None.

    """

    h = (a-b)/n_p
    nodo_x = [ci_1]
    nodo_y = [ci_2]
    nodo_t = [a]
    
    for z in range(n_p):
        
        #Pendientes para primera ecuacion con f_1
        k1 = f_1(nodo_t[z],nodo_x,nodo_y[z])
        k2 = f_1(nodo_t[z]+(h/2), nodo_x[z]+k1*(h/2),nodo_y[z]+(h/2))
        k3 = f_1(nodo_t[z]+(h/2), nodo_x[z]+k2*(h/2),nodo_y[z]+(h/2))
        k4 = f_1(nodo_t[z]+h, nodo_x[z]+k3*h,nodo_y[z]+h)
        
     
        #Pendientes segunda ecuacion f_2
        m1 = f_2(nodo_t[z],nodo_x[z],nodo_y[z])
        m2 = f_2(nodo_t[z]+(h/2),nodo_x[z]+(h/2), nodo_y[z]+m1*(h/2))
        m3 = f_2(nodo_t[z]+(h/2),nodo_x[z]+(h/2), nodo_y[z]+m2*(h/2))
        m4 = f_2(nodo_t[z]+h,nodo_x[z]+h, nodo_y[z]+m3*h)
        
       
        #Actualizando nodos
        nodo_x.append(nodo_x[z]+(h/6)*(k1+2*k2+2*k3+k4))
        nodo_y.append(nodo_y[z]+(h/6)*(m1+2*m2+2*m3+m4))
        nodo_t.append(nodo_t[z]+h)
        
    return nodo_x,nodo_y,nodo_t
