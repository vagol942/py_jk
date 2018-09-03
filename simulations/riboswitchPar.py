import simulations.master as ma
import numpy as np

riboswitchRun(num_cells):
    #RIBOSWITCH
    #Las siguientes ecuaciones diferenciales dictaminan la evolucion estocastica del sistema en cuestion.

    # (d/dt)Ta = KT + mum*R - mup*E*Ta - gammaT*Ta
    # (d/dt)P = KP*Ta - gammaP*P
    #(d/dt)R = mup*E*Ta - gammaR*R - mum*R

    #Definiciones de constantes
    KT = 20.0
    gammaT = 1/5.0
    gammaR = 1/5.0
    KP = 1.0
    gammaP = 1/30.0
    mum = 1/3.0
    mup = 1/20000.0
    #Molecula de entrada (nota: pensar como acoplar variaciones en el tiempo)
    E = 10000.0

    #Acciones de creacion o destruccion de Variables. Hay 7 eventos posibles
    actions = [[1.0, 0.0, 0.0], [1.0, 0.0, -1.0], [-1.0, 0.0, 1.0], [-1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, -1.0]]

    #Probabilidades para este sistema dadas las ecuaciones. Hay 7 eventos posibles.
    def darEventos(xs):
        return [KT, mum*xs[2], mup*E*xs[0], gammaT*xs[0], KP*xs[0], gammaP*xs[1], gammaR*xs[2]]


    #CAMBIAR AL GUSTO
    #Valores iniciales de cada compuesto
    initVals = [0.0, 0.0, 0.0]
    #horas
    hours = 10
    #resolucion minima en minutos, NO CAMBIAR POR EL MOMENTO.
    #nota: para hacer que esto funcione bien toca usar numeros imaginarios.
    dt = 1.0
    #Numero de celulas a simular
    N_cells = num_cells
    #nota: pensar como acoplar resolucion minima de moleculas "dm"

    #Se genera la simulacion para las condiciones dadas. Definiendo asi arreglos de tiempo y promedios estandarizados.
    #Junto con el cubo de informacion "many" que abarca toda la simulacion.
    T_s, s_vals, many = ma.sv_cells(N_cells, hours, dt, initVals, darEventos, actions)
    s_vals = np.matrix(s_vals)
    #Se definen las funciones de probabilidad no normalizadas
    densidad = ma.dens(many)

    return desindad
    
