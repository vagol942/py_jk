import math
import simulations.riboswitchPar as ribo
from multiprocessing import Process

CORES = 2
CELL_NUMBER = 1000

def main():
    p = Process(target=riboswitchRun, args(math.ceil(CELL_NUMBER/CORES)))
    p.start()
    np.savetxt('Riboswitch/RIBOSWITCH__T_s('+ str(N_cells)+'cells).txt', T_s, delimiter=',')
    np.savetxt('Riboswitch/RIBOSWITCH__dens_T('+ str(N_cells)+'cells).txt', densidad[0], delimiter=',')
    np.savetxt('Riboswitch/RIBOSWITCH__dens_P('+ str(N_cells)+'cells).txt', densidad[1], delimiter=',')
    np.savetxt('Riboswitch/RIBOSWITCH__dens_R('+ str(N_cells)+'cells).txt', densidad[2], delimiter=',')

main()
