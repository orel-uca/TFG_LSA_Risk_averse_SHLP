from gurobipy import Model, GRB, quicksum
import numpy as np
import os
import math 


def optimizareturnh(n, q, s, tiempo):

    modelo = 'Formulacionreturnh'
    print()
    print('========================================================')
    print(modelo)
    print('n= ', n)
    print('s= ', s)
    
    N = range(n)
    Q = range(q)
    Omega = range(s)

    dataP = np.zeros((n))
    dataW = np.zeros((n*n,4))

    # Nombre de los archivos (Ej: si n=25 y escenario=1 => flow_poisson50_apll25_1.txt)
    aux1 = np.loadtxt('./Datos_phub/' + str(n) + 'L.txt', max_rows=n)
    dataP[:] = aux1
    f = np.zeros((n))
    for i in N:
        f[i] = dataP[i]*0.5 #Costes fijos
        
    
    
            
    W = np.zeros((n, n))
    Ws = np.zeros((n,n,s))
    for scenarios in Omega:    
        aux1 = np.loadtxt('./Datos_phub/flow_poisson50_apll' + str(n) + '_' + str(scenarios+1) + '.txt', skiprows =1, max_rows=n*n)
        dataW[:, :] = aux1
        d = np.zeros((n, n))
        for i in range(n*n):
            Ws[int(dataW[i,0]-1), int(dataW[i,1]-1),scenarios] = dataW[i,2]*10 #Matriz de flujos
            d[int(dataW[i,0]-1), int(dataW[i,1]-1)] = dataW[i,3] #Matriz de distancias
       
    ##CALCULO LA MEDIA DE LOS FLUJOS PARA GENERAR UN ÚNICO ESCENARIO
    for i in N:
        for j in N:
            W[i,j]=sum(Ws[i,j,scenarios] for scenarios in Omega)/s

    
    alfa = [1, 0.9, 0.8, 0.7]
    beta = [0, 5, 15, 35]
    u = [50, 100, 200, GRB.INFINITY]



    
    #Flujo total con origen en i (todo lo que sale de i)
    O = np.zeros((n))
    for i in N:
        O[i]=W[i,:].sum()
    #Flujo total con destino en i (todo lo que entra en i)
    D = np.zeros((n))
    for i in N:
        D[i]=W[:,i].sum()

    # ------------------------------------------------------------------------
    #              CREAMOS LAS VARIABLES Y EL MODELO
    # ------------------------------------------------------------------------
    mo = Model("M1_incompleta")

    # ---------------------  VARIABLES ----------------------
    h = mo.addVars([(k) for k in N ], vtype=GRB.BINARY, name="h")
    x = mo.addVars([(i, k) for i in N for k in N], vtype=GRB.BINARY, name="x")
    z = mo.addVars([(k,m,q) for k in N for m in N for q in Q], vtype=GRB.BINARY, name="z")
    y = mo.addVars([(i, j, k, m, q) for i in N for j in N for k in N for m in N for q in Q], lb=0.0, name="y")


    # -----------------  FUNCION OBJETIVO  ------------------
    mo.setObjective((quicksum(f[k]*h[k] for k in N) + quicksum((O[i]+D[i])*d[i,k]*x[i,k] for i in N for k in N) + quicksum(d[k,m]*beta[q]*z[k,m,q] for k in N for m in N for q in Q) + quicksum(W[i,j] * d[k,m] * alfa[q] *y[i,j,k,m,q] for i in N for j in N for k in N for m in N for q in Q) ), GRB.MINIMIZE)

    # -------------------  RESTRICCIONES  -------------------
    mo.addConstrs((quicksum(x[i,k] for k in N) == 1 for i in N), name='R1')
    mo.addConstrs((x[i, k] <= h[k] for i in N for k in N), name='R2')
    mo.addConstrs((quicksum(y[i,j,k,m,q] for m in N for q in Q) == x[i,k] for i in N for j in N for k in N), name='R3')
    mo.addConstrs((quicksum(y[i,j,k,m,q] for k in N for q in Q) == x[j,m] for i in N for j in N for m in N), name='R4')
    mo.addConstrs((quicksum(z[k,m,q] for q in Q) >= h[k] + h[m] - 1 for k in N for m in N), name='R5')
    mo.addConstrs((quicksum(W[i,j] * y[i,j,k,m,q] for i in N for j in N) <= u[q] * z[k,m,q] for k in N for m in N for q in Q), name='R6')
    
    

    # ------------------------------------------------------------------------
    #               MODIFICAMOS PARÁMETROS Y RESOLVEMOS
    # ------------------------------------------------------------------------

    # Parámetros de Gurobi
    # mo.Params.NodeLimit = 1.0    # 1
    # mo.Params.PreCrush = 1       # 2 -
    #mo.Params.IntParam_Cuts = 0  # 3
    # mo.Params.VarBranch = 1      # 4 -
    # mo.Params.MIPFocus = 3       # 5 -
    #mo.Params.Heuristics = 0       # 6 -
    #mo.Params.Cuts = 0
    #mo.Params.Presolve = 0
    #mo.Params.PreCrush = 1

    mo.Params.timelimit = tiempo
    #mo.Params.Threads = 1
    mo.Params.MIPGap = 1e-6


    # ------- Resolvemos ---------
    mo.optimize()
    btime = mo.Runtime  # Tiempo total que ha tardado

    # ------------------------------------------------------------------------
    #                    SALIDAS QUE GUARDAMOS EN ARCHIVO
    # ------------------------------------------------------------------------

    #if not os.path.exists('Salidas'):
        #os.makedirs('Salidas')

    #fichero = 'Salidas/resultados.txt'
    #outfile = open(fichero, 'a')

    #if mo.status == GRB.Status.OPTIMAL:
        #print('{:10} {:4}  {:3}   {:8.1f}    {:7.1f}   {:5.2f}  {:7d}'.format(modelo, n, s, mo.ObjVal, btime, mo.MIPGap, int(mo.NodeCount)), file=outfile)
        #for k in N:
            #if h[k].X >0.9:
                #print( '{:3} '.format(k), file=outfile)
    #else:
        #if mo.SolCount > 0:
            #print('{:10} {:4}  {:3}   {:8.1f}    {:7.1f}   {:5.2f}  {:7d}'.format(modelo, n, s, mo.ObjVal, btime, mo.MIPGap, int(mo.NodeCount)), file=outfile)
            #if h[k].X >0.9:
                #print( '{:3} '.format(k), file=outfile)
        #else:
             #print('{:10} {:4}  {:3}   {:8.1f}    {:7.1f}   {:5.2f}  {:7d}'.format(modelo, n, s, mo.ObjVal, btime, mo.MIPGap, int(mo.NodeCount)), file=outfile)          
    #print('{:3} '.format( ' ' ), file=outfile)           
            
    #outfile.close()

    sol = np.zeros((n))
    for i in N:
        if h[i].X >0.9:
            sol[i]=1
        else: sol[i] = 0
    
    print(h[1].X)

    print(sol)
    return(sol)
    
    #for k in N:
    #Para escribir por pantalla los hubs que abrimos
        #print(k, ' = ', h[k].X)