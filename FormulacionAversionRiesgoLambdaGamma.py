from gurobipy import Model, GRB, quicksum
import numpy as np
import os
import math 


def optimizaLauAversionRiesgolambdagamma(n, q, s, lmbda, gamma, tiempo):

    modelo = 'FormulacionLauAversionRiesgolambdagamma'
    print()
    print('========================================================')
    print(modelo)
    print('lambda = ', lmbda)
    print('gamma = ', gamma)
    
    N = range(n)
    Q = range(q)
    Omega = range(s)
    
    #lmbda=0.5 #Cuanto mas cerca de 1 esté lambda más importancia se le da a la adversión al riesgo.
    #gamma=0.8

    dataP = np.zeros((n))
    dataW = np.zeros((n*n,4))
    
    p = np.zeros(s)
    for scenarios in Omega:
        p[scenarios]=1/s

    # Nombre de los archivos (Ej: si n=25 y escenario=1 => flow_poisson50_apll25_1.txt)
    aux1 = np.loadtxt('./Datos_phub/' + str(n) + 'L.txt', max_rows=n)
    dataP[:] = aux1
    f = np.zeros((n))
    for i in N:
        f[i] = dataP[i]*0.5 #Costes fijos
        


        
    W = np.zeros((n, n, s))
    for scenarios in Omega:    
        aux1 = np.loadtxt('./Datos_phub/flow_poisson50_apll' + str(n) + '_' + str(scenarios+1) + '.txt', skiprows =1, max_rows=n*n)
        dataW[:, :] = aux1
        d = np.zeros((n, n))
        for i in range(n*n):
            W[int(dataW[i,0]-1), int(dataW[i,1]-1),scenarios] = dataW[i,2]*10 #Matriz de flujos
            d[int(dataW[i,0]-1), int(dataW[i,1]-1)] = dataW[i,3] #Matriz de distancias

    
    alfa = [1, 0.9, 0.8, 0.7]
    beta = [0, 5, 15, 35]
    u = [50, 100, 200, GRB.INFINITY]



    
    O = np.zeros((n,s))
    #Flujo total con origen en i (todo lo que sale de i)
    for scenarios in Omega:
        for i in N:
            O[i, scenarios]=W[i,:,scenarios].sum()
    
    D = np.zeros((n,s))
    #Flujo total con desstino en i (todo lo que entra en i)
    for scenarios in Omega:
        for i in N:
            D[i, scenarios]=W[:,i,scenarios].sum()

    # ------------------------------------------------------------------------
    #              CREAMOS LAS VARIABLES Y EL MODELO
    # ------------------------------------------------------------------------
    mo = Model("M1_incompleta")

    # ---------------------  VARIABLES ----------------------
    
    h = mo.addVars([(k) for k in N ], vtype=GRB.BINARY, name="h")
    z = mo.addVars([(i, k, scenarios) for i in N for k in N for scenarios in Omega], vtype=GRB.BINARY, name="z")
    y = mo.addVars([(k,l,q, scenarios) for k in N for l in N for q in Q for scenarios in Omega], vtype=GRB.BINARY, name="y")
    x = mo.addVars([(i, j, k, l, q, scenarios) for i in N for j in N for k in N for l in N for q in Q for scenarios in Omega], lb=0.0, name="x")
    eta = mo.addVar(vtype=GRB.CONTINUOUS, name="eta")
    v = mo.addVars([(scenarios) for scenarios in Omega], lb = 0.0, name = "v")
    
    


    # -----------------  FUNCION OBJETIVO  ------------------
    mo.setObjective((quicksum(f[k]*h[k] for k in N) + lmbda*(eta+1/(1-gamma)*quicksum(p[scenarios]*v[scenarios] for scenarios in Omega)) + (1-lmbda)*quicksum(p[scenarios] *(quicksum((O[i, scenarios]+D[i, scenarios])*d[i,k]*z[i,k,scenarios] for i in N for k in N) + quicksum(d[k,l]*beta[q]*y[k,l,q,scenarios] for k in N for l in N for q in Q) + quicksum(W[i,j,scenarios]*d[k,l]*alfa[q]*x[i,j,k,l,q,scenarios] for i in N for j in N for k in N for l in N for q in Q )) for scenarios in Omega)), GRB.MINIMIZE)

    # -------------------  RESTRICCIONES  -------------------
    mo.addConstrs((v[scenarios]>= quicksum((O[i,scenarios]+D[i,scenarios])*d[i,k]*z[i,k,scenarios] for i in N for k in N)+quicksum(d[k,l]*beta[q]*y[k,l,q,scenarios] for k in N for l in N for q in Q)+ quicksum(W[i,j,scenarios]*d[k,l]*alfa[q]*x[i,j,k,l,q,scenarios]for i in N for j in N for k in N for l in N for q in Q)-eta for scenarios in Omega), name="R0")
    mo.addConstrs((quicksum(z[i,k,scenarios] for k in N) == 1 for i in N for scenarios in Omega), name='R1')
    mo.addConstrs((z[i, k, scenarios] <= h[k] for i in N for k in N for scenarios in Omega), name='R2')
    mo.addConstrs((quicksum(x[i,j,k,l,q, scenarios] for l in N for q in Q) == z[i,k,scenarios] for i in N for j in N for k in N for scenarios in Omega), name='R3')
    mo.addConstrs((quicksum(x[i,j,k,l,q, scenarios] for k in N for q in Q) == z[j,l,scenarios] for i in N for j in N for l in N for scenarios in Omega), name='R4')
    mo.addConstrs((quicksum(y[k,l,q,scenarios] for q in Q) >= h[k] + h[l] - 1 for k in N for l in N for scenarios in Omega), name='R5')
    mo.addConstrs((quicksum(W[i,j,scenarios] * x[i,j,k,l,q,scenarios] for i in N for j in N) <= u[q] * y[k,l,q,scenarios] for k in N for l in N for q in Q for scenarios in Omega), name='R6')
    
    

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

    if not os.path.exists('Salidas'):
        os.makedirs('Salidas')

    fichero = 'Salidas/resultados.txt'
    outfile = open(fichero, 'a')

    if mo.status == GRB.Status.OPTIMAL:
        print('{:10} {:4}  {:3}   {:8.1f}    {:7.1f}   {:5.2f}  {:7d}'.format(modelo, lmbda, gamma, mo.ObjVal, btime, mo.MIPGap, int(mo.NodeCount)), file=outfile)
        for k in N:
            if h[k].X >0.9:
                print( '{:3} '.format(k), file=outfile)
    else:
        if mo.SolCount > 0:
            print('{:10} {:4}  {:3}   {:8.1f}    {:7.1f}   {:5.2f}  {:7d}'.format(modelo, lmbda, gamma, mo.ObjVal, btime, mo.MIPGap, int(mo.NodeCount)), file=outfile)
            for k in N:
                if h[k].X >0.9:
                    print( '{:3} '.format(k), file=outfile)
        else:
             print('{:10} {:4}  {:3}   {:8.1f}    {:7.1f}   {:5.2f}  {:7d}'.format(modelo, lmbda, gamma, mo.ObjVal, btime, mo.MIPGap, int(mo.NodeCount)), file=outfile)          
    print('{:3} '.format( ' ' ), file=outfile)           


    outfile.close()
    
    #Para escribir por pantalla los hubs que abrimos
    #for k in N:
        #print(k, ' = ', h[k].X)
