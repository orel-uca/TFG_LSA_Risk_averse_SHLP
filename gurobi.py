import gurobipy as gp
from gurobipy import GRB
# Crear modelo
model = gp.Model("ejemplo")
# Variables
x = model.addVar(vtype=GRB.BINARY, name="x")
y = model.addVar(vtype=GRB.BINARY, name="y")
# Restricción
model.addConstr(x + y <= 1, "r1")
# Función objetivo
model.setObjective(x + 2 * y, GRB.MAXIMIZE)
# Resolver
model.optimize()
# Mostrar solución
for v in model.getVars():
    print(f'{v.varName}: {v.x}')
print(f'Objetivo: {model.objVal}')