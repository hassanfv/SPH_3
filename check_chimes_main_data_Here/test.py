from assimulo.solvers import CVode
from assimulo.problem import Explicit_Problem

def model(t, y):
    return np.array([-1000 * y[0] + 3000 - 2000 * t])

y0 = [0]
t0 = 0
problem = Explicit_Problem(model, y0, t0)
simulator = CVode(problem)
t, y = simulator.simulate(1)

