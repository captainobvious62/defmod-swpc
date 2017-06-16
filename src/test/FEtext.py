#!/usr/bin/python
from dolfin import *

#define mesh and function space
mesh=Interval(20, 0, 1)
V = FunctionSpace(mesh, "Lagrange", 1)

#define left boundary
class BoundaryLeft(SubDomain):
    def inside(self, x, on_boundary):
        return on_boundary and \
            (near(x[0], 0.))

#define right boundary
class BoundaryRight(SubDomain):
    def inside(self, x, on_boundary):
        return on_boundary and \
            (near(x[0], 1.))

#define inputs
u_0 = Constant(0.)
u_1 = Constant(1.)

#initialize subdomains
boundaryLeft = BoundaryLeft()
boundaryRight = BoundaryRight()

#apply essential boundary conditions
bcs = [DirichletBC(V, u_0, boundaryLeft),
       DirichletBC(V, u_1, boundaryRight)
      ]

#define functions
u = TrialFunction(V)
v = TestFunction(V)
f = Constant(0.)
g = Constant(0.)
nu = Constant(1.)
ac = Constant(1.)
 
#define problem
a = nu*inner(grad(u), grad(v))*dx+ac*inner(grad(u), v)*dx
L = f*v*dx + g*v*ds

#solve problem
u = Function(V)
solve(a == L, u, bcs)
 
#plot solution
plot(u, interactive=True)
