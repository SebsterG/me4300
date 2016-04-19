from dolfin import *
from mshr import *



#mesh = UnitCircleMesh(100)
"""domain = Ellipse(Point(0,0),a,b,density)
mesh = generate_mesh(domain,density)
"""
dens = 100
c1 = Circle(Point(0,0),5,dens)
c2 = Circle(Point(0,0),1,dens)
domain = c1-c2
mesh  = generate_mesh(domain,dens)



def boundary(x,on_boundary):
	return on_boundary

V = FunctionSpace(mesh, "Lagrange",1)

u = TrialFunction(V)
v = TestFunction(V)

nos = Expression("0")

bcs = DirichletBC(V,nos,boundary)



u_exact = Expression("-(1./4)*(25-(x[0]*x[0]+x[1]*x[1])+24*log(5/sqrt(x[0]*x[0]+x[1]*x[1]))/log(1./5))")
a = inner(nabla_grad(u), nabla_grad(v))*dx
f = Constant(-1.0)
L = f*v*dx
u = Function(V)
solve(a==L, u, bcs)
plot(u)
interactive()
u_e = project(u_exact, V, bcs=bcs)
u_error = errornorm (u_e, u, degree_rise =0)
print u_error 
plot(u_e, title="Exact velocity for concentric annulus")
interactive()






