from dolfin import *


mesh = UnitSquareMesh(100,100)
V = VectorFunctionSpace(mesh,"CG", 2)
Q = FunctionSpace(mesh,"CG", 1)
VQ = V * Q # The Mixed space , alternative writing :
#VQ = M i x e d F u n c t i o n S p a c e ([V , Q])
u , p = TrialFunctions ( VQ )
v , q = TestFunctions ( VQ )

class Top(SubDomain):
	def inside(self, x, on_boundary):
		return near(x[1], 1.0)

class Nos(SubDomain):
	def inside(self,x,on_boundary):
		return on_boundary and not near(x[1],1.0)

top = Top()
nos = Nos()
bound = FacetFunction("size_t", mesh)
bound.set_all(0)
top.mark(bound, 1)
nos.mark(bound, 2)
#plot(bound); interactive()

bc1 = DirichletBC(VQ.sub(0), (0,0), nos)
bc2 = DirichletBC(VQ.sub(0), (1.0,0), top)
bcs = [bc1,bc2]
mu = 100.0
F = mu * inner(grad(v), grad(u))* dx - inner(div(v) ,p)* dx -inner(q, div(u))* dx
up_ = Function ( VQ )
solve(lhs(F) == rhs(F), up_, bcs)
u_, p_= up_.split()
#plot(u_);interactive()

bc01 = DirichletBC(V, (0, 0), nos)
bc02 = DirichletBC(V, (1.0, 0), top)
uu = project(u_, V, bcs=[bc01,bc02])
#plot(uu)
#interactive()
psi = Function(Q)
p = TrialFunction(Q)
q = TestFunction(Q)
solve(inner(grad(p), grad(q))*dx == inner(curl(u_), q)*dx, psi, bcs=[DirichletBC(Q, 0, "on_boundary")])
pa = psi.vector().array().argmin()
# argsort 
xx = interpolate(Expression("x[0]"), Q)
yy = interpolate(Expression("x[1]"), Q)
xm = xx.vector()[pa]
ym = yy.vector()[pa]
plot(psi);interactive()
print "Center eddy: x: %.4f, y: %.4f " %(xm, ym)
f = File ("cavity_psi.pvd")
f << psi



