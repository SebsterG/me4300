from dolfin import *


mesh = Mesh("step.xml")
f = File ("mesh_step.pvd")
f << mesh
#plot(mesh,interactive=True)
V = VectorFunctionSpace(mesh,"CG", 2)
Q = FunctionSpace(mesh,"CG", 1)
VQ = V * Q # The Mixed space , alternative writing :
#VQ = M i x e d F u n c t i o n S p a c e ([V , Q])
u , p = TrialFunctions ( VQ )
v , q = TestFunctions ( VQ )

class Top(SubDomain):
	def inside(self, x, on_boundary):
		return near(x[1], 0.5)
class Bottom(SubDomain):
	def inside(self,x,on_boundary):
		return on_boundary and not (near(x[1], 0.5) or near(x[0],-0.5) or near(x[0],0.5))
class Left(SubDomain):
	def inside(self,x,on_boundary):
		return near(x[0],-0.5)
class Right(SubDomain):
	def inside(selv,x,on_boundary):
		return near(x[0],0.5)


top=Top()
nos = Bottom()
left = Left()
right = Right()
bound = FacetFunction("size_t",mesh)
bound.set_all(0)
left.mark(bound,1)
right.mark(bound,2)
nos.mark(bound,3)
#ds = ds[bound]
ds = Measure("ds",subdomain_data=bound)
#plot(bound);interactive();




#plot(bound); interactive()
u0 = -1.0

bc1 = DirichletBC(VQ.sub(0), (0,0), nos)
bc2 = DirichletBC(VQ.sub(0), (u0,0), top)
bc3 = DirichletBC(VQ.sub(1), 0.0, left)
bc4 = DirichletBC(VQ.sub(1), 0.0, right)
bcs = [bc1,bc3,bc4,bc2]
mu = 100.0
F = mu * inner(grad(v), grad(u))* dx - inner(div(v) ,p)* dx - inner(q, div(u))* dx
up_ = Function(VQ)
solve(lhs(F) == rhs(F), up_, bcs)
u_, p_= up_.split()
n = FacetNormal(mesh)
#plot(p_);interactive()
############## FLUX INTEGRAL

Integral_left = assemble(dot(u_,n)*ds(1))
Integral_right = assemble(dot(u_,n)*ds(2))
print "Flux left: ",Integral_left
print "Flux right: ",Integral_right

#plot(u_,interactive=True ,range_min=0.0, range_max = 0.05)

############## Normal stress

tau = -p_*Identity(2) 
print "Normal stress: " , assemble(dot(dot(tau,n),n)*ds(3))      





Q2 = FunctionSpace(mesh, "CG", 1)
curl_u = project(curl(u_), Q2)
#plot(curl_u,interactive=True)

psi = Function(Q2)
p = TrialFunction(Q2)
q = TestFunction(Q2)
grad_psi = as_vector((-u_[1],u_[0]))
G = inner(grad(p), grad(q))*dx - q*dot(grad_psi,n)*ds - inner(curl(u_),q)*dx

solve(lhs(G)==rhs(G), psi)

pmin = psi.vector().array().argmin()
pmax = psi.vector().array().argmax()
# argsort 
xx = interpolate(Expression("x[0]"), Q2)
yy = interpolate(Expression("x[1]"), Q2)
xmin = xx.vector()[pmin]
ymin = yy.vector()[pmin]

xmax = xx.vector()[pmax]
ymax = yy.vector()[pmax]

if u0>0:
	print "Center eddy: x: %.4f, y: %.4f " %(xmin, ymin)
if u0<0:
	print "Center eddy: x: %.4f, y: %.4f " %(xmax, ymax)

f = File ("psi_1.pvd")
f << psi
#plot(psi); interactive();







