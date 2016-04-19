from dolfin import *

mesh = Mesh("Cylinder_in_square.xml")
V = VectorFunctionSpace(mesh,"CG", 2)
Q = FunctionSpace(mesh,"CG", 1)
VQ = V * Q # The Mixed space , alternative writing :
#VQ = M i x e d F u n c t i o n S p a c e ([V , Q])
u , p = TrialFunctions(VQ)
v , q = TestFunctions(VQ)

class Circle(SubDomain):
	def inside(self, x, on_boundary):
		return on_boundary and not (near(x[0],0) or near(x[0],6) or near(x[1],-1) or near(x[1],1))

class Right(SubDomain):
	def inside(self,x,on_boundary):
		return near(x[0], 0)

class Left(SubDomain):
	def inside(self,x,on_boundary):
		return near(x[0], 6)
class Nos(SubDomain):
	def inside(self,x,on_boundary):
		return near(x[1], 1) or near(x[1],-1)


circle = Circle()
left = Left()
right = Right()
nos = Nos()
bound = FacetFunction("size_t", mesh)
bound.set_all(0)
nos.mark(bound, 3)
circle.mark(bound, 1)
left.mark(bound,4)
right.mark(bound,2)
#plot(bound); interactive()
inlet = Expression(("1-x[1]*x[1]", "0"))
v_theta = Expression(("-4*x[1]","4*(x[0]-1)"))

bc1 = DirichletBC(VQ.sub(0), inlet , left)
bc2 = DirichletBC(VQ.sub(0), (0, 0), nos)
bc3 = DirichletBC(VQ.sub(0), v_theta , circle)

bcs = [bc1,bc2,bc3]


mu = 1.0
rho = 1.0
up_ = Function(VQ)
u_ , p_ = split(up_)
#F = -mu * inner(grad(v), grad(u))* dx + inner(div(v) ,p)* dx - inner(q, div(u))* dx
F = -mu * inner(grad(v), grad(u_))* dx + inner(div(v) ,p_)* dx - rho*inner(dot(u_,grad(u_)),v)*dx - inner(q, div(u_))* dx


solve(F==0, up_, bcs,
      solver_parameters={'newton_solver':{'maximum_iterations': 10,
                                          'error_on_nonconvergence': False}})
                                          #'linear_solver': 'LU'}})

plot(u_, title='Velocity') ; interactive()
plot(p_, title='Pressure')

R = VectorFunctionSpace(mesh, 'R', 0)
c = TestFunction(R)
tau = -p_*Identity(2)+mu*(grad(u_)+grad(u_).T)
n = FacetNormal(mesh)
ds = ds[bound]
forces = -assemble(dot(dot(tau, n), c)*ds(1))
print "Drag = {}, Lift = {}".format(*forces)
