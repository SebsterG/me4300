from dolfin import *
import numpy as np
#set_log_active(False)

mesh = Mesh("von_karman_street.xml")
V = VectorFunctionSpace(mesh,"CG", 2)
Q = FunctionSpace(mesh,"CG", 1)
VQ = V*Q
#VQ = V * Q # The Mixed space , alternative writing :
#VQ = M i x e d F u n c t i o n S p a c e ([V , Q])
#u , p = TrialFunctions(VQ)


class Circle(SubDomain):
	def inside(self, x, on_boundary):
		return on_boundary and not (near(x[0],0) or near(x[0],2.2) or near(x[1],0) or near(x[1],0.41))

class Right(SubDomain):
	def inside(self,x,on_boundary):
		return near(x[0], 2.2)

class Left(SubDomain):
	def inside(self,x,on_boundary):
		return near(x[0], 0)
class Nos(SubDomain):
	def inside(self,x,on_boundary):
		return near(x[1], 0) or near(x[1],0.41)


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
Um = 0.3
H = 0.41
inlet = Expression(("4*Um*x[1]*(H-x[1])/(H*H)", "0"),Um = Um,H=H)
#U_mean = project(Expression(("4*Um*(H/2.0)*(H/2.0)/(H*H)","0"),Um=Um,H=H),V)
U_mean = 2.0*Um/3.0
v_theta = Expression(("0","0"))

bc1 = DirichletBC(VQ.sub(0), inlet , left)
bc2 = DirichletBC(VQ.sub(0), (0, 0), nos)
bc3 = DirichletBC(VQ.sub(0), v_theta , circle)
#bc4 = DirichletBC(Q, 0.0,right)
bcs = [bc1,bc2,bc3]


up = Function(VQ)
u, p = split(up)
vq = TestFunction(VQ)
v, q = split(vq)

nu = 0.001
nu = Constant(nu)

a = -inner(dot(grad(u),u), v)*dx +  nu*inner(grad(u), grad(v))*dx - dot(p,div(v))*dx
L = dot(q,div(u))*dx
F = a-L

#ufile = File("results/velocity.pvd")
#pfile = File("results/pressure.pvd")

solve(F==0,up,bcs)
u_,p_ = up.split(True)
#plot(u_,rescale = True)
#plot(p_,rescale = True)
interactive()
R = VectorFunctionSpace(mesh, 'R', 0)
c = TestFunction(R)
tau = -p_*Identity(2)+nu*(grad(u_)+grad(u_).T)
n = FacetNormal(mesh)
ds = ds[bound]
forces = -assemble(dot(dot(tau, n), c)*ds(1)).array()*2/(U_mean**2*0.1)
print "Drag = {}, Lift = {}".format(*forces)
p1 = p_.compute_vertex_values()
x0 = np.where(mesh.coordinates()[:,0]==0.15)
x1 = np.where(mesh.coordinates()[:,0]==0.25)
diff_p = p1[x0[0]]-p1[x1[0]]
print diff_p
#print pressure(array([0.15,0.20]))-pressure(array([0.25,0.20]))

#print "diff p",(p_.vector().array([0.15,0.20])-p_.vector().array([0.25,0.20]))

#interactive()"""
