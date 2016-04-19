from dolfin import *

L = 4 
mesh = IntervalMesh(100, 0, L)

V = FunctionSpace(mesh,"CG", 2)

def left(x,on_boundary):
	return near(x[0],0) and on_boundary
def right(x,on_boundary):
	return near(x[0],L) and on_boundary

VV = V*V

bc1 = DirichletBC(VV.sub(1), 0, left)
bc2 = DirichletBC(VV.sub(0), 0, left)
bc3 = DirichletBC(VV.sub(1), 1,right)
bcs = [bc1,bc2,bc3]

fh = TrialFunction(VV)
f,h = split(fh)

vf, vh = TestFunction(VV)


f_h = interpolate(Expression(("x[0]","x[0]")),VV)
f_, h_ = split(f_h)

F1 = -inner(grad(h),grad(vh))*dx + f_*h.dx(0)*vh*dx + (1-h*h_)*vh*dx
F2 = h*vf*dx - f.dx(0)*vf*dx
F = F1+F2 
eps = 1.0
tol = 1.0E-13
maxiter = 50
f_h_1 = Function(VV)
iter = 0

while eps > tol and iter < maxiter:
	
	iter += 1
	#solve(F==0,f_h_,bcs)
	solve(lhs(F)==rhs(F), f_h, bcs) 
	eps = errornorm(f_h, f_h_1,norm_type="l2",degree_rise=3)
	#plot(h);interactive()
	print "iteration: %d: norm: %.4f" %(iter,eps)
	f_h_1.assign(f_h)
plot(f_,interactive=True)
plot(h_,interactive=True)

