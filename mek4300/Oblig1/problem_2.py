from dolfin import *

L = 4 
for dens in [100]:
	mesh = IntervalMesh(dens, 0, L)

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
	vf, vh = TestFunctions(VV)

	f_h_ = Function(VV)
	f_,h_ = split(f_h_)
	#h_ = interpolate(Expression("0"),V)

	F1 = -inner(grad(h_),grad(vh))*dx + f_*h_.dx(0)*vh*dx + (1-h_*h_)*vh*dx
	F2 = h_*vf*dx - f_.dx(0)*vf*dx        
	F = F1+F2 

	#J = derivative(F, f_h_, fh)
	solve(F==0,f_h_,bcs)	
	plot(f_); interactive();
	"""
	problem = NonlinearVariationalProblem(F, f_h_, bcs, J)
	solver  = NonlinearVariationalSolver(problem)
	prm = solver.parameters
	prm['newton_solver']['absolute_tolerance'] = 1E-13
	prm['newton_solver']['relative_tolerance'] = 1E-13
	prm['newton_solver']['maximum_iterations'] = 50
	prm['newton_solver']['relaxation_parameter'] = 1.0
	solver.solve()
	plot(f_); interactive();
	"""
