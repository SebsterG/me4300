from dolfin import *
L = 4
mesh = IntervalMesh(50, 0, L)
V = FunctionSpace(mesh, 'CG', 1)
VV = V*V #Making mixed function space
vf, vh = TestFunctions(VV) #different testfunc for both equations
fh = TrialFunction(VV) #To the picard
f, h = split(fh)

def newton():
	#Mind Constant((0,0)) due to mixed function space
	bc0 = DirichletBC(VV, Constant((0,0)), "std::abs(x[0]) < 1e-10")

	#Fordi jeg ikke tar med sub sa gjelder bc for begge
	bc1 = DirichletBC(VV.sub(1), Constant(1), "std::abs(x[0]-%d) < 1e-10" % L)
	#Mind sub(1) cause this only yields for F' (H)

	#First "guess" of the solution
	fh_ = interpolate(Expression(("x[0]","x[0]")), VV)
	f_, h_ = split(fh_)

	F = h_*vf*dx -f_.dx(0)*vf*dx - inner(grad(h_),grad(vh))*dx + f_*h_.dx(0)*vh*dx + vh*dx - h_**2*vh*dx

	solve(F==0, fh_, [bc0, bc1])


	viz = plot(f_)
	viz.write_png('146newton1')
	viz = plot(h_)
	viz.write_png('146newton2')
	interactive()

def picard():
	#Mind Constant((0,0)) due to mixed function space
	bc0 = DirichletBC(VV, Constant((0,0)), "std::abs(x[0]) < 1e-10")

	#Mind sub(1) cause this only yields for F' (H)
	bc1 = DirichletBC(VV.sub(1), Constant(1), "std::abs(x[0]-%d) < 1e-10" % L)

	#Previous known functionvalue
	fh_ = interpolate(Expression(("x[0]","x[0]")), VV)
	f_, h_ = split(fh_)

	F = h*vf*dx - f.dx(0)*vf*dx - inner(grad(h),grad(vh))*dx + f_*h.dx(0)*vh*dx + vh*dx - h*h_*vh*dx
	#h*h_ is the linearsation of

	k = 0
	fh_1 = Function(VV) #The new unknown function
	error = 1

	while k < 100 and error > 1e-12:
		solve(lhs(F) == rhs(F), fh_, [bc0, bc1])
		#solve godtar ikke vf*dx arg..
		error = errornorm(fh_, fh_1)
		fh_1.assign(fh_) #updates the unknown for next it.
		print "Error = ", k, error
		k += 1

	viz = plot(f_)
	viz.write_png('146picard1')
	viz = plot(h_)
	viz.write_png('146picard2')
	interactive()

newton()
picard()
