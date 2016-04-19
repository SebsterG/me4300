from dolfin import *
from mshr import *
import numpy as np

for density in[20,30,40]:
	a = 2
	b = 3
	problem = "Ecentric"
	if problem == "Ellipse":
		domain = Ellipse(Point(0,0),a,b)
		mesh = generate_mesh(domain,density)
		#plot(mesh)
		u_exact = Expression(" -(1./2.)*(a*a*b*b/(a*a+b*b))* ( 1 - ((x[0]*x[0]) /(a*a) )-((x[1]*x[1])/(b*b)) ) ",a=a, b=b)
		Q_exact = -(pi/4)*(a*a*a*b*b*b/(a*a+b*b) )
	if  problem == "Triangle":
		pts = [Point(-sqrt(3),0),Point(sqrt(3)),Point(0,3)]
		a = 2*sqrt(3)
		domain = Polygon(pts)
		mesh = generate_mesh(domain,density)
		#plot(mesh);interactive()
		u_exact = Expression("-1/(2*sqrt(3)*a)*(x[1]-0.5*a*sqrt(3))*(3*x[0]*x[0]-x[1]*x[1])",a=a)
		Q_exact = -a**4*sqrt(3)/320 
	if problem == "Ecentric":
		a = 5
		b = 1
		c1 = Circle(Point(0,0),a, density)
		c2 = Circle(Point(0.5,0.0),b, density)
		domain = c1-c2
		mesh  = generate_mesh(domain,density)
		c = 0.5 
		F = (a**2-b**2+c**2)/(2*c)
		M = sqrt(F**2-a**2)
		alpha = 0.5*ln((F+M)/(F-M))
		beta = 0.5*ln((F-c+M)/(F-c-M))
		Q_exact = (-pi/8.0) * (a**4-b**4 - 4*c**2*M**2/(beta-alpha) - 8*c**2*M**2* \
		np.sum([(i*exp(-(beta+alpha))/np.sinh(i*beta-i*alpha)) for i in range(1,1000)]) ) 

	def boundary(x,on_boundary):
		return on_boundary

	V = FunctionSpace(mesh, "CG",2)

	u = TrialFunction(V)
	v = TestFunction(V)

	nos = Expression("0")

	bcs = DirichletBC(V,nos,boundary)



	a = inner(nabla_grad(u), nabla_grad(v))*dx
	f = Constant(-1.0)
	L = f*v*dx
	u = Function(V)
	solve(a==L, u, bcs)

	Q = assemble(u*dx)
	#print Q
	#plot(u)
	#nteractive()
	#u_e = project(u_exact, V, bcs=bcs)
	#u_error = errornorm (u_exact, u, degree_rise = 0)
	print "mesh size: ", mesh.hmin()
	#print "u -  Error: ", u_error 
	print "Q - Error: ", abs(sqrt(Q_exact**2) - sqrt(Q**2))
	#plot(u_e, title="Exact velocity for concentric annulus")
	#interactive()

