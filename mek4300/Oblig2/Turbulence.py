from dolfin import *
import numpy as np
import matplotlib.pyplot as plt
N = 1001
"""
mesh = IntervalMesh(N, 0, 1)
#plot(mesh);interactive()
x = mesh.coordinates()
x[:,0] = -np.arctan(pi*(x[:,0]))/np.arctan(pi)"""

mesh = UnitIntervalMesh(N)
V = FunctionSpace(mesh,"CG",1)
v = TestFunction(V)


class Nos(SubDomain):
    def inside(self, x, on_boundary):
        return near(x[0],0)

nos = Nos()
bound = FacetFunction("size_t",mesh)
bound.set_all(0)
nos.mark(bound,1)
#plot(bound,interactive=True)
bc1 = DirichletBC(V,0,nos)
bcs = [bc1]

k = 0.41
A = 26
B = 5.5
v_star = 0.05
nu = 0.05/1000.0

y_pluss = Expression("x[0]*v_star/nu",v_star=v_star,nu=nu)
l = Expression("k*x[0]*(1-exp(-y_pluss/A))",k=k,y_pluss = y_pluss, A = A)

nu = Constant(nu); v_star = Constant(v_star) ;k = Constant(k); A = Constant(A)
u_0 = Expression("0")
u = project(u_0,V)
#plot(u);interactive()
nu_t = l*l*abs(u.dx(0))
dpdx = Expression("0.05*0.05")




a = -nu*u.dx(0)*v.dx(0)*dx + dpdx*v*dx - nu_t*u.dx(0)*v.dx(0)*dx


solve(a==0,u,bcs)

print "nut", project(nu_t,V).vector().array()
plt.plot(project(nu_t,V).vector().array())
plt.show()

k = 0.41
A = 26
B = 5.16
v_star = 0.05
nu = 0.05/1000.0

y1 = project(Expression("x[0]"), V).vector().array()
value1 = (np.abs(y1-0.005)).argmin()

u_arr = u.vector().array()[value1:][::-1]
y_arr = y1[value1:][::-1]
#print y_arr
#print u_arr
plt.plot(y_arr,u_arr/v_star,label="Numerical")
plt.plot(y_arr*1000, label = "Analytical")
axes = plt.gca()
legend = axes.legend(loc='lower right', shadow=True)
plt.show()



value2 = (np.abs(y1-0.030)).argmin()
u_arr1 = u.vector().array()[:value2][::-1]
y_arr1 = y1[:value2][::-1]
u_exact = (1.0/k)*np.log(y_arr1*1000) + B
#print u_arr1
#print y_arr1
plt.plot(y_arr1,u_exact, label="Analytical")
plt.plot(y_arr1,u_arr1/v_star,label = "Numerical")
axes = plt.gca()
legend = axes.legend(loc='lower right', shadow=True)
plt.show()
