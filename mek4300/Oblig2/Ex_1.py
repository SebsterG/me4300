from dolfin import *
import matplotlib.pyplot as plt
#parameters["reorder_dofs_serial"] = False
beta = [-0.1]#1.0,0.3,0.0,-0.1,-0.18,-0.198838]
e_1 = Expression(("0","0"))
e_2 =  Expression(("1","1")); e_3 =  Expression(("10","10"))

guess = [e_1,e_2]
for i in range(len(beta)):
    for j in range(len(guess)):
        print "beta: ", beta[i]
        L = 6
        dens = 100
        mesh = IntervalMesh(dens,0,L)

        V = FunctionSpace(mesh,"CG",1)
        VV = V*V

        def left(x,on_boundary):
        	return near(x[0],0) and on_boundary
        def right(x,on_boundary):
        	return near(x[0],L) and on_boundary

        bc1 = DirichletBC(VV.sub(1), 0, left)
        bc2 = DirichletBC(VV.sub(0), 0, left)
        bc3 = DirichletBC(VV.sub(1), 1, right)
        bcs = [bc1,bc2,bc3]

        #fh = TrialFunction(VV)
        #f, h = split(fh)
        vf,vh = TestFunction(VV)
        #f_h_ = Function(VV)
        #f_h_ = interpolate(Expression(("0","0")),VV)
        f_h_ = project(guess[j],VV)
        f_,h_ = split(f_h_)
        beta_1 = Constant(beta[i])

        F1 = -inner(grad(h_),grad(vh))*dx + f_*h_.dx(0)*vh*dx +beta_1*(1-h_*h_)*vh*dx
        F2 = h_*vf*dx - f_.dx(0)*vf*dx
        F = F1+F2
        solve(F==0,f_h_,bcs)
        f,h = split(f_h_)
        h_derivative = project(h.dx(0),V)
        h = project(h,V)
        plt.plot(mesh.coordinates(), h.vector().array()[::-1],label=("guess: %.f beta: %.1f"%(j+1,beta[i])))
        #plt.plot(mesh.coordinates(), h_derivative.vector().array()[::-1],label=("beta: %.4f"%beta[i]))
        #plt.plot(h_plot)
        axes = plt.gca()
        #axes.set_xlim([0,6])
        #axes.set_ylim([0,1.0])
        legend = axes.legend(loc='lower right', shadow=True)
        #plot(h_)
plt.grid()
plt.show()
#interactive()
