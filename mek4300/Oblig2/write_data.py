import numpy as np
import matplotlib.pyplot as plt


lift = np.loadtxt(open("left.txt"))
drag = np.loadtxt(open("drag.txt"))
time = np.loadtxt(open("time.txt"))
pressure = np.loadtxt(open("pressure.txt"))

"""print "lift_max: ",max(lift)
print "drag_max: ",max(drag)
print "pressure_max: ",max(pressure)"""
Um = 1.5
U_mean = 2.0*Um/3.0
new_lift = lift[:]*2/(U_mean**2*0.1)
new_drag = drag[:]*2/(U_mean**2*0.1)

print "lift:" ,np.max(new_lift)
print "drag:",np.max(new_drag)
print "pressure:",np.max(pressure)



plt.plot(time,new_lift, label =("lift"))
plt.plot(time,new_drag, label =("drag"))
plt.plot(time,pressure, label =("pressure"))
plt.axis([0, 10, -1.2, 5])
axes = plt.gca()
legend = axes.legend(loc='upper right', shadow=True)
plt.show()
