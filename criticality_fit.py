import numpy as np
import matplotlib.pyplot as plt

Ts = np.linspace(2.2625, 2.278, 4)
Ts = list(Ts)
c1 = []
c2 = []
L = [8, 16, 32, 64]
L_reci = [1.0 / 16, 1.0 / 32 , 1.0 / 64]
cross_p = []
Tc = 2.269185

U2_8 = [1.067959, 1.069669, 1.071369, 1.072887]
U2_16 = [1.067631, 1.070811, 1.073274, 1.076795]
U2_32 = [1.061331, 1.068016, 1.076723, 1.084832]
U2_64 = [1.054256, 1.068769, 1.084052, 1.106944]
U2 = [U2_8, U2_16, U2_32, U2_64]

def lin_reg(x, y):
    '''using linear regression to get the coefficient of the linear function by least square method'''
    n = len(x)
    sum_x, sum_y, sum_x_y, sum_x_sq = 0.0, 0.0, 0.0, 0.0
    for i in range(n):
        sum_x += x[i]
        sum_y += y[i]
        sum_x_sq += x[i] ** 2
        sum_x_y += x[i] * y[i]
    a = (n * sum_x_y - sum_x * sum_y) / (n * sum_x_sq - sum_x ** 2)
    b = (sum_x_sq * sum_y - sum_x * sum_x_y) / (n * sum_x_sq - sum_x ** 2)
    return a, b

def get_cross():
    '''get the linear function and cross point between L and L / 2'''
    for i in range(len(L)):
        a, b = lin_reg(Ts, U2[i])
        c1.append(a)
        c2.append(b)
        if b >= 0:
            print("L = %d, y = %.7f x + %.7f" %(L[i], a, b))
        else:
            print("L = %d, y = %.7f x - %.7f" %(L[i], a, abs(b)))
    for i in range(len(c1) - 1):
        cross_p.append((c2[i + 1] -c2[i]) / (c1[i] - c1[i + 1]))
        
get_cross()

def get_crit():
    '''get the critical point of Ising model by linear regression'''
    a, b = lin_reg(L_reci, cross_p)
    print("critical temperature is %.6f" % b)
    return a, b

a, b = get_crit()  

'''plot the Binder Cumulant dependence of temperature for different sizes of the system'''
T = np.linspace(2.2625, 2.278)
for i in range(len(L)):
    Bin_cum = c1[i] * T + c2[i]
    plt.plot(T, Bin_cum, label = 'L = %d'%L[i])
    plt.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=3,
               ncol=3, mode='expand', borderaxespad=0.)
    plt.scatter(Ts, U2[i])
    plt.xlabel("Temperature T / J", fontsize=20)
    plt.ylabel("Binder Cumulant U2 ", fontsize=20)
plt.show()

'''plot finite size scaling of the Binder cumulant crossing points'''
T_diff = []
L_reci_space = np.linspace(1.0 / 16.0, 1.0 / 256)
for i in range(len(L) - 1):
    T_diff.append(cross_p[i] - Tc)
T_diff_space = a * L_reci_space + b - Tc
plt.scatter(L_reci, T_diff)
plt.plot(L_reci_space, T_diff_space)
plt.xlabel("1 / L", fontsize=20)
plt.ylabel("T* - Tc ", fontsize=20)
plt.axis([0.0, 0.07, -0.006, 0.001])
plt.show()
