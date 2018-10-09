import numpy as np
import matplotlib.pyplot as plt

'''define variables'''
H = 0 #external magnetic field
n = 8
epoch = 2 ** 12 
steps = 2 ** 10
T = np.linspace(1.2,3.8,2 ** 6)
T = list(T)
Energy = [] 
Magnetization = []
SpecificHeat = []
Susceptibility = []
    
def deltaE(s, H, J, Sn):
    '''Use PCB(Peridical Boundary Condition) to compute the change of energe when flip one random spin'''
    dE = 2 * s * (J * Sn + H)
    return dE

def init_lattice(n, H = 0, J = 1):

    '''Create a nxn lattice with random spin configuration'''
    lattice = np.random.choice([1, -1], size=(n, n))
    return lattice

def Ising(lattice, n, steps, T, H = 0, J = 1):
    '''Ising model'''
    spins = 0
    E_prime = 0
    for k in range(steps):
        i, j = np.random.randint(n), np.random.randint(n)
        Sn = lattice[(i + 1) % n, j] + lattice[(i - 1) % n, j] + lattice[i, (j + 1) % n] + lattice[i, (j - 1) % n]
        dE = deltaE(lattice[i][j], H, J, Sn)
        if dE < 0 or np.random.random() < np.exp(-dE / T):
            lattice[i, j] *= -1
    for i in range(n):
        for j in range(n):
            spins += lattice[i][j]
            Sn = lattice[(i + 1) % n, j] + lattice[(i - 1) % n, j] + lattice[i, (j + 1) % n] + lattice[i, (j - 1) % n]
            E_prime += -1 * lattice[i][j] * Sn 
    E_prime /= 4
    return spins, E_prime
        
c1 = 1.0 / (epoch * n * n)
c2 = 1.0 / (epoch * epoch * n * n)

''' compute variables '''
for t in T:
    c = 0
    x = 0
    M = 0
    E = 0
    init_config = init_lattice(n)
    for i in range(epoch):
        spins, E_prime = Ising(init_config, n, steps, t)
        M += spins
        E += E_prime 
        c += E_prime ** 2
        x += spins ** 2
    C = (c1 * c - c2 * E * E) / (t * t)
    X = (c1 * x - c2 * M * M) / t
    Energy.append(E * c1)
    Magnetization.append(M * c1)
    SpecificHeat.append(C)
    Susceptibility.append(X / c2)

''' plot the calculated values '''
Magnetization = np.array(Magnetization)

f = plt.figure(figsize=(18, 10))  
 
sp =  f.add_subplot(2, 2, 1 )
plt.plot(T, Energy, 'o', color="red")
plt.xlabel("Temperature (T)", fontsize=20)
plt.ylabel("Energy ", fontsize=20)
 
sp =  f.add_subplot(2, 2, 2 )
plt.plot(T, abs(Magnetization), 'o', color="blue")
plt.xlabel("Temperature (T)", fontsize=20)
plt.ylabel("Magnetization ", fontsize=20)

sp =  f.add_subplot(2, 2, 3 )
plt.plot(T, SpecificHeat, 'o', color="red")
plt.xlabel("Temperature (T)", fontsize=20)
plt.ylabel("Specific Heat ", fontsize=20)
 
sp =  f.add_subplot(2, 2, 4 )
plt.plot(T, Susceptibility, 'o', color="blue")
plt.xlabel("Temperature (T)", fontsize=20)
plt.ylabel("Susceptibility", fontsize=20)
