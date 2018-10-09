import math
import random
import numpy as np
import matplotlib.pyplot as plt

size = 8
steps = 80
epoch = 2 ** 14

Ts = np.linspace(2.2625,2.278,4)
Ts = list(Ts)
Energy = [] 
Magnetization = []
SpecificHeat = []
Susceptibility = []
Binder_cumulant = []

class Config:
    '''initiate the configuration, give the lattice site and return the cluster'''
    def __init__(self, size):
        self.size = size
        self.init_config()
        self.Kb = 1

    def init_config(self):
        self.config = np.random.choice([1, -1], size=(self.size, self.size))

    def get_neighbours(self, i, j):
        return [((i + 1) % self.size, j), (i, (j - 1) % self.size), ((i - 1) % self.size, j), (i, (j + 1) % self.size)]
    
    def get_spin(self, i, j):
        return self.config[i][j]

    def flip_spin(self, i, j):
        self.config[i][j] *= -1

    def get_cluster(self, i, j, T):
        pocket = []
        cluster = []
        rej_prob = 1 - math.exp(-2.0 / (self.Kb * T))
        spin = self.get_spin(i, j)
        cluster.append((i, j))
        pocket.append((i, j))
        while pocket != []:
            site = random.choice(pocket)
            for k in self.get_neighbours(site[0], site[1]):
                if self.get_spin(k[0], k[1]) == spin and k not in cluster and random.uniform(0.0, 1.0) < rej_prob:
                    cluster.append((k[0], k[1]))
                    pocket.append((k[0], k[1]))
            pocket.remove(site)
        return cluster

class Ising_model:
    '''Ising model with Wolff cluster updates algorithm, and compute the variables'''
    def __init__(self, size, T):
        self.config = Config(size)
        self.size = size
        self.T = T

    def flip_cluster(self, cluster):
        for site in cluster:
            self.config.flip_spin(site[0], site[1])

    def random_pick(self):
        return np.random.randint(self.size), np.random.randint(self.size)

    def Wolff(self):
        site = self.random_pick()
        cluster = self.config.get_cluster(site[0], site[1], self.T)
        self.flip_cluster(cluster)
            
    def calc_Mag(self):
        M = self.config.config.sum()
        m_abs = abs(M) / (self.size ** 2)
        m_sq = (M / (self.size ** 2)) ** 2
        return m_abs, m_sq
    
    def calc_U2(self, m_abs_exp, m_sq_exp):
        U2 = m_sq_exp / m_abs_exp ** 2
        return U2
        
    def calc_Susceptibility(self, m_abs_exp, m_sq_exp):
        X = (1.0 / self.T) * self.size ** 2 * (m_sq_exp - m_abs_exp ** 2)
        return X
    
    def simulate(self, steps, epoch):
        m_abs_sum = 0.0
        m_sq_sum = 0.0
        for i in range(steps):
            self.Wolff()
        for j in range(epoch):
            self.Wolff()
            m_abs, m_sq = self.calc_Mag()
            m_abs_sum += m_abs
            m_sq_sum += m_sq
        m_abs_exp = m_abs_sum / epoch
        m_sq_exp = m_sq_sum / epoch
        U2 = self.calc_U2(m_abs_exp, m_sq_exp)
        X = self.calc_Susceptibility(m_abs_exp, m_sq_exp)
        return U2, X, m_abs_exp

'''implementing the Ising model'''        
for T in Ts:
    ising = Ising_model(size, T)
    U2, X, m = ising.simulate(steps, epoch)
    Binder_cumulant.append(U2)
    Susceptibility.append(X)
    Magnetization.append(m)
    print('T = %f, U2 = %f' % (T, U2))

''' plot the calculated values '''
f = plt.figure(figsize=(18, 10))

sp =  f.add_subplot(2, 2, 1 )
plt.plot(Ts, Magnetization, 'o', color="red")
plt.xlabel("Temperature (T)", fontsize=20)
plt.ylabel("Magnetization ", fontsize=20)
 
sp =  f.add_subplot(2, 2, 2 )
plt.plot(Ts, Binder_cumulant, 'o', color="blue")
plt.xlabel("Temperature (T)", fontsize=20)
plt.ylabel("Binder cumulant ", fontsize=20)

sp =  f.add_subplot(2, 2, 3 )
plt.plot(Ts, Susceptibility, 'o', color="red")
plt.xlabel("Temperature (T)", fontsize=20)
plt.ylabel("Susceptibility ", fontsize=20)

    