import math
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.mlab as mlab
import matplotlib.gridspec as gridspec

class Person:
    def __init__(self, isIIS = False):
        
        if isIIS == False:
            self.a = 0.02
            self.b = -0.1
            self.c = -65
            self.d = 6
        else:
            self.a = -0.02
            self.b = -1
            self.c = -60
            self.d = 8

def count_snn_response(n, input):
    network_size = len(n)
    v = []
    u = []
    sim_time = 200
    dt = 1

    T = math.ceil(sim_time/dt)
    ii = np.zeros(T)

    for i in range(network_size):
        init_v = np.zeros(T)
        init_u = np.zeros(T)
        init_v[0] = -70          # Resting potential
        init_u[0] = -14          # Steady state
        v.append(init_v)
        u.append(init_u)

    time_range = range(0, T-1)

    for i in time_range:
        
        



def network_mapping_routing(input, wout, Cref):
    sampling = len(input[0])
    patterns = len(input)

    L = [0 for i in range(sampling)]
    avg = [0 for i in range(sampling)]
    w1 = [0 for i in range(sampling)]
    w2 = [0 for i in range(sampling)]
    n = [Person() for i in range(sampling)]
    C = [0 for i in range(sampling)]
    dt1, dt2 = 0, 0

    for i in range(sampling):
        for j in range(patterns):
            L[i] += input[i][j]
        
        avg[i] = L[i]/patterns
        #  w1(i) = f(avg(i))

        w1[i] = 2*(25 - avg[i])
        #  w1[i] = 1.4*(26 - avg[i])

        if w1[i] < 0:
            n[i] = Person(isIIS=True)
        # else:
        #     n[i] = Person()

        w2[i] = wout
        C[i] = Cref
    
    #  [dt1, dt2] = mode(SNN(IN)) - moda (najczęściej występująca wartość) z 2 pierwszych czasów impulsów z odpowiedzi pozytywnej (po dobraniu wzorca)

    return w1, w2, C, n, dt1, dt2



# Initialize parameters
dt = 0.5
d = 8
a = 0.02 
c = -65
b = 0.2
T = math.ceil(1000/dt)
ii = np.zeros(T)

# Reserve memory
v1 = np.zeros(T)
u1 = np.zeros(T)
v1[0] = -70          # Resting potential
u1[0] = -14          # Steady state

v2 = np.zeros(T)
u2 = np.zeros(T)
v2[0] = -70
u2[0] = -14

v3 = np.zeros(T)
u3 = np.zeros(T)
v3[0] = -70
u3[0] = -14

R1 = 7
R2 = 10

time_range = range(0, T-1)

# Loop over time 1:T-1
for t in time_range:

    # Get input
    if t*dt > 500:
        Iapp = 15
    else:
        Iapp = 10
   
    if (t*dt > 250 and t*dt < 500) or (t*dt > 750 and t*dt < 1000):
        Iapp2 = 15
    else:
        Iapp2 = 10
    
    # Neuron 1
    if v1[t] < 35:
        # Update ODE
        dv = (0.04*v1[t] + 5)*v1[t] + 140 - u1[t]
        v1[t + 1] = v1[t] + (dv + Iapp)*dt
        du = a*(b*v1[t] - u1[t])
        u1[t + 1] = u1[t] + dt*du
    else:
        # Spike
        v1[t] = 35
        v1[t + 1] = c
        u1[t + 1] = u1[t] + d
    
    # Neuron 2
    if v2[t] < 35:
        # Update ODE
        dv2 = (0.04*v2[t] + 5)*v2[t] + 140 - u2[t]
        
        v2[t + 1] = v2[t] + (dv2 + Iapp2)*dt
        du2 = a*(b*v2[t] - u2[t])
        u2[t + 1] = u2[t] + dt*du2
    else:
        # Spike
        v2[t] = 35
        v2[t + 1] = c
        u2[t + 1] = u2[t] + d
    
    I3 = (v1[t]-v3[t])/R1 + (v2[t]-v3[t])/R2
    ii[t] = I3
     
    # Neuron 3
    if v3[t] < 35:
        # Update ODE
        dv3 = (0.04*v3[t] + 5)*v3[t] + 140 - u3[t]
        
        v3[t + 1] = v3[t] + (dv3 + I3)*dt
        du3 = a*(b*v3[t] - u3[t])
        u3[t + 1] = u3[t] + dt*du3
    else:
        # Spike
        v3[t] = 35
        v3[t+1] = c
        u3[t + 1] = u3[t] + d


# Plot voltage trace

fig, (ax0, ax1, ax2, ax3) = plt.subplots(4, 1)

t_arranged = np.arange(0, 1000, dt)

ax0.plot(t_arranged, v1)
ax1.plot(t_arranged, v2)
ax2.plot(t_arranged, v3)
ax3.plot(ii)
plt.show()

# figure(1)
# plot((0:T-1)*dt,v,'g')
# figure(2)
# plot((0:T-1)*dt,v2,'b')

# xlabel('Time[ms]')
# ylabel('Membrane voltage[mV]')
# figure(3)
# plot(ii,'r')

# figure(4)
# plot((0:T-1)*dt,v3,'b')