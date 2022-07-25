import math
import numpy as np
import matplotlib.pyplot as plt

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

vesicle_fusion = [1, 1, 2, 2, 3, 2, 2, 1, 3, 4, 4, 3, 2, 2, 1, 2, 4, 4, 3, 4, 3, 2, 2, 1, 2, 2, 3, 3, 4, 5, 3, 4, 3, 4, 4, 3, 2, 2, 1, 2, 3, 5, 4, 4, 3, 4, 2, 3, 2, 3, 2, 2, 3, 3, 4, 4, 3, 2, 2, 3, 4, 5, 5, 4, 4, 3, 4, 6, 8, 15, 28, 30, 30, 31, 30, 32, 33, 31, 34, 35, 37, 33, 34, 35, 36, 37, 35, 34, 36, 37, 38, 38, 36, 35, 36, 38, 40, 42, 45, 42, 48, 50, 48, 46, 47, 48, 49, 52, 56, 58, 60, 62, 66, 75, 90, 115, 150, 195, 240, 275, 300, 305, 298, 290, 270, 250, 230, 200, 175, 160, 140, 125, 115, 105, 100, 90, 82, 75, 65, 68, 63, 60, 55, 50, 47, 45, 46, 45, 42, 40, 37, 34, 35, 32, 31, 28, 25, 26, 27, 25, 22, 20, 21, 20, 15, 17, 19, 16, 18, 17, 15, 13, 12, 14, 11, 10, 8, 6, 5, 3, 4]

#skalowanie
max_value = 80 #pA

skala = max_value/max(vesicle_fusion)
#vesicle_fusion = vesicle_fusion*skala
for i in range(0, len(vesicle_fusion)):
    vesicle_fusion[i] *= skala

#print(vesicle_fusion)
probkowanie = 20
k = 0
kk = 0
iter = len(vesicle_fusion)
original_data = []

while True:
    k = int(k+iter/probkowanie)
    if k < iter:
        original_data.append(vesicle_fusion[k] * 1.2)   # Scaling
    else:
        break
    
#print(len(original_data))

merged_input = []
merged_input.append(original_data)

for i in range(1, 40):
    L = []
    for j in range(0, len(original_data)):
        if j < 16 or (j > 25 and j < 29):
            L.append(original_data[j] + random() * 13)
        else:
            L.append(original_data[j] + random() * 10)

    merged_input.append(L)

w1, w2, C, n, dt1, dt2, dt1_mode, dt2_mode = network_mapping_routing(merged_input, w_out, 0)

print(dt1)
print(dt2)

wrong_input = []

for i in range(0, 40):
    N = []
    for j in range(0, len(original_data)):
        N.append((max(original_data) * random())/0.6)
    
    wrong_input.append(N)

dt1, dt2 = count_snn_response(n, wrong_input, 200, 150, w1, w2)

print(dt1)
print(dt2)