import math
from random import random
import numpy as np
from statistics import mode
import matplotlib.pyplot as plt
#import matplotlib.mlab as mlab
#import matplotlib.gridspec as gridspec

class Neuron:
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



w_out = 0.00225             # TODO: Evaluate if it's correct

def synapse_dt(w, vin, vout, vinp, voutp):
    c_i = 0.0083

    I_out_s = w * (vin - vout) + ((vin - vout) - (vinp - voutp)) * c_i

    return I_out_s


def count_snn_response(n, input, sim_time, start_time, w_in, w_out):
    network_size = len(n)
    dt = 1      # ms
    T = math.ceil(sim_time/dt)
    time_range = range(0, T-2)

    dt1_list = []
    dt2_list = []
    n_out = Neuron()

    for k in range(0, len(input)):
        v = []
        u = []
        Iapp = []

        v_out = np.zeros(T)
        u_out = np.zeros(T)
        I_out = np.zeros(T)

        for i in range(0, network_size):
            init_v = np.zeros(T)
            init_u = np.zeros(T)
            init_Iapp = np.zeros(T)
            init_v[0] = -70          # Resting potential
            init_u[0] = -14          # Steady state

            if i == 0:
                v_out[i] = -70       # Resting potential of the output
                u_out[i] = -14       # Steady state

            v.append(init_v)
            u.append(init_u)
            Iapp.append(init_Iapp)

        dt1 = 0
        dt2 = 0
        switch = 0          # 0 - no spikes, 1 - after 1st spike, 2 - after 2nd spike

        for i in time_range:
            # Setting the input current values
            if i < 3:
                for j in range(0, network_size):
                    Iapp[j][i] = 10
            elif i < start_time:
                for j in range(0, network_size):
                    Iapp[j][i] = 0
            else:
                for j in range(0, network_size):
                    Iapp[j][i] = input[k][j] * w_in[j]

            # counting the neuron response
            for j in range(0, network_size):
                if v[j][i] < 35:
                    # Update ODE
                    dv = (0.04 * v[j][i] + 5) * v[j][i] + 140 - u[j][i]
                    v[j][i + 1] = v[j][i] + (dv + Iapp[j][i]) * dt
                    du = n[j].a * (n[j].b * v[j][i] - u[j][i])
                    u[j][i + 1] = u[j][i] + dt * du
                else:
                    # Spike
                    v[j][i] = 35
                    v[j][i + 1] = n[j].c
                    u[j][i + 1] = u[j][i] + n[j].d

                if i > 0:
                    I_out[i] += synapse_dt(w_out[j], v[j][i], v_out[i], v[j][i - 1], v_out[i - 1])
                else:
                    I_out[i] += synapse_dt(w_out[j], v[j][i], v_out[i], 0, 0)

            # Count output
            if v_out[i] < 35:
                # Update ODE
                dv = (0.04 * v_out[i] + 5) * v_out[i] + 140 - u_out[i]
                v_out[i + 1] = v_out[i] + (dv + I_out[i]) * dt
                du = n_out.a * (n_out.b * v_out[i] - u_out[i])
                u_out[i + 1] = u_out[i] + dt * du
            else:
                # Spike
                v_out[i] = 35
                v_out[i + 1] = n_out.c
                u_out[i + 1] = u_out[i] + n_out.d

                if switch == 0:
                    dt1 = i - start_time
                    dt1_list.append(dt1)
                    switch += 1
                elif switch == 1:
                    dt2 = i - (dt1 + start_time)
                    dt2_list.append(dt2)
                    switch += 1

    plt.plot(v_out)
    plt.show()

    return dt1_list, dt2_list      # Get the output (times between 2 first spikes)



def network_mapping_routing(input, wout, Cref):
    patterns = len(input)
    sampling = len(input[0])

    L = [0 for i in range(sampling)]
    avg = [0 for i in range(sampling)]
    w1 = [0 for i in range(sampling)]
    w2 = [0 for i in range(sampling)]
    n = [Neuron() for i in range(sampling)]
    C = [0 for i in range(sampling)]
    dt1, dt2 = 0, 0

    for i in range(0, sampling):
        for j in range(0, patterns):
            L[i] += input[j][i]
        
        avg[i] = L[i]/patterns
        #  w1(i) = f(avg(i))

        w1[i] = 2*(25 - avg[i])
        #  w1[i] = 1.4*(26 - avg[i])

        if w1[i] < 0:
            n[i] = Neuron(isIIS=True)
        # else:
        #     n[i] = Neuron()

        w2[i] = wout
        C[i] = Cref
    
    #  [dt1, dt2] = mode(SNN(IN)) - moda (najczęściej występująca wartość) z 2 pierwszych czasów impulsów z odpowiedzi pozytywnej (po dobraniu wzorca)
    dt1, dt2 = count_snn_response(n, input, 200, 150, w1, w2)
    dt1_mode = mode(dt1)
    dt2_mode = mode(dt2)

    return w1, w2, C, n, dt1, dt2, dt1_mode, dt2_mode

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
