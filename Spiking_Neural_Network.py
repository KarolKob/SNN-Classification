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

    #plt.plot(v_out)
    #plt.show()

    return dt1_list, dt2_list      # Get the output (times between 2 first spikes)

def trained_networks_responses(n_A, w_in_A, w_out_A, n_N, w_in_N, w_out_N, n_X, w_in_X, w_out_X, input, sim_time, start_time):
    dt1_A, dt2_A = count_snn_response(n_A, input, sim_time, start_time, w_in_A, w_out_A)
    dt1_N, dt2_N = count_snn_response(n_N, input, sim_time, start_time, w_in_N, w_out_N)
    dt1_X, dt2_X = count_snn_response(n_X, input, sim_time, start_time, w_in_X, w_out_X)
    return [[dt1_A, dt2_A], [dt1_N, dt2_N], [dt1_X, dt2_X]]

def network_mapping_routing(input, wout, Cref, d, k):
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

        w1[i] = d*(k - avg[i])
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

# Read all the input data
with open("Letter_A.txt", "r+") as letter_A:
    A_lines = letter_A.readlines()

vals_A = []
for i in range(0, len(A_lines)):
    vals_A.append(float(A_lines[i]))

with open("Letter_N.txt", "r+") as letter_N:
    N_lines = letter_N.readlines()

vals_N = []
for i in range(0, len(N_lines)):
    vals_N.append(float(N_lines[i]))

with open("Letter_X.txt", "r+") as letter_X:
    X_lines = letter_X.readlines()

vals_X = []
for i in range(0, len(X_lines)):
    vals_X.append(float(X_lines[i]))

with open("Square.txt", "r+") as square:
    square_lines = square.readlines()

vals_square = []
for i in range(0, len(square_lines)):
    vals_square.append(float(square_lines[i]))

merged_input_A = []
merged_input_A.append(vals_A)

for i in range(1, 40):
    L = []
    for j in range(0, len(vals_A)):
        if j < 16 or (j > 25 and j < 29):
            L.append(vals_A[j] + random() * 0.6)
        else:
            L.append(vals_A[j] + random() * 0.4)

    merged_input_A.append(L)

#skalowanie
max_value = 80 #pA

skala = max_value/max(map(max, merged_input_A))
for i in range(0, len(merged_input_A)):
    for j in range(0, len(merged_input_A[0])):
        merged_input_A[i][j] *= skala

merged_input_N = []
merged_input_N.append(vals_N)

for i in range(1, 40):
    L = []
    for j in range(0, len(vals_N)):
        if j < 16 or (j > 25 and j < 29):
            L.append(vals_N[j] + random() * 0.7)
        else:
            L.append(vals_N[j] + random() * 0.5)

    merged_input_N.append(L)

#skalowanie
max_value = 80 #pA

skala = max_value/max(map(max, merged_input_N))
for i in range(0, len(merged_input_N)):
    for j in range(0, len(merged_input_N[0])):
        merged_input_N[i][j] *= skala

merged_input_X = []
merged_input_X.append(vals_X)

for i in range(1, 40):
    L = []
    for j in range(0, len(vals_X)):
        if j < 16 or (j > 25 and j < 29):
            L.append(vals_X[j] + random() * 1.3)
        else:
            L.append(vals_X[j] + random() * 1.1)

    merged_input_X.append(L)

#skalowanie
max_value = 80 #pA

skala = max_value/max(map(max, merged_input_X))
for i in range(0, len(merged_input_X)):
    for j in range(0, len(merged_input_X[0])):
        merged_input_X[i][j] *= skala

merged_input_sq = []
merged_input_sq.append(vals_square)

for i in range(1, 40):
    L = []
    for j in range(0, len(vals_square)):
        if j < 16 or (j > 25 and j < 29):
            L.append(vals_square[j] + random() * 10)
        else:
            L.append(vals_square[j] + random() * 8)

    merged_input_sq.append(L)

#skalowanie
max_value = 80 #pA

skala = max_value/max(map(max, merged_input_sq))
for i in range(0, len(merged_input_sq)):
    for j in range(0, len(merged_input_sq[0])):
        merged_input_sq[i][j] *= skala

w1_A, w2_A, C_A, n_A, dt1_A, dt2_A, dt1_mode_A, dt2_mode_A = network_mapping_routing(merged_input_A, w_out, 0, 2.1, 63.5)
w1_N, w2_N, C_N, n_N, dt1_N, dt2_N, dt1_mode_N, dt2_mode_N = network_mapping_routing(merged_input_N, w_out, 0, 0.68, 65)
w1_X, w2_X, C_X, n_X, dt1_X, dt2_X, dt1_mode_X, dt2_mode_X = network_mapping_routing(merged_input_X, w_out, 0, 0.96, 60.9)
print("A")
print(dt1_A)
print(dt2_A)

print("N")
print(dt1_N)
print(dt2_N)

print("X")
print(dt1_X)
print(dt2_X)

dt1, dt2 = count_snn_response(n_A, merged_input_sq, 200, 150, w1_A, w2_A)

print("SNN_A -> Square")
print(dt1)
print(dt2)
print(mode(dt1))
print(mode(dt2))

dt1, dt2 = count_snn_response(n_N, merged_input_sq, 200, 150, w1_N, w2_N)

print("SNN_N -> Square")
print(dt1)
print(dt2)
print(mode(dt1))
print(mode(dt2))

dt1, dt2 = count_snn_response(n_X, merged_input_sq, 200, 150, w1_X, w2_X)

print("SNN_X -> Square")
print(dt1)
print(dt2)
print(mode(dt1))
print(mode(dt2))

dts = trained_networks_responses(n_A, w1_A, w2_A, n_N, w1_N, w2_N, n_X, w1_X, w2_X, merged_input_A, 200, 150)
print("A")
print(dts[0])
print(dts[1])
print(dts[2])

dts = trained_networks_responses(n_A, w1_A, w2_A, n_N, w1_N, w2_N, n_X, w1_X, w2_X, merged_input_N, 200, 150)
print("A")
print(dts[0])
print(dts[1])
print(dts[2])

dts = trained_networks_responses(n_A, w1_A, w2_A, n_N, w1_N, w2_N, n_X, w1_X, w2_X, merged_input_X, 200, 150)
print("A")
print(dts[0])
print(dts[1])
print(dts[2])