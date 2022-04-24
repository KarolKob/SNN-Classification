import math
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.mlab as mlab
import matplotlib.gridspec as gridspec

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


def count_snn_response(n, input, patterns, start):
    network_size = len(n)
    v = []
    u = []
    Iapp = []
    sim_time = patterns
    dt = 1

    T = math.ceil(sim_time/dt)
    ii = np.zeros(T)

    v_out = np.zeros(T)
    u_out = np.zeros(T)
    I_out = np.zeros(T)
    n_out = Neuron()

    for i in range(network_size):
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

    time_range = range(0, T-1)

    for i in time_range:
        # Setting the input current values
        if i < 3:
            for j in range(network_size):
                Iapp[j][i] = 10
        elif i < start:
            for j in range(network_size):
                Iapp[j][i] = 0
        else:
            for j in range(network_size):
                Iapp[j][i] = input[j][i]

        # counting the neuron response
        for j in range(network_size):
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
                I_out[j] += synapse_dt(w_out, v[j][i], v_out[i], v[j][i - 1], v_out[i - 1])
            else:
                I_out[j] += synapse_dt(w_out, v[j][i], v_out[i], 0, 0)

        # TODO: Count output
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



def network_mapping_routing(input, wout, Cref):
    sampling = len(input[0])
    patterns = len(input)

    L = [0 for i in range(sampling)]
    avg = [0 for i in range(sampling)]
    w1 = [0 for i in range(sampling)]
    w2 = [0 for i in range(sampling)]
    n = [Neuron() for i in range(sampling)]
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
            n[i] = Neuron(isIIS=True)
        # else:
        #     n[i] = Neuron()

        w2[i] = wout
        C[i] = Cref
    
    #  [dt1, dt2] = mode(SNN(IN)) - moda (najczęściej występująca wartość) z 2 pierwszych czasów impulsów z odpowiedzi pozytywnej (po dobraniu wzorca)

    return w1, w2, C, n, dt1, dt2

