import numpy as np
import matplotlib.pyplot as plt

DC = np.loadtxt('DC_recorder.txt', delimiter=',', unpack=True)
t, i, V = np.loadtxt('spike_recorder.txt', delimiter=',', unpack=True)
t01, i01, V01 = np.loadtxt('spike_recorder_01.txt', delimiter=',', unpack=True)

F=[]
F01=[]

print("DC dt=0.1 dt=1")
for i_ in range(16):
    F.append(np.sum(i==i_+1)/10)
    F01.append(np.sum(i01==i_+1)/10)
    print(DC[i_], F01[-1], F[-1])



