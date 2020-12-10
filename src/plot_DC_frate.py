import numpy as np
import matplotlib.pyplot as plt

t, i, V = np.loadtxt('spike_recorder.txt', delimiter=',', unpack=True)
t01, i01, V01 = np.loadtxt('spike_recorder_01.txt', delimiter=',', unpack=True)

I = np.arange(300,9010,100)
F = np.zeros(np.shape(I))
F01 = np.zeros(np.shape(I))

for idx,i_ in enumerate(I):
    F[idx] = np.sum(i==(i_-300+100)/100)/10
    F01[idx] = np.sum(i01==(i_-300+100)/100)/10

plt.ion()

plt.figure(1)
plt.plot(t01*0.1,i01*100+200,' .',t,i*100+200,' .')
plt.xlabel('time (ms)')
plt.ylabel('DC current (pA)')
plt.title('DC current: 300 ~ 9000 (pA), 100 (pA) steps')
plt.legend({'dt=0.1 (ms)','dt=1.0 (ms)'})
plt.slim([200,250])
plt.savefig('DC_raster_200_250.png')
#plt.show()

plt.figure(2)
plt.plot(I,F01,' .',I,F,' .')
plt.xlabel('I (pA)')
plt.ylabel('F (H\)')
plt.title('DC current: 300 ~ 9000 (pA), 100 (pA) steps')
plt.legend({'dt=0.1 (ms)','dt=1.0 (ms)'})
plt.savefig('DC_FI.png')
plt.show()

