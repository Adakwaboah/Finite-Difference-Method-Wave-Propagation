# Written by Akwasi Darkwah Akwaboah
# Description: Finite Difference Method (FDM) to simulate tissue potential
# Date: 11/26/2019

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors

# Step 1: Initialize the grid numbers for the solution regions
delta_t = 0.01
#delta_t = 0.01 * 1.50

delta_h = 0.2
#delta_h = 0.2 * 1.5

time_limit = 27
t = np.arange(0.0, time_limit * delta_t, delta_t)
x = np.arange(0.0, 3 + delta_h, delta_h)
y = np.arange(0.0, 3 + delta_h, delta_h)

R = 1
Cm = 1

r = ((R / Cm) * (delta_t)) / (delta_h ** 2)
print(r)

#mesh grid node initialization
V = np.ones((len(x), len(y), time_limit), dtype='float64') * -80  #initial resting potential
Anode = 100
Cathode = -100

# Initialize boundary conditions for each time iteration
# Electrode Stimulation
V[-2:, -2:,:] = Anode  # top right corner
V[:2, -2:, :] = Cathode  # bottom left corner

# #impulse
# V[-2:, -2:, 0] = Anode  # top right corner
# V[:2, -2:, 0] = Cathode  # bottom left corner

V[:, 0, :] = V[:, 1, :]  # left - no flux
V[len(x) - 1, :, :] = V[len(x) - 2, :, :]  # top - no flux
V[0, :, :] = V[1, :, :]  # bottom - no flux
V[:, len(x) - 1, :] = (4 / 3) * V[:, len(x) - 2, :]  # right flux

# Step 2: approximate derivative | finite difference - written by Akwasi
for k in np.arange(0, len(t) - 1):  # time
    for j in np.arange(1, len(y) - 1):  # do not include boundaries
        for i in np.arange(1, len(x) - 1):
            V[i, j, k + 1] = r * (V[j, i + 1, k] + V[j, i - 1, k] + V[j + 1, i, k] + V[j - 1, i, k]) - (
                        (4 * r / R) + delta_t - 1) * V[j, i, k]
    # ensure boundary conditions are maintained through the time course
    V[:, len(x) - 1, :] = (4 / 3) * V[:, len(x) - 2, :]  # right
    V[1, :, :] = V[0, :, :]  # bottom
    V[:, 1, :] = V[:, 0, :]  # left
    V[len(x) - 1, :, :] = V[len(x) - 2, :, :]  # top

    V[-2:, -2:, :] = Anode  # top right corner
    V[:2, -2:, :] = Cathode  # bottom left

    # # impulse propagation
    # V[-2:, -2:, 0] = Anode  # top right corner
    # V[:2, -2:, 0] = Cathode  # bottom left corner

# V = np.flip(V, axis = 0)
for time in np.arange(0, 1, 5):
    fig = plt.figure(figsize=(15, 15))
    fig.suptitle('$\Delta t = 0.1ms, \Delta h = \Delta x = \Delta y = 0.2cm$')
    ax = fig.add_subplot(321)
    im = plt.imshow(V[:, :, time], origin='lower')
    plt.title('$t = 0ms$')
    im.axes.get_xaxis().set_visible(False)
    im.axes.get_yaxis().set_visible(False)
    fig.colorbar(im, label='Potential(mV)', fraction=0.046, pad=0.04)
    #fig.tight_layout()

    ax = fig.add_subplot(322)
    im = plt.imshow(V[:, :, 2], origin='lower')
    plt.title('$t = 2\Delta t = 0.2ms$')
    im.axes.get_xaxis().set_visible(False) #Akwasi
    im.axes.get_yaxis().set_visible(False)
    fig.colorbar(im, label='Potential(mV)', fraction=0.046, pad=0.04)

    ax = fig.add_subplot(323)
    im = plt.imshow(V[:, :, time + 5], origin='lower')
    plt.title('$t = 5\Delta t = 0.5ms, $')
    im.axes.get_xaxis().set_visible(False)
    im.axes.get_yaxis().set_visible(False)
    fig.colorbar(im, label='Potential(mV)', fraction=0.046, pad=0.04)
    #fig.tight_layout() #Akwasi

    ax = fig.add_subplot(324)
    im = plt.imshow(V[:, :, time + 10], origin='lower')
    plt.title('$t = 10\Delta t = 1ms$')
    im.axes.get_xaxis().set_visible(False)
    im.axes.get_yaxis().set_visible(False)
    fig.colorbar(im, label='Potential(mV)', fraction=0.046, pad=0.04)
    #fig.tight_layout()

    ax = fig.add_subplot(325)
    im = plt.imshow(V[:, :, time + 18], origin='lower')
    plt.title('$t = 18\Delta t = 1.8ms$')
    im.axes.get_xaxis().set_visible(False)
    im.axes.get_yaxis().set_visible(False)
    fig.colorbar(im, label='Potential(mV)', fraction=0.046, pad=0.04)
    #fig.tight_layout()

    ax = fig.add_subplot(326)
    im = plt.imshow(V[:, :, time + 25], origin='lower')
    plt.title('$t = 25\Delta t = 2.5ms$')
    im.axes.get_xaxis().set_visible(False)
    im.axes.get_yaxis().set_visible(False)
    fig.colorbar(im, label='Potential(mV)', fraction=0.046)
    #fig.tight_layout() #Akwasi

plt.show()