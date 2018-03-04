import numpy as np
import matplotlib.pyplot as plt
nus = np.array([1.233* i for i in range(21)])

sides = [
    0, 0.303, 0.978, 1.836, 2.78,
    3.78, 4.80, 5.80, 6.78, 7.74
]

g_nu = np.array(sides + [8.62] + sides[::-1])
print(nus)
print(g_nu)

plt.plot(nus, g_nu)
