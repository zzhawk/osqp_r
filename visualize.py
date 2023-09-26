import matplotlib.pyplot as plt
import numpy as np

spinlines = open("ans.csv", "r")

s = []
ref = []
upper = []
lower = []
plan = []


for ln in spinlines:
    spl = ln.split(',')
    s.append(float(spl[0]))
    ref.append(float(spl[1]))
    upper.append(float(spl[2]))
    lower.append(float(spl[3]))
    plan.append(float(spl[4]))

plt.plot(np.asarray(s), np.asarray(ref), "green")
plt.plot(np.asarray(s), np.asarray(upper), "red")
plt.plot(np.asarray(s), np.asarray(lower), "red")
plt.plot(np.asarray(s), np.asarray(plan), "blue")

plt.legend()

plt.show()

spinlines.close()
