import numpy as np

nRun=5

sum_x, sum_y, sum_z = None, None, None

for i in range(1,nRun):
    data = np.genfromtxt('../T300/runs/run0/data/'+str(i)+'/G_dprime_T300.data')
    x, y = data[:, 0], data[:, 1]

    if sum_x is None:
        sum_x, sum_y  = np.zeros_like(x), np.zeros_like(y)

    sum_x += x
    sum_y += y

# Compute averages
avg_x = sum_x / (nRun-1)
avg_y = sum_y / (nRun-1)

# Save to a new file
np.savetxt("G_dp_ave.dat", np.column_stack([avg_x, avg_y]),
           header="Average X  Average Y", fmt="%.6f")

print("Averaged coordinates saved to 'G_dp_ave.dat'")
