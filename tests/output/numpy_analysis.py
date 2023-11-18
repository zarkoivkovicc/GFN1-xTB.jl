import numpy as np

my = np.loadtxt("methanol_H_my")
ref = np.loadtxt("methanol_H_ref")
diff = abs(my - ref)
print(diff  > 1e-6)
print(diff)
