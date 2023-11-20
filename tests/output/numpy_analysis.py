import numpy as np

my = np.loadtxt("co/co_Hp_my")
ref = np.loadtxt("co/co_Hp_ref")
diff = abs(my - ref)
print(diff  > 1e-6)
print(diff)
