import numpy as np
import matplotlib.pyplot as plt 
import seaborn as sns
sns.set_style('white')
import warnings
warnings.filterwarnings("ignore")

# import data
fname = "kepler.dat" 
results = np.loadtxt(open(fname), delimiter = ' ')

# label data
x = np.asarray([var[0] for var in results])
y = np.asarray([var[1] for var in results]) 


# x vs y
fig = plt.figure()
plt.title("Kepler")
plt.xlabel("x")
plt.ylabel("y")
plt.plot(x, y, c = 'red', label='Kepler Data')
plt.legend()

plt.axis('scaled')
plt.gca().set_ylim(0.0, 0.0025)

plt.savefig("Kepler.pdf")
plt.show()