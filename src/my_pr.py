import matplotlib.pyplot as plt
import numpy as np

pathExy = '/home/alex/src/diplom_FDTD/results/mur_Exy.txt'
pathExz = '/home/alex/src/diplom_FDTD/results/mur_Exz.txt'


a = 79
b = 51
r1 = 75
r2 = 77
r = range(r1,r2)
Exy = np.loadtxt(pathExy)
Exz = np.loadtxt(pathExz)
data_collected = Exy+ Exz
stepSize = 0.01

#Generated vertices
positions = []
array_dna = []

# for r in range(15, 45):

t = 0
i = 0
while t <  np.pi*2:
    positions.append(np.array([r * np.cos(t) + a, r * np.sin(t) + b]).astype(np.int64).T)
    i+=1
    t += stepSize
    # array_dna.append(data_collected[positions[-1][0]][positions[-1][1]])
    summ = 0
    for vect in positions[-1]:
        summ += data_collected[vect[0]][vect[1]]
    array_dna.append(summ/(r2-r1))

# plt.plot(range(len(array_dna)),array_dna)  
# naming the x axis  

theta = np.arange(0, np.pi*2, 0.01)

fig, ax = plt.subplots(subplot_kw={'projection': 'polar'})


ax.plot(theta, array_dna)
# ax.set_rticks([0.5, 1, 1.5, 2])  # Less radial ticks
# ax.set_rlabel_position(-22.5)  # Move radial labels away from plotted line
ax.grid(True)

ax.set_title("A line plot on a polar axis", va='bottom')
plt.show()
# plt.xlabel('x - axis')  
# # naming the y axis  
# plt.ylabel('y - axis')  
    
# # giving a title to my graph  
# plt.title('DNA r = %d'%r)  
    
# function to show the plot  
  

# plt.imshow((np.abs(data_collected)))
# plt.show()