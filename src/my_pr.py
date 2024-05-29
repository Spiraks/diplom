import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

pathExy = '/home/alex/my_diplom/mur_Exy.txt'
pathExz = '/home/alex/my_diplom/mur_Exz.txt'
pathEyz = '/home/alex/my_diplom/mur_Eyz.txt'
pathEyx = '/home/alex/my_diplom/mur_Eyx.txt'

MPA = pd.read_csv('/home/alex/my_diplom/MPA - DirTotal.csv')
dp1 = np.asarray(MPA[['dB(DirTotal) [] - Freq=\'2.442GHz\' Phi=\'0deg\'']])
dp2 = np.asarray(MPA[['dB(DirTotal) [] - Freq=\'2.442GHz\' Phi=\'90.0000000000002deg\'']])



a = 51
b = 124
r1 = 0
r2 = 77
r = range(r1,r2)
Exy = np.loadtxt(pathExy)
Exz = np.loadtxt(pathExz)
Eyz = np.loadtxt(pathEyz)
Eyx = np.loadtxt(pathEyx)
# антенна1 Exz + Exy
data_collected = Exz + Exy
# data_collected = Eyz + Eyx

stepSize = 1
print(len(Exz))

#Generated vertices
positions = []

# for r in range(15, 45):
deg_start = 270
deg_end = 360
theta = np.arange(deg_start, deg_end, stepSize)
array_dna = np.zeros(len(theta))

def point_pos(x0, y0, d, _theta):
    _theta_rad = np.radians(_theta)
    return x0 + d*np.cos(_theta_rad), y0 + d*np.sin(_theta_rad)



t = deg_start
i = 0
while t <  deg_end:
    positions.append(np.array(point_pos(a,b,r,t)).astype(np.int64).T)
    t += stepSize
    # array_dna.append(data_collected[positions[-1][0]][positions[-1][1]])
    summ = 0
    for vect in positions[-1]:
        summ += data_collected[vect[0]][vect[1]]
    array_dna[i] = (summ/(r2-r1))
    i+=1


 

# theta = np.arange(deg_start, deg_end, 0.01)
fig, ax = plt.subplots()
# ax.plot(theta, array_dna)
# ax.grid(True)
# ax.set_title("A line plot on a polar axis", va='bottom')
def f_norm(x_):
    return (x_ - min(x_)) / (max(x_) - min(x_))
ax.plot(range(len(dp1)),f_norm(dp1),label=r"Эталонные данные",linestyle=":",linewidth=2)
ax.plot(range(len(array_dna)),f_norm(10*np.log10(array_dna)),label=r"Результаты измерения",linewidth=2)
ax.grid()
ax.set_xticks([0,10,20,30,40,50,60,70,80,90])
ax.set_xticklabels(["0", "10", "20", "30", "40", "50", "60", "70", "80", "90"])
ax.set(xlabel="Градусы", ylabel="!!!!!!!!!!")
ax.legend()
plt.show()

  

# plt.imshow((np.abs(data_collected)))
# plt.show()