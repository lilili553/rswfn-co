# 此脚本用于使用四次多项式拟合α-(节点、权重)数据

import numpy as np
import matplotlib.pyplot as plt

alps=np.loadtxt("alps.txt")
# roots=np.loadtxt("RYS_ROOTS_10.txt",delimiter=",")
# print(alps.shape)
# print(roots.shape)
# Np=4 # 使用几阶多项式拟合
# for i in range(5):
#     xs=alps
#     ys=roots[i,:]
#     cs=np.polyfit(xs,ys,Np) # 拟合
#     vs=cs[0]*xs+cs[1]
#     for i in range(2,Np+1):
#         vs=vs*xs+cs[i]
#     delta=np.linalg.norm(ys-vs)
#     print(delta)
#     plt.scatter(xs,ys)
#     plt.plot(xs,vs)
# plt.show()

plt.plot(alps)
plt.show()