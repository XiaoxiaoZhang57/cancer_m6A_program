from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import numpy as np
 
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

xeq = [1,2,3,4,5,6,7,8,9,10]
yeq = [89,75,24,27,35,60,178,284,184,287]
xmore = [1,2,3,4,5,6,7,8,9,10]
ymore = [109227,74651,84964,95406,105247,117455,131335,154240,192634,212293]
xless = [1,2,3,4,5,6,7,8,9,10]
yless = [76435,49006,48913,51884,53508,57045,63181,68981,88284,114844]
xlessnon = [1,2,3,4,5,6,7,8,9,10]
ylessnon = [1436094,773752,754488,807731,905858,996065,1076966,1203042,1629947,1828276]
xmorenon = [1,2,3,4,5,6,7,8,9,10]
ymorenon = [2109860,1255949,1360701,1503965,1721912,1954157,2179628,2505107,3245393,3198744]
#ax.bar(xeq, yeq, zs=0, zdir='y',alpha=0.8)
ax.bar(xless, yless, zs=1, zdir='y',alpha=0.8)
ax.bar(xmore, ymore, zs=2, zdir='y',alpha=0.8)
#ax.bar(xlessnon, ylessnon, zs=1.3, zdir='y',alpha=0.9)
#ax.bar(xmorenon, ymorenon, zs=2.3, zdir='y',alpha=0.9)
print(xs,ys) 
ax.set_xlabel('site')
ax.set_ylabel('ratio vs 1')
ax.set_zlabel('number')
 
#plt.savefig('./cancer+more+less.pdf')
plt.show()
