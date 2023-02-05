#Plotting a hydrofoil's profile from its geometry points
#Author: Danick Lamoureux
#LM2 project under Frédérick Gosselin's supervision
#Date: 2022-10-14

import matplotlib.pyplot as plt
import numpy as np
from scipy.interpolate import CubicSpline

profile_filename = r'TestFiles\Francis99\F0.dat'

profile = open(profile_filename, 'r')
raw_coords = profile.readlines()
npoints = len(raw_coords)
if npoints < 3:
    print("Not enough points, please try again")
profile_coords = np.zeros([npoints, 2])
for i in range(npoints):
    line_coords = raw_coords[i].split(" ")
    line_coords = [item for item in line_coords if item != '']
    profile_coords[i, 0] = line_coords[0]
    profile_coords[i, 1] = line_coords[1]

profile_coords = np.array([profile_coords[:,0]/max(profile_coords[:,0]), profile_coords[:,1]]).T
theta = np.arctan2(profile_coords[:,1],profile_coords[:,0]-0.5)
for i in range(len(theta)):
    if theta[i] < 0 and i>len(theta)/4:
        theta[i] += 2*np.pi
print(profile_coords)
cs = CubicSpline(theta, profile_coords, bc_type='natural')
theta_s = np.linspace(min(theta),max(theta),1000)

plt.figure('Profile')
plt.plot(profile_coords[[-1,0],0], profile_coords[[-1,0],1], color='black')
plt.plot(profile_coords[:,0], profile_coords[:,1], color='black')
plt.plot([profile_coords[-1,0], profile_coords[0,0]], [profile_coords[-1,1], profile_coords[0,1]], c='k')
#plt.plot(cs(theta_s)[:,0], cs(theta_s)[:,1], color='black')
plt.xlim(min(profile_coords[:,0])-0.1,max(profile_coords[:,0])+0.1)
plt.ylim(min(profile_coords[:, 1])-0.1, max(profile_coords[:, 1])+0.1)
plt.gca().set_aspect('equal', adjustable='box')

plt.show()