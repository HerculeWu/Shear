import numpy as np
from scipy.spatial import cKDTree
import sys
rt = 30.0
lensdens = 0.1
sourcedens = 5.0
radius = 1*60.0
args = sys.argv
if '-h' in args or len(args) != 4:
    print("Usage:")
    print("python3 generator.py [radius of the survey (deg)] [truncation radius (arcmin)] [output file name]")
    print("generate random source and halo positon")
    print("calculate relative position in polar coordinate")
    print("select only halos with in the truncation radius")
    exit(0)
radius = float(args[1])*60.0
rt = float(args[2])
# generate random points
# units in arcmin
Nsource = int(np.pi*radius*sourcedens)
Nhalo = int(np.pi*radius*lensdens)
r = np.sqrt(np.random.rand(Nsource))*radius + rt
phi = np.random.rand(Nsource)*2.*np.pi

spos = np.zeros((Nsource, 2))
spos[:, 0] = r*np.cos(phi)
spos[:, 1] = r*np.sin(phi)

r = np.sqrt(np.random.rand(Nhalo))*radius
phi = np.random.rand(Nhalo)*2.*np.pi

hpos = np.zeros((Nhalo, 2))
hpos[:, 0] = r*np.cos(phi)
hpos[:, 1] = r*np.sin(phi)

alphas = np.random.rand(Nhalo)*2.*np.pi
sid = np.arange(0, Nsource, 1).astype(int)

# print(Nsource, Nhalo)
# get neigboure list with in truncation radius
sourceTree = cKDTree(hpos)
res = sourceTree.query_ball_point(spos, rt)
Nhalos = np.array([len(n) for n in res])
Nhalomax = np.max(Nhalos)
# get relative position in polar coordinate (theta, phi)
# write to file
out = open(args[3], 'w+')
ns = 0
for i, neig in enumerate(res):
    hs = np.array(neig)
    if hs.size > 0:
        ns += 1
out.write('{0:d} {1:d}\n'.format(ns, Nhalomax))
for i, neig in enumerate(res):
    hs = np.array(neig)
    if hs.size > 0:
        outstr = ''
        rposx = spos[i, 0] - hpos[neig, 0]
        rposy = spos[i, 1] - hpos[neig, 1]
        cpos = rposx+(1.j*rposy)
        theta = np.abs(cpos)
        phi = np.angle(cpos)
        outstr += '{0:d} {1:d}'.format(sid[i], Nhalos[i])
        for j in range(len(theta)):
            outstr += ' {0:.15g} {1:.15g} {2:.15g}'.format(theta[j], phi[j], alphas[hs[j]])
        outstr += ' 0 0 0'*(Nhalomax-len(theta))
        out.write(outstr+'\n')
out.close()
# plot a view
# import matplotlib.pyplot as plt
# plt.figure(figsize=(10, 10))
# plt.scatter(spos[:, 0], spos[:, 1], s=1, label='source')
# plt.scatter(hpos[:, 0], hpos[:, 1], s=10, label='lens')
# plt.xlabel('x [arcmin]')
# plt.ylabel('y [arcmin]')
# plt.legend()
# plt.show()