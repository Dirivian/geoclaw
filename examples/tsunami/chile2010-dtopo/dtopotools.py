from numpy import linspace, diff
import numpy as np
from pylab import find

class DTopo(object):

    def __init__(self, dtopo_params={}):
        self.dtopo_params = dtopo_params
        self.subfaults = []
        self.dz_list = []
        self.times = []
        self.x = None
        self.y = None


def read_dtopo(fname, dtopotype):
    if dtopotype==3:
        fid = open(fname)
        mx = int(fid.readline().split()[0])
        my = int(fid.readline().split()[0])
        mt = int(fid.readline().split()[0])
        xlower = float(fid.readline().split()[0])
        ylower = float(fid.readline().split()[0])
        t0 = float(fid.readline().split()[0])
        dx = float(fid.readline().split()[0])
        dy = float(fid.readline().split()[0])
        dt = float(fid.readline().split()[0])
        fid.close()

        xupper = xlower + (mx-1)*dx
        yupper = ylower + (my-1)*dy
        x=linspace(xlower,xupper,mx)
        y=linspace(ylower,yupper,my)
        times = linspace(t0, t0+(mt-1)*dt, mt)

        dZvals = np.loadtxt(fname, skiprows=9)
        dz_list = []
        for k,t in enumerate(times):
            dZk = np.reshape(dZvals[k*my:(k+1)*my, :], (my,mx))
            dZk = np.flipud(dZk)
            dz_list.append(dZk)
            
        dtopo = DTopo()
        dtopo.mx = mx
        dtopo.my = my
        dtopo.x = x
        dtopo.y = y
        dtopo.times = times
        dtopo.dz_list = dz_list
    elif dtopotype==1:
        txydz = np.loadtxt(fname)
        tdiff = np.diff(txydz[:,0])
        ind_tdiff = [0] + list(find(tdiff>0))
        xlist = txydz[0:ind_tdiff[1]+1,1]
        ylist = txydz[0:ind_tdiff[1]+1,2]
        ind_xdiff = [0] + list(find(diff(xlist)<0))
        mx = ind_xdiff[1] + 1
        my = int(len(xlist) / mx)
        x = xlist[0:ind_xdiff[1]+1]
        y = np.flipud(ylist[ind_xdiff])
        times = txydz[ind_tdiff,0]

        dz_list = []
        for k,t in enumerate(times):
            dZk = np.reshape(txydz[k*mx*my:(k+1)*mx*my, 3], (my,mx))
            dZk = np.flipud(dZk)
            dz_list.append(dZk)
            
        dtopo = DTopo()
        dtopo.mx = mx
        dtopo.my = my
        dtopo.x = x
        dtopo.y = y
        dtopo.times = times
        dtopo.dz_list = dz_list

    else:
        raise Exception("*** Unrecognized dtopotype: %s" % dtopotype)

    return dtopo
        
