"""
Plot transect of seafloor deformation and surface at various times,
as recorded by using the "fixed grid" output feature.  Set the fg so that 
it captures transect by setting ypoints=1, for example in setrun.py:
    # == setfixedgrids.data values ==
    geodata.fixedgrids = []
    # for fixed grids append lines of the form
    # [t1,t2,noutput,x1,x2,y1,y2,xpoints,ypoints,ioutarrivaltimes,ioutsurfacemax]
    geodata.fixedgrids.append([0., 1.4, 11, 140., 145., 37., 37., 500, 1, 0, 0])
"""

from pylab import *
import glob

fgfiles = glob.glob('_output/fort.fg*')
fgnum = len(fgfiles)
if fgnum==0:
    print "*** No fg files found in _output"

figure(23)
clf()
ylim(-2,3)
title('surface')

figure(24)
clf()
ylim(-2,3)
title('seafloor deformation')

fg0 = loadtxt('_output/fort.fg01_0001', skiprows=9)

plots1 = []
plots2 = []
for i in range(1,fgnum+1):
    fgfile = '_output/fort.fg01_%s' % str(i).zfill(4)
    line = open(fgfile).readline()
    print line
    time = line.split()[0]
    print "Time = ", time
    fg = loadtxt(fgfile, skiprows=9)

    figure(23)
    plots1 = plots1 + plot(fg[:,0]+fg[:,3], 'r') 
    for p in plots1[:-1]: 
        p.set_color('b')
    show()

    figure(24)
    plots2 = plots2 + plot(fg[:,3]-fg0[:,3], 'r') 
    for p in plots2[:-1]: 
        p.set_color('b')
    show()
    ans = raw_input('Hit return ')
    if ans=='q': break


# Add final deformation from dtopo file to plot:
# Not working since plots above don't include x, would need to 
# construct x from mx, xlow, xhi in fg files.

if 0:
    y_dtopo_slice = 37.0

    try:
        import dtopotools as D
        dtopo = D.read_dtopo('disp_source.dtopo3',3)
        dz = dtopo.dz_list[-1]
        xt = dtopo.x
        yt = dtopo.y
        j_dtopo_slice = min(find(yt>y_dtopo_slice))
    except:
        print "*** could not read dtopo file!"

    figure(24)
    plot(xt,dz[j_dtopo_slice,:],'g',label='final deformation from dtopo')
