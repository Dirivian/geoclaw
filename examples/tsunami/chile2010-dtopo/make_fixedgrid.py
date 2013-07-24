"""
Create fixed grid output file in new style.
"""

from numpy import linspace,sin,cos,meshgrid


def makefgrid():
    num_fgrids = 1
    xlower = -110
    xupper = -70
    ylower = -50
    yupper = -10
    dx = 0.5
    dy = dx
    mx = int(round((xupper-xlower)/dx)) + 1  
    my = int(round((yupper-ylower)/dy)) + 1  
    print "Creating fixed grid: %s by %s " % (mx,my)
    if abs((mx-1)*dx + xlower - xupper) > 1e-6:
        print "Warning: abs((mx-1)*dx + xlower - xupper) = ", \
              abs((mx-1)*dx + xlower - xupper)
    if abs((my-1)*dy + ylower - yupper) > 1e-6:
        print "Warning: abs((my-1)*dy + ylower - yupper) = ", \
              abs((my-1)*dy + ylower - yupper)
    tstart_max =   0.      # when to start monitoring max values
    tend_max = 1.e10       # when to stop monitoring max values
    dt_for_max = 60.       # target time increment between updating max values
    minlevel_for_max = 3   # which levels to monitor max on

    # Not yet used...
    tstart_output = 15.    # when to start output on fgrid
    tend_output = 15.      # when to stop output on fgrid
    dt_for_output = 15.    # target time increment between output
    minlevel_for_output = 3   # which levels to output fixed grid on


    x = linspace(xlower, xupper, mx)
    y = linspace(ylower, yupper, my)
    npts = mx*my

    fname = 'setfixedgrids2.data'
    fid = open(fname,'w')
    fid.write("%s               # num_fgrids\n" % num_fgrids)
    fid.write("%g  %g             # tstart_max, tend_max\n" \
               % (tstart_max, tend_max))
    fid.write("%g                 # dt_for_max\n" % dt_for_max)
    fid.write("%i                 # minlevel_for_max\n" % minlevel_for_max)
    fid.write("%g  %g             # tstart_output, tend_output\n" \
                % ( tstart_output, tend_output))
    fid.write("%g                 # dt_for_output\n" % dt_for_output)
    fid.write("%i                 # minlevel_for_output\n" % minlevel_for_output)
    fid.write("%g  %g  %g                  # npts, mx, my \n" % (npts,mx,my))



    for j in range(my):
        for i in range(mx):
            fid.write("%20.10e %20.10e\n" % (x[i],y[j]))


    print "Created file ", fname
    fid.close()

if __name__ == "__main__":
    makefgrid()


