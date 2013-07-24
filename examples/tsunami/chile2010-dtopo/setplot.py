
""" 
Set up the plot figures, axes, and items to be done for each frame.

This module is imported by the plotting routines and then the
function setplot is called to set the plot parameters.
    
""" 

import numpy as np
import matplotlib.pyplot as plt

from clawpack.geoclaw import topotools
from clawpack.clawutil import clawdata

from pylab import find


# For slice of dtopo:

y_dtopo_slice = -35.0

try:
    import dtopotools as D
    dtopo = D.read_dtopo('usgs100227.tt1',1)
    dz = dtopo.dz_list[-1]
    xt = dtopo.x
    yt = dtopo.y
    j_dtopo_slice = min(find(yt>y_dtopo_slice))
except:
    print "*** could not read dtopo file!"



#--------------------------
def setplot(plotdata):
#--------------------------
    
    """ 
    Specify what is to be plotted at each frame.
    Input:  plotdata, an instance of pyclaw.plotters.data.ClawPlotData.
    Output: a modified version of plotdata.
    
    """ 


    from clawpack.visclaw import colormaps, geoplot
    from numpy import linspace

    plotdata.clearfigures()  # clear any old figures,axes,items data


    # To plot gauge locations on pcolor or contour plot, use this as
    # an afteraxis function:

    def addgauges(current_data):
        from clawpack.visclaw import gaugetools
        gaugetools.plot_gauge_locations(current_data.plotdata, \
             gaugenos='all', format_string='ko', add_labels=True)
    
    # ========================================================================
    #  Water helper functions
    # ========================================================================
    def b(cd):
        return cd.q[3,:,:] - cd.q[0,:,:]
        
    def extract_eta(h,eta,DRY_TOL=10**-3):
        index = np.nonzero((np.abs(h) < DRY_TOL) + (h == np.nan))
        eta[index[0],index[1]] = np.nan
        return eta
    
    def extract_velocity(h,hu,DRY_TOL=10**-8):
        u = np.zeros(hu.shape)
        index = np.nonzero((np.abs(h) > DRY_TOL) * (h != np.nan))
        u[index[0],index[1]] = hu[index[0],index[1]] / h[index[0],index[1]]
        return u
    
    def eta(cd):
        return extract_eta(cd.q[0,:,:],cd.q[3,:,:])
        
    def water_u(cd):
        return extract_velocity(cd.q[0,:,:],cd.q[1,:,:])
        
    def water_v(cd):
        return extract_velocity(cd.q[0,:,:],cd.q[2,:,:])
        
    def water_speed(current_data):
        u = water_u(current_data)
        v = water_v(current_data)
            
        return np.sqrt(u**2+v**2)


    #-----------------------------------------
    # Figure for surface
    #-----------------------------------------
    plotfigure = plotdata.new_plotfigure(name='Surface', figno=0)

    # Set up for axes in this figure:
    plotaxes = plotfigure.new_plotaxes('pcolor')
    plotaxes.title = 'Surface'
    plotaxes.scaled = True

    def fixup(current_data):
        import pylab
        addgauges(current_data)
        t = current_data.t
        t = t / 3600.  # hours
        pylab.title('Surface at %4.2f hours' % t, fontsize=20)
        pylab.xticks(fontsize=15)
        pylab.yticks(fontsize=15)
    plotaxes.afteraxes = fixup

    # Water
    plotitem = plotaxes.new_plotitem(plot_type='2d_pcolor')
    # plotitem.plot_var = geoplot.surface
    plotitem.plot_var = geoplot.surface_or_depth
    plotitem.pcolor_cmap = geoplot.tsunami_colormap
    plotitem.pcolor_cmin = -0.2e0
    plotitem.pcolor_cmax = 0.2e0
    plotitem.add_colorbar = True
    plotitem.amr_celledges_show = [0,0,0]
    plotitem.patchedges_show = 1

    # Land
    plotitem = plotaxes.new_plotitem(plot_type='2d_pcolor')
    plotitem.plot_var = geoplot.land
    plotitem.pcolor_cmap = geoplot.land_colors
    plotitem.pcolor_cmin = 0.0
    plotitem.pcolor_cmax = 100.0
    plotitem.add_colorbar = False
    plotitem.amr_celledges_show = [0,0,0]
    plotitem.patchedges_show = 1
    plotaxes.xlimits = [-120,-60]
    plotaxes.ylimits = [-60,0]

    # add contour lines of bathy if desired:
    plotitem = plotaxes.new_plotitem(plot_type='2d_contour')
    plotitem.show = False
    plotitem.plot_var = geoplot.topo
    plotitem.contour_levels = linspace(-3000,-3000,1)
    plotitem.amr_contour_colors = ['y']  # color on each level
    plotitem.kwargs = {'linestyles':'solid','linewidths':2}
    plotitem.amr_contour_show = [1,0,0]  
    plotitem.celledges_show = 0
    plotitem.patchedges_show = 0

    #-----------------------------------------
    # Figure for velocities
    #-----------------------------------------
    plotfigure = plotdata.new_plotfigure(name='Speeds', figno=1)
    plotfigure.show = False

    # Set up for axes in this figure:
    plotaxes = plotfigure.new_plotaxes('speeds')
    plotaxes.title = 'Speeds'
    plotaxes.scaled = True

    def fixup(current_data):
        import pylab
        addgauges(current_data)
        t = current_data.t
        t = t / 3600.  # hours
        pylab.title('Speeds at %4.2f hours' % t, fontsize=20)
        pylab.xticks(fontsize=15)
        pylab.yticks(fontsize=15)
    plotaxes.afteraxes = fixup

    # Speed
    plotitem = plotaxes.new_plotitem(plot_type='2d_pcolor')
    plotitem.plot_var = water_speed
    plotitem.pcolor_cmap = plt.get_cmap('PuBu')
    plotitem.pcolor_cmin = 0.0
    plotitem.pcolor_cmax = 0.01
    plotitem.add_colorbar = True
    plotitem.amr_celledges_show = [0,0,0]
    plotitem.patchedges_show = 1

    # Land
    plotitem = plotaxes.new_plotitem(plot_type='2d_pcolor')
    plotitem.plot_var = geoplot.land
    plotitem.pcolor_cmap = geoplot.land_colors
    plotitem.pcolor_cmin = 0.0
    plotitem.pcolor_cmax = 100.0
    plotitem.add_colorbar = False
    plotitem.amr_celledges_show = [1,1,0]
    plotitem.patchedges_show = 1
    plotaxes.xlimits = [-120,-60]
    plotaxes.ylimits = [-60,0]


    #-----------------------------------------
    # Figures for gauges
    #-----------------------------------------
    plotfigure = plotdata.new_plotfigure(name='Surface & topo', figno=300, \
                    type='each_gauge')
    plotfigure.clf_each_gauge = True

    # Set up for axes in this figure:
    plotaxes = plotfigure.new_plotaxes()
    plotaxes.xlimits = 'auto'
    plotaxes.ylimits = 'auto'
    plotaxes.title = 'Surface'

    # Plot surface as blue curve:
    plotitem = plotaxes.new_plotitem(plot_type='1d_plot')
    plotitem.plot_var = 3
    plotitem.plotstyle = 'b-'

    # Plot topo as green curve:
    plotitem = plotaxes.new_plotitem(plot_type='1d_plot')
    plotitem.show = False

    def gaugetopo(current_data):
        q = current_data.q
        h = q[0,:]
        eta = q[3,:]
        topo = eta - h
        return topo
        
    # plotitem.plot_var = gaugetopo
    # plotitem.plotstyle = 'g-'

    def add_zeroline(current_data):
        from pylab import plot, legend, xticks, floor
        t = current_data.t
        #legend(('surface','topography'),loc='lower left')
        plot(t, 0*t, 'k')
        n = int(floor(t.max()/3600.) + 2)
        xticks([3600*i for i in range(n)])

    plotaxes.afteraxes = add_zeroline

    #-----------------------------------------
    # Figure for line plot
    #-----------------------------------------
    plotfigure = plotdata.new_plotfigure(name='line', figno=22)
    #plotfigure.show = False

    def slice(current_data):
        """
        Return a slice of x,eta,B at a fixed y value y_dtopo_slice.
        Does not interpolate, uses minimum j for which y(j) > y_dtopo_slice.
        Return arrays of nan's if this grid doesn't include the y value.
        """
        from pylab import where,find,nan,array
        x = current_data.x[:,0]
        y = current_data.y[0,:]
        #print "+++ y.min() = ",y.min()
        if (y.min() <= y_dtopo_slice):
            jfind = find(y > y_dtopo_slice)
            if len(jfind) > 0:
                j = min(jfind)
                #print '+++ level, j, ymin, ymax: ', \
                #        current_data.level,j,y.min(),y.max()
                q = current_data.q
                h = q[0,:,j]
                eta = q[3,:,j]
                eta = where(h>0, eta, 0.)
                B = eta - h
                return x,eta,B
            else:
                return array([nan]),array([nan]),array([nan])
        else:
            return array([nan]),array([nan]),array([nan])

    def eta_slice(current_data):
        x,eta,B = slice(current_data)
        return x,eta

    def B_slice(current_data):
        x,eta,B = slice(current_data)
        return x,B

    def plot_dtopo_slice(current_data):
        from pylab import find,figure,subplot,plot,legend,xlim
        #subplot(211)
        plot(xt,dz[j_dtopo_slice,:],'r',label='from dtopo')
        xlim(-82,-72)
        legend()

    # Set up for axes in this figure:
    plotaxes = plotfigure.new_plotaxes('line')
    #plotaxes.axescmd = 'subplot(211)'
    plotaxes.afteraxes = plot_dtopo_slice
    plotaxes.title = 'Surface transect at y = %s' % y_dtopo_slice
    plotaxes.xlimits = [-120,-60]


    # Water
    plotitem = plotaxes.new_plotitem(plot_type='1d_from_2d_data')
    #plotitem.show = False
    plotitem.map_2d_to_1d = eta_slice
    plotitem.color = 'b'
    plotitem.kwargs = {'linewidth':2}

 


    #-----------------------------------------
    
    # Parameters used only when creating html and/or latex hardcopy
    # e.g., via pyclaw.plotters.frametools.printframes:

    plotdata.printfigs = True                # print figures
    plotdata.print_format = 'png'            # file format
    plotdata.print_framenos = 'all'          # list of frames to print
    plotdata.print_gaugenos = 'all'          # list of gauges to print
    plotdata.print_fignos = 'all'            # list of figures to print
    plotdata.html = True                     # create html files of plots?
    plotdata.html_homelink = '../README.html'   # pointer for top of index
    plotdata.latex = True                    # create latex file of plots?
    plotdata.latex_figsperline = 2           # layout of plots
    plotdata.latex_framesperline = 1         # layout of plots
    plotdata.latex_makepdf = False           # also run pdflatex?

    return plotdata

