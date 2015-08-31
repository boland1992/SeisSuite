# -*- coding: utf-8 -*-
"""
Created on Mon Aug  3 08:54:34 2015

@author: boland
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from descartes.patch import PolygonPatch

class Density:
    """
    Class defined to perform to density field operations e.g. 2d histogram, 
    gradient fields and averages and standard deviations thereof. 
    """
    def __init__(self, paths=None, nbins=200):
        # initialise path points for shape (2,N) to calculate point density
        self.paths = paths
        # initialise the number of bins per axis for 2d histogram
        self.nbins = nbins
        # initialise the density calculation
        self.H = 0.
        self.H_masked = 0.
        # initialise the density gradient field calculation
        self.grad = 0.
        self.grad_masked = 0.
        # initialise x and y tridiagonal coordinate matrices
        self.xedges = 0.
        self.yedges = 0.
        
    def hist2d(self, paths=None):
        """
        Function that calculates the 2D histogram and generates H, xedges, 
        and yedges. 
        """
        if paths is None:
            paths = self.paths
        self.H, self.xedges, self.yedges = np.histogram2d(paths[:,0],
                                           paths[:,1],
                                           bins=self.nbins)
        return self.H, self.xedges, self.yedges
            
    def hgrad(self, H=None):
        """
        Function that calculates the 2D histogram and generates H, xedges, 
        and yedges. 
        """
        #, xedges=None, yedges=None
        if H is None: 
            H, xedges, yedges = self.hist2d()
            
        self.grad = np.abs(np.asarray(np.gradient(H)[0]))
        
        return self.grad
            
    def transform_h(self, H=None):
        """
        Function that rotates, flips and masks the H density field 
        in order for it to be plotted etc.
        """
        if H is None:
            if all( [self.H != 0., self.xedges != 0, self.yedges != 0] ):
                H, xedges, yedges = self.H, self.xedges, self.yedges
            else:
                H, xedges, yedges = self.hist2d()
            
        H = np.rot90(H)
        H = np.flipud(H)
        self.H_masked = np.ma.masked_where(H==0,H)  
        return self.H_masked
        
    def transform_grad(self, grad=None):
        """
        Function that rotates, flips and masks the H density gradient field 
        in order for it to be plotted etc.
        """
        if grad is None: 
            grad, xedges, yedges = self.hgrad()
        grad = np.rot90(grad)
        grad = np.flipud(grad)
        self.grad_masked = np.ma.masked_where(grad==0,grad)
        return self.grad_masked
    

    def plot_lims(self, paths=None):
        if paths is None:
            try:
                lons, lats = self.paths[:,0], self.paths[:,1]
            except Exception as error:
                raise error
        else:
            try:
                lons, lats = paths[:,0], paths[:,1]
            except Exception as error:
                raise error
            
        return np.min(lons), np.max(lons), np.min(lats), np.max(lats)
    
    def select_points(self, perc=0.1, high=None):
        """
        Function that returns the lat-lon coordinates of points below a 
        certain density. This is taken by H < perc*np.average(H) OR if
        high=True or is not None, then high density points are chosen from
        H >perc*np.average(H). Perc=0.1 by default.
        """
        
        if high is None:
            search = np.where(self.H<perc*np.average(self.H))
        else:
            search = np.where(self.H>perc*np.average(self.H))
         
        xmin, xmax =  np.min(self.xedges), np.max(self.xedges)
        ymin, ymax =  np.min(self.yedges), np.max(self.yedges) 

        Hdensx, Hdensy =  search[1], search[0]    
        Hdensx = (xmax-xmin)/(self.nbins) * Hdensx + xmin
        Hdensy = (ymax-ymin)/(self.nbins) * Hdensy + ymin
        return np.column_stack((Hdensx, Hdensy))

    
    def plot_field(self, grad=False, SHAPE=None, swell=0.05):
        
        lonmin, lonmax, latmin, latmax = self.plot_lims()
        
        fig = plt.figure(figsize=(15,10), dpi=100)
        
        plt.xlabel('longitude (degrees)')
        plt.ylabel('latitude (degrees)')
        plt.xlim(lonmin-swell*abs(lonmax-lonmin),\
                 lonmax+swell*abs(lonmax-lonmin))
        plt.ylim(latmin-swell*abs(latmax-latmin),\
                 latmax+swell*abs(latmax-latmin))
    
        if not grad:
            plt.title("Path Density Distribution")
            if self.H_masked is not 0.:
                H_masked = self.H_masked
            else:                    
                H_masked = self.transform_h()
            plt.pcolor(self.xedges, self.yedges, H_masked, norm=LogNorm(\
            vmin=np.min(H_masked), vmax=np.max(H_masked)), cmap='rainbow',\
            alpha=0.6, zorder = 3)
            col = plt.colorbar()
            col.ax.set_ylabel('Points Per Bin')
        elif grad:
            plt.title("Gradient Path Density Distribution")
            
            if self.grad_masked is not 0.:
                grad_masked = self.grad_masked
            else:           
                raise Exception("grad_masked has not yet been defined. please\
                run the necessary functions e.g. transform_grad before plotting")
                
            plt.pcolor(self.xedges, self.yedges, grad_masked, norm=LogNorm(\
            vmin=np.min(grad), vmax=np.max(grad)), cmap='rainbow',\
            alpha=0.6, zorder = 3)
            col = plt.colorbar()
            col.ax.set_ylabel('Gradient Points Per Bin')      
        else:
            raise Exception("Either you have not chosen to plot anything OR\n\
                   both H and grad are inputed and the function doesn't\
                   know what to do.")

        if SHAPE is not None:
            patch = PolygonPatch(SHAPE, facecolor='white',\
            edgecolor='k', zorder=1)
            ax = fig.add_subplot(111)
            ax.add_patch(patch)
            ax.set_xlim(lonmin-0.05*abs(lonmax-lonmin), \
                        lonmax+0.05*abs(lonmax-lonmin))
            ax.set_ylim(latmin-0.05*abs(latmax-latmin), \
                        latmax+0.05*abs(latmax-latmin))
                            
        #plt.scatter(new_coords[:,0], new_coords[:,1],c='r', s=30)
        fig.savefig("plot_density.png")
        fig.clf()