import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from scipy.interpolate import interpn
from math import sqrt,pow,sin,cos,tan,atan2,asin,pi
from .BasicFig import BasicFig


class QualFig(BasicFig):
    """
    Inherit the plot format defined in the __init__ of BasicFig
    
    利用LES输出结果进行定性的分析，主要包括三个面的云图、风矢量图以及流线图等等
    """
    def __init__(self):
        super().__init__()  #继承LESfigure类的rcParams

    def VerCutPlane(self,
                    ax,
                    p,
                    field1,
                    yloc,
                    xlim=[],
                    zlim=[],
                    WTxc=[],
                    WTzc=[],
                    cmap="jet",
                    norm=None):
        """
        ax is the representation of fig
        
        p is the instance of Params
        
        WindStreamLine, WindVector are logical variable
        
        xlim, zlim(list) and yind are integral
        
        WTxc and WTzc(list) are nondimensional coordinates(using rotor diameter)
        
        cmap='jet' is alterable
        """
        zmax=np.shape(field1)[2]
        X = p.X23[xlim[0]:xlim[1], zlim[0]:zlim[1]]
        Z = p.Y23_uv[xlim[0]:xlim[1], zlim[0]:zlim[1]]
        points=np.zeros((p.nx,zmax,3))
        for i in range(p.nx):
            points[i,:,0]=p.x[i]
        points[:,:,1]=yloc
        for k in range(zmax):
            points[:,k,2]=p.z_uv[k]
        u=interpn([p.x,p.y,p.z_uv[0:zmax]],field1,points,bounds_error=False)
        self.cs=ax.contourf(X,
                    Z,
                    u[xlim[0]:xlim[1], zlim[0]:zlim[1]],
                    100,
                    cmap=cmap,
                    norm=norm)
        # ax.set_xlabel(r"$x/D$")
        # ax.set_ylabel(r"$z/D$")
        #TODO improve this
        if WTxc is not None:
            for i in range(len(WTxc)):
                self.VerWindTurbine(ax, WTxc[i], WTzc[i])
        # if WindVector:
        #     self.field2 = field2[xlim[0]:xlim[1], yind, zlim[0]:zlim[1]]
        #     self.field3 = field3[xlim[0]:xlim[1], yind, zlim[0]:zlim[1]]
        #     self.WindVector(ax, X, Z)
        # if WindStreamLine:
        #     X = p.X1[xlim[0]:xlim[1]]
        #     Z = p.Z1_uv[zlim[0]:zlim[1]]
        #     self.field2 = field2[xlim[0]:xlim[1], yind, zlim[0]:zlim[1]]
        #     self.field3 = field3[xlim[0]:xlim[1], yind, zlim[0]:zlim[1]]
        #     self.WindStreamLine(ax, X, Z)

    def HorCutPlane(self,
                    ax,
                    p,
                    field1,
                    zloc,
                    xlim=[],
                    ylim=[],
                    WTxc=[],
                    WTyc=[],
                    cmap="jet",
                    norm=None):
        """
        ax is the representation of fig
        
        p is the instance of Params
        
        WindStreamLine, WindVector are logical variable
        
        xlim, zlim(list) and yind are integral
        
        WTxc and WTzc(list) are nondimensional coordinates(using rotor diameter)
        
        cmap='jet' is alterable
        """
        zmax=np.shape(field1)[2]
        X = p.X21[xlim[0]:xlim[1], ylim[0]:ylim[1]]
        Y = p.Y21[xlim[0]:xlim[1], ylim[0]:ylim[1]]
        points=np.zeros((p.nx,p.ny,3))
        for i in range(p.nx):
            points[i,:,0]=p.x[i]
        for j in range(p.ny):
            points[:,j,1]=p.y[j]
        points[:,:,2]=zloc
        u=interpn([p.x,p.y,p.z_uv[0:zmax]],field1,points,bounds_error=False)
        self.cs=ax.contourf(X,
                    Y,
                    u[xlim[0]:xlim[1], ylim[0]:ylim[1]],
                    100,
                    cmap=cmap,
                    norm=norm)
        # ax.set_xlabel(r"$x/D$")
        # ax.set_ylabel(r"$y/D$")
        if WTxc is not None:
            for i in range(len(WTxc)):
                self.HorWindTurbine(ax, WTxc[i], WTyc[i])
        # if WindVector:
        #     self.field2 = field2[xlim[0]:xlim[1], ylim[0]:ylim[1], zind]
        #     self.field3 = field3[xlim[0]:xlim[1], ylim[0]:ylim[1], zind]
        #     self.WindVector(ax, X, Y)
        # if WindStreamLine:
        #     X = p.X1[xlim[0]:xlim[1]]
        #     Y = p.Y1[ylim[0]:ylim[1]]
        #     self.field2 = field2[xlim[0]:xlim[1], ylim[0]:ylim[1], zind]
        #     self.field3 = field3[xlim[0]:xlim[1], ylim[0]:ylim[1], zind]
        #     self.WindStreamLine(ax, X, Y)

    def CrossCutPlane(self,
                      ax,
                      p,
                      field1,
                      xloc,
                      ylim=[],
                      zlim=[],
                      WTyc=[],
                      WTzc=[],
                      cmap="coolwarm",
                      norm=None):
        """
        ax is the representation of fig
        
        p is the instance of Params
        
        WindStreamLine, WindVector are logical variables
        
        ylim, zlim(list) and xind are integral
        
        WTyc and WTzc(list) are nondimensional coordinates(using rotor diameter)
        
        cmap='jet' is alterable
        """
        zmax=np.shape(field1)[2]
        Y = p.X22[ylim[0]:ylim[1], zlim[0]:zlim[1]]
        Z = p.Y22_uv[ylim[0]:ylim[1], zlim[0]:zlim[1]]
        points=np.zeros((p.ny,zmax,3))
        points[:,:,0]=xloc
        for j in range(p.ny):
            points[j,:,1]=p.y[j]
        for k in range(zmax):
            points[:,k,2]=p.z_uv[k]
        u=interpn([p.x,p.y,p.z_uv[0:zmax]],field1,points,bounds_error=False)
        self.cs=ax.contourf(Y,
                    Z,
                    u[ylim[0]:ylim[1], zlim[0]:zlim[1]],
                    100,
                    cmap=cmap,
                    norm=norm)
        # ax.set_xlabel(r"$y/D$")
        # ax.set_ylabel(r"$z/D$")
        if WTyc is not None:
            for i in range(len(WTyc)):
                self.CrossWindTurbine(ax, WTyc[i], WTzc[i])
        # if WindVector:
        #     self.field2 = field2[xind, ylim[0]:ylim[1], zlim[0]:zlim[1]]
        #     self.field3 = field3[xind, ylim[0]:ylim[1], zlim[0]:zlim[1]]
        #     self.WindVector(ax, Y, Z)
        # if WindStreamLine:
        #     Y = p.Y1[ylim[0]:ylim[1]]
        #     Z = p.Z1_uv[zlim[0]:zlim[1]]
        #     self.field2 = field2[xind, ylim[0]:ylim[1], zlim[0]:zlim[1]]
        #     self.field3 = field3[xind, ylim[0]:ylim[1], zlim[0]:zlim[1]]
        #     self.WindStreamLine(ax, Y, Z)

    # def WindStreamLine(self, ax, a, b):
    #     """
    #     a, b are 1D arrays corresponding to the region where you are interested
    #     in
        
    #     streamplot(x, y, u, v, density=1, linewidth=None, color=None, cmap=None,
    #     norm=None, arrowsize=1, arrowstyle='-|>', minlength=0.1, transform=None, 
    #     zorder=None, start_points=None, maxlength=4.0, integration_direction='both', *, 
    #     data=None),
        
    #     x, y : 1D arrays
        
    #     An evenly spaced grid.
        
    #     u, v : 2D arrays
        
    #     x and y-velocities. The number of rows and columns must match the 
    #     ength of y and x, respectively.
    #     """
    #     u2 = np.transpose(self.field2)
    #     u3 = np.transpose(self.field3)
    #     self.strm = ax.streamplot(a,
    #                               b,
    #                               u2,
    #                               u3,
    #                               density=10,
    #                               color=np.hypot(u2, u3),
    #                               cmap='winter')

    # def WindVector(self, ax, a, b):
    #     """
    #     a, b are 2D arrays corresponding to the region where you are interested
    #     in
    #     """
    #     u2 = self.field2
    #     u3 = self.field3
    #     ax.quiver(a, b, u2, u3, units='xy', scale_units='xy')

    def VerWindTurbine(self, ax, xc, zc):
        #nacelle
        x = [xc, xc + 0.3]
        z = [zc - 0.075, zc - 0.075]
        ax.plot(x, z, 'k', linewidth=0.5)
        x = [xc, xc + 0.3]
        z = [zc + 0.075, zc + 0.075]
        ax.plot(x, z, 'k', linewidth=0.5)
        x = [xc + 0.3, xc + 0.3]
        z = [zc - 0.075, zc + 0.075]
        ax.plot(x, z, 'k', linewidth=0.5)
        theta = np.arange(np.pi / 2, np.pi * 3 / 2, 0.01)
        x = np.zeros(theta.size)
        z = np.zeros(theta.size)
        for i in range(theta.size):
            x[i] = 0.075 * np.cos(theta[i]) + xc
            z[i] = 0.075 * np.sin(theta[i]) + zc
        ax.plot(x, z, 'k', linewidth=0.5)
        ax.fill(x, z, 'w')
        ax.fill([xc, xc, xc + 0.3, xc + 0.3],
                [zc - 0.075, zc + 0.075, zc + 0.075, zc - 0.075], 'w')
        #tower
        x = [xc+0.1,xc+0.1]
        z = [0,zc-0.075]
        ax.plot(x,z,'k',linewidth=0.5)
        x = [xc+0.1,xc+0.2]
        z = [zc-0.075,zc-0.075]
        ax.plot(x,z,'k',linewidth=0.5)
        x = [xc+0.2,xc+0.2]
        z = [zc-0.075, 0]
        ax.plot(x,z,'k',linewidth=0.5)
        ax.fill([xc + 0.1, xc + 0.2, xc + 0.2, xc + 0.1],
                 [0, 0, zc - 0.075, zc - 0.075], 'w')
        #blade1
        x = [xc- 0.075*sin(pi/6), xc- 0.075*sin(pi/6)]
        z = [zc - 0.075*cos(pi/6), zc - 0.5 - 0.075]
        ax.plot(x, z, 'k', linewidth=0.5)
        x = [xc- 0.075*sin(pi/6), xc- 0.5*0.075*sin(pi/6)]
        z = [zc - 0.5 - 0.075, zc - 0.5 - 0.075]
        ax.plot(x, z, 'k', linewidth=0.5)
        x = [xc - 0.5*0.075*sin(pi/6), xc + 0.075]
        z = [zc - 0.5 - 0.075, zc - 0.15 -0.075]
        ax.plot(x, z, 'k', linewidth=0.5)
        x = [xc + 0.075, xc + 0.075*sin(pi/6)]
        z = [zc - 0.15 - 0.075, zc - 0.05 -0.075] 
        ax.plot(x, z, 'k', linewidth=0.5)
        x = [xc + 0.075 * sin(pi/6), xc + 0.075 * sin(pi/6)]
        z = [zc - 0.05 -0.075, zc - 0.075]
        ax.plot(x, z, 'k', linewidth=0.5)
        ax.fill([xc- 0.075*sin(pi/6), xc- 0.075*sin(pi/6), xc- 0.5*0.075*sin(pi/6), 
                 xc + 0.075, xc + 0.075 * sin(pi/6), xc + 0.075 * sin(pi/6)],
                [zc - 0.075*cos(pi/6), zc - 0.5 - 0.075, zc - 0.5 - 0.075, 
                 zc - 0.15 - 0.075, zc - 0.05 -0.075, zc - 0.075],
                'w')
        #blade2
        x = [xc - 0.075*sin(pi/6), xc - 0.075*sin(pi/6)]
        z = [zc - 0.075*cos(pi/6), zc + 0.5 + 0.075]
        ax.plot(x, z, 'k', linewidth=0.5)
        x = [xc- 0.075*sin(pi/6), xc - 0.5*0.075*sin(pi/6)]
        z = [zc + 0.5 + 0.075, zc + 0.5 + 0.075]
        ax.plot(x, z, 'k', linewidth=0.5)
        x = [xc - 0.5*0.075*sin(pi/6), xc + 0.075]
        z = [zc + 0.5 + 0.075, zc + 0.15 + 0.075]
        ax.plot(x, z, 'k', linewidth=0.5)
        x = [xc + 0.075, xc + 0.075*sin(pi/6)]
        z = [zc + 0.15 + 0.075, zc + 0.05 + 0.075] 
        ax.plot(x, z, 'k', linewidth=0.5)
        x = [xc + 0.075 * sin(pi/6), xc + 0.075 * sin(pi/6)]
        z = [zc + 0.05 + 0.075, zc + 0.075]
        ax.plot(x, z, 'k', linewidth=0.5)
        ax.fill([xc - 0.075*sin(pi/6), xc - 0.075*sin(pi/6), xc- 0.5*0.075*sin(pi/6),
                 xc + 0.075, xc + 0.075 * sin(pi/6), xc + 0.075 * sin(pi/6)],
                [zc + 0.075*cos(pi/6), zc + 0.5 + 0.075, zc + 0.5 + 0.075,
                 zc + 0.15 + 0.075, zc + 0.05 + 0.075, zc + 0.075],
                'w')

    def HorWindTurbine(self, ax, xc, yc):
        #nacelle
        x = [xc, xc + 0.3]
        z = [yc - 0.075, yc - 0.075]
        ax.plot(x, z, 'k', linewidth=0.5)
        x = [xc, xc + 0.3]
        z = [yc + 0.075, yc + 0.075]
        ax.plot(x, z, 'k', linewidth=0.5)
        x = [xc + 0.3, xc + 0.3]
        z = [yc - 0.075, yc + 0.075]
        ax.plot(x, z, 'k', linewidth=0.5)
        theta = np.arange(np.pi / 2, np.pi * 3 / 2, 0.01)
        x = np.zeros(theta.size)
        z = np.zeros(theta.size)
        for i in range(theta.size):
            x[i] = 0.075 * np.cos(theta[i]) + xc
            z[i] = 0.075 * np.sin(theta[i]) + yc
        ax.plot(x, z, 'k', linewidth=0.5)
        ax.fill(x, z, 'w')
        ax.fill([xc, xc, xc + 0.3, xc + 0.3],
                [yc - 0.075, yc + 0.075, yc + 0.075, yc - 0.075], 'w')
        
        #blade1
        x = [xc- 0.075*sin(pi/6), xc- 0.075*sin(pi/6)]
        z = [yc - 0.075*cos(pi/6), yc - 0.5 - 0.075]
        ax.plot(x, z, 'k', linewidth=0.5)
        x = [xc- 0.075*sin(pi/6), xc- 0.5*0.075*sin(pi/6)]
        z = [yc - 0.5 - 0.075, yc - 0.5 - 0.075]
        ax.plot(x, z, 'k', linewidth=0.5)
        x = [xc - 0.5*0.075*sin(pi/6), xc + 0.075]
        z = [yc - 0.5 - 0.075, yc - 0.15 -0.075]
        ax.plot(x, z, 'k', linewidth=0.5)
        x = [xc + 0.075, xc + 0.075*sin(pi/6)]
        z = [yc - 0.15 - 0.075, yc - 0.05 -0.075] 
        ax.plot(x, z, 'k', linewidth=0.5)
        x = [xc + 0.075 * sin(pi/6), xc + 0.075 * sin(pi/6)]
        z = [yc - 0.05 -0.075, yc - 0.075]
        ax.plot(x, z, 'k', linewidth=0.5)
        ax.fill([xc- 0.075*sin(pi/6), xc- 0.075*sin(pi/6), xc- 0.5*0.075*sin(pi/6), 
                 xc + 0.075, xc + 0.075 * sin(pi/6), xc + 0.075 * sin(pi/6)],
                [yc - 0.075*cos(pi/6), yc - 0.5 - 0.075, yc - 0.5 - 0.075, 
                 yc - 0.15 - 0.075, yc - 0.05 -0.075, yc - 0.075],
                'w')
        #blade2
        x = [xc - 0.075*sin(pi/6), xc - 0.075*sin(pi/6)]
        z = [yc - 0.075*cos(pi/6), yc + 0.5 + 0.075]
        ax.plot(x, z, 'k', linewidth=0.5)
        x = [xc- 0.075*sin(pi/6), xc - 0.5*0.075*sin(pi/6)]
        z = [yc + 0.5 + 0.075, yc + 0.5 + 0.075]
        ax.plot(x, z, 'k', linewidth=0.5)
        x = [xc - 0.5*0.075*sin(pi/6), xc + 0.075]
        z = [yc + 0.5 + 0.075, yc + 0.15 + 0.075]
        ax.plot(x, z, 'k', linewidth=0.5)
        x = [xc + 0.075, xc + 0.075*sin(pi/6)]
        z = [yc + 0.15 + 0.075, yc + 0.05 + 0.075] 
        ax.plot(x, z, 'k', linewidth=0.5)
        x = [xc + 0.075 * sin(pi/6), xc + 0.075 * sin(pi/6)]
        z = [yc + 0.05 + 0.075, yc + 0.075]
        ax.plot(x, z, 'k', linewidth=0.5)
        ax.fill([xc - 0.075*sin(pi/6), xc - 0.075*sin(pi/6), xc- 0.5*0.075*sin(pi/6),
                 xc + 0.075, xc + 0.075 * sin(pi/6), xc + 0.075 * sin(pi/6)],
                [yc + 0.075*cos(pi/6), yc + 0.5 + 0.075, yc + 0.5 + 0.075,
                 yc + 0.15 + 0.075, yc + 0.05 + 0.075, yc + 0.075],
                'w')

    def CrossWindTurbine(self, ax, yc, zc):
        theta = np.arange(0, 2 * np.pi, 0.01)
        y = np.zeros(theta.size)
        z = np.zeros(theta.size)
        for i in range(theta.size):
            y[i] = 0.5 * np.cos(theta[i]) + yc
            z[i] = 0.5 * np.sin(theta[i]) + zc
        ax.plot(y, z, 'k--')
        #self.threedimturbine(yc, zc)
    def rotateblade(self, points, center, angle):
        for i in points:
            r=sqrt(pow(i[0]-center[0],2)+pow(i[1]-center[1],2))
            theta=atan2(i[1]-center[1],i[0]-center[0])
            theta=theta+angle
            i[0]=center[0]+r*cos(theta)
            i[1]=center[1]+r*sin(theta)
        return points
    
    def threedimturbine(self, xloc, yloc):
        hhub=yloc
        lblade=0.5
        rhub=0.075
        wtower=0.1
        #tower line
        p2=[xloc+wtower/2,yloc]
        p3=[xloc+wtower/2,yloc+hhub-rhub*cos(asin(wtower/2/rhub))]
        p0=[xloc-wtower/2,yloc+hhub-rhub*cos(asin(wtower/2/rhub))]
        p1=[xloc-wtower/2,yloc]
        p=[p0,p1,p2,p3]
        for i in range(4):
            if i<3:
                lx=[p[i][0],p[i+1][0]];ly=[p[i][1],p[i+1][1]]
            else:
                #lx=[p[i][0],p[0][0]];ly=[p[i][1],p[0][1]]
                pass
            plt.plot(lx,ly,'k',linewidth=0.5)
        x=np.zeros(4)
        y=np.zeros(4)
        for i in range(4):
            x[i]=p[i][0]
            y[i]=p[i][1]
        plt.fill(x, y, 'w')
        #nacelle
        theta = np.arange(0, 2*pi, 0.01)
        x = np.zeros(theta.size)
        y = np.zeros(theta.size)
        for i in range(theta.size):
            x[i] = rhub * np.cos(theta[i]) +xloc
            y[i] = rhub * np.sin(theta[i]) +yloc+hhub
        plt.plot(x, y, 'k',linewidth=0.5)
        plt.fill(x, y, 'w')
        #blade 1
        p0=[xloc-wtower*sin(pi/6),yloc+hhub+wtower*cos(pi/6)]
        p1=[xloc-wtower*sin(pi/6),yloc+hhub+wtower*cos(pi/6)+wtower]
        p2=[xloc-wtower*sin(pi/6)-wtower,yloc+hhub+wtower*cos(pi/6)+wtower+3/2*wtower]
        p3=[xloc-wtower/2*sin(pi/6),yloc+hhub+lblade+wtower/2*cos(pi/6)]
        p4=[xloc+wtower/2*sin(pi/6),yloc+hhub+lblade+wtower/2*cos(pi/6)]
        p5=[xloc+wtower*sin(pi/6),yloc+hhub+wtower*cos(pi/6)]
        b1=[p0,p1,p2,p3,p4,p5]
        for i in range(6):
            if i<5:
                lx=[b1[i][0],b1[i+1][0]];ly=[b1[i][1],b1[i+1][1]]
            else:
                #lx=[b1[i][0],b1[0][0]];ly=[b1[i][1],b1[0][1]]
                pass
            plt.plot(lx,ly,'k',linewidth=0.5)
        x=np.zeros(6)
        y=np.zeros(6)
        for i in range(6):
            x[i]=b1[i][0]
            y[i]=b1[i][1]
        plt.fill(x, y, 'w')
        #blade 2
        b2 = self.rotateblade(b1, [xloc,yloc+hhub], 2*pi/3)
        for i in range(6):
            if i<5:
                lx=[b2[i][0],b2[i+1][0]];ly=[b2[i][1],b2[i+1][1]]
            else:
                pass
            #lx=[b2[i][0],b2[0][0]];ly=[b2[i][1],b2[0][1]]
            plt.plot(lx,ly,'k',linewidth=0.5)
        x=np.zeros(6)
        y=np.zeros(6)
        for i in range(6):
            x[i]=b2[i][0]
            y[i]=b2[i][1]
        plt.fill(x, y, 'w')
        #blade 3
        b3 = self.rotateblade(b2, [xloc,yloc+hhub], 2*pi/3)
        for i in range(6):
            if i<5:
                lx=[b3[i][0],b3[i+1][0]];ly=[b3[i][1],b3[i+1][1]]
            else:
                #lx=[b3[i][0],b3[0][0]];ly=[b3[i][1],b3[0][1]]
                pass
            plt.plot(lx,ly,'k',linewidth=0.5)
        x=np.zeros(6)
        y=np.zeros(6)
        for i in range(6):
            x[i]=b3[i][0]
            y[i]=b3[i][1]
        plt.fill(x, y, 'w')
    
class CmapNorm:
    #??? how to make it come true
    #TODO basing the number of visualization fields yield one colorbar automatically
    def __init__(self, cases, plane, index, ran1, ran2):
        '''
        create a user defined colorbar for normalize
        cases is a list consists of different scenes
        plane is one of "ver", "hor" or "cross"
        index, ran1 and ran2 are even permutation 
        '''
        self.code = ""
        if plane == "ver":
            ran = "[ran2[0]:ran2[1],index,ran1[0]:ran1[1]],"
        elif plane == "hor":
            ran = "[ran1[0]:ran1[1],ran2[0]:ran2[1],index],"
        else:
            ran = "[index,ran1[0]:ran1[1],ran2[0]:ran2[1]],"
        for i in cases:
            self.code += i + ran
        print(self.code)

    def create(self):
        self.vmin = min(np.concatenate(eval(self.code)))
        self.vmax = max(np.concatenate(eval(self.code)))
        self.norm = matplotlib.colors.Normalize(self.vmin, self.vmax)


#    def VerCmapNorm(self,cases,yind,xlim=[],zlim=[]):
#        for i in len(cases):
#            codemin = "np.min(cases[i][xlim[0]:xlim[1],yind,zlim[0]:zlim[1]]),"
#
#        vmin=min(np.min(a[xlim[0]:xlim[1],yind,zlim[0]:zlim[1]]),
#                 np.min(b[xlim[0]:xlim[1],yind,zlim[0]:zlim[1]]))
#        vmax=max(np.max(a[xlim[0]:xlim[1],yind,zlim[0]:zlim[1]]),
#                          np.max(b[xlim[0]:xlim[1],yind,zlim[0]:zlim[1]]))
#
#    def HorCmapNorm(self,a,b,zind,xlim=[],ylim=[]):
#        vmin=min(np.min(a[xlim[0]:xlim[1],ylim[0]:ylim[1],zind]),
#                 np.min(b[xlim[0]:xlim[1],ylim[0]:ylim[1],zind]))
#        vmax=max(np.max(a[xlim[0]:xlim[1],ylim[0]:ylim[1],zind]),
#                           np.max(b[xlim[0]:xlim[1],ylim[0]:ylim[1],zind]))
#        self.norm = matplotlib.colors.Normalize(vmin, vmax)
#    def CrossCmapNorm(self,a,b,xind,ylim=[],zlim=[]):
#        vmin=min(np.min(a[xind,ylim[0]:ylim[1],zlim[0]:zlim[1]]),
#                 np.min(b[xind,ylim[0]:ylim[1],zlim[0]:zlim[1]]))
#        vmax=max(np.max(a[xind,ylim[0]:ylim[1],zlim[0]:zlim[1]]),
#                 np.max(b[xind,ylim[0]:ylim[1],zlim[0]:zlim[1]]))
#        self.norm = matplotlib.colors.Normalize(vmin, vmax)