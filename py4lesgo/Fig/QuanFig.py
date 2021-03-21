# coding: utf-8
import numpy as np
from scipy import interpolate
from scipy.interpolate import interpn
from .BasicFig import BasicFig


class QuanFig(BasicFig):
    """
    """
    def __init__(self):
        """
    Inherit the plot format defined in the __init__ of BasiFig
    """
        super().__init__()
    def casefig(self, p, axs, tavg, hub_height, ylim, dim, case, ABL, theta=True):
        """
        p is the instance of Params
        
        axs is the representation of fig
        
        tavg is the instance of Tavg
        
        hub_height is a parameter using diameter as nondimensional length
        
        ylim is the range you want to look
        
        dim is the instance of DataProcess
        """
        zmax=np.shape(tavg.u)[2]
        if "CBL" in case:
            color='red'
        elif "NBL-H" in case:
            color='green'
        elif "NBL-L" in case:
            color='purple'
        elif "SBL" in case:
            color='blue'
        axs[0,0].plot(tavg.uMean * dim.u_star, p.Z1_uv[0:zmax], color, label=case)
        axs[0,0].plot([0, 15], [hub_height - 0.5, hub_height - 0.5], 'k--')
        axs[0,0].plot([0, 15], [hub_height, hub_height], 'k--')
        axs[0,0].plot([0, 15], [hub_height + 0.5, hub_height + 0.5], 'k--')
        axs[0,0].set(xlim=[0,15])
        
        
        
        if theta:
            if ABL=="CBL":
                thetal=290
                thetar=315
                TIr=0.2
                TIl=0.0
            else:
                thetal=285
                thetar=295
                TIr=0.20
                TIl=0.0
            
            axs[0,1].plot(tavg.thetaMean * 300, p.Z1_uv[0:zmax], color, label=case)
            axs[0,1].plot([thetal, thetar], [hub_height - 0.5, hub_height - 0.5], 'k--')
            axs[0,1].plot([thetal, thetar], [hub_height, hub_height], 'k--')
            axs[0,1].plot([thetal, thetar], [hub_height + 0.5, hub_height + 0.5], 'k--')
            axs[0,1].set(xlim=[thetal, thetar])
            axs[0,1].set(xlabel='$\Theta(K)$')
        else:
            if ABL == "CBL":
                uwl=-0.5
                TIr=0.2
                TIl=0.0
                
            else:
                uwl=-0.3
                TIr=0.20
                TIl=0.0
            axs[0,1].plot((tavg.txzMean + tavg.uwMean) * dim.u_star**2,
                      p.Z1_uv[0:zmax], color, label=case)
            axs[0,1].plot([uwl, 0], [hub_height - 0.5, hub_height - 0.5], 'k--')
            axs[0,1].plot([uwl, 0], [hub_height, hub_height], 'k--')
            axs[0,1].plot([uwl, 0], [hub_height + 0.5, hub_height + 0.5], 'k--')
            axs[0,1].set(xlim=[uwl, 0])
            axs[0,1].set(xlabel='kinematic shear stress($m^2s^{-2}$)')
        
        axs[0,2].plot(tavg.TIMean / dim.U_hub, p.Z1_uv[0:zmax], color, label=case)
        axs[0,2].plot([TIl, TIr], [hub_height - 0.5, hub_height - 0.5], 'k--')
        axs[0,2].plot([TIl, TIr], [hub_height, hub_height], 'k--')
        axs[0,2].plot([TIl, TIr], [hub_height + 0.5, hub_height + 0.5], 'k--')
        axs[0,2].set(xlim=[TIl, TIr])
        
        axs[1,0].plot(np.sqrt(tavg.uuMean) / dim.U_hub, p.Z1_uv[0:zmax], color, label=case)
        axs[1,1].plot(np.sqrt(tavg.vvMean) / dim.U_hub, p.Z1_uv[0:zmax], color, label=case)
        axs[1,2].plot(np.sqrt(tavg.wwMean) / dim.U_hub, p.Z1_w[0:zmax], color, label=case)
        for j in range(3):
            axs[1,j].plot([TIl, TIr], [hub_height - 0.5, hub_height - 0.5], 'k--')
            axs[1,j].plot([TIl, TIr], [hub_height, hub_height], 'k--')
            axs[1,j].plot([TIl, TIr], [hub_height + 0.5, hub_height + 0.5], 'k--')
            axs[1,j].set(xlim=[TIl,TIr])
        for i in range(2):
            for j in range(3):
                axs[i,j].set(ylim=ylim)
        axs[0,0].set(xlabel='$U(m/s)$',ylabel="$z/D$")
        
        axs[0,2].set(xlabel=r'$TI$')
        axs[1,0].set(xlabel='$I_u$',ylabel="$z/D$")
        axs[1,1].set(xlabel='$I_v$')
        axs[1,2].set(xlabel='$I_w$')
        # axs[1,0].set(xlabel='Streamwise turbulence intensity'+r'$I_u$',ylabel=r"$z/D$")
        # axs[1,1].set(xlabel='Spanwise turbulence intensity'+r'$I_v$')
        # axs[1,2].set(xlabel='Vertical turbulence intensity'+r'$I_w$')

        
    def Profile(self, ax, p, data, loc='', label='', linestyle='-', marker=None):
        if loc == 'x':
            if marker is not None:
                ax.plot(data, p.X1, linestyle, marker=marker,markersize=4, label=label)
            else:
                ax.plot(data, p.X1, linestyle, label=label)
            #ax.legend()
        elif loc == 'y':
            if marker is not None:
                ax.plot(data, p.Y1, linestyle, marker=marker,markersize=4, label=label)
            else:
                ax.plot(data, p.Y1, linestyle, label=label)
            #ax.legend()
        elif loc == 'z':
            if marker is not None:
                ax.plot(data, p.Z1_uv[0:len(data)], linestyle, marker=marker,markersize=4,label=label)
            else:
                ax.plot(data, p.Z1_uv[0:len(data)], linestyle, label=label)
            #ax.legend()
        else:
            print("Please input 'x', 'y' or 'z'!")

    def Profiles(self, ax, p, data, ind, loc='', xloc=[], turb_geo=[],scatter=True):
        """
        ax is the representation of the fig
        
        p is the instance of Params
        
        data is what you are willing to draw
        
        ind is one index corresponding to the loc
        
        loc is a character which represents the plane of the profiles
        
        xloc is a list which consists of the positions you are interested in
        """
        #TODO找最近的索引值，避免插值
        zmax=np.shape(data)[2]
        marker=['o','s','v','3','*','h','+','D',]
        if loc == 'y':
            points=np.zeros((len(xloc),p.ny,3))
            for m in range(len(xloc)):
                points[m,:,0]=xloc[m]
            for n in range(p.ny):
                points[:,n,1]=p.y[n]
            points[:,:,2]=ind
            u=interpn([p.x,p.y,p.z_uv[0:zmax]],data,points,bounds_error=False)
            for i in range(len(xloc)):
                label = "H:x="+str(int(round((xloc[i]-turb_geo[0])/turb_geo[1]))) + 'D'
                if scatter==True:
                    self.Profile(ax, p, u[i, :], loc=loc, 
                                 label=label,marker=marker[i])
                else:
                    self.Profile(ax, p, u[i, :], loc=loc, 
                                 label=label)
        if loc == 'z':
            points=np.zeros((len(xloc),zmax,3))
            for m in range(len(xloc)):
                points[m,:,0]=xloc[m]
            points[:,:,1]=ind
            for n in range(zmax):
                points[:,n,2]=p.z_uv[n]
            u=interpn([p.x,p.y,p.z_uv[0:zmax]],data,points,bounds_error=False)
            for i in range(len(xloc)):
                label = "V:x="+str(int(round((xloc[i]-turb_geo[0])/turb_geo[1]))) + 'D'
                if scatter==True:
                    self.Profile(ax, p, u[i, :], loc=loc, 
                                 label=label, linestyle='--',marker=marker[i])
                else:
                    self.Profile(ax, p, u[i, :], loc=loc, 
                                 label=label,linestyle='--')

    def VelSelfSim(self,
                   ax,
                   p,
                   Nvd,
                   xloc=[],
                   ylim=[],
                   zlim=[],
                   line=False,
                   yloc=None,
                   zloc=None,
                   half=True):
        zmax=np.shape(Nvd)[2]
        marker=['o','s','v','3','*','h','+','D',]
        # index = [0] * len(xloc)
        # for i in range(len(xloc)):
        #     index[i] = (np.abs(xloc[i] - p.x)).argmin()
        Nvdmax = np.zeros(len(xloc))
        rw = np.zeros(len(xloc))
        if half:
            xs = np.arange(0, 3, 0.01)
            ys = np.zeros(xs.size)
            for i in range(xs.size):
                ys[i] = 1 / 2 * np.exp(-np.log(2) * (xs[i]**2 - 1))
            ax.plot(xs, ys, 'k', label='Gaussian')
            if yloc is not None:
                labelv = ['Gaussian'] + [
                    'V:x=' + str(int(round((xloc[i]-0.4)/0.08))) + 'D'
                    for i in range(len(xloc))
                ]
                points=np.zeros((len(xloc),p.nz_tot,3))
                for m in range(len(xloc)):
                    points[m,:,0]=xloc[m]
                points[:,:,1]=yloc
                for n in range(zmax):
                    points[:,n,2]=p.z_uv[n]
                u=interpn([p.x,p.y,p.z_uv[0:zmax]],Nvd,points,bounds_error=False)
                for i in range(len(xloc)):
                    Nvdmax[i] = np.max(u[i, zlim[0]:zlim[1]])
                    f = interpolate.interp1d(u[i, zlim[0]:zlim[1]],
                                             p.Z1_uv[zlim[0]:zlim[1]])
                    rw[i] = abs(f(Nvdmax[i] / 2))
                    ax.plot(abs(p.Z1_uv[zlim[0]:zlim[1]] - 0.875) /
                            (rw[i] - 0.875),
                            u[i, zlim[0]:zlim[1]] / Nvdmax[i],
                            marker=marker[i],markersize=4,
                            label=labelv[i + 1])
                #ax.legend(labelv,loc="upper right")
            if zloc is not None:
                labelh = ['Gaussian'] + [
                    'H:x=' + str(int(round((xloc[i]-0.4)/0.08))) + 'D'
                    for i in range(len(xloc))
                ]
                points=np.zeros((len(xloc),p.ny,3))
                for m in range(len(xloc)):
                    points[m,:,0]=xloc[m]
                for n in range(p.ny):
                    points[:,n,1]=p.y[n]
                points[:,:,2]=zloc
                u=interpn([p.x,p.y,p.z_uv[0:zmax]],Nvd,points,bounds_error=False)
                for i in range(len(xloc)):
                    Nvdmax[i] = np.max(u[i, ylim[0]:ylim[1]])
                    f = interpolate.interp1d(u[i, ylim[0]:ylim[1]],
                                             p.Y1[ylim[0]:ylim[1]])
                    rw[i] = abs(f(Nvdmax[i] / 2))
                    ax.plot(abs(p.Y1[ylim[0]:ylim[1]]) / rw[i],
                            u[i, ylim[0]:ylim[1]] / Nvdmax[i],
                            marker=marker[i],markersize=4,
                            label=labelh[i + 1])
                #ax.legend(labelh,loc="upper right")
        else:
            xs = np.arange(-3, 3, 0.01)
            ys = np.zeros(xs.size)
            for i in range(xs.size):
                ys[i] = 1 / 2 * np.exp(-np.log(2) * (xs[i]**2 - 1))
            ax.plot(xs, ys, 'k', label='Gaussian')
            if yloc is not None:
                labelv = ['Gaussian'] + [
                    'V:x=' + str(int(round((xloc[i]-0.4)/0.08))) + 'D'
                    for i in range(len(xloc))
                ]
                points=np.zeros((len(xloc),zmax,3))
                for m in range(len(xloc)):
                    points[m,:,0]=xloc[m]
                points[:,:,1]=yloc
                for n in range(zmax):
                    points[:,n,2]=p.z_uv[n]
                u=interpn([p.x,p.y,p.z_uv[0:zmax]],Nvd,points,bounds_error=False)
                for i in range(len(xloc)):
                    a = u[i, zlim[0]:zlim[1]].reshape(zlim[1]-zlim[0])
                    Nvdmax[i] = np.max(a)
                    z = np.where(a == np.max(a))
                    print(z)
                    f = interpolate.interp1d(u[i,z[0][0]:zlim[1]],
                                             p.Z1_uv[z[0][0]:zlim[1]])
                    rw[i] = abs(f(Nvdmax[i] / 2))
                    print(p.Z1_uv[z])
                    if line:
                        ax.plot((p.Z1_uv[zlim[0]:zlim[1]] - p.Z1_uv[z])/
                                (rw[i] - p.Z1_uv[z]),
                                u[i, zlim[0]:zlim[1]] / Nvdmax[i],
                                '--',
                                marker=marker[i],markersize=4,
                                label=labelv[i + 1])
                    else:
                        ax.scatter((p.Z1_uv[zlim[0]:zlim[1]] - p.Z1_uv[z])/
                                   (rw[i] - p.Z1_uv[z]),
                                   u[i, zlim[0]:zlim[1]] / Nvdmax[i],
                                   marker=marker[i],
                                   label=labelv[i + 1])
                print(rw)
                #ax.legend(labelv,loc="upper right")
            if zloc is not None:
                labelh = ['Gaussian'] + [
                    'H:x=' + str(int(round((xloc[i]-0.4)/0.08))) + 'D'
                    for i in range(len(xloc))
                ]
                points=np.zeros((len(xloc),p.ny,3))
                for m in range(len(xloc)):
                    points[m,:,0]=xloc[m]
                for n in range(p.ny):
                    points[:,n,1]=p.y[n]
                points[:,:,2]=zloc
                u=interpn([p.x,p.y,p.z_uv[0:zmax]],Nvd,points,bounds_error=False)
                for i in range(len(xloc)):
                    a = u[i, ylim[0]:ylim[1]].reshape(ylim[1]-ylim[0])
                    Nvdmax[i] = np.max(a)
                    y = np.where(a == np.max(a))
                    print(y)
                    f = interpolate.interp1d(u[i, 48:ylim[1]],
                                             p.Y1[48:ylim[1]])
                    rw[i] = abs(f(Nvdmax[i] / 2))
                    if line:
                        ax.plot((p.Y1[ylim[0]:ylim[1]] - p.Y1[y[0][0]]) /
                                (rw[i] - p.Y1[y[0][0]]),
                                a / Nvdmax[i],
                                marker=marker[i],markersize=4,
                                label=labelh[i + 1])
                    else:
                        ax.scatter((p.Y1[ylim[0]:ylim[1]] - p.Y1[y[0][0]]) /
                                   (rw[i] - p.Y1[y[0][0]]),
                                   a / Nvdmax[i],
                                   marker=marker[i],
                                   label=labelh[i + 1])
                print(rw)
                #ax.legend(labelh,loc="upper right")
                
