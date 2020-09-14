# Copyright 2019 NREL

# Licensed under the Apache License, Version 2.0 (the "License"); you may not use
# this file except in compliance with the License. You may obtain a copy of the
# License at http://www.apache.org/licenses/LICENSE-2.0

# Unless required by applicable law or agreed to in writing, software distributed
# under the License is distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR
# CONDITIONS OF ANY KIND, either express or implied. See the License for the
# specific language governing permissions and limitations under the License.

import numpy as np
from scipy import interpolate
from .LoadLESOutput import Snap
from .PartialInterp import interp_w_uv,partialz_uv_w


class ABLProperties:
    """
    A class consists of useful data for precursor domain comparison
    
    Hhub is nondimensional length using z_i as characteristic length
    
    Uhub is the dimensional velocity
    
    ABL is a string which imply the stability of ABL
    
    U_hub,u_star,<U>,<w'\theta'>,<\theta>,<u'w'>,<v'w'>,TIst,TIsp,TIv
    """
    def __init__(self):
        pass
    def extract(self, p, Tavg, Hhub, Uhub, ABL):
        self.p = p
        zmax=np.shape(Tavg.u)[2]
        f = interpolate.interp1d(p.z_uv[0:zmax], Tavg.uMean)
        self.U_hub = f(Hhub)
        ustar = np.arange(0, 1, 0.001)
        error = np.zeros_like(ustar)
        for i in range(len(ustar)):
            error[i] = Uhub - self.U_hub * ustar[i]
        temp = np.min(abs(error))
        error = list(abs(error))
        self.u_star = ustar[error.index(temp)]
        
        self.u = Tavg.uMean * self.u_star
        self.wpthetap = (Tavg.wpthetapMean + Tavg.sgst3Mean) * self.u_star * p.T_scale
        self.theta = Tavg.thetaMean * p.T_scale
        self.upwp = (Tavg.txzMean + Tavg.uwMean) * self.u_star**2
        self.vpwp = (Tavg.tyzMean + Tavg.vwMean) * self.u_star**2
        self.TIst = np.sqrt(Tavg.uuMean) / self.U_hub
        self.TIsp = np.sqrt(Tavg.vvMean) / self.U_hub
        self.TIv = np.sqrt(Tavg.wwMean) / self.U_hub
        self.TIv_uv = np.zeros(zmax)
        for j in range(zmax-1):
            self.TIv_uv[j] = (self.TIv[j]+self.TIv[j+1])/2
        f = interpolate.interp1d(p.z_uv[0:zmax], self.TIst)
        self.Iu=f(Hhub)
        f = interpolate.interp1d(p.z_uv[0:zmax], self.TIsp)
        self.Iv=f(Hhub)
        f = interpolate.interp1d(p.z_uv[0:zmax], self.TIv_uv)
        self.Iw=f(Hhub)
        f = interpolate.interp1d(p.z_uv[0:zmax], np.sqrt(1/3*(self.TIst**2+self.TIsp**2+self.TIv_uv**2)))
        self.TI = f(Hhub)
        
        L = np.zeros((p.nx,p.ny))
        if Tavg.sgst3[:,:,0].any()!=0:
            L = -self.u_star**2/(p.vonk*p.g*Tavg.sgst3[:,:,0])
        self.L = np.mean(L)
        dudz_w = partialz_uv_w(p, Tavg.u)  # w grid
        # use wall model to calculate dudz, dvdz
        ustar = np.zeros((p.nx, p.ny))
        u_avg = np.zeros((p.nx, p.ny))
        vonk = p.vonk
        z0 = p.zo
        k = 0  # wall model
        Psi = np.zeros((p.nx,p.ny))
        Phi = np.ones((p.nx,p.ny))
        L = np.zeros((p.nx,p.ny))
        ghat = 9.81 * p.z_i / p.u_star**2
        if Tavg.sgst3[:,:,0].any() != 0:
            L = Tavg.theta[:,:,0]/(vonk*ghat*Tavg.sgst3[:,:,0])
        if ABL == "CBL":
            x = np.power(1-16*0.5*p.dz/L,0.25)
            Psi = 2*np.log(1/2*(1+x))+np.log(1/2*(1+x**2))-2*np.arctan(x)+np.pi/2
            Phi = x**(-1)
        elif ABL == "SBL":
            Psi = -5*0.5*p.dz/L
            Phi = 1 + 5*0.5*p.dz/L
        # theoretical velocity in the first uv grid
        demo = np.log(0.5 * p.dz / z0) - Psi
        u_avg[:, :] = np.sqrt(Tavg.u[:, :, k]**2 + Tavg.v[:, :, k]**2)
        ustar[:, :] = u_avg[:, :] * vonk / demo
        # w grid
        dudz_w[:, :, k] = ustar[:, :] / (
        0.5 * p.dz * vonk) * Tavg.u[:, :, k] / u_avg[:, :]*Phi
        # uv grid
        dudz = interp_w_uv(p, dudz_w)
        self.Rif=np.zeros(zmax)
        for i in range(zmax):
            self.Rif[i] = p.g*self.wpthetap[i]/self.theta[i]/\
                (self.upwp[i]*np.mean(dudz[:,:,i]*self.u_star/p.z_i))
        self.Rit=np.zeros(zmax)
        dthetadz_w = partialz_uv_w(p, Tavg.theta)
        dthetadz = interp_w_uv(p, dthetadz_w)
        for i in range(zmax):
            self.Rit[i] = p.g/self.theta[i]*np.mean(dthetadz[:,:,i])*p.T_scale/p.z_i/\
                (np.mean(dudz[:,:,i])*self.u_star/p.z_i)**2
        self.Ribb=p.g/self.theta[8]*(self.theta[8]-self.theta[3])*(p.z_uv[8]-p.z_uv[3])*p.z_i/\
            (self.u[8]-self.u[3])**2
        self.Ribt=p.g/self.theta[8]*(self.theta[13]-self.theta[8])*(p.z_uv[13]-p.z_uv[8])*p.z_i/\
            (self.u[13]-self.u[8])**2
        

    def calnvd(self, blue, red):
            """
            blue is time-averaged streamwise velocity in blue domain
            
            red is horizonal and time-averaged streamwise velocity in red domain
            
            how to define nvd
            """
            temp = np.zeros_like(blue)
            for k in range(np.shape(blue)[2]):
                #temp[:, :, k] = red[k]
                temp[:, :, k] = np.mean(blue[0:8, :, k])
            self.nvd = (temp - blue) / self.U_hub
        
def ToVTK(p, u, v, w, filename):
    """
    Save FlowData Object to vtk
    u,v,w are three dimensional flow feild

    Args:
        filename (str): Write-to path for vtk file
    """
    u = np.reshape(u, p.nx * p.ny * p.nz_tot)
    v = np.reshape(v, p.nx * p.ny * p.nz_tot)
    w = np.reshape(w, p.nx * p.ny * p.nz_tot)
    velfield = Vec3(u, v, w)

    n_points = p.nx * p.ny * p.nz_tot
    vtk_file = Output(filename)
    vtk_file.write_line('# vtk DataFile Version 3.0')
    vtk_file.write_line('array.VInstant')
    vtk_file.write_line('ASCII')
    vtk_file.write_line('DATASET STRUCTURED_POINTS')
    vtk_file.write_line('DIMENSIONS {}'.format(
        Vec3(p.nx, p.ny, p.nz_tot)))
    vtk_file.write_line('ORIGIN {}'.format(Vec3(0, 0, 0)))
    vtk_file.write_line('SPACING {}'.format(
        Vec3(p.dx, p.dy, p.dz)))
    vtk_file.write_line('POINT_DATA {}'.format(n_points))
    vtk_file.write_line('FIELD attributes 1')
    vtk_file.write_line('VInstant 3 {} float'.format(n_points))
    for u, v, w in zip(velfield.x1, velfield.x2, velfield.x3):
        vtk_file.write_line('{}'.format(Vec3(u, v, w)))

def WakeMeandering(p, coord, filepath, red, ABL, xloc=[], yr=[], zr=[]):
    """
    this function use method presented in "Inﬂuence of the Coriolis force 
    on the structure and evolution of wind turbine wakes" to quantify the
    intensity in different directions 

    red is horzional and time-averaged velocity of red domain
    
    xloc is the sequene of x locations(lesgo coordinate)
    
    yr, zr are ranges which are in nondimensional forms(lesgo 
    coordinate)
    """
    snap_time = np.arange(p.domain_nstart, p.domain_nend, p.domain_nskip)
    temp = [0] * len(p.x)
    index = [0] * len(xloc)
    for i in range(len(xloc)):
        for j in range(len(p.x)):
            temp[j] = abs(xloc[i] - p.x[j])
            for j in range(len(p.x)):
                if temp[j] == np.min(temp):
                    index[i] = j
    ycstd, zcstd = [np.zeros(len(xloc)), np.zeros(len(xloc))]

    zc = np.zeros((len(xloc), len(snap_time)))
    yc = np.zeros((len(xloc), len(snap_time)))

    u_in = ABL.U_hub * np.ones((p.ny, p.nz_tot))
    for k in range(p.nz_tot):
        u_in[:, k] = red[k]

    for i in range(len(snap_time)):
        snaptime = snap_time[i]
        c = Snap(p, coord, snaptime, filepath, Scalar=False)
        for j in range(len(xloc)):
            utemp = c.ui[index[j], :, :]
            yc[j][i] = wakecenteryc(p, utemp - u_in, yr, zr)
            zc[j][i] = wakecenterzc(p, utemp - u_in, yr, zr)

    for i in range(len(xloc)):
        ycstd[i] = np.std(yc[i, :])
        zcstd[i] = np.std(zc[i, :])
    return ycstd, zcstd

def wakecenteryc(p, x, yr, zr):
    """
    calculate y coordinate of the wake center defined as the momentum 
    deficit center
    
    x should be an instance of NVD, Y2 is p.X22
    
    yr, zr are ranges of coordinates(lesgo coordinate), respectively
    """
    ys = int(yr[0] / p.dy)
    ye = int(yr[1] / p.dy)
    zs = int(zr[0] / p.dz)
    ze = int(zr[1] / p.dz)
    no = np.sum(np.power(x[ys:ye, zs:ze], 2) * p.X22[ys:ye, zs:ze])
    deno = np.sum(np.power(x[ys:ye, zs:ze], 2))
    yc = no / deno
    return yc

def wakecenterzc(p, x, yr, zr):
    """
    calculate z coordinate of the wake center defined as the momentum
    deficit center
    
    x should be an instance of NVD, Z2_uv is p.Y22_uv
    
    yr,zr are ranges of coordinates(lesgo coordinate), repectively
    """
    ys = int(yr[0] / p.dy)
    ye = int(yr[1] / p.dy)
    zs = int(zr[0] / p.dz)
    ze = int(zr[1] / p.dz)
    no = np.sum(np.power(x[ys:ye, zs:ze], 2) * p.Y22_uv[ys:ye, zs:ze])
    deno = np.sum(np.power(x[ys:ye, zs:ze], 2))
    zc = no / deno
    return zc   
