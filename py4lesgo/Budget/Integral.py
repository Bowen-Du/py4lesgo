import numpy as np


class Integral:
    """
    this moudle is based on "Effects of a three-dimensional hill on the wake 
    characteristics of a model wind turbine" and "Distribution of mean kinetic 
    energy aint an isolated wind turbine and a characteristic wind turbine of
    a very large wind farm"
    """
    def __init__(self, p, xr=[], yr=[], zr=[]):
        xs = int(xr[0] / p.dx)
        xe = int(xr[1] / p.dx) +1
        ys = int(yr[0] / p.dy)
        ye = int(yr[1] / p.dy) +1
        zs = int(zr[0] / p.dz)
        ze = int(zr[1] / p.dz) +1
        #alphayz
        self.alphayz=np.zeros((ye-ys,ze-zs))
        for i in range(1,ye-ys-1):
            for j in range(1,ze-zs-1):
                self.alphayz[i,j]=1
        i=0
        j=0
        self.alphayz[i,j]=(yr[0]-(ys+1)*p.dy)*(zr[0]-(zs+1)*p.dz)/p.dy/p.dz
        j=0
        for i in range(1,ye-ys-1):
            self.alphayz[i,j]=-(zr[0]-(zs+1)*p.dz)/p.dz
        i=ye-ys-1
        j=0
        self.alphayz[i,j]=(yr[1]-(ye-1)*p.dy)*-(zr[0]-(zs+1)*p.dz)/p.dy/p.dz
        i=ye-ys-1
        for j in range(1,ze-zs-1):
            self.alphayz[i,j]=(yr[1]-(ye-1)*p.dy)/p.dy
        i=ye-ys-1
        j=ze-zs-1
        self.alphayz[i,j]=(yr[1]-(ye-1)*p.dy)*(zr[1]-(ze-1)*p.dz)/p.dy/p.dz
        j=ze-zs-1
        for i in range(1,ye-ys-1):
            self.alphayz[i,j]=(zr[1]-(ze-1)*p.dz)/p.dz
        i=0
        j=ze-zs-1
        self.alphayz[i,j]=-(yr[0]-(ys+1)*p.dy)*(zr[1]-(ze-1)*p.dz)/p.dy/p.dz
        i=0
        for j in range(1,ze-zs-1):
            self.alphayz[i,j]=-(yr[0]-(ys+1)*p.dy)/p.dy
        #alphaxz
        self.alphaxz=np.zeros((xe-xs,ze-zs))
        for i in range(1,xe-xs-1):
            for j in range(1,ze-zs-1):
                self.alphaxz[i,j]=1
        i=0
        j=0
        self.alphaxz[i,j]=(xr[0]-(xs+1)*p.dx)*(zr[0]-(zs+1)*p.dz)/p.dx/p.dz
        j=0
        for i in range(1,xe-xs-1):
            self.alphaxz[i,j]=-(zr[0]-(zs+1)*p.dz)/p.dz
        i=xe-xs-1
        j=0
        self.alphaxz[i,j]=(xr[1]-(xe-1)*p.dx)*-(zr[0]-(zs+1)*p.dz)/p.dx/p.dz
        i=xe-xs-1
        for j in range(1,ze-zs-1):
            self.alphaxz[i,j]=(xr[1]-(xe-1)*p.dx)/p.dx
        i=xe-xs-1
        j=ze-zs-1
        self.alphaxz[i,j]=(xr[1]-(xe-1)*p.dx)*(zr[1]-(ze-1)*p.dz)/p.dx/p.dz
        j=ze-zs-1
        for i in range(1,xe-xs-1):
            self.alphaxz[i,j]=(zr[1]-(ze-1)*p.dz)/p.dz
        i=0
        j=ze-zs-1
        self.alphaxz[i,j]=-(xr[0]-(xs+1)*p.dx)*(zr[1]-(ze-1)*p.dz)/p.dx/p.dz
        i=0
        for j in range(1,ze-zs-1):
            self.alphaxz[i,j]=-(xr[0]-(xs+1)*p.dx)/p.dx
        #alphaxy
        self.alphaxy=np.zeros((xe-xs,ye-ys))
        for i in range(1,xe-xs-1):
            for j in range(1,ye-ys-1):
                self.alphaxy[i,j]=1
        i=0
        j=0
        self.alphaxy[i,j]=(xr[0]-(xs+1)*p.dx)*(yr[0]-(ys+1)*p.dy)/p.dx/p.dy
        j=0
        for i in range(1,xe-xs-1):
            self.alphaxy[i,j]=-(yr[0]-(ys+1)*p.dy)/p.dy
        i=xe-xs-1
        j=0
        self.alphaxy[i,j]=(xr[1]-(xe-1)*p.dx)*-(yr[0]-(ys+1)*p.dy)/p.dx/p.dy
        i=xe-xs-1
        for j in range(1,ye-ys-1):
            self.alphaxy[i,j]=(xr[1]-(xe-1)*p.dx)/p.dx
        i=xe-xs-1
        j=ye-ys-1
        self.alphaxy[i,j]=(xr[1]-(xe-1)*p.dx)*(yr[1]-(ye-1)*p.dy)/p.dx/p.dy
        j=ye-ys-1
        for i in range(1,xe-xs-1):
            self.alphaxy[i,j]=(yr[1]-(ye-1)*p.dy)/p.dy
        i=0
        j=ye-ys-1
        self.alphaxy[i,j]=-(xr[0]-(xs+1)*p.dx)*(yr[1]-(ye-1)*p.dy)/p.dx/p.dy
        i=0
        for j in range(1,ye-ys-1):
            self.alphaxy[i,j]=-(xr[0]-(xs+1)*p.dx)/p.dx
        #alpha
        self.alpha=np.zeros((xe-xs,ye-ys,ze-zs))
        for i,j,k in zip(range(1,xe-xs-1),range(1,ye-ys-1),range(1,ze-zs-1)):
            self.alpha[i,j,k]=1
        self.alpha[0,:,:]=((xs+1)*p.dx-xr[0])*self.alphayz/p.dx
        self.alpha[xe-xs-1,:,:]=(xr[1]-(xe-1)*p.dx)*self.alphayz/p.dx
        self.alpha[1:xe-xs-1,0,:]=((ys+1)*p.dy-yr[0])*self.alphaxz[1:xe-xs-1,:]/p.dy
        self.alpha[1:xe-xs-1,ye-ys-1,:]=(yr[1]-(ye-1)*p.dx)*self.alphaxz[1:xe-xs-1,:]/p.dy
        self.alpha[1:xe-xs-1,1:ye-ys-1,1]=((zs+1)*p.dz-zr[0])*self.alphaxy[1:xe-xs-1,1:ye-ys-1]/p.dz
        self.alpha[1:xe-xs-1,1:ye-ys-1,ze-zs-1]=(zr[1]-(ze-1)*p.dz)*self.alphaxy[1:xe-xs-1,1:ye-ys-1]/p.dz
        #TODO 在__init__中设置转化条件通过整合以减少重复代码，FLORIS.Tools
   
    def tripleintegral(self, p, MKE, xr=[], yr=[], zr=[]):
        """
        xr, yr, zr are nondimensional integral range
        
        """
        xs = int(xr[0] / p.dx)
        xe = int(xr[1] / p.dx) +1
        ys = int(yr[0] / p.dy)
        ye = int(yr[1] / p.dy) +1
        zs = int(zr[0] / p.dz)
        ze = int(zr[1] / p.dz) +1
        
        
        self.MC = np.sum(MKE.mc[xs:xe, ys:ye, zs:ze]*self.alpha * p.dx * p.dy * p.dz)
        self.PT = np.sum(MKE.pt[xs:xe, ys:ye, zs:ze]*self.alpha * p.dx * p.dy * p.dz)
        self.TC = np.sum(MKE.tc[xs:xe, ys:ye, zs:ze]*self.alpha * p.dx * p.dy * p.dz)
        self.DF = np.sum(MKE.df[xs:xe, ys:ye, zs:ze]*self.alpha * p.dx * p.dy * p.dz)
        self.TP = np.sum(MKE.tp[xs:xe, ys:ye, zs:ze]*self.alpha * p.dx * p.dy * p.dz)
        self.DP = np.sum(MKE.dp[xs:xe, ys:ye, zs:ze]*self.alpha * p.dx * p.dy * p.dz)
        if MKE.wt is not None:
            self.WT = np.sum(MKE.wt[xs:xe, ys:ye, zs:ze]*self.alpha * p.dx * p.dy * p.dz)
        if MKE.ga is not None:
            self.GA = np.sum(MKE.ga[xs:xe, ys:ye, zs:ze]*self.alpha * p.dx * p.dy * p.dz)
        if MKE.fp is not None:
            self.FP = np.sum(MKE.fp[xs:xe, ys:ye, zs:ze]*self.alpha * p.dx * p.dy * p.dz)

    def doubleintegralx(self, p, MKE, x, yr=[], zr=[], TC=False, MC=True):
        """
        yr, zr are nondimensional integral range, and x is nondimensional 
        streamwise coordinate, this function can be used to interpret the streamwise
        development of wind turbine's wake region, TC(turbulence convection) 
        includes the different contribution by three directions
        """
        i = round(x / p.dx)
        ys = int(yr[0] / p.dy)
        ye = int(yr[1] / p.dy)+1
        zs = int(zr[0] / p.dz)
        ze = int(zr[1] / p.dz)+1
        
        self.MC_x = (np.sum(MKE.mc[i, ys:ye, zs:ze]*self.alphayz)
                     +np.sum(MKE.mc[i+1, ys:ye, zs:ze]*self.alphayz))/2 * p.dy * p.dz
        self.PT_x = (np.sum(MKE.pt[i, ys:ye, zs:ze]*self.alphayz)
                     +np.sum(MKE.pt[i+1, ys:ye, zs:ze]*self.alphayz))/2 * p.dy * p.dz
        self.TC_x = (np.sum(MKE.tc[i, ys:ye, zs:ze]*self.alphayz)
                     +np.sum(MKE.tc[i+1, ys:ye, zs:ze]*self.alphayz))/2 * p.dy * p.dz
        self.DF_x = (np.sum(MKE.df[i, ys:ye, zs:ze]*self.alphayz)
                     +np.sum(MKE.df[i+1, ys:ye, zs:ze]*self.alphayz))/2 * p.dy * p.dz
        self.TP_x = (np.sum(MKE.tp[i, ys:ye, zs:ze]*self.alphayz)
                     +np.sum(MKE.tp[i+1, ys:ye, zs:ze]*self.alphayz))/2 * p.dy * p.dz
        self.DP_x = (np.sum(MKE.dp[i, ys:ye, zs:ze]*self.alphayz)
                     +np.sum(MKE.dp[i+1, ys:ye, zs:ze]*self.alphayz))/2 * p.dy * p.dz
        if MKE.wt is not None:
            self.WT_x = (np.sum(MKE.wt[i, ys:ye, zs:ze]*self.alphayz)
                         +np.sum(MKE.wt[i+1, ys:ye, zs:ze]*self.alphayz))/2 * p.dy * p.dz
        if MKE.ga is not None:
            self.GA_x = (np.sum(MKE.ga[i, ys:ye, zs:ze]*self.alphayz)
                         +np.sum(MKE.ga[i+1, ys:ye, zs:ze]*self.alphayz))/2 * p.dy * p.dz
        if MKE.fp is not None:
            self.FP_x = (np.sum(MKE.fp[i, ys:ye, zs:ze]*self.alphayz)
                         +np.sum(MKE.fp[i+1, ys:ye, zs:ze]*self.alphayz))/2 * p.dy * p.dz
        if TC:
            self.TCx = (np.sum(MKE.tcx[i, ys:ye, zs:ze]*self.alphayz)
                        +np.sum(MKE.tcx[i+1, ys:ye, zs:ze]*self.alphayz))/2 * p.dy * p.dz
            self.TCy = (np.sum(MKE.tcy[i, ys:ye, zs:ze]*self.alphayz)
                        +np.sum(MKE.tcy[i+1, ys:ye, zs:ze]*self.alphayz))/2 * p.dy * p.dz
            self.TCz = (np.sum(MKE.tcz[i, ys:ye, zs:ze]*self.alphayz)
                        +np.sum(MKE.tcz[i+1, ys:ye, zs:ze]*self.alphayz))/2 * p.dy * p.dz
            self.TCy3 = (np.sum(MKE.tcf34[i, ys, zs:ze]*self.alphayz[1,:])
                         +np.sum(MKE.tcf34[i+1, ys, zs:ze]*self.alphayz[1,:]))/2 * p.dz
            self.TCy4 = (np.sum(MKE.tcf34[i, ye, zs:ze]*self.alphayz[1,:])
                         +np.sum(MKE.tcf34[i+1, ye, zs:ze]*self.alphayz[1,:]))/2 * p.dz
            self.TCz5 = (np.sum(MKE.tcf56[i, ys:ye, zs]*self.alphayz[:,1])
                         +np.sum(MKE.tcf56[i+1, ys:ye, zs]*self.alphayz[:,1]))/2 * p.dy
            self.TCz6 = (np.sum(MKE.tcf56[i, ys:ye, ze]*self.alphayz[:,1])
                         +np.sum(MKE.tcf56[i+1, ys:ye, ze]*self.alphayz[:,1]))/2 * p.dy
        if MC:
            self.Mc = (np.sum(MKE.mcf12[i, ys:ye, zs:ze]*self.alphayz)
                       +np.sum(MKE.mcf12[i+1, ys:ye, zs:ze]*self.alphayz))/2 * p.dy * p.dz
            self.MKE = (np.sum(MKE.u2[i, ys:ye, zs:ze]*self.alphayz)
                       +np.sum(MKE.u2[i+1, ys:ye, zs:ze]*self.alphayz))/2 * p.dy * p.dz

    def doubleintegraly(self, p, MKE, y, xr=[], zr=[]):
        """
        xr, zr are nondimensional integral range, and y is nondimensional 
        lateral coordinate
        """
        xs = int(xr[0] / p.dx)
        xe = int(xr[1] / p.dx)+1
        j = round(y / p.dy)
        zs = int(zr[0] / p.dz)
        ze = int(zr[1] / p.dz)+1
        
        self.MC_y = (np.sum(MKE.mc[xs:xe, j, zs:ze]*self.alphaxz)
                     +np.sum(MKE.mc[xs:xe, j+1, zs:ze]*self.alphaxz))/2 * p.dx * p.dz
        self.PT_y = (np.sum(MKE.pt[xs:xe, j, zs:ze]*self.alphaxz)
                     +np.sum(MKE.pt[xs:xe, j+1, zs:ze]*self.alphaxz))/2 * p.dx * p.dz
        self.TC_y = (np.sum(MKE.tc[xs:xe, j, zs:ze]*self.alphaxz)
                     +np.sum(MKE.tc[xs:xe, j+1, zs:ze]*self.alphaxz))/2 * p.dx * p.dz
        self.DF_y = (np.sum(MKE.df[xs:xe, j, zs:ze]*self.alphaxz)
                     +np.sum(MKE.df[xs:xe, j+1, zs:ze]*self.alphaxz))/2 * p.dx * p.dz
        self.TP_y = (np.sum(MKE.tp[xs:xe, j, zs:ze]*self.alphaxz)
                     +np.sum(MKE.tp[xs:xe, j+1, zs:ze]*self.alphaxz))/2 * p.dx * p.dz
        self.DP_y = (np.sum(MKE.dp[xs:xe, j, zs:ze]*self.alphaxz)
                     +np.sum(MKE.dp[xs:xe, j+1, zs:ze]*self.alphaxz))/2 * p.dx * p.dz
        if MKE.wt is not None:
            self.WT_y = (np.sum(MKE.wt[xs:xe, j, zs:ze]*self.alphaxz)
                         +np.sum(MKE.wt[xs:xe, j+1, zs:ze]*self.alphaxz))/2 * p.dx * p.dz
        if MKE.ga is not None:
            self.GA_y = (np.sum(MKE.ga[xs:xe, j, zs:ze]*self.alphaxz)
                         +np.sum(MKE.ga[xs:xe, j+1, zs:ze]*self.alphaxz))/2 * p.dx * p.dz
        if MKE.fp is not None:
            self.FP_y = (np.sum(MKE.fp[xs:xe, j, zs:ze]*self.alphaxz)
                         +np.sum(MKE.fp[xs:xe, j+1, zs:ze]*self.alphaxz))/2 * p.dx * p.dz

    def doubleintegralz(self, p, MKE, z, xr=[], yr=[]):
        """
        xr, zr are nondimensional integral range, and z is nondimensional 
        vertical coordinate
        """
        xs = int(xr[0] / p.dx)
        xe = int(xr[1] / p.dx)+1
        ys = int(yr[0] / p.dy)
        ye = int(yr[1] / p.dy)+1
        k = round(z / p.dz)
        
        self.MC_z = (np.sum(MKE.mc[xs:xe, ys:ye, k]*self.alphaxy)
                     +np.sum(MKE.mc[xs:xe, ys:ye, k+1]*self.alphaxy))/2 * p.dy * p.dx
        self.PT_z = (np.sum(MKE.pt[xs:xe, ys:ye, k]*self.alphaxy)
                     +np.sum(MKE.pt[xs:xe, ys:ye, k+1]*self.alphaxy))/2 * p.dy * p.dx
        self.TC_z = (np.sum(MKE.tc[xs:xe, ys:ye, k]*self.alphaxy)
                     +np.sum(MKE.tc[xs:xe, ys:ye, k+1]*self.alphaxy))/2 * p.dy * p.dx
        self.DF_z = (np.sum(MKE.df[xs:xe, ys:ye, k]*self.alphaxy)
                     +np.sum(MKE.df[xs:xe, ys:ye, k+1]*self.alphaxy))/2 * p.dy * p.dx
        self.TP_z = (np.sum(MKE.tp[xs:xe, ys:ye, k]*self.alphaxy)
                     +np.sum(MKE.tp[xs:xe, ys:ye, k+1]*self.alphaxy))/2 * p.dy * p.dx
        self.DP_z = (np.sum(MKE.dp[xs:xe, ys:ye, k]*self.alphaxy)
                     +np.sum(MKE.dp[xs:xe, ys:ye, k+1]*self.alphaxy))/2 * p.dy * p.dx
        if MKE.wt is not None:
            self.WT_z = (np.sum(MKE.wt[xs:xe, ys:ye, k]*self.alphaxy)
                         +np.sum(MKE.wt[xs:xe, ys:ye, k+1]*self.alphaxy))/2 * p.dy * p.dx
        if MKE.ga is not None:
            self.GA_z = (np.sum(MKE.ga[xs:xe, ys:ye, k]*self.alphaxy)
                         +np.sum(MKE.ga[xs:xe, ys:ye, k+1]*self.alphaxy))/2 * p.dy * p.dx
        if MKE.fp is not None:
            self.FP_z = (np.sum(MKE.fp[xs:xe, ys:ye, k]*self.alphaxy)
                         +np.sum(MKE.fp[xs:xe, ys:ye, k+1]*self.alphaxy))/2 * p.dy * p.dx

    def doubleintegral1(self, p, mcfx, xr=[], yr=[], zr=[]):
        """
        Quantify x-direction flux
        """
        i1 = round(xr[0] / p.dx)
        i2 = round(xr[1] / p.dx)
        ys = int(yr[0] / p.dy)
        ye = int(yr[1] / p.dy)+1
        zs = int(zr[0] / p.dz)
        ze = int(zr[1] / p.dz)+1
        
        self.fl1 = (np.sum(-mcfx[i1, ys:ye, zs:ze]*self.alphayz)
                    +np.sum(-mcfx[i1+1, ys:ye, zs:ze]*self.alphayz))/2 * p.dy * p.dz
        self.fl2 = (np.sum(mcfx[i2, ys:ye, zs:ze]*self.alphayz)
                    +np.sum(mcfx[i1+1, ys:ye, zs:ze]*self.alphayz))/2 * p.dy * p.dz

    def doubleintegral2(self, p, mcfy, xr=[], yr=[], zr=[]):
        """
        Quantify y-direction flux
        """
        xs = int(xr[0] / p.dx)
        xe = int(xr[1] / p.dx)+1
        j1 = round(yr[0] / p.dy)
        j2 = round(yr[1] / p.dy)
        zs = int(zr[0] / p.dz)
        ze = int(zr[1] / p.dz)+1
        
        self.fl3 = (np.sum(-mcfy[xs:xe, j1, zs:ze]*self.alphaxz)
                    +np.sum(-mcfy[xs:xe, j1+1, zs:ze]*self.alphaxz))/2 * p.dx * p.dz
        self.fl4 = (np.sum(mcfy[xs:xe, j2, zs:ze]*self.alphaxz)
                    +np.sum(mcfy[xs:xe, j1+1, zs:ze]*self.alphaxz))/2 * p.dx * p.dz

    def doubleintegral3(self, p, mcfz, xr=[], yr=[], zr=[]):
        """
        Quantify z-direction flux
        """
        xs = int(xr[0] / p.dx)
        xe = int(xr[1] / p.dx)+1
        ys = int(yr[0] / p.dy)
        ye = int(yr[1] / p.dy)+1
        k1 = round(zr[0] / p.dz)
        k2 = round(zr[1] / p.dz)
        
        self.fl5 = (np.sum(-mcfz[xs:xe, ys:ye, k1]*self.alphaxy)
                    +np.sum(-mcfz[xs:xe, ys:ye, k1+1]*self.alphaxy))/2* p.dy * p.dx
        self.fl6 = (np.sum(mcfz[xs:xe, ys:ye, k2]*self.alphaxy)
                    +np.sum(mcfz[xs:xe, ys:ye, k2+1]*self.alphaxy))/2* p.dy * p.dx