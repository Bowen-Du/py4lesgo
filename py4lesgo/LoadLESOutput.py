import numpy as np
from .PartialInterp import interp_w_uv,interp_uv_w

class Tavg:
    """
    A class which consists of time-averaged flow statistics
    
    p is an instance of Params
    
    coord can determine array size you postprocess to save memory size
    
    Scalar is a logical variable which implies whether the code uses scalar module
    
    filepath is the same as Params class arguments
    """
    def __init__(self, p, coord, filepath, Scalar=True, Budget=False, modify=True):
        #velocity
        self.getVelocity(p, coord, filepath)
        #square of velocity
        self.getSqOVel(p, coord, filepath)
        #Reynold stress
        self.getReyStress(p, coord, filepath)
        if modify == True:
            self.uw = self.uw2-interp_uv_w(p, self.u)*self.w
        self.ww_uv = interp_w_uv(p, self.ww)
        self.TI = np.sqrt(1 / 3 * (self.uu + self.vv + self.ww_uv))
        #sgs stress
        self.getSgsStress(p, coord, filepath)
        #square of Smagorinsky coefficient
        self.getCs_opt2(p, coord, filepath)
        #pressure*
        self.getPre(p, coord, filepath)
        #product of pressure and velocity(pressure)
        self.getPV(p, coord, filepath)
        #velocity(pressure)-pressure covariance
        self.getPpVp(p, coord, filepath)
        #build a linear-momentum-transporting velocity field
        self.buildum()
        #build the kinetic energy transport velocity field
        self.builduk()
        #wind turbine induced force
        if "blue" in filepath:
            self.getForce(p, coord, filepath)
            #self.getCt(p, filepath)
            if "simulation" in filepath:
                self.getCt(p, filepath)
        if Budget:
            self.getppSijp(p, coord, filepath)
            self.getuipujpukp(p, coord, filepath)
            self.getuptaup(p, coord, filepath)
            self.gettaupduidxjp(p, coord, filepath)
            if "blue" in filepath:
                self.getfipujp(p, coord, filepath)
        if Scalar:
            #potential temperature
            self.getScalar(p, coord, filepath)
            #unresolved heat flux
            self.getSgsTheta(p, coord, filepath)
            #product of potential temperature and velocity(temperature,pressure)
            self.getScalar2(p, coord, filepath)
            #velocity(temperature,pressure)-temperature covariance
            self.getScalarRS(p, coord, filepath)
        #Horizontal-plane mean
        self.getMean(p, coord, Scalar)
        


    def getVelocity(self, p, coord, filepath):
        self.u = np.zeros((p.nx,p.ny,p.zmax_buf[coord-1]))
        self.v = np.zeros((p.nx,p.ny,p.zmax_buf[coord-1]))
        self.w_uv = np.zeros((p.nx,p.ny,p.zmax_buf[coord-1]))
        self.w = np.zeros((p.nx,p.ny,p.zmax_buf[coord-1]))
        temp=[self.u,self.v,self.w_uv]
        for n in range(coord):
            zmin = p.zmin_buf[n]
            zmax = p.zmax_buf[n]
            filename = filepath + r"\output\veluv_avg.c" + str(n) + ".bin"
            dummy = np.fromfile(filename)
            Size = p.nx * p.ny * p.nz2
            for i in range(len(temp)):
                temp[i][:, :, zmin:zmax] = dummy[Size * i:Size * (i+1)].\
                    reshape((p.nx, p.ny, p.nz2),order="F")
            filename = filepath + r"\output\velw_avg.c" + str(n) + ".bin"
            dummy = np.fromfile(filename)
            self.w[:, :, zmin:zmax] = dummy.reshape((p.nx, p.ny, p.nz2),order="F")
        

    def getSqOVel(self, p, coord, filepath):
        self.u2 = np.zeros((p.nx,p.ny,p.zmax_buf[coord-1]))
        self.v2 = np.zeros((p.nx,p.ny,p.zmax_buf[coord-1]))
        self.w2 = np.zeros((p.nx,p.ny,p.zmax_buf[coord-1]))
        self.uw2 = np.zeros((p.nx,p.ny,p.zmax_buf[coord-1]))
        self.vw2 = np.zeros((p.nx,p.ny,p.zmax_buf[coord-1]))
        self.uv2 = np.zeros((p.nx,p.ny,p.zmax_buf[coord-1]))
        temp=[self.u2,self.v2,self.w2,self.uw2,self.vw2,self.uv2]
        for n in range(coord):
            zmin = p.zmin_buf[n]
            zmax = p.zmax_buf[n]
            filename = filepath + r"\output\vel2_avg.c" + str(n) + ".bin"
            dummy = np.fromfile(filename)
            Size = p.nx * p.ny * p.nz2
            for i in range(len(temp)):
                temp[i][:, :, zmin:zmax] = dummy[Size * i:Size * (i+1)].reshape((p.nx, p.ny, p.nz2),order="F")
        

    def getReyStress(self, p, coord, filepath):
        self.uu = np.zeros((p.nx,p.ny,p.zmax_buf[coord-1]))
        self.vv = np.zeros((p.nx,p.ny,p.zmax_buf[coord-1]))
        self.ww = np.zeros((p.nx,p.ny,p.zmax_buf[coord-1]))
        self.uw = np.zeros((p.nx,p.ny,p.zmax_buf[coord-1]))
        self.vw = np.zeros((p.nx,p.ny,p.zmax_buf[coord-1]))
        self.uv = np.zeros((p.nx,p.ny,p.zmax_buf[coord-1]))
        temp=[self.uu,self.vv,self.ww,self.uw,self.vw,self.uv]
        for n in range(coord):
            zmin = p.zmin_buf[n]
            zmax = p.zmax_buf[n]
            filename = filepath + r"\output\rs.c" + str(n) + ".bin"
            dummy = np.fromfile(filename)
            Size = p.nx * p.ny * p.nz2
            for i in range(len(temp)):
                temp[i][:, :, zmin:zmax] = dummy[Size * i:Size * (i+1)].reshape((p.nx, p.ny, p.nz2),order="F")
        

    def getSgsStress(self, p, coord, filepath):
        self.txx = np.zeros((p.nx,p.ny,p.zmax_buf[coord-1]))
        self.txy = np.zeros((p.nx,p.ny,p.zmax_buf[coord-1]))
        self.tyy = np.zeros((p.nx,p.ny,p.zmax_buf[coord-1]))
        self.txz = np.zeros((p.nx,p.ny,p.zmax_buf[coord-1]))
        self.tyz = np.zeros((p.nx,p.ny,p.zmax_buf[coord-1]))
        self.tzz = np.zeros((p.nx,p.ny,p.zmax_buf[coord-1]))
        temp=[self.txx,self.txy,self.tyy,self.txz,self.tyz,self.tzz]
        for n in range(coord):
            zmin = p.zmin_buf[n]
            zmax = p.zmax_buf[n]
            filename = filepath + r"\output\tau_avg.c" + str(n) + ".bin"
            dummy = np.fromfile(filename)
            Size = p.nx * p.ny * p.nz2
            for i in range(len(temp)):
                temp[i][:, :, zmin:zmax] = dummy[Size * i:Size * (i+1)].reshape((p.nx, p.ny, p.nz2),order="F")
        

    def getCs_opt2(self, p, coord, filepath):
        self.cs_opt2 = np.zeros((p.nx,p.ny,p.zmax_buf[coord-1]))
        for n in range(coord):
            zmin = p.zmin_buf[n]
            zmax = p.zmax_buf[n]
            filename = filepath + r"\output\cs_opt2.c" + str(n) + ".bin"
            dummy = np.fromfile(filename)
            self.cs_opt2[:, :, zmin:zmax] = dummy.reshape((p.nx, p.ny, p.nz2),order="F")
        

    def getPre(self, p, coord, filepath):
        self.pre = np.zeros((p.nx,p.ny,p.zmax_buf[coord-1]))
        for n in range(coord):
            zmin = p.zmin_buf[n]
            zmax = p.zmax_buf[n]
            #filename = filepath + r"\output\pres_avg.c" + str(n) + ".bin"
            filename = filepath + r"\output\pre_uv_avg.c" + str(n) + ".bin"
            dummy = np.fromfile(filename)
            self.pre[:, :, zmin:zmax] = dummy.reshape((p.nx, p.ny, p.nz2),order="F")
        

    def getPV(self, p, coord, filepath):
        self.p2 = np.zeros((p.nx,p.ny,p.zmax_buf[coord-1]))
        self.pu = np.zeros((p.nx,p.ny,p.zmax_buf[coord-1]))
        self.pv = np.zeros((p.nx,p.ny,p.zmax_buf[coord-1]))
        self.pw = np.zeros((p.nx,p.ny,p.zmax_buf[coord-1]))
        temp=[self.p2,self.pu,self.pv,self.pw]
        for n in range(coord):
            zmin = p.zmin_buf[n]
            zmax = p.zmax_buf[n]
            filename = filepath + r"\output\ps2_uv_avg.c" + str(n) + ".bin"
            dummy = np.fromfile(filename)
            Size = p.nx * p.ny * p.nz2
            for i in range(len(temp)):
                temp[i][:, :, zmin:zmax] = dummy[Size * i:Size * (i+1)].reshape((p.nx, p.ny, p.nz2),order="F")
        

    def getPpVp(self, p, coord, filepath):
        self.pp2 = np.zeros((p.nx,p.ny,p.zmax_buf[coord-1]))
        self.ppup = np.zeros((p.nx,p.ny,p.zmax_buf[coord-1]))
        self.ppvp = np.zeros((p.nx,p.ny,p.zmax_buf[coord-1]))
        self.ppwp = np.zeros((p.nx,p.ny,p.zmax_buf[coord-1]))
        temp=[self.pp2,self.ppup,self.ppvp,self.ppwp]
        for n in range(coord):
            zmin = p.zmin_buf[n]
            zmax = p.zmax_buf[n]
            filename = filepath + r"\output\ps_uv_avg.c" + str(n) + ".bin"
            dummy = np.fromfile(filename)
            Size = p.nx * p.ny * p.nz2
            for i in range(len(temp)):
                temp[i][:, :, zmin:zmax] = dummy[Size * i:Size * (i+1)].reshape((p.nx, p.ny, p.nz2),order="F")
        

    def getForce(self, p, coord, filepath):
        self.fx = np.zeros((p.nx,p.ny,p.zmax_buf[coord-1]))
        self.fy = np.zeros((p.nx,p.ny,p.zmax_buf[coord-1]))
        self.fz = np.zeros((p.nx,p.ny,p.zmax_buf[coord-1]))
        temp=[self.fx,self.fy,self.fz]
        for n in range(coord):
            zmin = p.zmin_buf[n]
            zmax = p.zmax_buf[n]
            filename = filepath + r"\output\force_avg.c" + str(n) + ".bin"
            dummy = np.fromfile(filename)
            Size = p.nx * p.ny * p.nz2
            for i in range(len(temp)):
                temp[i][:, :, zmin:zmax] = dummy[Size * i:Size * (i+1)].reshape((p.nx, p.ny, p.nz2),order="F")
        

    def getScalar(self, p, coord, filepath):
        self.theta = np.zeros((p.nx,p.ny,p.zmax_buf[coord-1]))
        for n in range(coord):
            zmin = p.zmin_buf[n]
            zmax = p.zmax_buf[n]
            filename = filepath + r"\output\scalar_avg.c" + str(n) + ".bin"
            dummy = np.fromfile(filename)
            self.theta[:, :, zmin:zmax] = dummy.reshape((p.nx, p.ny, p.nz2),order="F")
        

    def getScalar2(self, p, coord, filepath):
        self.theta2 = np.zeros((p.nx,p.ny,p.zmax_buf[coord-1]))
        self.utheta = np.zeros((p.nx,p.ny,p.zmax_buf[coord-1]))
        self.vtheta = np.zeros((p.nx,p.ny,p.zmax_buf[coord-1]))
        self.wtheta = np.zeros((p.nx,p.ny,p.zmax_buf[coord-1]))
        self.ptheta = np.zeros((p.nx,p.ny,p.zmax_buf[coord-1]))
        temp=[self.theta2,self.utheta,self.vtheta,self.wtheta,self.ptheta]
        for n in range(coord):
            zmin = p.zmin_buf[n]
            zmax = p.zmax_buf[n]
            filename = filepath + r"\output\scalar2_avg.c" + str(n) + ".bin"
            dummy = np.fromfile(filename)
            Size = p.nx * p.ny * p.nz2
            for i in range(len(temp)):
                temp[i][:, :, zmin:zmax] = dummy[Size * i:Size * (i+1)].reshape((p.nx, p.ny, p.nz2),order="F")
        

    def getScalarRS(self, p, coord, filepath):
        self.thetap2 = np.zeros((p.nx,p.ny,p.zmax_buf[coord-1]))
        self.upthetap = np.zeros((p.nx,p.ny,p.zmax_buf[coord-1]))
        self.vpthetap = np.zeros((p.nx,p.ny,p.zmax_buf[coord-1]))
        self.wpthetap = np.zeros((p.nx,p.ny,p.zmax_buf[coord-1]))
        self.ppthetap = np.zeros((p.nx,p.ny,p.zmax_buf[coord-1]))
        temp=[self.thetap2,self.upthetap,self.vpthetap,self.wpthetap,self.ppthetap]
        for n in range(coord):
            zmin = p.zmin_buf[n]
            zmax = p.zmax_buf[n]
            filename = filepath + r"\output\scalar_rs.c" + str(n) + ".bin"
            dummy = np.fromfile(filename)
            Size = p.nx * p.ny * p.nz2
            for i in range(len(temp)):
                temp[i][:, :, zmin:zmax] = dummy[Size * i:Size * (i+1)].reshape((p.nx, p.ny, p.nz2),order="F")
        

    def getSgsTheta(self, p, coord, filepath):
        self.sgst3 = np.zeros((p.nx,p.ny,p.zmax_buf[coord-1]))
        for n in range(coord):
            zmin = p.zmin_buf[n]
            zmax = p.zmax_buf[n]
            filename = filepath + r"\output\scalar_sgs.c" + str(n) + ".bin"
            dummy = np.fromfile(filename)
            self.sgst3[:, :, zmin:zmax] = dummy.reshape((p.nx, p.ny, p.nz2),order="F")
        
    def buildum(self):
        self.um = self.u + (self.uu-self.txx)/self.u
        self.vm = self.v + (self.uv-self.txy)/self.u
        self.wm = self.w_uv + (self.uw-self.txz)/self.u
    def builduk(self):
        K=1/2*(self.u**2+self.v**2+self.w_uv**2)
        self.uk = self.u + ((self.uu-self.txx)*self.u+(self.uv-self.txy)*self.v+
                            (self.uw-self.txz)*self.w_uv)/K
        self.vk = self.v + ((self.uv-self.txy)*self.u+(self.vv-self.tyy)*self.v+
                            (self.vw-self.tyz)*self.w_uv)/K
        self.wk = self.w_uv + ((self.uw-self.txz)*self.u+(self.vw-self.tyz)*self.v+
                            (self.ww-self.tzz)*self.w_uv)/K
    def getppSijp(self, p, coord, filepath):
        self.ppS11p = np.zeros((p.nx,p.ny,p.zmax_buf[coord-1]))
        self.ppS12p = np.zeros((p.nx,p.ny,p.zmax_buf[coord-1]))
        self.ppS13p = np.zeros((p.nx,p.ny,p.zmax_buf[coord-1]))
        self.ppS22p = np.zeros((p.nx,p.ny,p.zmax_buf[coord-1]))
        self.ppS23p = np.zeros((p.nx,p.ny,p.zmax_buf[coord-1]))
        self.ppS33p = np.zeros((p.nx,p.ny,p.zmax_buf[coord-1]))
        temp=[self.ppS11p,self.ppS12p,self.ppS13p,
                          self.ppS22p,self.ppS23p,
                                      self.ppS33p]
        for n in range(coord):
            zmin = p.zmin_buf[n]
            zmax = p.zmax_buf[n]
            filename = filepath + r"\output\ppSijp.c" + str(n) + ".bin"
            dummy = np.fromfile(filename)
            Size = p.nx * p.ny * p.nz2
            for i in range(len(temp)):
                temp[i][:, :, zmin:zmax] = dummy[Size * i:Size * (i+1)].reshape((p.nx, p.ny, p.nz2),order="F")
        
        
    def getuipujpukp(self, p, coord, filepath):
        self.upupup = np.zeros((p.nx,p.ny,p.zmax_buf[coord-1]))
        self.vpvpvp = np.zeros((p.nx,p.ny,p.zmax_buf[coord-1]))
        self.wpwpwp = np.zeros((p.nx,p.ny,p.zmax_buf[coord-1]))
        self.upupvp = np.zeros((p.nx,p.ny,p.zmax_buf[coord-1]))
        self.upupwp = np.zeros((p.nx,p.ny,p.zmax_buf[coord-1]))
        self.vpvpup = np.zeros((p.nx,p.ny,p.zmax_buf[coord-1]))
        self.vpvpwp = np.zeros((p.nx,p.ny,p.zmax_buf[coord-1]))
        self.wpwpup = np.zeros((p.nx,p.ny,p.zmax_buf[coord-1]))
        self.wpwpvp = np.zeros((p.nx,p.ny,p.zmax_buf[coord-1]))
        self.upvpwp = np.zeros((p.nx,p.ny,p.zmax_buf[coord-1]))
        temp=[self.upupup,self.vpvpvp,self.wpwpwp,
              self.upupvp,self.upupwp,self.vpvpup,self.vpvpwp,self.wpwpup,self.wpwpvp,
              self.upvpwp]
        for n in range(coord):
            zmin = p.zmin_buf[n]
            zmax = p.zmax_buf[n]
            filename = filepath + r"\output\uipujpukp.c" + str(n) + ".bin"
            dummy = np.fromfile(filename)
            Size = p.nx * p.ny * p.nz2
            for i in range(len(temp)):
                temp[i][:, :, zmin:zmax] = dummy[Size * i:Size * (i+1)].reshape((p.nx, p.ny, p.nz2),order="F")
        
        
    def getuptaup(self, p, coord, filepath):
        self.uptau11p = np.zeros((p.nx,p.ny,p.zmax_buf[coord-1]))
        self.uptau12p = np.zeros((p.nx,p.ny,p.zmax_buf[coord-1]))
        self.uptau13p = np.zeros((p.nx,p.ny,p.zmax_buf[coord-1]))
        self.uptau22p = np.zeros((p.nx,p.ny,p.zmax_buf[coord-1]))
        self.uptau23p = np.zeros((p.nx,p.ny,p.zmax_buf[coord-1]))
        self.uptau33p = np.zeros((p.nx,p.ny,p.zmax_buf[coord-1]))
        self.vptau11p = np.zeros((p.nx,p.ny,p.zmax_buf[coord-1]))
        self.vptau12p = np.zeros((p.nx,p.ny,p.zmax_buf[coord-1]))
        self.vptau13p = np.zeros((p.nx,p.ny,p.zmax_buf[coord-1]))
        self.vptau22p = np.zeros((p.nx,p.ny,p.zmax_buf[coord-1]))
        self.vptau23p = np.zeros((p.nx,p.ny,p.zmax_buf[coord-1]))
        self.vptau33p = np.zeros((p.nx,p.ny,p.zmax_buf[coord-1]))
        self.wptau11p = np.zeros((p.nx,p.ny,p.zmax_buf[coord-1]))
        self.wptau12p = np.zeros((p.nx,p.ny,p.zmax_buf[coord-1]))
        self.wptau13p = np.zeros((p.nx,p.ny,p.zmax_buf[coord-1]))
        self.wptau22p = np.zeros((p.nx,p.ny,p.zmax_buf[coord-1]))
        self.wptau23p = np.zeros((p.nx,p.ny,p.zmax_buf[coord-1]))
        self.wptau33p = np.zeros((p.nx,p.ny,p.zmax_buf[coord-1]))
        temp=[self.uptau11p,self.uptau12p,self.uptau13p,
                            self.uptau22p,self.uptau23p,
                                          self.uptau33p,
              self.vptau11p,self.vptau12p,self.vptau13p,
                            self.vptau22p,self.vptau23p,
                                          self.vptau33p,
              self.wptau11p,self.wptau12p,self.wptau13p,
                            self.wptau22p,self.wptau23p,
                                          self.wptau33p]
        for n in range(coord):
            zmin = p.zmin_buf[n]
            zmax = p.zmax_buf[n]
            filename = filepath + r"\output\uptaup.c" + str(n) + ".bin"
            dummy = np.fromfile(filename)
            Size = p.nx * p.ny * p.nz2
            for i in range(len(temp)):
                temp[i][:, :, zmin:zmax] = dummy[Size * i:Size * (i+1)].reshape((p.nx, p.ny, p.nz2),order="F")
        
    
    def gettaupduidxjp(self, p, coord, filepath):
        self.tau11pdudxp = np.zeros((p.nx,p.ny,p.zmax_buf[coord-1]))
        self.tau11pdvdxp = np.zeros((p.nx,p.ny,p.zmax_buf[coord-1]))
        self.tau11pdwdxp = np.zeros((p.nx,p.ny,p.zmax_buf[coord-1]))
        self.tau12pdudyp = np.zeros((p.nx,p.ny,p.zmax_buf[coord-1]))
        self.tau12pdvdyp = np.zeros((p.nx,p.ny,p.zmax_buf[coord-1]))
        self.tau12pdwdyp = np.zeros((p.nx,p.ny,p.zmax_buf[coord-1]))
        self.tau13pdudzp = np.zeros((p.nx,p.ny,p.zmax_buf[coord-1]))
        self.tau13pdvdzp = np.zeros((p.nx,p.ny,p.zmax_buf[coord-1]))
        self.tau13pdwdzp = np.zeros((p.nx,p.ny,p.zmax_buf[coord-1]))
        self.tau21pdudxp = np.zeros((p.nx,p.ny,p.zmax_buf[coord-1]))
        self.tau21pdvdxp = np.zeros((p.nx,p.ny,p.zmax_buf[coord-1]))
        self.tau21pdwdxp = np.zeros((p.nx,p.ny,p.zmax_buf[coord-1]))
        self.tau22pdudyp = np.zeros((p.nx,p.ny,p.zmax_buf[coord-1]))
        self.tau22pdvdyp = np.zeros((p.nx,p.ny,p.zmax_buf[coord-1]))
        self.tau22pdwdyp = np.zeros((p.nx,p.ny,p.zmax_buf[coord-1]))
        self.tau23pdudzp = np.zeros((p.nx,p.ny,p.zmax_buf[coord-1]))
        self.tau23pdvdzp = np.zeros((p.nx,p.ny,p.zmax_buf[coord-1]))
        self.tau23pdwdzp = np.zeros((p.nx,p.ny,p.zmax_buf[coord-1]))
        self.tau31pdudxp = np.zeros((p.nx,p.ny,p.zmax_buf[coord-1]))
        self.tau31pdvdxp = np.zeros((p.nx,p.ny,p.zmax_buf[coord-1]))
        self.tau31pdwdxp = np.zeros((p.nx,p.ny,p.zmax_buf[coord-1]))
        self.tau32pdudyp = np.zeros((p.nx,p.ny,p.zmax_buf[coord-1]))
        self.tau32pdvdyp = np.zeros((p.nx,p.ny,p.zmax_buf[coord-1]))
        self.tau32pdwdyp = np.zeros((p.nx,p.ny,p.zmax_buf[coord-1]))
        self.tau33pdudzp = np.zeros((p.nx,p.ny,p.zmax_buf[coord-1]))
        self.tau33pdvdzp = np.zeros((p.nx,p.ny,p.zmax_buf[coord-1]))
        self.tau33pdwdzp = np.zeros((p.nx,p.ny,p.zmax_buf[coord-1]))
        temp = [self.tau11pdudxp, self.tau11pdvdxp, self.tau11pdwdxp,
                    self.tau12pdudyp, self.tau12pdvdyp, self.tau12pdwdyp,
                    self.tau13pdudzp, self.tau13pdvdzp, self.tau13pdwdzp,
                    self.tau21pdudxp, self.tau21pdvdxp, self.tau21pdwdxp,
                    self.tau22pdudyp, self.tau22pdvdyp, self.tau22pdwdyp,
                    self.tau23pdudzp, self.tau23pdwdzp, self.tau23pdwdzp,
                    self.tau31pdudxp, self.tau31pdvdxp, self.tau31pdwdxp,
                    self.tau31pdwdxp, self.tau32pdudyp, self.tau32pdwdyp,
                    self.tau33pdudzp, self.tau33pdvdzp, self.tau33pdwdzp]
        for n in range(coord):
            zmin = p.zmin_buf[n]
            zmax = p.zmax_buf[n]
            filename = filepath + r"\output\taupduidxjp.c" + str(n) + ".bin"
            dummy = np.fromfile(filename)
            Size = p.nx * p.ny * p.nz2
            for i in range(len(temp)):
                temp[i][:, :, zmin:zmax] = dummy[Size * i:Size * (i+1)].reshape((p.nx, p.ny, p.nz2),order="F")
        
        
    def getfipujp(self, p, coord, filepath):
        self.fxpup = np.zeros((p.nx,p.ny,p.zmax_buf[coord-1]))
        self.fypup = np.zeros((p.nx,p.ny,p.zmax_buf[coord-1]))
        self.fzpup = np.zeros((p.nx,p.ny,p.zmax_buf[coord-1]))
        self.fxpvp = np.zeros((p.nx,p.ny,p.zmax_buf[coord-1]))
        self.fypvp = np.zeros((p.nx,p.ny,p.zmax_buf[coord-1]))
        self.fzpvp = np.zeros((p.nx,p.ny,p.zmax_buf[coord-1]))
        self.fxpwp = np.zeros((p.nx,p.ny,p.zmax_buf[coord-1]))
        self.fypwp = np.zeros((p.nx,p.ny,p.zmax_buf[coord-1]))
        self.fzpwp = np.zeros((p.nx,p.ny,p.zmax_buf[coord-1]))
        temp=[self.fxpup,self.fypup,self.fzpup,
              self.fxpvp,self.fypvp,self.fzpvp,
              self.fxpwp,self.fypwp,self.fzpwp]
        for n in range(coord):
            zmin = p.zmin_buf[n]
            zmax = p.zmax_buf[n]
            filename = filepath + r"\output\fipujp.c" + str(n) + ".bin"
            dummy = np.fromfile(filename)
            Size = p.nx * p.ny * p.nz2
            for i in range(len(temp)):
                temp[i][:, :, zmin:zmax] = dummy[Size * i:Size * (i+1)].reshape((p.nx, p.ny, p.nz2),order="F")
        
        
    def getMean(self, p, coord, Scalar):
        self.uMean = np.mean(self.u,axis=(0,1))
        self.vMean = np.mean(self.v,axis=(0,1))
        self.wMean = np.mean(self.w,axis=(0,1))

        self.uuMean = np.mean(self.uu,axis=(0,1))
        self.vvMean = np.mean(self.vv,axis=(0,1))
        self.wwMean = np.mean(self.ww,axis=(0,1))
        self.uwMean = np.mean(self.uw,axis=(0,1))
        self.vwMean = np.mean(self.vw,axis=(0,1))
        self.uvMean = np.mean(self.uv,axis=(0,1))
        self.TIMean = np.mean(self.TI,axis=(0,1))
        
        self.txzMean = np.mean(self.txz,axis=(0,1))
        self.tyzMean = np.mean(self.tyz,axis=(0,1))
            
        if Scalar:
            self.thetaMean = np.mean(self.theta,axis=(0,1))
            self.uthetaMean = np.mean(self.utheta,axis=(0,1))
            self.vthetaMean = np.mean(self.vtheta,axis=(0,1))
            self.wpthetapMean = np.mean(self.wpthetap,axis=(0,1))
            self.pthetaMean = np.mean(self.ptheta,axis=(0,1))
            self.sgst3Mean = np.mean(self.sgst3,axis=(0,1))
        else:
            self.thetaMean = np.ones(p.zmax_buf[coord-1])
    def getCt(self, p, filepath):
        dummy=np.loadtxt(filepath+"\\turbine\\turbine_1.dat")
        self.tdim=dummy[:,0]
        self.ucenter=dummy[:,1]
        self.vcenter=dummy[:,2]
        self.wcenter=dummy[:,3]
        self.u_d=dummy[:,4]
        self.u_d_T=dummy[:,5]
        self.theta1=dummy[:,6]
        self.theta2=dummy[:,7]
        self.Ctprime=dummy[:,8]

class Snap:
    def __init__(self, p, coord, snap_time, filepath, Scalar=False):
        self.getVelocity(p, coord, snap_time, filepath)
        if Scalar:
            self.getPre(p, coord, snap_time, filepath)
            self.getScalar(p, coord, snap_time, filepath)

    def getVelocity(self, p, coord, snap_time, filepath): 
        self.u = np.zeros((p.nx,p.ny,p.zmax_buf[coord-1]))
        self.v = np.zeros((p.nx,p.ny,p.zmax_buf[coord-1]))
        self.w_uv = np.zeros((p.nx,p.ny,p.zmax_buf[coord-1]))
        self.mag = np.zeros((p.nx, p.ny, p.zmax_buf[coord-1]))
        for n in range(coord):
            zmin = p.zmin_buf[n]
            zmax = p.zmax_buf[n]
            filename = filepath + r"\output\vel." + str(
                snap_time) + ".c" + str(n) + ".bin"
            dummy = np.fromfile(filename)
            Size = p.nx * p.ny * p.nz2
            self.u[:, :, zmin:zmax] = dummy[0:Size].reshape((p.nx, p.ny, p.nz2),order="F")
            self.v[:, :, zmin:zmax] = dummy[Size:Size * 2].reshape(
                p.nz2, p.ny, p.nx)
            self.w_uv[:, :, zmin:zmax] = dummy[Size * 2:Size * 3].reshape(
                p.nz2, p.ny, p.nx)
        self.mag = np.sqrt(self.v**2 + self.w_uv**2)

    def getCs_opt2(self, p, coord, snap_time, filepath):
        self.cs_opt2 = np.zeros((p.nx,p.ny,p.zmax_buf[coord-1]))
        for n in range(coord):
            zmin = p.zmin_buf[n]
            zmax = p.zmax_buf[n]
            filename = filepath + r"\output\cs_opt2." + str(
                snap_time) + ".c" + str(n) + ".bin"
            dummy = np.fromfile(filename)
            self.cs_opt2[:, :, zmin:zmax] = dummy.reshape((p.nx, p.ny, p.nz2),order="F")
        

    def getPre(self, p, coord, snap_time, filepath):
        self.pre = np.zeros((p.nx,p.ny,p.zmax_buf[coord-1]))
        for n in range(coord):
            zmin = p.zmin_buf[n]
            zmax = p.zmax_buf[n]
            filename = filepath + r"\output\pre." + str(
                snap_time) + ".c" + str(n) + ".bin"
            dummy = np.fromfile(filename)
            self.pre[:, :, zmin:zmax] = dummy.reshape((p.nx, p.ny, p.nz2),order="F")
        

    def getScalar(self, p, coord, snap_time, filepath):
        self.theta = np.zeros((p.nx,p.ny,p.zmax_buf[coord-1]))
        for n in range(coord):
            zmin = p.zmin_buf[n]
            zmax = p.zmax_buf[n]
            filename = filepath + r"\output\scalar." + str(
                snap_time) + ".c" + str(n) + ".bin"
            dummy = np.fromfile(filename)
            self.theta[:, :, zmin:zmax] = dummy.reshape((p.nx, p.ny, p.nz2),order="F")
        


class Snap_X:
    def __init__(self, p, coord, snap_time, filepath, xloc, Scalar=False):
        if Scalar:
            self.getVelocity(p, coord, snap_time, filepath, xloc)
            self.getPre(p, coord, snap_time, filepath, xloc)
            self.getScalar(p, coord, snap_time, filepath, xloc)
        else:
            self.getVelocity(p, coord, snap_time, filepath, xloc)
            self.getPre(p, coord, snap_time, filepath, xloc)
    def getVelocity(self, p, coord, snap_time, filepath, xloc):
        self.ux = np.zeros((p.ny, p.zmax_buf[coord-1]))
        self.vx = np.zeros((p.ny, p.zmax_buf[coord-1]))
        self.wx = np.zeros((p.ny, p.zmax_buf[coord-1]))
        for n in range(coord):
            filename = filepath + r"\output\vel.x-" + "{loc:.5f}".format(
                loc=xloc) + "." + str(snap_time) + ".c" + str(n) + ".bin"
            dummy = np.fromfile(filename)
            zmin = p.zmin_buf[n]
            zmax = p.zmax_buf[n]
            self.ux[:,zmin:zmax] = dummy[0:p.ny * p.nz2].reshape((p.ny,p.nz2),order="F")
            self.vx[:,zmin:zmax] = dummy[p.ny * p.nz2:p.ny * p.nz2 *
                                        2].reshape((p.ny,p.nz2),order="F")
            self.wx[:,zmin:zmax] = dummy[p.ny * p.nz2 * 2:p.ny * p.nz2 *
                                        3].reshape((p.ny,p.nz2),order="F")
        
    def getPre(self, p, coord, snap_time, filepath, xloc):
        self.prex = np.zeros((p.ny,p.zmax_buf[coord-1]))
        for n in range(coord):
            zmin = p.zmin_buf[n]
            zmax = p.zmax_buf[n]
            filename = filepath + r"\output\pre.x-" + "{loc:.5f}".format(
                loc=xloc) + "." + str(snap_time) + ".c" + str(n) + ".bin"
            dummy = np.fromfile(filename)
            self.prex[:, zmin:zmax] = dummy.reshape((p.ny,p.nz2),order="F")
            
    def getScalar(self, p, coord, snap_time, filepath, xloc):
        self.thetax = np.zeros((p.ny, p.zmax_buf[coord-1]))
        for n in range(coord):
            zmin = p.zmin_buf[n]
            zmax = p.zmax_buf[n]
            filename = filepath + r"\output\scalar.x-" + "{loc:.5f}".format(
                loc=xloc) + "." + str(snap_time) + ".c" + str(n) + ".bin"
            dummy = np.fromfile(filename)
            self.thetax[:, zmin:zmax] = dummy.reshape((p.ny,p.nz2),order="F")
    


class Snap_Y:
    def __init__(self, p, coord, snap_time, filepath, yloc, Scalar=False):
        if Scalar:
            self.getVelocity(p, coord, snap_time, filepath, yloc)
            self.getPre(p, coord, snap_time, filepath, yloc)
            self.getScalar(p, coord, snap_time, filepath, yloc)
        else:
            self.getVelocity(p, coord, snap_time, filepath, yloc)
            self.getPre(p, coord, snap_time, filepath, yloc)
    def getVelocity(self, p, coord, snap_time, filepath, yloc):
        self.uy = np.zeros((p.nx, p.zmax_buf[coord-1]))
        self.vy = np.zeros((p.nx, p.zmax_buf[coord-1]))
        self.wy = np.zeros((p.nx, p.zmax_buf[coord-1]))
        for n in range(coord):
            filename = filepath + r"\output\vel.y-" + "{loc:.5f}".format(
                loc=yloc) + "." + str(snap_time) + ".c" + str(n) + ".bin"
            dummy = np.fromfile(filename)
            zmin = p.zmin_buf[n]
            zmax = p.zmax_buf[n]
            self.uy[:,zmin:zmax] = dummy[0:p.nx * p.nz2].reshape((p.nx,p.nz2), order="F")
            self.vy[:,zmin:zmax] = dummy[p.nx * p.nz2:p.nx * p.nz2 *
                                        2].reshape((p.nx,p.nz2), order="F")
            self.wy[:,zmin:zmax] = dummy[p.nx * p.nz2 * 2:p.nx * p.nz2 *
                                        3].reshape((p.nx,p.nz2), order="F")
        
    def getPre(self, p, coord, snap_time, filepath, yloc):
        self.prey = np.zeros((p.nx, p.zmax_buf[coord-1]))
        for n in range(coord):
            zmin = p.zmin_buf[n]
            zmax = p.zmax_buf[n]
            filename = filepath + r"\output\pre.y-" + "{loc:.5f}".format(
                loc=yloc) + "." + str(snap_time) + ".c" + str(n) + ".bin"
            dummy = np.fromfile(filename)
            self.prey[:, zmin:zmax] = dummy.reshape((p.nx,p.nz2), order="F")
        
    def getScalar(self, p, coord, snap_time, filepath, yloc):
        self.thetay = np.zeros((p.nx, p.zmax_buf[coord-1]))
        for n in range(coord):
            zmin = p.zmin_buf[n]
            zmax = p.zmax_buf[n]
            filename = filepath + r"\output\scalar.y-" + "{loc:.5f}".format(
                loc=yloc) + "." + str(snap_time) + ".c" + str(n) + ".bin"
            dummy = np.fromfile(filename)
            self.thetay[:, zmin:zmax] = dummy.reshape((p.nx,p.nz2), order="F")



class Snap_Z:
    def __init__(self, p, snap_time, filepath, zloc, Scalar=False):
        if Scalar:
            self.getVelocity(p, snap_time, filepath, zloc)
            self.getPre(p, snap_time, filepath, zloc)
            self.getScalar(p, snap_time, filepath, zloc)
        else:
            self.getVelocity(p, snap_time, filepath, zloc)
            self.getPre(p, snap_time, filepath, zloc)
    def getVelocity(self, p, snap_time, filepath, zloc):
        self.uz = np.zeros((p.nx, p.ny))
        self.vz = np.zeros((p.nx, p.ny))
        self.wz = np.zeros((p.nx, p.ny))
        for n in range(1):
            filename = filepath + r"\output\vel.z-" + "{loc:.5f}".format(
                loc=zloc) + "." + str(snap_time) + ".c" + str(n) + ".bin"
            dummy = np.fromfile(filename)
            self.uz = dummy[0:p.nx * p.ny].reshape((p.nx, p.ny), order="F")
            self.vz = dummy[p.nx * p.ny:p.nx * p.ny * 2].reshape((p.nx, p.ny), order="F")
            self.wz = dummy[p.nx * p.ny * 2:p.nx * p.ny * 3].reshape((p.nx, p.ny), order="F")
    def getPre(self, p, snap_time, filepath, zloc):
        self.prez = np.zeros((p.nx, p.ny))
        for n in range(1):
            filename = filepath + r"\output\pre.z-" + "{loc:.5f}".format(
                loc=zloc) + "." + str(snap_time) + ".c" + str(n) + ".bin"
            dummy = np.fromfile(filename)
            self.prez = dummy.reshape((p.nx, p.ny), order="F")
    def getScalar(self, p, snap_time, filepath, zloc):
        self.thetaz = np.zeros((p.nx, p.ny))
        for n in range(1):
            filename = filepath + r"\output\scalar.z-" + "{loc:.5f}".format(
                loc=zloc) + "." + str(snap_time) + ".c" + str(n) + ".bin"
            dummy = np.fromfile(filename)
            self.thetaz = dummy.reshape((p.nx, p.ny), order="F")
        

class Point:
    def __init__(self, filepath, xloc, yloc, zloc, Scalar=False):
        if Scalar:
            self.getVelocity(filepath, xloc, yloc, zloc)
            self.getPre(filepath, xloc, yloc, zloc)
            self.getScalar(filepath, xloc, yloc, zloc)
        else:
            self.getVelocity(filepath, xloc, yloc, zloc)
            self.getPre(filepath, xloc, yloc, zloc)
    def getVelocity(self, filepath, xloc, yloc, zloc):
        f = open(filepath + r"\output\vel.x-" + "{loc:.5f}".format(loc=xloc) +
                  ".y-" + "{loc:.5f}".format(loc=yloc) + ".z-" +
                  "{loc:.5f}".format(loc=zloc) + ".dat")
        data = []
        for lines in f.readlines():
            data.append(list(map(float, lines.split())))
        f.close()
        dummy = np.array(data)
        self.t = dummy[:, 0]
        self.u = dummy[:, 1]
        self.v = dummy[:, 2]
        self.w = dummy[:, 3]
    def getPre(self, filepath, xloc, yloc, zloc):
        f = open(filepath + r"\output\pre.x-" + "{loc:.5f}".format(loc=xloc) +
                  ".y-" + "{loc:.5f}".format(loc=yloc) + ".z-" +
                  "{loc:.5f}".format(loc=zloc) + ".dat")
        data = []
        for lines in f.readlines():
            data.append(list(map(float, lines.split())))
        f.close()
        dummy = np.array(data)
        self.t = dummy[:, 0]
        self.pre = dummy[:, 1]
    def getScalar(self, filepath, xloc, yloc, zloc):
        f = open(filepath + r"\output\scalar.x-" + "{loc:.5f}".format(loc=xloc) +
                  ".y-" + "{loc:.5f}".format(loc=yloc) + ".z-" +
                  "{loc:.5f}".format(loc=zloc) + ".dat")
        data = []
        for lines in f.readlines():
            data.append(list(map(float, lines.split())))
        f.close()
        dummy = np.array(data)
        self.t = dummy[:, 0]
        self.scalar = dummy[:, 1]
