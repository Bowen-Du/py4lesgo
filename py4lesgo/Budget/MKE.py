import numpy as np
from scipy.interpolate import interpn
def partialx(p, x):
    y = np.zeros_like(x)

    i = 0
    y[i, :, :] = (x[i + 1, :, :] - x[p.nx - 1, :, :]) / 2 / p.dx
    y[1:p.nx - 1, :, :] = (x[2:p.nx, :, :] - x[0:p.nx - 2, :, :]) / 2 / p.dx
    i = p.nx - 1
    y[i, :, :] = (x[0, :, :] - x[i - 1, :, :]) / 2 / p.dx

    return y


def partialy(p, x):
    y = np.zeros_like(x)

    j = 0
    y[:, j, :] = (x[:, j + 1, :] - x[:, p.ny - 1, :]) / 2 / p.dy
    y[:, 1: p.ny - 1, :] = (x[:, 2: p.ny, :] - x[:, 0: p.ny - 2, :]) / 2 / p.dy
    j = p.ny - 1
    y[:, j, :] = (x[:, 0, :] - x[:, j - 1, :]) / 2 / p.dy

    return y


def partialz_uv_w(p, x):
    y = np.zeros_like(x)

    k = 0
    y[:, :, k] = x[:, :, k] * 2 / p.dz  #see it in the MKE.py, need log law
    y[:, :, 1:np.shape(x)[2] - 1] = (x[:, :, 1: np.shape(x)[2] - 1] - x[:, :, 0: np.shape(x)[2] - 2]) / p.dz
    k = np.shape(x)[2] - 1
    y[:, :, k] = 0

    return y


def partialz_w_uv(p, x):
    y = np.zeros_like(x)

    y[:, :, 0:np.shape(x)[2] - 1] = (x[:, :, 1:np.shape(x)[2]] - x[:, :, 0:np.shape(x)[2] - 1]) / p.dz
    k = np.shape(x)[2] - 1
    y[:, :, k] = 0

    return y


def interp_w_uv(p, x):
    y = np.zeros_like(x)
    
    y[:, :, 0:np.shape(x)[2] - 1] = (x[:, :, 1:np.shape(x)[2]] + x[:, :, 0:np.shape(x)[2] - 1]) / 2
    k = np.shape(x)[2] - 1
    y[:, :, k] = x[:, :, k]

    return y


def interp_uv_w(p, x):
    y = np.zeros_like(x)

    k = 0
    y[:, :, k] = x[:, :, k]
    y[:, :, 1:np.shape(x)[2]] = (x[:, :, 1:np.shape(x)[2]] + x[:, :, 0:np.shape(x)[2]-1]) / 2
    
    return y

class MKE:
    """
    calculate the mean kinetic energy transport equation terms
    """
    def __init__(self, p, avg, Scalar, domain, ABL):
        """
        avg is an instance of class Tavg
        
        ABL is string imply which atmospheric stability
        
        Calculate MKE(the square of time-averaged velocity) and partial
        derivative terms, take care of the grid plane where variables are
        stored
        """
        # Calculate MKE
        u2 = avg.u * avg.u / 2  # uv grid
        v2 = avg.v * avg.v / 2  # uv grid
        w2_uv = avg.w_uv * avg.w_uv / 2  # uv grid
        w2 = avg.w * avg.w / 2  # w grid
        # calculate the partial derivative of u2
        self.du2dx=partialx(p, u2) # uv grid
        self.du2dy=partialy(p, u2) # uv grid
        du2dz_w=partialz_uv_w(p, u2) # w grid
        self.du2dz=interp_w_uv(p, du2dz_w) # uv grid
        # calculate the partial derivative of v2
        self.dv2dx=partialx(p, v2) # uv grid
        self.dv2dy=partialy(p, v2) # uv grid
        dv2dz_w=partialz_uv_w(p, v2) # w grid
        self.dv2dz=interp_w_uv(p, dv2dz_w) # uv grid
        # calculate the partial derivative of w2
        self.dw2dx=partialx(p, w2_uv) # uv grid
        self.dw2dy=partialy(p, w2_uv) # uv grid
        self.dw2dz=partialz_w_uv(p, w2) # uv grid
        # interpolate Reynold Stress variables in w grid plane to uv grid
        # plane
        self.uw_uv = interp_w_uv(p, avg.uw)  # uv grid
        self.vw_uv = interp_w_uv(p, avg.vw)  # uv grid
        self.ww_uv = interp_w_uv(p, avg.ww)  # uv grid
        # interpolate SGS Stress variables in w grid plane to uv grid plane
        self.txz_uv = interp_w_uv(p, avg.txz)  # uv grid
        self.tyz_uv = interp_w_uv(p, avg.tyz)  # uv grid
        # calculate the partial derivative of u
        self.dudx = partialx(p, avg.u)  # uv grid
        self.dudy = partialy(p, avg.u)  # uv grid
        self.dudz_w = partialz_uv_w(p, avg.u)  # w grid
        # calculate the partial derivative of v
        self.dvdx = partialx(p, avg.v)  # uv grid
        self.dvdy = partialy(p, avg.v)  # uv grid
        self.dvdz_w = partialz_uv_w(p, avg.v)  # w grid
        # calculate the partial derivative of w
        self.dwdx = partialx(p, avg.w_uv)  # uv grid
        self.dwdy = partialy(p, avg.w_uv)  # uv grid
        self.dwdz = partialz_w_uv(p, avg.w)  # uv grid
        # use wall model to calculate dudz, dvdz
        ustar = np.zeros((p.nx, p.ny))
        u_avg = np.zeros((p.nx, p.ny))
        vonk = p.vonk
        z0 = p.zo
        k = 0  # wall model
        Psi = np.zeros((p.nx,p.ny))
        L = np.zeros((p.nx,p.ny))
        ghat = 9.81 * p.z_i / p.u_star**2
        if avg.sgst3[:,:,0].any() != 0:
            L = avg.theta[:,:,0]/(vonk*ghat*avg.sgst3[:,:,0])
        if ABL == "CBL":
            x = np.power(1-16*0.5*p.dz/L,0.25)
            Psi = 2*np.log(1/2*(1+x))+np.log(1/2*(1+x**2))-2*np.arctan(x)+np.pi/2
        elif ABL == "SBL":
            Psi = -5*0.5*p.dz/L
        # theoretical velocity in the first uv grid
        demo = np.log(0.5 * p.dz / z0) - Psi
        u_avg[:, :] = np.sqrt(avg.u[:, :, k]**2 + avg.v[:, :, k]**2)
        ustar[:, :] = u_avg[:, :] * vonk / demo
        # w grid
        self.dudz_w[:, :, k] = ustar[:, :] / (
            0.5 * p.dz * vonk) * avg.u[:, :, k] / u_avg[:, :]
        # uv grid
        self.dudz = interp_w_uv(p, self.dudz_w)
        # w grid
        self.dvdz_w[:, :, k] = ustar[:, :] / (
            0.5 * p.dz * vonk) * avg.v[:, :, k] / u_avg[:, :]
        # uv grid
        self.dvdz = interp_w_uv(p, self.dvdz_w)

        # call different functions to calculate corresponding terms in MKE
        # transport equation
        self.MKE_MC(avg)
        self.MKE_TC(p, avg)
        self.MKE_TP(avg)
        self.MKE_DP(avg)
        self.MKE_DF(p, avg)
        self.MKE_PT(p, avg)
        if Scalar:
            self.MKE_GA(p, avg)
        else:
            self.ga = np.zeros_like(avg.u)
        if domain == "blue":
            self.MKE_WT(avg)
            self.fp = np.zeros_like(avg.u)
        elif domain == "red":
            self.wt = np.zeros_like(avg.u)
            self.fp = avg.u * p.mean_p_force
        # this term must include the residual, so it may be different with pt
        # term which is calculated directly
        self.pteq = -(self.mc + self.tc + self.df + self.tp + self.dp +
                      self.wt + self.ga + self.fp)

    def MKE_MC(self, avg):
        """
        calculate the mean convection term
        """
        # numerator, all terms on uv grid
        mc1 = (avg.u * self.du2dx + avg.v * self.du2dy + avg.w_uv * self.du2dz)
        mc2 = (avg.u * self.dv2dx + avg.v * self.dv2dy + avg.w_uv * self.dv2dz)
        mc3 = (avg.u * self.dw2dx + avg.v * self.dw2dy + avg.w_uv * self.dw2dz)
        # denominator, all terms on uv grid
        mcx = avg.u * (self.du2dx + self.dv2dx + self.dw2dx)
        mcy = avg.v * (self.du2dy + self.dv2dy + self.dw2dy)
        mcz = avg.w * (self.du2dz + self.dv2dz + self.dw2dz)
        # also can write self.mc=-(mcx+mcy+mcz)
        self.mc = -(mc1 + mc2 + mc3)

        u2 = avg.u * avg.u / 2  # uv grid
        v2 = avg.v * avg.v / 2  # uv grid
        w2_uv = avg.w_uv * avg.w_uv / 2  # uv grid
        w2 = avg.w * avg.w / 2  # w grid
        # prepare for mc flux terms
        self.u2 = u2
        self.mcf12 = avg.u * (u2 + v2 + w2_uv)  # uv grid
        self.mcf34 = avg.v * (u2 + v2 + w2_uv)  # uv grid
        self.mcf56 = avg.w_uv * (u2 + v2 + w2_uv)  # uv grid

    def MKE_TC(self, p, avg):
        """
        calculate the turbulence convection term
        """
        # all terms on uv grid plane
        tc11 = partialx(p, avg.uu * avg.u)  # uv grid
        tc12 = partialy(p, avg.uv * avg.u)  # uv grid
        tc13_w = partialz_uv_w(p, self.uw_uv * avg.u)  # w grid
        tc13 = interp_w_uv(p, tc13_w)  # uv grid
        # all terms on uv grid plane
        tc21 = partialx(p, avg.uv * avg.v)  # uv grid
        tc22 = partialy(p, avg.vv * avg.v)  # uv grid
        tc23_w = partialz_uv_w(p, self.vw_uv * avg.v)  # w grid
        tc23 = interp_w_uv(p, tc23_w)  # uv grid
        # all terms on uv grid plane
        tc31 = partialx(p, self.uw_uv * avg.w_uv)  # uv grid
        tc32 = partialy(p, self.vw_uv * avg.w_uv)  # uv grid
        tc33 = partialz_w_uv(p, avg.ww * avg.w)  # uv grid
        # numerator,
        tc1 = (tc11 + tc12 + tc13)
        tc2 = (tc21 + tc22 + tc23)
        tc3 = (tc31 + tc32 + tc33)
        # denominator, these will be used in interpreting the contribution of
        # tc term in different directions
        self.tcx = (tc11 + tc21 + tc31)
        self.tcy = (tc12 + tc22 + tc32)
        self.tcz = (tc13 + tc23 + tc33)
        # also can write self.tc=tcx+tcy+tcz
        self.tc = -(tc1 + tc2 + tc3)

        # prepare for tc flux term
        # uv grid
        self.tcf12 = avg.uu * avg.u + avg.uv * avg.v + self.uw_uv * avg.w_uv
        # uv grid
        self.tcf34 = avg.uv * avg.u + avg.vv * avg.v + self.vw_uv * avg.w_uv
        # uv grid
        self.tcf56 = self.uw_uv * avg.u + self.vw_uv * avg.v
        + interp_w_uv(p, avg.ww * avg.w)

    def MKE_DF(self, p, avg):
        """
        calculate the diffusion term
        """
        # avoid incorrect values in the boundary
        avg.txx[:, :, -1] = 0
        avg.txy[:, :, -1] = 0
        avg.tyy[:, :, -1] = 0
        avg.tzz[:, :, -1] = 0
        # all terms on uv grid plane
        df11 = partialx(p, avg.u * avg.txx)  # uv grid
        df12 = partialy(p, avg.u * avg.txy)  # uv grid
        df13_w = partialz_uv_w(p, avg.u * self.txz_uv)  # w grid
        df13 = interp_w_uv(p, df13_w)  # uv grid
        # all terms on uv grid plane
        df21 = partialx(p, avg.v * avg.txy)  # uv grid
        df22 = partialy(p, avg.v * avg.tyy)  # uv grid
        df23_w = partialz_uv_w(p, avg.v * self.tyz_uv)  # w grid
        df23 = interp_w_uv(p, df23_w)  # uv grid
        # all terms on uv grid plane
        df31 = partialx(p, avg.w_uv * self.txz_uv)  # uv grid
        df32 = partialy(p, avg.w_uv * self.tyz_uv)  # uv grid
        df33_w = partialz_uv_w(p, avg.w_uv * avg.tzz)  # w grid
        df33 = interp_w_uv(p, df33_w)  #uv grid
        # sum up according to numerator
        df1 = df11 + df12 + df13
        df2 = df21 + df22 + df23
        df3 = df31 + df32 + df33
        # sum up according to denominator
        dfx = df11 + df21 + df31
        dfy = df12 + df22 + df32
        dfz = df13 + df23 + df33
        # also can write self.df=-(dfx+dfy+dfz)
        self.df = -(df1 + df2 + df3)
        # prepare for df flux term
        self.dff12 = avg.txx * avg.u + avg.txy * avg.v + self.txz_uv * avg.w_uv  #uv grid
        self.dff34 = avg.txy * avg.u + avg.tyy * avg.v + self.tyz_uv * avg.w_uv  #uv grid
        self.dff56 = self.txz_uv * avg.u + self.tyz_uv * avg.v + interp_w_uv(
            p, avg.tzz * avg.w)  # uv grid

    def MKE_TP(self, avg):
        """
        calculate the turbulence production term, take care of the wall model 
        because this term is associated with dudz, dvdz directly
        """
        #all terms on uv grid plane
        tp11 = avg.uu * self.dudx  #uv grid
        tp12 = avg.uv * self.dudy  #uv grid
        tp13 = self.uw_uv * self.dudz  #uv grid
        #all terms on uv grid plane
        tp21 = avg.uv * self.dvdx  #uv grid
        tp22 = avg.vv * self.dvdy  #uv grid
        tp23 = self.vw_uv * self.dvdz  #uv grid
        #all terms on uv grid plane
        tp31 = self.uw_uv * self.dwdx  #uv grid
        tp32 = self.vw_uv * self.dwdy  #uv grid
        tp33 = self.ww_uv * self.dwdz  #uv grid
        #sum up according to numerator
        tp1 = tp11 + tp12 + tp13
        tp2 = tp21 + tp22 + tp23
        tp3 = tp31 + tp32 + tp33
        #sum up according to denominator
        tpx = tp11 + tp21 + tp31
        tpy = tp12 + tp22 + tp32
        tpz = tp13 + tp23 + tp33
        #also can write self.tp=tpx+tpy+tpz
        self.tp = tp1 + tp2 + tp3

    def MKE_DP(self, avg):
        """
        calculate dissipation term using dudz, dvdz which already take wall model
        into account during the calculation of tp term
        """
        #uv grid
        dp11 = avg.txx * self.dudx  #uv grid
        dp12 = avg.txy * self.dudy  #uv grid
        dp13 = self.txz_uv * self.dudz  #uv grid
        #uv grid
        dp21 = avg.txy * self.dvdx  #uv grid
        dp22 = avg.tyy * self.dvdy  #uv grid
        dp23 = self.tyz_uv * self.dvdz  #uv grid
        #uv grid
        dp31 = self.txz_uv * self.dwdx  #uv grid
        dp32 = self.tyz_uv * self.dwdy  #uv grid
        dp33 = avg.tzz * self.dwdz  #uv grid
        #sum up according to numerator
        dp1 = dp11 + dp12 + dp13
        dp2 = dp21 + dp22 + dp23
        dp3 = dp31 + dp32 + dp33
        #sum up according to denominator
        dpx = dp11 + dp21 + dp31
        dpy = dp12 + dp22 + dp32
        dpz = dp13 + dp23 + dp33
        #also can write self.dp=dpx+dpy+dpz
        self.dp = (dp1 + dp2 + dp3)

    def MKE_GA(self, p, avg):
        #???need examination
        # calculate the Ga term which is related to the buoyancy effect
        thetaMean = avg.thetaMean  #uv grid
        # calculate the nondimensional gravity accelaration which will appear
        # in the Navier-Stokes equation
        ghat = 9.81 * p.z_i / p.u_star**2
        self.ga = np.zeros_like(avg.u)
        for k in range(np.shape(avg.u)[2]):
            self.ga[:, :, k] = ghat * (avg.theta[:, :, k] / thetaMean[k] -
                                       1) * avg.w_uv[:, :, k]  # uv grid

    def MKE_PT(self, p, avg):
        """
        calculate pressure-related term
        """
        ptx = avg.u * partialx(p, avg.pre)  #uv grid
        pty = avg.v * partialy(p, avg.pre)  #uv grid
        ptz_w = avg.w * partialz_uv_w(p, avg.pre)  #w grid
        ptz = interp_w_uv(p, ptz_w)  #uv grid
        self.pt = -(ptx + pty + ptz)

    def MKE_WT(self, avg):
        """
        calculate the WT(wind turbine) term
        """
        wtx = avg.fx * avg.u  #uv grid
        wty = avg.fy * avg.v  #uv grid
        wtz = avg.fz * avg.w_uv  #uv grid
        self.wt = wtx + wty + wtz
    def fdmtofvm(self, p, coord):
        points=np.zeros((p.nx,p.ny,p.zmax_buf[coord-1],3))
        for i in range(p.nx):
            for j in range(p.ny):
                for k in range(p.zmax_buf[coord-1]):
                    points[i,j,k,0]=p.x[i]+p.dx/2
                    points[i,j,k,1]=p.y[j]+p.dy/2
                    points[i,j,k,2]=p.z_uv[k]+p.dz/2
        for term in [self.mc,self.pt,self.tc,self.df,self.tp,self.dp,self.wt,self.ga,self.fp,
                      self.tcx,self.tcy,self.tcz,self.tcf12,self.tcf34,self.tcf56,
                      self.mcf12,self.mcf34,self.mcf56]:
            term=interpn([p.x,p.y,p.z_uv],term,points,bounds_error=False)
            

class Resolved_TKE:
    """
    calculatute different terms in turbulent kinetic energy transport equation
    """
    def __init__(self, p, avg, snap, mke, Scalar, m, n):
        """
        avg is an object which consists of time-averaged flow properties
        snap is an object which consists of instantaneous flow field
        mke is an object which consists of direvatives which can be used to 
        calculate mean strain-rate tensor
        m, n are indexes of Reynold Stress component or "Tke" 
        -i represents instantaneous which is related to snap field
        -p represents fluctuation relative to mean flow field
        """
        #according to the relationship between m and n to determine the function
        #that will be called
        #calculate different partial derivative for different velocity components
        duidx = partialx(p, snap.ui)  #uv grid
        duidy = partialy(p, snap.ui)  #uv grid
        duidz_w = partialz_uv_w(p, snap.ui)  #w grid
        dvidx = partialx(p, snap.vi)  #uv grid
        dvidy = partialy(p, snap.vi)  #uv grid
        dvidz_w = partialz_uv_w(p, snap.vi)  #w grid
        dwidx = partialx(p, snap.w_uvi)  #uv grid
        dwidy = partialy(p, snap.w_uvi)  #uv grid
        dwidz = interp_w_uv(p,
                                          partialz_uv_w(
                                              p, snap.w_uvi))  #uv grid
        #use wall model to determine duidz and dvidz
        ustar = np.zeros((p.nx, p.ny), dtype=np.float32)
        u_avg = np.zeros((p.nx, p.ny), dtype=np.float32)
        #calculate instantaneous strain-rate tensor
        vonk = p.vonk
        z0 = p.zo
        k = 0
##        #theoretical velocity in the first uv grid
##        demo = np.log(0.5 * p.dz / z0)
##        u_avg[:, :] = np.sqrt(snap.ui[:, :, k]**2 + snap.vi[:, :, k]**2)
##        ustar[:, :] = u_avg[:, :] * vonk / demo
##        #w grid
        duidz_w[:, :, k] = mke.dudz_w[:,:,k]
        #uv grid
        duidz = interp_w_uv(p, duidz_w)
##        #w grid
##        dvidz_w[:, :, k] = ustar[:, :] / (
##            0.5 * p.dz * vonk) * snap.vi[:, :, k] / u_avg[:, :]
        dvidz_w[:, :, k] = mke.dvdz_w[:,:,k]
        #uv grid
        dvidz = interp_w_uv(p, dvidz_w)
        #resolved-scale strain-rate tensor, all on uv grid plane
        S_11i = 1 / 2 * (duidx + duidx)
        S_12i = 1 / 2 * (duidy + dvidx)
        S_13i = 1 / 2 * (duidz + dwidx)
        S_22i = 1 / 2 * (dvidy + dvidy)
        S_23i = 1 / 2 * (dvidz + dwidy)
        S_33i = 1 / 2 * (dwidz + dwidz)
        S_21i = S_12i
        S_31i = S_13i
        S_32i = S_23i
        # calculate mean strain-rate tensor according to mke object
        S_11 = 1 / 2 * (mke.dudx + mke.dudx)
        S_12 = 1 / 2 * (mke.dudy + mke.dvdx)
        S_13 = 1 / 2 * (mke.dudz + mke.dwdx)
        S_22 = 1 / 2 * (mke.dvdy + mke.dvdy)
        S_23 = 1 / 2 * (mke.dvdz + mke.dwdy)
        S_33 = 1 / 2 * (mke.dwdz + mke.dwdz)
        
        #calculate the magnitude  of S
        S_total = np.zeros_like(avg.u)
        for i in range(1, 4):
            for j in range(1, 4):
                S_total = S_total + 2 * np.power(
                    eval('S_' + str(i) + str(j) + 'i'), 2)
        S_M = np.sqrt(S_total)
        #calculate the fliter width
        Delta = np.power(p.dx * p.dy * p.dz, 1 / 3)
        if Scalar:
            self.tau_11i = -2 * Delta**2 * snap.cs_opt2i * S_M * S_11i
            self.tau_12i = -2 * Delta**2 * snap.cs_opt2i * S_M * S_12i
            self.tau_13i = -2 * Delta**2 * snap.cs_opt2i * S_M * S_13i
            self.tau_22i = -2 * Delta**2 * snap.cs_opt2i * S_M * S_22i
            self.tau_23i = -2 * Delta**2 * snap.cs_opt2i * S_M * S_23i
            self.tau_33i = -2 * Delta**2 * snap.cs_opt2i * S_M * S_33i
            self.tau_21i = self.tau_12i
            self.tau_31i = self.tau_13i
            self.tau_32i = self.tau_23i
        else:
            self.tau_11i = -2 * Delta**2 * interp_w_uv(p, avg.cs_opt2) * S_M * S_11i
            self.tau_12i = -2 * Delta**2 * interp_w_uv(p, avg.cs_opt2) * S_M * S_12i
            self.tau_13i = -2 * Delta**2 * interp_w_uv(p, avg.cs_opt2) * S_M * S_13i
            self.tau_22i = -2 * Delta**2 * interp_w_uv(p, avg.cs_opt2) * S_M * S_22i
            self.tau_23i = -2 * Delta**2 * interp_w_uv(p, avg.cs_opt2) * S_M * S_23i
            self.tau_33i = -2 * Delta**2 * interp_w_uv(p, avg.cs_opt2) * S_M * S_33i
        
        
        #calculate fluctuation variable fields
        self.up = snap.ui - avg.u  #uv grid
        self.vp = snap.vi - avg.v  #uv grid
        self.w_uvp = snap.w_uvi - avg.w_uv  #uv grid
        self.tke = 1 / 2 * (self.up**2 + self.vp**2 + self.w_uvp**2)  #uv grid
        self.prep = snap.prei - avg.pre  #uv gird
        
        self.tau_11p = self.tau_11i - avg.txx
        self.tau_12p = self.tau_12i - avg.txy
        self.tau_13p = self.tau_13i - interp_w_uv(p, avg.txz)
        self.tau_22p = self.tau_22i - avg.tyy
        self.tau_23p = self.tau_23i - interp_w_uv(p, avg.tyz)
        self.tau_33p = self.tau_33i - avg.tzz
        self.tau_21p = self.tau_12p
        self.tau_31p = self.tau_13p
        self.tau_32p = self.tau_23p
        
        self.S_11p = S_11i - S_11
        self.S_12p = S_12i - S_12
        self.S_13p = S_13i - S_13
        self.S_22p = S_22i - S_22
        self.S_23p = S_23i - S_23
        self.S_33p = S_33i - S_33
        self.S_21p = self.S_12p
        self.S_31p = self.S_13p
        self.S_32p = self.S_23p
        if Scalar:
            self.thetap = snap.thetai - avg.theta  #uv grid
        #build a dictionary to maxmium the ability of function eval
        #use number indexes, circulation and eval to shorten the length of code
        self.m1 = avg.u
        self.m2 = avg.v
        self.m3 = avg.w_uv
        self.p1 = self.up
        self.p2 = self.vp
        self.p3 = self.w_uvp
        items = [("m1", self.m1), ("m2", self.m2), ("m3", self.m3),
                 ("p1", self.p1), ("p2", self.p2), ("p3", self.p3),
                 ("tau_11p", self.tau_11i), ("tau_12p", self.tau_12i),
                 ("tau_13p", self.tau_13i), ("tau_21p", self.tau_21i),
                 ("tau_22p", self.tau_22i), ("tau_23p", self.tau_23i),
                 ("tau_31p", self.tau_31i), ("tau_32p", self.tau_32i),
                 ("tau_33p", self.tau_33i), ("S_11p", self.S_11p),
                 ("S_12p", self.S_21p), ("S_13p", self.S_31p),
                 ("S_21p", self.S_21p), ("S_22p", self.S_22p),
                 ("S_23p", self.S_32p), ("S_31p", self.S_31p),
                 ("S_32p", self.S_32p), ("S_33p", self.S_33p)]
        d = dict(items)
        if m == n == "Tke":
            self.adv(p, avg)
            self.t_p(p, avg)
            self.t_t(p, avg)
            self.t_sgs(p, avg)
            self.p_s(p, avg, mke)
            self.p_t(p, avg)
            if Scalar:
                self.p_theta(p, avg)
            else:
                self.P_theta = np.zeros_like(avg.u)
            self.Epsilon(p, avg, d)
        else:
            self.c_mn(p, avg, m, n, d)
            self.p_mn(p, avg, m, n, d)
            self.phi_mn(p, avg, m, n, d)
            self.d_mn(p, avg, m, n, d)
            self.Epsilon_mn(p, avg, m, n, d)
            self.pt_mn(p, avg, m, n)
            if Scalar:
                self.Theta_mn(p, avg, m, n, d)
            else:
                self.theta_mn = np.zeros_like(avg.u)
    #rewrite function eval to use self.eval
    def eval(self, str, d):
        return d[str]
    

    #TKE transport equation

    def adv(self, p, avg):
        '''
        convection term
        '''
        adv1 = avg.u * partialx(p, self.tke)  #uv grid
        adv2 = avg.v * partialy(p, self.tke)  #uv grid
        adv3_w = avg.w * partialz_uv_w(p, self.tke)  #w grid
        adv3 = interp_w_uv(p, adv3_w)  #uv grid
        self.Adv = -(adv1 + adv2 + adv3)  #uv grid

    def t_t(self, p, avg):
        """
        T_t represent the transport of the resolved-scale TKE by the fluctuations of velocity
        """
        t_t1 = partialx(p, self.up * self.tke)  #uv grid
        t_t2 = partialy(p, self.vp * self.tke)  #uv grid
        t_t3_w = partialz_uv_w(p, self.w_uvp * self.tke)  #w grid
        t_t3 = interp_w_uv(p, t_t3_w)  #uv grid
        self.T_t = -(t_t1 + t_t2 + t_t3)  #uv grid

    def t_p(self, p, avg):
        '''
        T_p represent the transport of the resolved-scale TKE by the fluctuations pressure
        '''
        t_p1 = partialx(p, self.up * self.prep)  #uv grid
        t_p2 = partialy(p, self.vp * self.prep)  #uv grid
        t_p3_w = partialz_uv_w(p, 
                                             self.w_uvp * self.prep)  #w grid
        t_p3 = interp_w_uv(p, t_p3_w)  #uv grid
        self.T_p = -(t_p1 + t_p2 + t_p3)  #uv grid

    def t_sgs(self, p, avg):
        """
        T_sgs represent the transport of the resolved-scale TKE by the fluctuations of SGS stresses
        """

        upper1 = self.up * self.tau_11p + self.vp * self.tau_21p +\
        self.w_uvp * self.tau_31p #uv grid
        upper2 = self.up * self.tau_12p  + self.vp * self.tau_22p +\
        self.w_uvp * self.tau_32p #uv grid
        upper3 = self.up * self.tau_13p + self.vp * self.tau_23p +\
        self.w_uvp * self.tau_33p #uv grid

        t_sgs1 = partialx(p, upper1)  #uv grid
        t_sgs2 = partialy(p, upper2)  #uv grid
        t_sgs3_w = partialz_uv_w(p, upper3)  #w grid
        t_sgs3 = interp_w_uv(p, t_sgs3_w)  #uv grid
        self.T_sgs = -(t_sgs1 + t_sgs2 + t_sgs3)  #uv grid

    def p_s(self, p, avg, mke):
        """
        P_s represents the conversion of mean kinetic energy to resolved-scale TKE
        """
        p_s1 = self.up * self.up * mke.dudx + self.vp * self.up * mke.dvdx + self.w_uvp * self.up * mke.dwdx  #uv grid
        p_s2 = self.up * self.vp * mke.dudy + self.vp * self.vp * mke.dvdy + self.w_uvp * self.vp * mke.dwdy  #uv grid
        p_s3 = self.up * self.w_uvp * mke.dudz + self.vp * self.w_uvp * mke.dvdz + self.w_uvp * self.w_uvp * mke.dwdz  #uv grid
        #p_s3=interp_w_uv(p,p_s3_w) #uv grid
        self.P_s = -(p_s1 + p_s2 + p_s3)

    def p_theta(self, p, avg):
        """
        P_theta represents the buoyancy production or destruction depending on the sign of the vertical heat flux
        """
        
        thetaMean = avg.thetaMean  #uv grid
        ghat = 9.81 * p.z_i / p.u_star**2
        self.P_theta = np.zeros_like(avg.u)
        for k in range(np.shape(avg.u)[2]):
            self.P_theta[:, :,
                         k] = ghat * self.thetap[:, :, k] * self.w_uvp[:, :, 
                                                k] / thetaMean[k]
        self.P_theta[:, :, 0] = +ghat * (self.thetap[:, :, 0] * self.w_uvp[:, :, 
                                                0]+avg.sgst3[:, :, 0]) / thetaMean[0]
    def p_t(self, p, avg):
        """
        P_t represents the rate of the work against the turbine-induced forces
        since there is no instantaneous wind turbine-induced force fields
        we let this term is zero
        """
        self.P_t = np.zeros_like(avg.u)

    def Epsilon(self, p, avg, d):
        """
        epsilon is the SGS dissipation rate which represents the energy transfer 
        of resolved-scale TKE to the subgrid scales
        """
        epsi = np.zeros_like(avg.u)
        for i in range(1, 4):
            for j in range(1, 4):
                epsi += self.eval(
                    'tau_' + str(i) + str(j) + 'p', d) * self.eval(
                        'S_' + str(i) + str(j) + 'p', d)
        self.epsilon = +epsi

#Reynold stress transport equation

    def c_mn(self, p, avg, m, n, d):
        """
        convection term
        """
        C_mn1 = avg.u * partialx(
            p,
            self.eval('p' + str(m), d) * self.eval('p' + str(n), d))  #uv grid
        C_mn2 = avg.v * partialy(
            p,
            self.eval('p' + str(m), d) * self.eval('p' + str(n), d))  #uv grid
        C_mn3_w = avg.w * partialz_uv_w(
            p,
            self.eval('p' + str(m), d) * self.eval('p' + str(n), d))  #w grid
        C_mn3 = interp_w_uv(p, C_mn3_w)  #uv grid
        self.C_mn = -(C_mn1 + C_mn2 + C_mn3)  #uv grid

    def p_mn(self, p, avg, m, n, d):
        """
        production term by shear
        """
        P_mn1=self.eval('p'+str(n),d)*self.up*partialx(p,self.eval('m'+str(m),d))+\
        self.eval('p'+str(m),d)*self.up*partialx(p,self.eval('m'+str(n),d)) #uv grid
        P_mn2=self.eval('p'+str(n),d)*self.vp*partialy(p,self.eval('m'+str(m),d))+\
        self.eval('p'+str(m),d)*self.vp*partialy(p,self.eval('m'+str(n),d)) #uv grid
        P_mn3=self.eval('p'+str(n),d)*self.w_uvp*interp_w_uv(p,partialz_uv_w(p,self.eval('m'+str(m),d)))+\
        self.eval('p'+str(m),d)*self.w_uvp*interp_w_uv(p,partialz_uv_w(p,self.eval('m'+str(n),d))) #uv grid
        self.P_mn = -(P_mn1 + P_mn2 + P_mn3)

    def Theta_mn(self, p, avg, m, n, d):
        '''
        production term by buoyancy effect
        '''
        #term theta_mn
        delta = np.array([[0, 0, 0, 0], [0, 1, 0, 0], [0, 0, 1, 0],
                          [0, 0, 0, 1]])
       
        thetaMean = avg.thetaMean  #uv grid
        ghat = 9.81 * p.z_i / p.u_star**2
        self.theta_mn = np.zeros_like(avg.u)
        for k in range(np.shape(avg.u)[2]):
            self.theta_mn[:, :, k] = +ghat * (delta[m, 3] * self.thetap[:, :, k] * self.eval('p' + str(n), d)[:, :, k]
            + delta[n, 3] * self.thetap[:, :, k] * self.eval('p' + str(m), d)[:, :, k]) / thetaMean[k]
        if m == 3 and n == 3:
            self.theta_mn[:, :, 0] = +ghat * avg.sgst3[:, :, 0] / thetaMean[0]

    def phi_mn(self, p, avg, m, n, d):
        """
        product of pressure fluctation and strain-rate fluctuation
        """
        self.Phi_mn = +self.prep * 2 * self.eval('S_' + str(m) + str(n) + 'p',
                                                 d)

    def d_mn(self, p, avg, m, n, d):
        """
        diffusion term which has a divergence form
        """
        delta = np.array([[0, 0, 0, 0], [0, 1, 0, 0], [0, 0, 1, 0],
                          [0, 0, 0, 1]])
        D_mn1=partialx(p,self.prep*self.eval('p'+str(m),d)*delta[n,1]+\
                                     self.prep*self.eval('p'+str(n),d)*delta[m,1]+\
                                     self.eval('p'+str(m),d)*self.eval('p'+str(n),d)*self.up+\
                                     self.eval('p'+str(m),d)*self.eval('tau_'+str(n)+str(1)+'p',d)+\
                                     self.eval('p'+str(n),d)*self.eval('tau_'+str(m)+str(1)+'p',d)) #uv grid
        D_mn2=partialy(p,self.prep*self.eval('p'+str(m),d)*delta[n,2]+\
                                     self.prep*self.eval('p'+str(n),d)*delta[m,2]+\
                                     self.eval('p'+str(m),d)*self.eval('p'+str(n),d)*self.vp+\
                                     self.eval('p'+str(m),d)*self.eval('tau_'+str(n)+str(2)+'p',d)+\
                                     self.eval('p'+str(n),d)*self.eval('tau_'+str(m)+str(2)+'p',d)) #uv grid
        D_mn3_w=partialz_uv_w(p,self.prep*self.eval('p'+str(m),d)*delta[n,3]+\
                                            self.prep*self.eval('p'+str(n),d)*delta[m,3]+\
                                            self.eval('p'+str(m),d)*self.eval('p'+str(n),d)*self.w_uvp+\
                                            self.eval('p'+str(m),d)*self.eval('tau_'+str(n)+str(3)+'p',d)+\
                                            self.eval('p'+str(n),d)*self.eval('tau_'+str(m)+str(3)+'p',d)) #w grid
        D_mn3 = interp_w_uv(p, D_mn3_w)  #uv grid
        self.D_mn = -(D_mn1 + D_mn2 + D_mn3)  #uv grid

    def Epsilon_mn(self, p, avg, m, n, d):
        """
        disspative term
        """
        epsilon_mn1=self.eval('tau_'+str(n)+str(1)+'p',d)*partialx(p,self.eval('p'+str(m),d))+\
        self.eval('tau_'+str(m)+str(1)+'p',d)*partialx(p,self.eval('p'+str(n),d))#uv grid
        epsilon_mn2=self.eval('tau_'+str(n)+str(2)+'p',d)*partialy(p,self.eval('p'+str(m),d))+\
        self.eval('tau_'+str(m)+str(2)+'p',d)*partialy(p,self.eval('p'+str(n),d))#uv grid
        pw1 = interp_uv_w(p, self.p1)  #w grid
        pw2 = interp_uv_w(p, self.p2)  #w grid
        pw3 = interp_uv_w(p, self.p3)  #w grid
        epsilon_mn3=self.eval('tau_'+str(n)+str(3)+'p',d)*partialz_w_uv(p, eval('pw'+str(m)))+\
        self.eval('tau_'+str(m)+str(3)+'p',d)*partialz_w_uv(p, eval('pw'+str(n)))#uv grid
        #epsilon_mn3 = interp_w_uv(p, epsilon_mn3_w)#uv grid
        self.epsilon_mn = (epsilon_mn1 + epsilon_mn2 + epsilon_mn3)

    def pt_mn(self, p, avg, m, n):
        """
        term relevant with turbine-induced forces
        """
        self.Pt_mn = np.zeros_like(avg.u)
        