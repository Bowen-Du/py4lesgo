# -*- coding: utf-8 -*-
"""
Created on Mon Feb 24 10:02:58 2020

@author: Du
"""
import numpy as np
from ..PartialInterp import partialx,partialy,partialz_uv_w,partialz_w_uv,interp_w_uv


class ReyStress:
    def __init__(self,Params,Tavg,ABLProperties,ABL,Domain,Dim=True):
        """
        Params,Tavg,ABLProperties are udf class,
        ABL and Domain are strings represent some details about LESGO
        """
        self.lesgotkebudget(Params,Tavg,ABL,Domain)
        if Dim:
            self.dimlesgotkebudget(Params,ABLProperties,Domain)
        else:
            self.nondimlesgotkebudget(ABLProperties,Domain)
    def lesgotkebudget(self,p,tavg,ABL,domain):
        # calculate the partial derivative of u
        dudx = partialx(p, tavg.u)  # uv grid
        dudy = partialy(p, tavg.u)  # uv grid
        dudz_w = partialz_uv_w(p, tavg.u)  # w grid
        # calculate the partial derivative of v
        dvdx = partialx(p, tavg.v)  # uv grid
        dvdy = partialy(p, tavg.v)  # uv grid
        dvdz_w = partialz_uv_w(p, tavg.v)  # w grid
        # calculate the partial derivative of w
        dwdx = partialx(p, tavg.w_uv)  # uv grid
        dwdy = partialy(p, tavg.w_uv)  # uv grid
        dwdz = partialz_w_uv(p, tavg.w)  # uv grid
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
        if ABL != "NBL":
            L = tavg.theta[:,:,0]/(vonk*ghat*tavg.sgst3[:,:,0])
        if ABL == "CBL":
            x = np.power(1-16*0.5*p.dz/L,0.25)
            Psi = 2*np.log(1/2*(1+x))+np.log(1/2*(1+x**2))-2*np.arctan(x)+np.pi/2
            Phi = x**(-1)
        elif ABL == "SBL":
            Psi = -5*0.5*p.dz/L
            Phi = 1 + 5*0.5*p.dz/L
        # theoretical velocity in the first uv grid
        demo = np.log(0.5 * p.dz / z0) - Psi
        u_avg[:, :] = np.sqrt(tavg.u[:, :, k]**2 + tavg.v[:, :, k]**2)
        ustar[:, :] = u_avg[:, :] * vonk / demo
        # w grid
        dudz_w[:, :, k] = ustar[:, :] / (
        0.5 * p.dz * vonk) * tavg.u[:, :, k] / u_avg[:, :]*Phi
        # uv grid
        dudz = interp_w_uv(p, dudz_w)
        # w grid
        dvdz_w[:, :, k] = ustar[:, :] / (
        0.5 * p.dz * vonk) * tavg.v[:, :, k] / u_avg[:, :]*Phi
        # uv grid
        dvdz = interp_w_uv(p, dvdz_w)
        duidxj=[[dudx,dudy,dudz],
                [dvdx,dvdy,dvdz],
                [dwdx,dwdy,dwdz]]
        uw_uv=interp_w_uv(p, tavg.uw)
        vw_uv=interp_w_uv(p, tavg.vw)
        rs=[[tavg.uu,tavg.uv,uw_uv],
            [tavg.uv,tavg.vv,vw_uv],
            [uw_uv  ,vw_uv  ,tavg.ww]]
        txz_uv=interp_w_uv(p, tavg.txz)
        tyz_uv=interp_w_uv(p, tavg.tyz)
        tau=[[tavg.txx,tavg.txy,txz_uv],
             [tavg.txy,tavg.tyy,tyz_uv],
             [txz_uv,  tyz_uv,  tavg.tzz]]
        delta =np.array([[1,0,0],[0,1,0],[0,0,1]])
        if domain == "blue":
            fpup = [[tavg.fxpup,tavg.fxpvp,tavg.fxpwp],
                    [tavg.fypup,tavg.fypvp,tavg.fypwp],
                    [tavg.fzpup,tavg.fzpvp,tavg.fzpwp]]
        #deltap =np.array([[0,0,1],[0,0,1],[0,0,2]])
        ppSijp = [[tavg.ppS11p,tavg.ppS12p,tavg.ppS13p],
                  [None       ,tavg.ppS22p,tavg.ppS23p],
                  [None       ,None       ,tavg.ppS33p]]
        ppup = [tavg.ppup,tavg.ppvp,tavg.ppwp]
        wpthetap_uv = interp_w_uv(p, tavg.wpthetap+tavg.sgst3)
        thetapup = [tavg.upthetap,tavg.vpthetap,wpthetap_uv]
        # uipujpukp = [[[tavg.upupup,tavg.upupvp,tavg.upupwp],
        #               [tavg.upupvp,tavg.vpvpup,tavg.upvpwp],
        #               [tavg.upupwp,tavg.upvpwp,tavg.wpwpup]],
        #              [[tavg.upupvp,tavg.vpvpup,tavg.upvpwp],
        #               [tavg.vpvpup,tavg.vpvpvp,tavg.vpvpwp],
        #               [tavg.upvpwp,tavg.vpvpwp,tavg.wpwpvp]],
        #              [[tavg.upupwp,tavg.upvpwp,tavg.wpwpup],
        #               [tavg.upvpwp,tavg.vpvpwp,tavg.wpwpvp],
        #               [tavg.wpwpup,tavg.wpwpvp,tavg.wpwpwp]]]
        uipujpukp = [[[tavg.upupup,tavg.upupvp,tavg.upupwp],
                      [None       ,tavg.vpvpup,tavg.upvpwp],
                      [None       ,None       ,tavg.wpwpup]],
                     [[tavg.upupvp,tavg.vpvpup,tavg.upvpwp],
                      [None       ,tavg.vpvpvp,tavg.vpvpwp],
                      [None       ,None       ,tavg.wpwpvp]],
                     [[tavg.upupwp,tavg.upvpwp,tavg.wpwpup],
                      [None       ,tavg.vpvpwp,tavg.wpwpvp],
                      [None       ,None       ,tavg.wpwpwp]]]
        uptaup = [[[tavg.uptau11p, tavg.uptau12p, tavg.uptau13p],
                   [tavg.uptau12p, tavg.uptau22p, tavg.uptau23p],
                   [tavg.uptau13p, tavg.uptau23p, tavg.uptau33p]],
                  [[tavg.vptau11p, tavg.vptau12p, tavg.vptau13p],
                   [tavg.vptau12p, tavg.vptau22p, tavg.vptau23p],
                   [tavg.vptau13p, tavg.vptau23p, tavg.vptau33p]],
                  [[tavg.wptau11p, tavg.wptau12p, tavg.wptau13p],
                   [tavg.wptau12p, tavg.wptau22p, tavg.wptau23p],
                   [tavg.wptau13p, tavg.wptau23p, tavg.wptau33p]]]
        # taupduidxjp = [[[tavg.tau11pdudxp, tavg.tau12pdudyp, tavg.tau13pdudzp],
        #                 [tavg.tau21pdudxp, tavg.tau22pdudyp, tavg.tau23pdudzp],
        #                 [tavg.tau31pdudxp, tavg.tau32pdudyp, tavg.tau33pdudzp]],
        #                [[tavg.tau11pdvdxp, tavg.tau12pdvdyp, tavg.tau13pdvdzp],
        #                 [tavg.tau21pdvdxp, tavg.tau22pdvdyp, tavg.tau23pdvdzp],
        #                 [tavg.tau31pdvdxp, tavg.tau32pdvdyp, tavg.tau33pdvdzp]],
        #                [[tavg.tau11pdwdxp, tavg.tau12pdwdyp, tavg.tau13pdwdzp],
        #                 [tavg.tau21pdwdxp, tavg.tau22pdwdyp, tavg.tau23pdwdzp],
        #                 [tavg.tau31pdwdxp, tavg.tau32pdwdyp, tavg.tau33pdwdzp]]]
        taupduidxjp = [[[tavg.tau11pdudxp, tavg.tau11pdvdxp, tavg.tau11pdwdxp],
                        [tavg.tau21pdudxp, tavg.tau21pdvdxp, tavg.tau21pdwdxp],
                        [tavg.tau31pdudxp, tavg.tau31pdvdxp, tavg.tau31pdwdxp]],
                       [[tavg.tau12pdudyp, tavg.tau12pdvdyp, tavg.tau12pdwdyp],
                        [tavg.tau22pdudyp, tavg.tau22pdvdyp, tavg.tau22pdwdyp],
                        [tavg.tau32pdudyp, tavg.tau32pdvdyp, tavg.tau32pdwdyp]],
                       [[tavg.tau13pdudzp, tavg.tau13pdvdzp, tavg.tau13pdwdzp],
                        [tavg.tau23pdudzp, tavg.tau23pdvdzp, tavg.tau23pdwdzp],
                        [tavg.tau33pdudzp, tavg.tau33pdvdzp, tavg.tau33pdwdzp]]]
        u=[tavg.u,tavg.v,tavg.w_uv]
        duudx=partialx(p, tavg.uu)
        duudy=partialy(p, tavg.uu)
        duudz=interp_w_uv(p, partialz_uv_w(p, tavg.uu))
        duvdx=partialx(p, tavg.uv)
        duvdy=partialy(p, tavg.uv)
        duvdz=interp_w_uv(p, partialz_uv_w(p, tavg.uv))
        dvvdx=partialx(p, tavg.vv)
        dvvdy=partialy(p, tavg.vv)
        dvvdz=interp_w_uv(p, partialz_uv_w(p, tavg.vv))
        dwwdx=partialx(p, tavg.ww)
        dwwdy=partialy(p, tavg.ww)
        dwwdz=interp_w_uv(p, partialz_uv_w(p, tavg.ww))
        duwdx=interp_w_uv(p, partialx(p, tavg.uw))
        duwdy=interp_w_uv(p, partialy(p, tavg.uw))
        duwdz=partialz_w_uv(p, tavg.uw)
        dvwdx=interp_w_uv(p, partialx(p, tavg.vw))
        dvwdy=interp_w_uv(p, partialy(p, tavg.vw))
        dvwdz=partialz_w_uv(p, tavg.vw)
        drsdxj=[[[duudx,duvdx,duwdx],
                 [None ,dvvdx,dvwdx],
                 [None ,None ,dwwdx]],
                [[duudy,duvdy,duwdy],
                 [None ,dvvdy,dvwdy],
                 [None ,None ,dwwdy]],
                [[duudz,duvdz,duwdz],
                 [None ,dvvdz,dvwdz],
                 [None ,None ,dwwdz]]]
        
        self.Epsilon=[[None,None,None],[None,None,None],[None,None,None]]
        self.D=[[None,None,None],[None,None,None],[None,None,None]]
        self.Phi=[[None,None,None],[None,None,None],[None,None,None]]
        self.Ptheta=[[None,None,None],[None,None,None],[None,None,None]]
        self.Ps=[[None,None,None],[None,None,None],[None,None,None]]
        self.Adv=[[None,None,None],[None,None,None],[None,None,None]]
        if domain == "blue":
            self.Pt=[[None,None,None],[None,None,None],[None,None,None]]
        self.R=[[None,None,None],[None,None,None],[None,None,None]]
        
        
        self.epsilon=[[None,None,None],[None,None,None],[None,None,None]]
        self.d=[[None,None,None],[None,None,None],[None,None,None]]
        self.phi=[[None,None,None],[None,None,None],[None,None,None]]
        self.ptheta=[[None,None,None],[None,None,None],[None,None,None]]
        self.ps=[[None,None,None],[None,None,None],[None,None,None]]
        self.adv=[[None,None,None],[None,None,None],[None,None,None]]
        if domain == "blue":
            self.pt=[[None,None,None],[None,None,None],[None,None,None]]
        self.r=[[None,None,None],[None,None,None],[None,None,None]]
        for i in range(3):
            for j in range(i,3):
                self.Epsilon[i][j] = taupduidxjp[0][i][j]+taupduidxjp[1][i][j]+taupduidxjp[2][i][j] +\
                                taupduidxjp[0][j][i]+taupduidxjp[1][j][i]+taupduidxjp[2][j][i]
                self.epsilon[i][j] = np.mean(self.Epsilon[i][j],axis=(0,1))
                temp=ppup[i]*delta[j,0]+ppup[j]*delta[i,0]+uipujpukp[0][i][j]+uptaup[i][j][0]+uptaup[j][i][0]
                D1=partialx(p, temp)
                temp=ppup[i]*delta[j,1]+ppup[j]*delta[i,1]+uipujpukp[1][i][j]+uptaup[i][j][1]+uptaup[j][i][1]
                D2=partialy(p, temp)
                temp=ppup[i]*delta[j,2]+ppup[j]*delta[i,2]+uipujpukp[2][i][j]+uptaup[i][j][2]+uptaup[j][i][2]
                D3=partialz_uv_w(p, temp)
                D3=interp_w_uv(p, D3)
                self.D[i][j]=-(D1+D2+D3)
                self.d[i][j]=np.mean(self.D[i][j],axis=(0,1))
                self.Phi[i][j]=2*ppSijp[i][j]
                self.phi[i][j]=np.mean(self.Phi[i][j],axis=(0,1))
                ghat = 9.81 * p.z_i / p.u_star**2
                self.Ptheta[i][j]=ghat*(delta[i,2]*thetapup[j]+delta[j,2]*thetapup[i])/tavg.thetaMean
                self.ptheta[i][j]=np.mean(self.Ptheta[i][j],axis=(0,1))
                # Ps[i][j] = -((rs[i][0]+tau[i][0])*duidxj[j][0]+(rs[i][1]+tau[i][1])*duidxj[j][1]+\
                #               (rs[i][2]+tau[i][2])*duidxj[j][2]+(rs[j][0]+tau[j][0])*duidxj[i][0]+\
                #               (rs[j][1]+tau[j][1])*duidxj[i][1]+(rs[j][2]+tau[j][2])*duidxj[i][2])
                self.Ps[i][j] = -(rs[i][0]*duidxj[j][0]+rs[i][1]*duidxj[j][1]+rs[i][2]*duidxj[j][2]+\
                             rs[j][0]*duidxj[i][0]+rs[j][1]*duidxj[i][1]+rs[j][2]*duidxj[i][2])
                self.ps[i][j]=np.mean(self.Ps[i][j],axis=(0,1))
                self.Adv[i][j]=-(u[0]*drsdxj[0][i][j]+u[1]*drsdxj[1][i][j]+u[2]*drsdxj[2][i][j])
                self.adv[i][j]=np.mean(self.Adv[i][j],axis=(0,1))
                if domain == "red":
                    self.R[i][j]=self.Epsilon[i][j]+self.D[i][j]+self.Phi[i][j]+\
                        self.Ptheta[i][j]+self.Ps[i][j]+self.Adv[i][j]
                    self.r[i][j]=self.ps[i][j]+self.ptheta[i][j]+self.phi[i][j]+\
                        self.d[i][j]+self.epsilon[i][j]+self.adv[i][j]
                else:
                    self.Pt[i][j]=fpup[i][j]+fpup[j][i]
                    self.pt[i][j]=np.mean(self.Pt[i][j],axis=(0,1))
                    self.R[i][j]=self.Epsilon[i][j]+self.D[i][j]+self.Phi[i][j]+\
                        self.Ptheta[i][j]+self.Ps[i][j]+self.Adv[i][j]+self.Pt[i][j]
                    self.r[i][j]=self.ps[i][j]+self.ptheta[i][j]+self.phi[i][j]+\
                        self.d[i][j]+self.epsilon[i][j]+self.adv[i][j]+self.pt[i][j]
    def dimlesgotkebudget(self,p,ABLProperties,Domain):
        """
        Using z_i and u_star to nondimensionalize different terms in ReyStress
        """
        for i in range(3):
            for j in range(i,3):
                self.Epsilon[i][j] = self.Epsilon[i][j] * ABLProperties.u_star**3/p.z_i
                self.epsilon[i][j] = self.epsilon[i][j] * ABLProperties.u_star**3/p.z_i
                self.Phi[i][j] = self.Phi[i][j] * ABLProperties.u_star**3/p.z_i
                self.phi[i][j] = self.phi[i][j] * ABLProperties.u_star**3/p.z_i
                self.D[i][j] = self.D[i][j] * ABLProperties.u_star**3/p.z_i
                self.d[i][j] = self.d[i][j] * ABLProperties.u_star**3/p.z_i
                self.Ptheta[i][j] = self.Ptheta[i][j] * ABLProperties.u_star**3/p.z_i
                self.ptheta[i][j] = self.ptheta[i][j] * ABLProperties.u_star**3/p.z_i
                self.Ps[i][j] = self.Ps[i][j] * ABLProperties.u_star**3/p.z_i
                self.ps[i][j] = self.ps[i][j] * ABLProperties.u_star**3/p.z_i
                self.Adv[i][j] = self.Adv[i][j] * ABLProperties.u_star**3/p.z_i
                self.adv[i][j] = self.adv[i][j] * ABLProperties.u_star**3/p.z_i
                if Domain=="red":
                    self.R[i][j] = self.R[i][j] * ABLProperties.u_star**3/p.z_i
                    self.r[i][j] = self.r[i][j] * ABLProperties.u_star**3/p.z_i
                else:
                    self.Pt[i][j] = self.Pt[i][j] * ABLProperties.u_star**3/p.z_i
                    self.pt[i][j] = self.pt[i][j] * ABLProperties.u_star**3/p.z_i
                    self.R[i][j] = self.R[i][j] * ABLProperties.u_star**3/p.z_i
                    self.r[i][j] = self.r[i][j] * ABLProperties.u_star**3/p.z_i
    def nondimlesgotkebudget(self,ABLProperties,Domain):
        """
        Using Uhub and D(!) to nondimensionalize different terms in ReyStress
        """
        for i in range(3):
            for j in range(i,3):
                self.Epsilon[i][j] = self.Epsilon[i][j] / 12.5 / ABLProperties.U_hub**3
                self.epsilon[i][j] = self.epsilon[i][j] / 12.5 / ABLProperties.U_hub**3
                self.Phi[i][j] = self.Phi[i][j] / 12.5 / ABLProperties.U_hub**3
                self.phi[i][j] = self.phi[i][j] / 12.5 / ABLProperties.U_hub**3
                self.D[i][j] = self.D[i][j] / 12.5 / ABLProperties.U_hub**3
                self.d[i][j] = self.d[i][j] / 12.5 / ABLProperties.U_hub**3
                self.Ptheta[i][j] = self.Ptheta[i][j] / 12.5 / ABLProperties.U_hub**3
                self.ptheta[i][j] = self.ptheta[i][j] / 12.5 / ABLProperties.U_hub**3
                self.Ps[i][j] = self.Ps[i][j] / 12.5 / ABLProperties.U_hub**3
                self.ps[i][j] = self.ps[i][j] / 12.5 / ABLProperties.U_hub**3
                self.Adv[i][j] = self.Adv[i][j] / 12.5 / ABLProperties.U_hub**3
                self.adv[i][j] = self.adv[i][j] / 12.5 / ABLProperties.U_hub**3
                if Domain=="red":
                    self.R[i][j] = self.R[i][j] / 12.5 / ABLProperties.U_hub**3
                    self.r[i][j] = self.r[i][j] / 12.5 / ABLProperties.U_hub**3
                else:
                    self.Pt[i][j] = self.Pt[i][j] / 12.5 / ABLProperties.U_hub**3
                    self.pt[i][j] = self.pt[i][j] / 12.5 / ABLProperties.U_hub**3
                    self.R[i][j] = self.R[i][j] / 12.5 / ABLProperties.U_hub**3
                    self.r[i][j] = self.r[i][j] / 12.5 / ABLProperties.U_hub**3
                
                
                