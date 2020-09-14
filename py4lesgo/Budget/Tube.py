# -*- coding: utf-8 -*-
"""
Created on Mon Apr 27 12:45:54 2020

@author: Du
"""
# xloc=[0.08, 0.24, 0.40, 0.56, 0.72, 0.88, 1.04, 1.20]
# label=["x/D=-4","x/D=-2","x/D=0","x/D=2","x/D=4","x/D=6","x/D=8","x/D=10","x/D=12","x/D=14"]
# fig,axs=plt.subplots(1,2,sharex=True,sharey=True)
# for i,tube in zip(range(2),[a,b]):
#     tube.crossection(xloc)
#     for j in range(8):
#         axs[i].plot(tube.tubecrossection[j,:,1],tube.tubecrossection[j,:,2],
#                     label=label[j])
#         axs[i].set(title=case[i])
#         axs[i].axis("image")
# axs[0].legend(bbox_to_anchor=(0, 1.2, 2, .102), loc='lower left',
#             ncol=9, mode="expand", borderaxespad=0.)
# plt.show()

# #hub height
# plt.figure()
# y=a.tubecrossection[7,:,1]
# z=a.tubecrossection[7,:,2]
# plt.plot(y,z)
# dis=np.zeros(101)
# for i in range(101):
#     dis[i]=np.sqrt((y[i]-p_SBL.y[48])**2+(z[i]-p_SBL.z_uv[8])**2)
# r=max(dis)
# theta = np.arange(0, 2 * np.pi, 0.01)
# y = np.zeros(theta.size)
# z = np.zeros(theta.size)
# for j in range(11):
#     for i in range(theta.size):
#         y[i] = r*j/10 * np.cos(theta[i]) + p_SBL.y[48]
#         z[i] = r*j/10 * np.sin(theta[i]) + p_SBL.z_uv[8]
#     plt.plot(y, z, 'k--')
# for i in range(11):
#     plt.plot([p_SBL.y[48],p_SBL.y[48]+r*np.sin(i*2*np.pi/10)],
#               [p_SBL.z_uv[8],p_SBL.z_uv[8]+r*np.cos(i*2*np.pi/10)],'k--')
# dy, dz = p_SBL.dy, p_SBL.dz
# Y, Z = np.mgrid[slice(0, 1.6, dy),
#                 slice(p_SBL.dz/2, 1+dz, dz)]
# # Z,Y = np.meshgrid(p_SBL.z_uv,p_SBL.y)
# plt.pcolor(Y[45:52,2:16],Z[45:52,2:16],tavg[1].u[50,45:52,2:16])
# plt.colorbar()
# plt.xlim([0.75,0.85])
# plt.ylim([0.0125,0.1375])
# plt.xticks(np.arange(0.75,0.85,p_SBL.dy))
# plt.yticks(np.arange(0.0125,0.1375,p_SBL.dz))
# plt.grid(linestyle="--")
# plt.axis("image")
# plt.show()
# #centroid
# plt.figure()
# y=a.tubecrossection[7,:,1]
# z=a.tubecrossection[7,:,2]
# plt.plot(y,z)
# yc=np.mean(y)
# zc=np.mean(z)
# dis=np.zeros(101)
# for i in range(101):
#     dis[i]=np.sqrt((y[i]-yc)**2+(z[i]-zc)**2)
# r=max(dis)
# theta = np.arange(0, 2 * np.pi, 0.01)
# y = np.zeros(theta.size)
# z = np.zeros(theta.size)
# for j in range(11):
#     for i in range(theta.size):
#         y[i] = r*j/10 * np.cos(theta[i]) + yc
#         z[i] = r*j/10 * np.sin(theta[i]) + zc
#     plt.plot(y, z, 'k--')
# for i in range(11):
#     plt.plot([yc,yc+r*np.sin(i*2*np.pi/10)],
#               [zc,zc+r*np.cos(i*2*np.pi/10)],'k--')
# dy, dz = p_SBL.dy, p_SBL.dz
# Y, Z = np.mgrid[slice(0, 1.6, dy),
#                 slice(p_SBL.dz/2, 1+dz, dz)]
# # Z,Y = np.meshgrid(p_SBL.z_uv,p_SBL.y)
# plt.pcolor(Y[45:52,2:16],Z[45:52,2:16],tavg[1].u[50,45:52,2:16])
# plt.colorbar()
# plt.xlim([0.75,0.85])
# plt.ylim([0.0125,0.1375])
# plt.xticks(np.arange(0.75,0.85,p_SBL.dy))
# plt.yticks(np.arange(0.0125,0.1375,p_SBL.dz))
# plt.grid(linestyle="--")
# plt.axis("image")
# plt.show()
import numpy as np
from scipy.interpolate import interpn
class Tube:
    def __init__(self,p,vectorfield,closedcurvepoints,steps,deltat=1E-3,scheme="RK4"):
        points=max(np.shape(closedcurvepoints))
        tubelines=np.zeros((points,steps,3))
        tubelines[:,0,:]=closedcurvepoints
        def dxdt(point):
            u=interpn([p.x,p.y,p.z_uv],vectorfield[0],point,bounds_error=False)
            return u
        def dydt(point):
            v=interpn([p.x,p.y,p.z_uv],vectorfield[1],point,bounds_error=False)
            return v
        def dzdt(point):
            w=interpn([p.x,p.y,p.z_uv],vectorfield[2],point,bounds_error=False)
            return w
        if scheme=="explicit Euler":
            for i in range(points):
                for j in range(steps):
                    for k,velocity in zip(range(3),[dxdt(tubelines[i,j,:]),
                                                    dydt(tubelines[i,j,:]),
                                                    dzdt(tubelines[i,j,:])]):
                        if j < steps-1:
                            tubelines[i,j+1,k]=tubelines[i,j,k]+velocity*deltat
        elif scheme=="SSPRK(3,3)": 
            for i in range(points):
                for j in range(steps):
                    #SSPRK(3,3)
                    #stage 1
                    Fu0=dxdt(tubelines[i,j,:])[0]
                    Fv0=dydt(tubelines[i,j,:])[0]
                    Fw0=dzdt(tubelines[i,j,:])[0]
                    x1=tubelines[i,j,0]+deltat*Fu0
                    y1=tubelines[i,j,1]+deltat*Fv0
                    z1=tubelines[i,j,2]+deltat*Fw0
                    #stage 2
                    Fu1=dxdt([x1,y1,z1])[0]
                    Fv1=dydt([x1,y1,z1])[0]
                    Fw1=dzdt([x1,y1,z1])[0]
                    x2=3/4*tubelines[i,j,0]+1/4*x1+1/4*deltat*Fu1
                    y2=3/4*tubelines[i,j,1]+1/4*y1+1/4*deltat*Fv1
                    z2=3/4*tubelines[i,j,2]+1/4*z1+1/4*deltat*Fw1
                    #stage 3
                    Fu2=dxdt([x2,y2,z2])[0]
                    Fv2=dydt([x2,y2,z2])[0]
                    Fw2=dzdt([x2,y2,z2])[0]
                    if j < steps-1:
                        tubelines[i,j+1,0]=1/3*tubelines[i,j,0]+2/3*x2+2/3*deltat*Fu2
                        tubelines[i,j+1,1]=1/3*tubelines[i,j,1]+2/3*y2+2/3*deltat*Fv2
                        tubelines[i,j+1,2]=1/3*tubelines[i,j,2]+2/3*z2+2/3*deltat*Fw2
        elif scheme=="RK4":
            for i in range(points):
                for j in range(steps):
                    x=tubelines[i,j,0]
                    y=tubelines[i,j,1]
                    z=tubelines[i,j,2]
                    
                    K0x=dxdt([x,y,z])[0]
                    K0y=dydt([x,y,z])[0]
                    K0z=dzdt([x,y,z])[0]
                    
                    K1x=dxdt([x+deltat/2*K0x,y+deltat/2*K0y,z+deltat/2*K0z])[0]
                    K1y=dydt([x+deltat/2*K0x,y+deltat/2*K0y,z+deltat/2*K0z])[0]
                    K1z=dzdt([x+deltat/2*K0x,y+deltat/2*K0y,z+deltat/2*K0z])[0]
                    
                    K2x=dxdt([x+deltat/2*K1x,y+deltat/2*K1y,z+deltat/2*K1z])[0]
                    K2y=dydt([x+deltat/2*K1x,y+deltat/2*K1y,z+deltat/2*K1z])[0]
                    K2z=dzdt([x+deltat/2*K1x,y+deltat/2*K1y,z+deltat/2*K1z])[0]
                    
                    K3x=dxdt([x+deltat*K2x,y+deltat*K2y,z+deltat*K2z])[0]
                    K3y=dydt([x+deltat*K2x,y+deltat*K2y,z+deltat*K2z])[0]
                    K3z=dzdt([x+deltat*K2x,y+deltat*K2y,z+deltat*K2z])[0]
                    if j < steps-1:
                        tubelines[i,j+1,0]=x+deltat/6*(K0x+2*K1x+2*K2x+K3x)
                        tubelines[i,j+1,1]=y+deltat/6*(K0y+2*K1y+2*K2y+K3y)
                        tubelines[i,j+1,2]=z+deltat/6*(K0z+2*K1z+2*K2z+K3z)
        self.tubelines=tubelines
        self.points=points
    def crossection(self,xloc):
        self.tubecrossection=np.zeros((len(xloc),self.points,3))
        for j in range(len(xloc)):
            for i in range(self.points):
                self.tubecrossection[j,i,0]=xloc[j]
                self.tubecrossection[j,i,1]=np.interp(xloc[j],self.tubelines[i,:,0],self.tubelines[i,:,1],)
                self.tubecrossection[j,i,2]=np.interp(xloc[j],self.tubelines[i,:,0],self.tubelines[i,:,2],)
    def CartoPolar(self,p,field,xloc):
        self.Polarfield=np.zeros((len(xloc),self.points,self.points))
        for i in range(len(xloc)):
            y=self.tubecrossection[i,:,1]
            z=self.tubecrossection[i,:,2]
            yc=np.mean(y)
            zc=np.mean(z)
            dis=np.zeros(self.points)
            for j in range(self.points):
                dis[j]=np.sqrt((y[j]-yc)**2+(z[j]-zc)**2)
            r=max(dis)
            theta = np.linspace(0, 2 * np.pi, self.points)
            polarpoints=np.zeros((self.points,self.points,3))
            for m in range(self.points):
                for n in range(theta.size):
                    polarpoints[m,n,0] = xloc[i]
                    polarpoints[m,n,1] = r*(m/100-1/200) * np.cos(theta[n]) + yc
                    polarpoints[m,n,2] = r*(m/100-1/200) * np.sin(theta[n]) + zc
            self.Polarfield[i,:,:]=interpn([p.x,p.y,p.z_uv],field,polarpoints)
    def proandarea(self,p,xloc):
        self.proportion=np.zeros((len(xloc),self.points,self.points))
        self.area=np.zeros((len(xloc),self.points,self.points))
        for i in range(len(xloc)):
            y=self.tubecrossection[i,:,1]
            z=self.tubecrossection[i,:,2]
            yc=np.mean(y)
            zc=np.mean(z)
            dis=np.zeros(self.points)
            midpointdis=np.zeros(self.points)
            for j in range(self.points):
                dis[j]=np.sqrt((y[j]-yc)**2+(z[j]-zc)**2)
            for j in range(self.points-1):
                midpointdis[j]=(dis[j]+dis[j+1])/2
            midpointdis[self.points-1]=(dis[0]+dis[self.points-1])/2
            theta = np.linspace(0, 2 * np.pi, self.points)
            r=max(dis)
            for m in range(self.points):
                for n in range(theta.size):
                    rc=r*m/100
                    larc=r*(m-1)/100
                    self.area[i,m,n]=np.pi*(rc**2-larc**2)/theta.size
                    if rc<=midpointdis[n]:
                        self.proportion[i,m,n]=1
                    elif rc>midpointdis[n]>=larc:
                        self.proportion[i,m,n]=(midpointdis[n]-larc)/(r/100)
                    else:
                        self.proportion[i,m,n]=0
    def tubecrossectioninte(self,p,MKE,xloc):
        num=len(xloc)
        self.Mc = np.zeros(num)
        self.Pt = np.zeros(num)
        self.Tc = np.zeros(num)
        self.Df = np.zeros(num)
        self.Tp = np.zeros(num)
        self.Dp = np.zeros(num)
        self.Ga = np.zeros(num)
        self.Wt = np.zeros(num)
        self.Fp = np.zeros(num)
        for term,inte in zip([MKE.mc,MKE.pt,MKE.tc,MKE.df,MKE.tp,MKE.dp,MKE.ga,MKE.wt,MKE.fp],
                             [self.Mc,self.Pt,self.Tc,self.Df,self.Tp,self.Dp,self.Ga,self.Wt,self.Fp]):
            self.CartoPolar(p,term,xloc)
            for i in range(num):
                inte[i]=np.sum(self.Polarfield[i,:,:]*self.proportion[i,:,:]*self.area[i,:,:])
            
            
            
            
            
            
            
            
            