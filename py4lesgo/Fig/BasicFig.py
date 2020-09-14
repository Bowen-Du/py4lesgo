import numpy as np
import matplotlib.pyplot as plt

class BasicFig:
    """
    set the format of figure globally which will be inheritted by QualFig and
    QuanFig
    
    include a function plot the basic properties of the flow like the 
    matlab-scripts
    
    Returns
    -------
    None.

    """
    def __init__(self):
        """
        set the format of figure globally, which can be inheritted by other 
        class such as QualFig and QuanFig

        Returns
        -------
        None.

        """
        
        # plt.rcParams['savefig.format'] = 'pdf'
        plt.rcParams['savefig.format']='svg'
        #most journals: 9 cm (or 3.5 inch) for single column width and 18.5 cm (or 7.3 inch) for double column width.
        plt.rcParams['figure.figsize'] = 8.27, 4.88
        plt.rcParams['figure.dpi'] = 240
        plt.rcParams['axes.titlesize']='large'
        plt.rcParams['axes.labelsize']='large'
        plt.rcParams['axes.titlepad']= 4.0
        plt.rcParams['axes.labelpad']= 1.0
        plt.rcParams['xtick.direction']='in'
        plt.rcParams['xtick.major.pad']=5.0
        # plt.rcParams['xtick.major.size']=3.0
        plt.rcParams['xtick.labelsize'] = 'large'
        plt.rcParams['ytick.direction']='in'
        # plt.rcParams['ytick.major.pad']=2.0
        # plt.rcParams['ytick.major.size']=3.0
        plt.rcParams['ytick.labelsize'] = 'large'
        plt.rcParams['text.usetex'] = True
        plt.rcParams['legend.fontsize'] = 'medium'
        plt.rcParams['figure.max_open_warning'] = 100
        # plt.rcParams['axes.labelsize'] = 20
        # plt.rcParams['axes.titlesize'] = 20
        # 
        # plt.rcParams['ytick.labelsize'] = 15
        # plt.rcParams['legend.fontsize'] = 15
        # plt.rcParams['lines.linewidth'] = 2
        # plt.rcParams['lines.markersize'] = 10
        # plt.rcParams['font.family'] = 'sans-serif'
        # plt.rcParams['font.sans-serif'] = 'DejaVu Sans'
        # plt.rcParams['font.family'] = 'serif'
        # plt.rcParams['font.serif'] = 'Times New Roman'
        # plt.rcdefaults()
    def basicplots(self,
                   p,
                   Tavg=None,
                   Snap=None,
                   Snap_X=None,
                   Snap_Y=None,
                   Snap_Z=None,
                   Point=None,
                   dns_profiles=False):
        """
        Plot the horizontally and time-averaged flow properties profiles
        
        Plot the countourf according to the logical parameters in arguments
        """
        #draw the basic figures
        #self.getMean(p, Tavg, avgScalar)
        zmax=np.shape(Tavg.u)[2]
        if dns_profiles:
            self.getdns()
        #avgScalar plays the role of Scalar
        # basic plots
        if p.avgVelocities and Tavg is not None:
            #plt.figure(figsize=[30,20])
            #specify the length and width of the figure, unit is inch, 1 inch=2.54cm.
            #dpi of my computer:13_U_U*7_U8, dpi is the number of pixel point per inch ?

            #compare LES streamwise velocity profile with that of D_US
            fig, axs = plt.subplots(2, 3, tight_layout=True)
            axs[0, 0].plot(p.z_uv[0:zmax], Tavg.uMean, 'b')
            if dns_profiles:
                axs[0, 0].plot(self.dns_z, self.dns_u, 'r')
            axs[0, 0].set(xlabel=r'$z/z_i$', ylabel=r'$\langle u \rangle$')
            

            #compare LES streamwise velocity profile with that of D_US and Log Law in a semilog coordinate
            kappa = p.vonk
            z0 = p.zo  #aerodynamic surface roughness, demensionaless
            loglaw = 1 / kappa * np.log(p.z_uv[0:zmax] / z0)  # rough wall
            p1, = axs[0, 1].semilogx(
                p.z_uv[0:zmax], loglaw, 'k')  # , is necessary for plot, semilog etc.
            p2 = axs[0, 1].scatter(p.z_uv[0:zmax],
                                   Tavg.uMean,
                                   marker='o',
                                   c='w',
                                   edgecolors='b')
            if dns_profiles:
                p3, = axs[0, 1].semilogx(self.dns_z, self.dns_u, 'r')
            axs[0, 1].set(xlabel=r'$z/z_i$', ylabel=r'$\langle \bar{\tilde{u}} \rangle$',
               xlim=[0.01, 1])
            if dns_profiles:
                axs[0, 1].legend([p1, p2, p3], ['Log Law', 'LES', 'DNS'])
            else:
                axs[0, 1].legend([p1, p2], ['Log Law', 'LES'])

            #compare LES's profiles of uw,txz,txz+uw with that of D_US
            p1 = axs[0, 2].scatter(p.z_w[0:zmax],
                                   -Tavg.uwMean,
                                   marker='o',
                                   c='w',
                                   edgecolors='b')
            p2 = axs[0, 2].scatter(p.z_w[0:zmax],
                                   -Tavg.txzMean,
                                   marker='o',
                                   c='w',
                                   edgecolors='c')
            p3 = axs[0, 2].scatter(p.z_w[0:zmax],
                                   -Tavg.txzMean - Tavg.uwMean,
                                   marker='o',
                                   c='w',
                                   edgecolors='k')
            if dns_profiles:
                p4, = axs[0, 2].plot(self.dns_z, -self.dns_uw, 'b')
                p5, = axs[0, 2].plot(self.dns_z, self.dns_tau, 'c')
                p6, = axs[0, 2].plot(self.dns_z, self.dns_tau - self.dns_uw,
                                     'k')  #what are these variables' meaning?
            axs[0, 2].plot(p.z_w[0:zmax], (1 - p.z_w[0:zmax]), 'g')
            #x+y=1
            axs[0, 2].set(xlabel=r'$z/z_i$', ylabel=r'$\langle -u^\prime w^\prime \rangle$')
              #use \ to acquire ' in the label text
            if dns_profiles:
                axs[0, 2].legend([p1, p2, p3, p4, p5, p6], [
                    "LES:"+r'$-\langle u^\prime w^\prime\rangle$',
                    "LES:"+r'$-\langle txz \rangle$',
                    "LES:"+r'$-(\langle txz \rangle + \langle u^\prime w^\prime \rangle)$',
                    "DNS:"+r'$-\langle u^\prime w^\prime\rangle$',
                    "DNS:"+r'$-\langle txz \rangle$',
                    "DNS:"+r'$-(\langle txz \rangle + \langle u^\prime w^\prime \rangle)$'
                ])
            else:
                axs[0, 2].legend(
                    [p1, p2, p3],
                    ["LES:"+r'$-\langle u^\prime w^\prime \rangle$',
                     "LES:"+r'$-\langle txz \rangle$',
                     "LES:"+r'$-(\langle txz \rangle + \langle u^\prime w^\prime \rangle)$'])

            #compare LES's profile of uu with that of D_US
            p1 = axs[1, 0].scatter(p.z_uv[0:zmax],
                                   Tavg.uuMean,
                                   marker='o',
                                   c='w',
                                   edgecolors='b')
            if dns_profiles:
                p2, = axs[1, 0].plot(self.dns_z, self.dns_uu, 'b')
            
            if dns_profiles:
                axs[1, 0].legend([p1, p2], ['LES', 'DNS'])
            else:
                axs[1, 0].legend([p1], ['LES'])
            axs[1, 0].set(xlabel=r'$z/z_i$', ylabel=r'$\langle u^\prime u^\prime \rangle$',
               ylim=[0, 8])

            #compare LES's profile of ww with that of D_US
            axs[1, 1].scatter(p.z_uv[0:zmax],
                              Tavg.vvMean,
                              marker='o',
                              c='w',
                              edgecolors='b')
            if dns_profiles:
                axs[1, 1].plot(self.dns_z, self.dns_vv, 'b')
            axs[1, 1].set(xlabel=r'$z/z_i$', ylabel=r'$\langle v^\prime v^\prime \rangle$',
               ylim=[0, 4])
            if dns_profiles:
                axs[1, 1].legend([p1, p2], ['LES', 'DNS'])
            else:
                axs[1, 1].legend([p1], ['LES'])

            #compare LES's profile of ww with that of DNS
            axs[1, 2].scatter(p.z_w[0:zmax],
                              Tavg.wwMean,
                              marker='o',
                              c='w',
                              edgecolors='b')
            if dns_profiles:
                axs[1, 2].plot(self.dns_z, self.dns_ww, 'b')
            axs[1, 2].set(xlabel=r'$z/z_i$', ylabel=r'$\langle w^\prime w^\prime \rangle$',
               ylim=[0, 2])
            if dns_profiles:
                axs[1, 2].legend([p1, p2], ['LES', 'DNS'])
            else:
                axs[1, 2].legend([p1], ['LES'])
            fig.show()

        if p.avgScalar and Tavg is not None:
            fig, ax = plt.subplots(tight_layout=True)
            ax.scatter(p.z_w[0:zmax],
                       Tavg.wpthetapMean,
                       marker='o',
                       c='w',
                       edgecolors='b')
            ax.scatter(p.z_w[0:zmax], Tavg.sgst3Mean, marker='o', c='w', edgecolors='c')
            ax.scatter(p.z_w[0:zmax],
                       Tavg.wpthetapMean + Tavg.sgst3Mean,
                       marker='o',
                       c='w',
                       edgecolors='k')
            ax.plot(p.z_w[0:zmax],p.wt_s*np.ones(zmax)/p.T_scale/p.u_star,'--')
            ax.set(xlabel=r'$z/z_i$', ylabel=r'$\langle  w^\prime\theta^\prime \rangle$',
                   ylim=[-0.001, 0.001])
            fig.show()

        if p.domain_snapshots and Snap is not None:
            fig, ax = plt.subplots(tight_layout=True)
            Y, X = np.meshgrid(p.y, p.x)
            mappable = ax.pcolor(X, Y, Snap.ui[:, :, 4])
            ax.set(xlabel=r'$x/z_i$', ylabel=r'$y/z_i$', title=r'Streamwise velocity at z = ' + str(p.z_uv[4]))
            ax.axis("image")
            fig.colorbar(mappable)
            fig.show()

        if p.x_snapshots and Snap_X is not None:
            fig, ax = plt.subplots(tight_layout=True)
            Z, Y = np.meshgrid(p.z_uv[0:zmax], p.y)
            mappable = ax.pcolor(Y, Z, Snap_X.ux)
            ax.set(xlabel=r'$y/z_i$', ylabel=r'$z/z_i$', title='Streamwise velocity at x = ' + str(p.nx_planes))
            ax.axis("image")
            fig.colorbar(mappable)
            fig.show()

        if p.y_snapshots and Snap_Y is not None:
            fig, ax = plt.subplots(tight_layout=True)
            Z, X = np.meshgrid(p.z_uv[0:zmax], p.x)
            mappable = ax.pcolor(X, Z, Snap_Y.uy)
            ax.set(xlabel=r'$x/z_i$', ylabel=r'$z/z_i$', title='Streamwise velocity at y = ' + str(p.ny_planes))
            ax.axis("image")
            fig.colorbar(mappable)
            fig.show()

        if p.z_snapshots and Snap_Z is not None:
            fig, ax = plt.subplots(tight_layout=True)
            Y, X = np.meshgrid(p.y, p.x)
            mappable = ax.pcolor(X, Y, Snap_Z.uz)
            ax.axis("image")
            ax.set(xlabel=r'$x/z_i$', ylabel=r'$y/z_i$', title='Streamwise velocity at z = ' + str(p.nz_planes))
            fig.colorbar(mappable)
            fig.show()

        if p.points and Point is not None:
            fig, ax = plt.subplots(tight_layout=True)
            ax.plot(Point.t, Point.up, label='u')
            ax.plot(Point.t, Point.vp, label='v')
            ax.plot(Point.t, Point.wp, label='w')
            ax.set(xlabel=r'$t$', ylabel=r'$u/u_*$', xlim=[min(Point.t), max(Point.t)])
            ax.legend()
            fig.show()

    def getdns(self):
        """
        load dns data as reference
        """
        filename = r"dns_profiles.txt"
        fp = open(filename)
        data = []
        for lines in fp.readlines():
            data.append(list(map(float, lines.split())))
            #use the fundamental string function 'append','split' to extract floating point number
        fp.close()
        dns_data = np.array(data)  #transfer list to array
        self.dns_z = dns_data[:, 0] / 1000  #z-plus -> z/h
        self.dns_u = dns_data[:, 1]  # u-plus
        self.dns_uw = dns_data[:, 2]
        self.dns_uu = dns_data[:, 3]
        self.dns_ww = dns_data[:, 4]
        self.dns_vv = dns_data[:, 5]
        self.dns_tau = dns_data[:, 7]
        self.dns_tot = dns_data[:, 8]
