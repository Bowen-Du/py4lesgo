import numpy as np
import re


class Params():
    def __init__(self, filepath, xc, yc, d, num_x, num_y):
        """
        all parametes are nondimensional
        
        filepath is the path of the data file, which must consists of 
        subfolder-output, do not include backslash
        
        xc and yc are locations of turbine, respectively, which are in 
        nondimensional form
        
        d is the diameter of turbine, using z_i as characteristic length scale.

        num_x and num_y are row numbers and column numbers of turbine, respectively
        """
        filename = filepath + r"\output\lesgo_param.out"
        f = open(filename)
        # read the content in lesgo_param.out by line
        EOF = False
        self.avgScalar = False
        while not EOF:
            lesgoParam = f.readline()
            linelist = lesgoParam.split(' ')
            if len(lesgoParam) == 0:
                EOF = True
            else:
                if linelist[0]=="nproc":
                    match = re.findall(r'\d{2}', lesgoParam)
                    self.nproc = int(match[0])
                if linelist[0]=="nx,":
                    match = re.findall(r'\d{1,3}', lesgoParam)
                    self.nx = int(match[0])
                    self.ny = int(match[1])
                    self.nz2 = int(match[2])
                    self.nz_tot = int(match[3])
                if linelist[0]=="z_i":
                    match = re.findall(r'\d\.\d{7}E\+\d{2}', lesgoParam)
                    self.z_i = float(match[0])
                if linelist[0]=="L_x,":
                    match = re.findall(r'\d\.\d{7}E\+\d{2}', lesgoParam)
                    self.L_x = float(match[0])
                    self.L_y = float(match[1])
                    self.L_z = float(match[2])
                if linelist[0]=="dx,":
                    match = re.findall(r'\d\.\d{7}E\-\d{2}', lesgoParam)
                    self.dx = float(match[0])
                    self.dy = float(match[1])
                    self.dz = float(match[2])
                if linelist[0]=="u_star":
                    match = re.findall(r'\d\.\d{7}E\+\d{2}', lesgoParam)
                    self.u_star = float(match[0])
                if linelist[0]=="vonk":
                    match = re.findall(r'\d\.\d{7}E\+\d{2}', lesgoParam)
                    self.vonk = float(match[0])
                if linelist[0]=="nu_molec":
                    match = re.findall(r'\d\.\d{7}E\-\d{2}', lesgoParam)
                    self.nu_molec = float(match[0])
                if linelist[0]=="nsteps":
                    match = re.findall(r'\d+', lesgoParam)
                    self.nsteps = int(match[0])
                if linelist[0]=="dt":
                    match = re.findall(r'\d\.\d{7}E\-\d{2}', lesgoParam)
                    if len(match)>0:
                        self.dt = float(match[0])
                    else:
                        self.dt = False
                if linelist[0]=="zo":
                    match = re.findall(r'\d\.\d{7}E\-\d{2}', lesgoParam)
                    self.zo = float(match[0])
                #Scalar block
                if linelist[0]=="wt_s":
                    self.avgScalar=True
                    match = re.findall(r'\d\.\d{7}E\+\d{2}|\d\.\d{7}E\-\d{2}|\-\d\.\d{7}E\-\d{2}', lesgoParam)
                    self.wt_s = float(match[0])
                if linelist[0]=="T_scale":
                    match = re.findall(r'\d\.\d{7}E\+\d{2}', lesgoParam)
                    self.T_scale = float(match[0])
                if linelist[0]=="gravity":
                    match = re.findall(r'\d\.\d{7}E\+\d{2}', lesgoParam)
                    self.g = float(match[0])
                # tavg block
                if linelist[0]=="tavg_calc":
                    if 'T' in lesgoParam:
                        self.avgVelocities=True
                        lesgoParam = f.readline()
                        match = re.findall(r'\d+', lesgoParam)
                        self.tavg_nstart = int(match[0])
                        self.tavg_nend = int(match[1])
                        self.tavg_nskip = int(match[2])
                    else:
                        self.avgVelocities=False
                # domain snap block
                if linelist[0]=="domain_calc":
                    if 'T' in lesgoParam:
                        self.domain_snapshots=True
                        lesgoParam = f.readline()
                        match = re.findall(r'\d+', lesgoParam)
                        self.domain_nstart = int(match[0])
                        self.domain_nend = int(match[1])
                        self.domain_nskip = int(match[2])
                    else:
                        self.domain_snapshots=False
                # point block
                if linelist[0]=="point_calc":
                    if 'T' in lesgoParam:
                        self.points=True
                        lesgoParam = f.readline()
                        match = re.findall(r'\d+', lesgoParam)
                        self.point_nstart = int(match[0])
                        self.point_nend = int(match[1])
                        self.point_nskip = int(match[2])
                        lesgoParam = f.readline()
                        match = re.findall(r'\d+', lesgoParam)
                        self.points = int(match[0])
                        self.point = np.zeros((self.points, 3))
                        for i in range(self.points):
                            lesgoParam = f.readline()
                            match = re.findall(r'\d\.\d{7}E[+,-]\d{2}', lesgoParam)
                            for j in range(3):
                                self.point[i,j] = float(match[j])
                    else:
                        self.points=False
                # xplane block
                if linelist[0]=="xplane_calc":
                    if 'T' in lesgoParam:
                        self.x_snapshots=True
                        lesgoParam = f.readline()
                        match = re.findall(r'\d+', lesgoParam)
                        self.xplane_nstart = int(match[0])
                        self.xplane_nend = int(match[1])
                        self.xplane_nskip = int(match[2])
                        lesgoParam = f.readline()
                        match = re.findall(r'\d+', lesgoParam)
                        self.nx_planes = int(match[0])
                        self.x_plane = np.zeros(self.nx_planes)
                        for i in range(self.nx_planes):
                            lesgoParam = f.readline()
                            match = re.findall(r'\d\.\d{7}E[+,-]\d{2}', lesgoParam)
                            self.x_plane[i] = float(match[i])
                    else:
                        self.x_snapshots=False
                # yplane block
                if linelist[0]=="yplane_calc":
                    if 'T' in lesgoParam:
                        self.y_snapshots=True
                        lesgoParam = f.readline()
                        match = re.findall(r'\d+', lesgoParam)
                        self.yplane_nstart = int(match[0])
                        self.yplane_nend = int(match[1])
                        self.yplane_nskip = int(match[2])
                        lesgoParam = f.readline()
                        match = re.findall(r'\d+', lesgoParam)
                        self.ny_planes = int(match[0])
                        self.y_plane = np.zeros(self.ny_planes)
                        for i in range(self.ny_planes):
                            lesgoParam = f.readline()
                            match = re.findall(r'\d\.\d{7}E[+,-]\d{2}', lesgoParam)
                            self.y_plane[i] = float(match[i])
                    else:
                        self.y_snapshots=False
                # zplane block
                if linelist[0]=="zplane_calc":
                    if 'T' in lesgoParam:
                        self.z_snapshots=True
                        lesgoParam = f.readline()
                        match = re.findall(r'\d+', lesgoParam)
                        self.zplane_nstart = int(match[0])
                        self.zplane_nend = int(match[1])
                        self.zplane_nskip = int(match[2])
                        lesgoParam = f.readline()
                        match = re.findall(r'\d+', lesgoParam)
                        self.nz_planes = int(match[0])
                        self.z_plane = np.zeros(self.nz_planes)
                        for i in range(self.nz_planes):
                            lesgoParam = f.readline()
                            match = re.findall(r'\d\.\d{7}E[+,-]\d{2}', lesgoParam)
                            self.z_plane[i] = float(match[i])
                    else:
                        self.z_snapshots=False
        f.close()
        self.num_x=num_x
        self.num_y=num_y
        #notice the difference between random.randint and nself.random.randint
        #random.randint can produce a random int in [a,b](only one)
        #nself.random.randint can produce a series of random int
        #for example,[a,b,c) can produce cth random numbers in [a,b)
        #notice the index different between matlab and python
        self.zmin_buf = np.array(np.zeros(self.nproc), dtype=np.int64)
        self.zmax_buf = np.array(np.zeros(self.nproc), dtype=np.int64)
        for dummynproc in range(self.nproc):
            self.zmin_buf[dummynproc] = dummynproc * (self.nz2 - 1)
            self.zmax_buf[dummynproc] = dummynproc * (self.nz2 - 1) + self.nz2
        self.x = np.arange(0, self.L_x, self.dx)
        self.y = np.arange(0, self.L_y, self.dy)
        # for avg vels, avg (rs's, dudz,dvdz,txz,tyz)
        self.z_w = np.arange(0, self.L_z + 0.001, self.dz)
        # for inst vels, avg (txx,txy,tzz,txy)
        self.z_uv = np.arange(self.dz / 2, self.L_z + self.dz * 1 / 2 + 0.001,
                              self.dz)

        self.getGrid(xc, yc, d)

    def getGrid(self, xc, yc, d):
        #获取计算域的网格坐标
        '''
        obtain the grid coordinates
        input:
            for single wind turbine:
                xc:turbine's x location (nondimensional) relative to origin domain
                yc:turbine's y location (nondimensional) relative to origin domain
            for wind farm:
                xc=0, zc=0 is OK
            d:turbine rotor's diameter (nondimensional)
        ouput:
            1D, 2D and 3D coordinates which are used to plot figures
        '''
        X, Y, Z_w = np.mgrid[self.x[0]:self.x[self.nx - 1] +
                             1E-6:self.dx, self.y[0]:self.y[self.ny - 1] +
                             1E-6:self.dy, self.z_w[0]:self.z_w[self.nz_tot -
                                                                1] +
                             1E-6:self.dz]
        X, Y, Z_uv = np.mgrid[self.x[0]:self.x[self.nx - 1] +
                              1E-6:self.dx, self.y[0]:self.y[self.ny - 1] +
                              1E-6:self.dy, self.
                              z_uv[0]:self.z_uv[self.nz_tot - 1] +
                              1E-6:self.dz]
        self.X = (X - xc) / d
        self.Y = (Y - yc) / d
        self.Z_uv = Z_uv / d
        self.Z_w = Z_w / d

        #2D array to draw contourf
        #Slice xy
        self.X21 = self.X[:, :, 0]
        self.Y21 = self.Y[:, :, 0]
        #Slice yz
        self.X22 = self.Y[0, :, :]
        self.Y22_uv = self.Z_uv[0, :, :]
        self.Y22_w = self.Z_w[0, :, :]
        #Slice xz
        self.X23 = self.X[:, 0, :]
        self.Y23_uv = self.Z_uv[:, 0, :]
        self.Y23_w = self.Z_w[:, 0, :]

        #1D array to profiles
        self.X1 = self.X[:, 0, 0]
        self.Y1 = self.Y[0, :, 0]
        self.Z1_w = self.Z_w[0, 0, :]
        self.Z1_uv = self.Z_uv[0, 0, :]