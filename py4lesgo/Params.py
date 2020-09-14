import numpy as np
import re


class Params():
    """
        all parametes are nondimensional
        
        filepath is the path of the data file, which must consists of 
        subfolder-output, do not include backslash
        
        xc and yc are locations of turbine, respectively, which are in 
        nondimensional form
        
        d is the diameter of turbine, using z_i as characteristic length scale.
    """
    #TODO build different arrays according to uv-plane or w-plane
    def __init__(self, filepath, xc, yc, d):
        
        filename = filepath + r"\output\lesgo_param.out"
        f = open(filename)
        # read the content in lesgo_param.out by line
        lesgoParam = f.readlines()  
        match = re.findall(r'\d{2}', lesgoParam[7])
        self.nproc = int(match[0])
        match = re.findall(r'\d{1,3}', lesgoParam[8])
        self.nx = int(match[0])
        self.ny = int(match[1])
        self.nz2 = int(match[2])
        self.nz_tot = int(match[3])
        match = re.findall(r'\d\.\d{7}E\+\d{2}', lesgoParam[9])
        self.z_i = float(match[0])
        match = re.findall(r'\d\.\d{7}E\+\d{2}', lesgoParam[10])
        self.L_x = float(match[0])
        self.L_y = float(match[1])
        self.L_z = float(match[2])
        match = re.findall(r'\d\.\d{7}E\-\d{2}', lesgoParam[11])
        self.dx = float(match[0])
        self.dy = float(match[1])
        self.dz = float(match[2])
        match = re.findall(r'\d\.\d{7}E\+\d{2}', lesgoParam[20])
        self.u_star = float(match[0])
        match = re.findall(r'\d\.\d{7}E\+\d{2}', lesgoParam[21])
        self.vonk = float(match[0])
        match = re.findall(r'\d\.\d{7}E\-\d{2}', lesgoParam[24])
        self.nu_molec = float(match[0])
        match = re.findall(r'\d+', lesgoParam[30])
        self.nsteps = int(match[0])
        match = re.findall(r'\d\.\d{7}E\-\d{2}', lesgoParam[45])
        self.zo = float(match[0])
        match = re.findall(r'\d\.\d{7}E\+\d{2}', lesgoParam[51])
        self.mean_p_force = float(match[0])
        match = re.findall(r'LITTLE_ENDIAN|BIG_ENDIAN', lesgoParam[59])
        self.write_endian = match[0]
        match = re.findall(r'LITTLE_ENDIAN|BIG_ENDIAN', lesgoParam[60])
        self.read_endian = match[0]
        if self.write_endian == 'BIG_ENDIAN':
            self.fmt = 's'
        else:
            self.fmt = 'a'
        if 'T' in lesgoParam[66]:
            self.avgVelocities=True
            match = re.findall(r'\d+', lesgoParam[67])
            self.tavg_nstart = int(match[0])
            self.tavg_nend = int(match[1])
            self.tavg_nskip = int(match[2])
        else:
            self.avgVelocities=False
        if 'T' in lesgoParam[68]:
            self.points=True
            match = re.findall(r'\d+', lesgoParam[69])
            self.point_nstart = int(match[0])
            self.point_nend = int(match[1])
            self.point_nskip = int(match[2])
            #TODO improve this
            match = re.findall(r'\d+', lesgoParam[70])
            self.points = int(match[0])
            self.point = np.zeros((self.points, 3))
            for i in range(self.points):
                match = re.findall(r'\d\.\d{7}E[+,-]\d{2}', lesgoParam[71+i])
                for j in range(3):
                    self.point[i,j] = float(match[j])
        else:
            self.points=False
            self.points=1
        if 'T' in lesgoParam[72+self.points-1]:
            self.domain_snapshots=True
            match = re.findall(r'\d+', lesgoParam[73+self.points-1])
            self.domain_nstart = int(match[0])
            self.domain_nend = int(match[1])
            self.domain_nskip = int(match[2])
        else:
            self.domain_snapshots=False
        if 'T' in lesgoParam[74+self.points-1]:
            self.x_snapshots=True
            match = re.findall(r'\d+', lesgoParam[75+self.points-1])
            self.xplane_nstart = int(match[0])
            self.xplane_nend = int(match[1])
            self.xplane_nskip = int(match[2])
            match = re.findall(r'\d+', lesgoParam[76+self.points-1])
            self.nx_planes = int(match[0])
            self.x_plane = np.zeros(self.nx_planes)
            for i in range(self.nx_planes):
                match = re.findall(r'\d\.\d{7}E[+,-]\d{2}', lesgoParam[76+self.points+i])
                self.x_plane[i] = float(match[i])
        else:
            self.x_snapshots=False
            self.nx_planes=3
        if 'T' in lesgoParam[76+self.points+self.nx_planes]:
            self.y_snapshots=True
            match = re.findall(r'\d+', lesgoParam[77+self.points+self.nx_planes])
            self.yplane_nstart = int(match[0])
            self.yplane_nend = int(match[1])
            self.yplane_nskip = int(match[2])
            match = re.findall(r'\d+', lesgoParam[78+self.points+self.nx_planes])
            self.ny_planes = int(match[0])
            self.y_plane = np.zeros(self.ny_planes)
            for i in range(self.nx_planes):
                match = re.findall(r'\d\.\d{7}E[+,-]\d{2}', lesgoParam[79+self.points+self.nx_planes+i])
                self.y_plane[i] = float(match[i])
        else:
            self.y_snapshots=False
            self.ny_planes=3
        if 'T' in lesgoParam[79+self.points+self.nx_planes+self.ny_planes]:
            self.z_snapshots=True
            match = re.findall(r'\d+', lesgoParam[80+self.points+self.nx_planes+self.ny_planes])
            self.zplane_nstart = int(match[0])
            self.zplane_nend = int(match[1])
            self.zplane_nskip = int(match[2])
            match = re.findall(r'\d+', lesgoParam[81+self.points+self.nx_planes+self.ny_planes])
            self.nz_planes = int(match[0])
            self.z_plane = np.zeros(self.nz_planes)
            for i in range(self.nz_planes):
                match = re.findall(r'\d\.\d{7}E[+,-]\d{2}', lesgoParam[82+self.points+self.nx_planes+self.ny_planes+i])
                self.z_plane[i] = float(match[i])
        else:
            self.z_snapshots=False
            self.nz_planes=3
        self.avgScalar=None
        if len(lesgoParam)>100:
            self.avgScalar=True
            match = re.findall(r'\d\.\d{7}E\+\d{2}|\d\.\d{7}E\-\d{2}|\-\d\.\d{7}E\-\d{2}', lesgoParam[len(lesgoParam)-9])
            self.wt_s = float(match[0])
            match = re.findall(r'\d\.\d{7}E\+\d{2}', lesgoParam[len(lesgoParam)-10])
            self.T_scale = float(match[0])
            match = re.findall(r'\d\.\d{7}E\+\d{2}', lesgoParam[len(lesgoParam)-35])
            self.g = float(match[0])
        f.close()
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
        acquire the grid coordinates
        input:
            xc:turbine's x location(nondimensional) relative to origin domain
            yc:turbine's y location(nondimensional) relative to origin domain
            d:turbine rotor's diameter(nondimensional)
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