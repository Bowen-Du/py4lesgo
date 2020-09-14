# py4lesgo

**py4lesgo**是基于python编写的用于后处理JHU的[LESGO](https://lesgo.me.jhu.edu)求解器输出数据的包。

安装方式

```
git clone https://github.com/Bowen-Du/py4lesgo.git
python setup.py install
```

## Params.py

该程序定义`class Params`，其中包含filepath+"\output\lesgo_param.out"中的参数(需根据LES版本进行对应修改)。

```python
from py4lesgo.Params import Params
filepath=r"..."
#x, y 分别是无量纲的风力机位置参数，充当参照物的作用, D是无量纲的风轮直径参数
p=Params(filepath,x,y,D)
```

## LoadLESOutput.py

该程序定义了`class Tavg, class Snap, class Snap_X, class Snap_Y, class Snap_Z, class Point`,不同的class可以将数据载入到程序的工作空间。

| class  | file name                                               | 'uv' grid                                                    | 'w'grid                                                      |
| ------ | ------------------------------------------------------- | ------------------------------------------------------------ | ------------------------------------------------------------ |
| Tavg   | veluv_avg.bin                                           | $\bar{u},\bar{v},\bar{w}$                                    |                                                              |
| Tavg   | velw_avg.bin                                            |                                                              | $\bar{w}$                                                    |
| Tavg   | vort_avg.bin                                            |                                                              | $\overline{\omega_x},\overline{\omega_y},\overline{\omega_z}$ |
| Tavg   | scalar_avg.bin                                          | $\bar{\theta}$                                               |                                                              |
| Tavg   | scalar2_avg.bin                                         | $\overline{u\theta},\overline{v\theta},\overline{w\theta},\overline{\theta\theta},\overline{\theta p^*_s}$ |                                                              |
| Tavg   | vel2_avg.bin                                            | $\overline{uu},\overline{uv},\overline{vv},\overline{uw},\overline{vw}$ | $\overline{ww}$                                              |
| Tavg   | force_avg.bin                                           | $\bar{f_x},\bar{f_y},\bar{f_z}$                              |                                                              |
| Tavg   | pres_avg.bin                                            | $\bar{p}^*_s=\bar{p}+1/3\overline{\tau_{kk}}$                |                                                              |
| Tavg   | ps_uv_avg                                               | $\overline{p_s^*u},\overline{p_s^*v},\overline{p_s^*w},\overline{p^*_sp^*_s}$ |                                                              |
| Tavg   | c2_opt2.bin                                             | $\bar{C_S^2}$                                                |                                                              |
| Tavg   | tau_avg.bin                                             | $\overline{\tau_{xx}},\overline{\tau_{xy}},\overline{\tau_{yy}},\overline{\tau_{zz}}$ | $\overline{\tau_{xz}},\overline{\tau_{yz}}$                  |
| Tavg   | rs.bin                                                  | $\overline{u^\prime u^\prime},\overline{u^\prime v^\prime}, \overline{v^\prime v^\prime},\overline{w^\prime w^\prime}$ | $\overline{u^\prime w^\prime},\overline{v^\prime w^\prime}$  |
| Tavg   | scalar_rs.bin                                           | $\overline{\theta^\prime u^\prime},\overline{\theta^\prime v^\prime},\overline{\theta^\prime w^\prime},\overline{\theta^\prime \theta^\prime},\overline{\theta^\prime p_s^{*\prime}}$ |                                                              |
| Tavg   | ps2_avg.bin                                             | $\overline{p_s^{*\prime} p_s^{*\prime}},\overline{p_s^{*\prime} u^\prime},\overline{p_s^{*\prime} v^\prime},\overline{p_s^{*\prime} w^\prime}$ |                                                              |
| Tavg   | ppSijp.bin                                              | $\overline{p_s^{*\prime} S_{11}^\prime},\overline{p_s^{*\prime} S_{12}^\prime},\overline{p_s^{*\prime} S_{13}^\prime},\overline{p_s^{*\prime} S_{22}^\prime},\overline{p_s^{*\prime} S_{23}^\prime},\overline{p_s^{*\prime} S_{33}^\prime},$ |                                                              |
| Tavg   | uipujpukp.bin                                           | $\overline{u^\prime u^\prime u^\prime},\overline{v^\prime v^\prime v^\prime},\overline{w^\prime w^\prime w^\prime},\overline{u^\prime u^\prime v^\prime},\overline{u^\prime u^\prime w^\prime},\overline{v^\prime v^\prime u^\prime},\overline{v^\prime v^\prime w^\prime},\\\overline{w^\prime w^\prime u^\prime},\overline{w^\prime w^\prime v^\prime},\overline{u^\prime v^\prime w^\prime}$ |                                                              |
| Tavg   | uptaup.bin                                              | $\overline{u^\prime\tau_{xx}^\prime},\overline{u^\prime\tau_{xy}^\prime},\overline{u^\prime\tau_{xz}^\prime},\overline{u^\prime\tau_{yy}^\prime},\overline{u^\prime\tau_{yz}^\prime},\overline{u^\prime\tau_{zz}^\prime},\\\overline{v^\prime\tau_{xx}^\prime},\overline{v^\prime\tau_{xy}^\prime},\overline{v^\prime\tau_{xz}^\prime},\overline{v^\prime\tau_{yy}^\prime},\overline{v^\prime\tau_{yz}^\prime},\overline{v^\prime\tau_{zz}^\prime},\\\overline{w^\prime\tau_{xx}^\prime},\overline{w^\prime\tau_{xy}^\prime},\overline{w^\prime\tau_{xz}^\prime},\overline{w^\prime\tau_{yy}^\prime},\overline{w^\prime\tau_{yz}^\prime},\overline{w^\prime\tau_{zz}^\prime},$ |                                                              |
| Tavg   | fipujp.bin                                              | $\overline{f_x^\prime u^\prime},\overline{f_x^\prime v^\prime},\overline{f_x^\prime w^\prime},\overline{f_y^\prime u^\prime},\overline{f_y^\prime v^\prime},\overline{f_y^\prime w^\prime},\overline{f_z^\prime u^\prime},\overline{f_z^\prime v^\prime},\overline{f_z^\prime w^\prime}$ |                                                              |
| Snap   | vel.&lt;timestep&gt;.bin                                | $u,v$                                                        | $w$                                                          |
| Snap   | scalar.&lt;timestep&gt;.bin                             | $\theta$                                                     |                                                              |
| Snap   | pres.&lt;timestep&gt;.bin                               | $p^*_s=p+1/3\tau_{kk}$                                       |                                                              |
| Snap   | vort.&lt;timestep&gt;.bin                               |                                                              | $\overline{\omega_x},\overline{\omega_y},\overline{\omega_z}$ |
| Snap_X | vel.x-&lt;xloc&gt;.&lt;timestep&gt;.bin                 | $u,v$                                                        | $w$                                                          |
| Snap_X | scalar.x-&lt;xloc&gt;.&lt;timestep&gt;.bin              | $\theta$                                                     |                                                              |
| Snap_X | pre.x-&lt;xloc&gt;.&lt;timestep&gt;.bin                 | $p^*_s=p+1/3\tau_{kk}$                                       |                                                              |
| Snap_Y | vel.y-&lt;yloc&gt;.&lt;timestep&gt;.bin                 | $u,v$                                                        | $w$                                                          |
| Snap_Y | scalar.y-&lt;yloc&gt;.&lt;timestep&gt;.bin              | $\theta$                                                     |                                                              |
| Snap_Y | pre.y-&lt;yloc&gt;.&lt;timestep&gt;.bin                 | $p^*_s=p+1/3\tau_{kk}$                                       |                                                              |
| Snap_Z | vel.z-&lt;zloc&gt;.&lt;timestep&gt;.bin                 | $u,v$                                                        | $w$                                                          |
| Snap_Z | scalar.z-&lt;zloc&gt;.&lt;timestep&gt;.bin              | $\theta$                                                     |                                                              |
| Snap_Z | pre.z-&lt;zloc&gt;.&lt;timestep&gt;.bin                 | $p^*_s=p+1/3\tau_{kk}$                                       |                                                              |
| Point  | vel.x-&lt;xloc&gt;.y-&lt;yloc&gt;.z-&lt;zloc&gt;.dat    | $u,v$                                                        | $w$                                                          |
| Point  | scalar.x-&lt;xloc&gt;.y-&lt;yloc&gt;.z-&lt;zloc&gt;.dat | $\theta$                                                     |                                                              |
| Point  | pre.x-&lt;xloc&gt;.y-&lt;yloc&gt;.z-&lt;zloc&gt;.dat    | $p^*_s=p+1/3\tau_{kk}$                                       |                                                              |

`Notice`:$\bar{}$ represents time-averaged. $S_{ij}=\frac{1}{2}(\frac{\partial u_i}{\partial x_j}+\frac{\partial u_j}{\partial x_i})$, $\bar{p}_s^*$ is modified($^*$) static($_s$) pressure.

```python
from py4lesgo.LoadLESOutput import Tavg,Snap,Snap_X, Snap_Y, Snap_Z, Point
#p是Params的一个实例，coord用来控制读入数据的大小，filepath是代表数据的存储路径
tavg=Tavg(p, coord, filepath, Scalar=False, Budget=False, modify=True)
snap=Snap(p, coord, snap_time, filepath, Scalar=False)
snap_x=Snap_X(p, coord, snap_time, filepath, xloc, Scalar=False)
snap_y=Snap_Y(p, coord, snap_time, filepath, yloc, Scalar=False)
snap_z=Snap_Z(p, snap_time, filepath, zloc, Scalar=False)
point_xyz=Point(filepath, xloc, yloc, zloc, Scalar=False)
#Scalar用于控制是否读入标量数据，Budget用于控制是否读入计算Reynolds Stress Budget的时均统计量ps2_avg.bin，ppSijp.bin,uipujpukp.bin,uptaup.bin,fipujp.bin,modify用于控制对rs.bin中数据的修正，通过看LESGO程序的版本或者风力机附近u'w'的分布来确定
```

## Budget

### MKE.py

MKE budget:
$$
\begin{equation}
\underbrace{\frac{\partial \bar{u}_i\bar{u}_i/2}{\partial t}}_{Storage\approx0}=\underbrace{-\bar{u}_j\frac{\partial \bar{u}_i\bar{u}_i/2}{\partial x_j}}_{mc}\underbrace{-\frac{\partial \bar{p}^{*}\bar{u}_j}{\partial x_j}}_{pt}\underbrace{-\frac{\partial \overline{u^\prime_iu^\prime_j}\bar{u}_i}{\partial x_j}}_{tc}\underbrace{+\frac{\partial \bar{\tau }_{ij}^d}{\partial x_j}}_{df}\underbrace{+\overline{u^\prime_iu^\prime_j}\frac{\partial \bar{u}_i}{\partial x_j}}_{tp}\underbrace{-\bar{\tau }_{ij}^d\frac{\partial \bar{u}_i}{\partial x_j}}_{dp}\underbrace{+\mathrm{g}\bar{u}_3(\frac{\bar{\theta}-\langle{\bar{\theta}}\rangle}{\langle{\bar{\theta}}\rangle})}_{\mathrm{g}a}\underbrace{+\frac{\bar{f}_i\bar{u}_i}{\rho}}_{wt}\underbrace{+F_P\bar{u}_1}_{fp}
\end{equation}
$$

基于LESGO输出数据计算上面MKE budget中的每一项(注意亚格子应力项)
### Integral.py

取一长方体控制体$(x_1\sim x_2, y_1\sim y_2, z_1\sim z_2)$，对控制体内MKE budget中的不同项的体积分以及面积分的大小进行分析
$$
Term=\int_{z_2}^{z_1} \int_{y_2}^{y_1} term \ \mathrm{d}y  \ \mathrm{d}z\\
\mathsf{TERM} =\int_{z_2}^{z_1} \int_{y_2}^{y_1} \int_{x_2}^{x_1} t \ \mathrm{d}x\ \mathrm{d}y  \ \mathrm{d}z\\
term\text{ can be }mc,pt,tc,df,tp,dp,\mathrm{g}a,wt,fp,\text{ and } Term\text{ is the corresponding
}Mc,Pt,Tc,Df,Tp,Dp,\mathrm{G}a,Wt,Fp.\\
\mathsf{TERM}\text{ represent }MC,PT,TC,DF,TP,DP,\mathrm{G}A,WT,FP
$$
**参考文献**：[Yang et al.](https://doi.org/10.1063/1.4907685), [Cortina et al.](https://doi.org/10.1103/PhysRevFluids.1.074402)

### Tube.py

利用平均速度场计算平均速度场的流管，通过对流管内部MKE budget中的不同项的流向变化规律进行分析

**参考文献**：[West et al.](https://doi.org/10.3390/en13051078), [Ge et al.](https://doi.org/10.1016/j.renene.2020.04.134)

### ReyStress.py

Reynolds stress budget:
$$
\begin{align}
\nonumber \frac{\partial \overline{{\widetilde{u}_i}^{'}{\widetilde{u}_j}^{'}}}{\partial t}=\underbrace{-\overline {{\widetilde u}}_k\frac{\partial \overline{{\widetilde{u}_i}^{'}{ {\widetilde u}_j^{'}}}}{\partial x_k}}_{C_{ij}}&\underbrace{-(\overline{{\widetilde{u}_j}^{'}{{\widetilde u}_k}^{'}}\frac{\partial {\overline {\widetilde u}_i}}{\partial x_k}+\overline{{\widetilde{u}_i}^{'}{{\widetilde u}_k}^{'}}\frac{\partial {\overline {\widetilde u}_j}}{\partial x_k})}_{P_{ij}}\underbrace{+\mathrm{g}\frac{(\delta_{i3}\overline{\widetilde{\theta}^{'}{\widetilde{u}_j}^{'}}+\delta_{j3}\overline{\widetilde{\theta}^{'}{\widetilde{u}_i}^{'}})}{\theta_0}}_{P_{\theta}}\underbrace{+\overline{{{\widetilde{p}}^{*'}}(\frac{\partial {\widetilde{u}_i}^{'}}{\partial x_j}+\frac{\partial {\widetilde{u}_j}^{'}}{\partial x_i})}}_{\Phi_{ij}}\\
  \nonumber &\underbrace{-\frac{\partial}{\partial x_k}(\overline{{{\widetilde{p}}^{*'}}{{\widetilde{u}_i}^{'}}\delta_{jk}+{{\widetilde{p}}^{*'}}{{\widetilde{u}_j}^{'}}\delta_{ik}+{{\widetilde{u}_i}^{'}}{{\widetilde{u}_j}^{'}}{{\widetilde{u}_k}^{'}}+\widetilde{u}_i^{'}\tau^{d'}_{jk}+\widetilde{u}_j^{'}\tau^{d'}_{ik}})}_{D_{ij}}
  \underbrace{+(\overline{\tau^{d'}_{jk}\frac{\partial {\widetilde{u}_i}^{'}}{\partial x_k}+\tau^{d'}_{ik}\frac{\partial {\widetilde{u}_j}^{'}}{\partial x_k}})}_{-\epsilon_{ij}}\\
  &\underbrace{-(\overline{\frac{{\widetilde{u}_j}^{'}f_i^{'}}{\rho}+\frac{{\widetilde{u}_i}^{'}f_j^{'}}{\rho}})}_{Pt_{ij}}
\end{align}
$$

基于LESGO输出数据计算Reynolds stress budget中的每一项(注意亚格子应力项，$\tau_{ij}^d$是亚格子应力的trace-free部分)


## PartialInterp.py

计算budget中的导数项

- uv-grid上物理量的处理

$$
uv-\text{grid}\longrightarrow uv-\text{grid}\\
\frac{\partial M^{uv}}{\partial x}|_i=\frac{M^{uv}_{i+1}-M^{uv}_{i-1}}{2\Delta x}\quad\text{partial x}\\
\frac{\partial M^{uv}}{\partial y}|_i=\frac{M^{uv}_{i+1}-M^{uv}_{i-1}}{2\Delta y}\quad\text{partial y}\\
uv-\text{grid}\longrightarrow w-\text{grid}\longrightarrow uv-\text{grid}\\
M^w_j=\frac{M^{uv}_{j-1}+M^{uv}_j}{2}\quad\text{interp_uv_w}\\
\frac{\partial M_w}{\partial z}|_j=\frac{M^w_{j+1}-M^w_j}{\Delta z}\quad\text{partialz_w_uv}
$$

- w-grid 上物理量的处理

$$
w-\text{grid}\longrightarrow uv-\text{grid}\\
\frac{\partial M_w}{\partial z}|_j=\frac{M^w_{j+1}-M^w_j}{\Delta z}\quad\text{partialz_w_uv}\\
w-\text{grid}\longrightarrow w-\text{grid}\longrightarrow uv-\text{grid}\\
\frac{\partial M^{w}}{\partial x}|_i=\frac{M^{w}_{i+1}-M^{w}_{i-1}}{2\Delta x}\quad\text{partial x}\\
\frac{\partial M^{w}}{\partial y}|_i=\frac{M^{w}_{i+1}-M^{w}_{i-1}}{2\Delta y}\quad\text{partial y}\\
M^{uv}_j=\frac{M^{w}_{j}+M^{w}_{j+1}}{2}\quad\text{interp_w_uv}\\
$$

## Differentiator

### explicit_differentiator.py

代码来源：[FloATPy](https://github.com/FPAL-Stanford-University/FloATPy)

class `ExplicitDifferentiator`

- ddx, ddy, ddz, d2dx2, d2dy2, d2dz2求给定数组(一维，二维或三维)的一阶以及二阶偏导数，可采用二阶，四阶或六阶精度中心差分格式，边界处可采用相同阶精度的单边差分格式
- gradient求给定数组(一维，二维或三维)的梯度ddx(一维)，(ddx, ddy)(二维)，(ddx, ddy, ddz)(三维)
- divergence求给定数组(一维，二维或三维)的散度
- curl求给定数组(二维或三维)的旋度
- laplacian求$\nabla^2$作用于给定数组(一维，二维或三维)的结果

可用来代替PartialInterp.py

### explicit_first.py

### explicit_second.py

## ABLProperties.py

- 计算尾流中心的函数`wakecenterzc`与`wakecenteryc`以及尾流摆振强度的函数`WakeMeandering`。
- class `ABLProperties`用来计算表征大气边界层状态的重要物理量($\langle\bar{u}\rangle,\langle\bar{\theta}\rangle,\langle\overline{w^\prime\theta^\prime}\rangle,\langle\overline{u^\prime w^\prime}\rangle,TI,I_u,I_v,I_w$)以及参数$L, Ri_f, Ri_b, Ri_t$

$$
L=-\frac{u_*^3\theta_0}{\kappa\mathrm{g}Q_0}\quad Ri_f=\frac{\frac{\mathrm{g}}{\theta_0}\overline{\theta^\prime w^\prime}}{\overline{u^\prime w^\prime}\frac{\partial U}{\partial z}}\quad Ri_t=\frac{\frac{\mathrm{g}}{\theta_0}\frac{\partial \theta}{\partial z}}{(\frac{\partial U}{\partial z})^2}\quad Ri_b=\frac{\frac{\mathrm{g}}{\theta_0}(\frac{\Delta \theta}{\Delta z})}{(\frac{\Delta U}{\Delta z})^2}
$$

- 数据输出为.vtk格式进行三维可视化（参考NREL的FLORIS）

##  Fig

### BasicFig.py

设置`plt.rcParams`控制输出图片的格式，输出基本图片

- Tavg
- Snap
- Snap_X
- Snap_Y
- Snap_Z
- Point

### QualFig.py

class `ABLProperties`

- VerCutPlane给出$x-z$ 平面某一位置处的contourf
- HorCutPlane给出$x-y$ 平面某一位置处的contourf
- CrossCutPlane给出$y-z$ 平面某一位置处的contourf
- VerWindTurbine给出风力机在$x-z$ 平面的侧视图
- HorWindTurbine给出风力机在$x-y$ 平面的俯视图
- CrossWindTurbine给出风力机在$y-z$ 平面的风轮轮廓

### QuanFig.py

class `QuanFig`

- casefig对比不同算例下计算域内流动的基本特征，风廓线，位温廓线，湍流剪切应力，总湍流强度，不同方向的湍流强度分量等
- Profiles画出风力机下游不同位置处(x)某一物理量随位置坐标(y or z)的变化规律
- VelSelfSim对风力机下游不同位置处(x)速度廓线的自相似特性进行验证（参考文献：[Bastankhah et al.](https://doi.org/10.1016/j.renene.2014.01.002))

