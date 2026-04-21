# 考虑不可穿透壁面约束速度散度计算
## 目的
对于近壁面粒子，计算速度散度时需要满足不可穿透壁面边界条件来防止壁面穿透现象。

## 不可穿透壁面边界条件
对于近壁面粒子，不可穿透壁面边界条件写为，
$$
\mathbf{n} \cdot \mathbf{u}_i = \mathbf{n} \cdot \mathbf{u}_{wall}
$$
其中$\mathbf{u}_i$和$\mathbf{u}_{wall}$分别为流体粒子和壁面粒子的临时速度大小，不可穿透壁面边界的物理意义是靠近壁面的流体粒子和壁面在法向上的相对速度为零，这就保证了粒子不会穿透壁面
其中壁面粒子的临时速度为
$$
u_{wall}^* = u_{wall} + \Delta t\mathbf{g}
$$
也即是，计算流体粒子速度散度的时候，壁面粒子相应地需要考虑重力作用，对于壁面静止的算例，壁面临时速度应该为$\Delta t\mathbf{g}$而不是0

## 具体实现
**1、计算基向量及邻域流体粒子约束**  
基向量$\mathbf{P}_{ij}$为，
$$
\mathbf{P}_{ij}=[\frac{x_{ij}}{r_s} \frac{y_{ij}}{r_s} \frac{x_{ij}^2}{r_s^2} \frac{y_{ij}^2}{r_s^2} \frac{x_{ij}y_{ij}}{r_s^2} ]^\top
$$
流体粒子邻域贡献$M_i$
$$
\mathbf{M}_i= \sum_{j \in fluid} w_{ij}P_{ij}P_{ij}^\top
$$

**2、计算壁面边界条件约束**  
边界条件约束项$\mathbf{L}_{wall}$
$$
L_i^{wall} = \sum w_{ij} (n_j n_j^\top) \otimes (P_{ij} P_{ij}^\top)
$$
也即是
$$
L_i^{wall} = \begin{bmatrix}
 \sum w_{ij} n_x^2 (P_{ij} P_{ij}^\top)   & \sum w_{ij} n_x n_y (P_{ij} P_{ij}^\top) \\ 
 \sum w_{ij} n_x n_y (P_{ij} P_{ij}^\top) & \sum w_{ij} n_y^2 (P_{ij} P_{ij}^\top) 
 \end{bmatrix}
$$
$\mathbf{L}_{wall}$最终是一个10*10的矩阵

**3、计算最终的矩矩阵**    
原始矩矩阵为，
$$
C_i = I_d \otimes M_i + \gamma L_{wall}
$$

求逆得到$C_i^{-1}$


**4、计算源项**  
对于源项，邻域中的流体粒子和壁面粒子分开计算
$$
f_{fluid} = \sum_{j \in fluid} w_{ij}(\mathbf{u}_j - \mathbf{u}_i) \otimes P_{ij}
$$

$$
f_{wall} = \sum_{j \in wall} w_{ij}(n_j \otimes P_{ij})n_j^\top(\mathbf{u}_{wall} - \mathbf{u}_i)
$$

$f_{fluid}$和$f_{wall}$均为10维向量，最终源项为$f_i=f_{fluid}+ \gamma f_{wall}$

**5、计算散度**

假设速度为$\mathbf{u}=[u,v]$
则速度散度为$\nabla \cdot \mathbf{u} = \frac{\partial u}{\partial x} + \frac{\partial v}{\partial y}$  
令$c_i = C_i^{-1}f_i$，最终散度为
$$
\nabla \cdot \mathbf{u} = \frac{1}{r_s}(c_1 + c_7)
$$