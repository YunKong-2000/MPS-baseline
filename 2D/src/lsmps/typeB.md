# lsmps type-B格式的实现
## 泰勒展开式
对于物理量$\phi$和i粒子以及与其相邻的粒子j，j粒子的泰勒展开式为
$$
\phi_j \simeq \phi_i + \nabla \phi_i \cdot \mathbf{r}_{ij} + \frac{1}{2} \nabla^2 \phi_i : (\mathbf{r}_{ij} \otimes \mathbf{r}_{ij})
$$
令$P_i(j) = \phi_i + \nabla \phi_i \cdot \mathbf{r}_{ij} + \frac{1}{2} \nabla^2 \phi_i : (\mathbf{r}_{ij} \otimes \mathbf{r}_{ij})$
则有
$$
P_i(j) = \phi_i + \phi_{x}x_{ij} + \phi_{y}y_{ij} + \frac{1}{2}\phi_{xx}x_{ij}^2 + \frac{1}{2}\phi_{yy}y_{ij}^2 + \phi_{xy}x_{ij}y_{ij}
$$
计算基函数，  
当j粒子为流体粒子时，  
$$
p_{ij} = 
[1,
\frac{x_{ij}}{r_s},
\frac{y_{ij}}{r_s},
\frac{x_{ij}^2}{2r_s^2},
\frac{y_{ij}^2}{2r_s^2},
\frac{x_{ij}y_{ij}}{r_s^2}]^T
$$
当j粒子为壁面粒子时，
$$
q_{ij} = 
[0,
n_x,
n_y,
\frac{x_{ij}n_x}{r_s},
\frac{y_{ij}n_y}{r_s},
\frac{x_{ij}n_y + y_{ij}n_x}{r_s}
]
$$
原始矩矩阵
$$
M_{i}^{fluid} = \sum_{j \in fluid} w_{ij}p_{ij}p_{ij}^T \\
M_{i}^{wall} = \sum_{j \in wall} w_{ij}q_{ij}q_{ij}^T
$$
lsmps矩阵
$$
C_{i} = (M_{i}^{fluid} + M_{i}^{wall})^{-1}
$$
压力梯度算子
$$
<\nabla P>_i= \frac{1}{r_s}\sum_{j \in fluid} w_{ij}P_j
\begin{bmatrix}
C_{i,2} \\
C_{i,3}
\end{bmatrix}
p_{ij}
+
\sum_{j \in wall} w_{ij}\rho \mathbf{g} \mathbf{n}_j
\begin{bmatrix}
C_{i,2} \\
C_{i,3}
\end{bmatrix}
q_{ij}
$$