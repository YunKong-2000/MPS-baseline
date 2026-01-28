# moment矩阵原有计算方法

## 离散方法
对于不需要加入壁面边界条件约束的算子，例如速度相关算子，使用常规周围粒子约束即可，梯度及拉普拉斯算子的离散仅考虑周围粒子的影响。 
基函数为，  
$\mathbf{P}_{ij} = [\frac{x}{r_e}, \frac{y}{r_e}, \frac{x^2}{r_e^{2}}, \frac{y^2}{r_e^{2}}, \frac{xy}{r_e^{2}}]$   
moment矩阵计算，  
$\mathbf{M_i}=\sum_{j \in neighbour} w_{ij} \mathbf{P}_{ij} \mathbf{P}_{ij}^{T}$  
梯度算子，  
$\nabla \phi = \frac{1}{r_e} \sum_{j \in neighbour} w_{ij}(\phi_j - \phi_i)
\begin{bmatrix}
M_{i,0} \\
M_{i,1} \\
\end{bmatrix}P_{ij}$   
散度算子，  
$\nabla \cdot \phi = \frac{1}{r_e} \sum_{j \in neighbour} w_{ij}(\phi_j - \phi_i) \cdot
\begin{bmatrix}
M_{i,0} \\
M_{i,1} \\
\end{bmatrix}P_{ij}$   
拉普拉斯算子，  
$\nabla^2 \phi = \frac{2}{r_e^2} \sum_{j \in neighbour} w_{ij}(\phi_j - \phi_i)
\begin{bmatrix}
M_{i,2} + M_{i,3}
\end{bmatrix}P_{ij}$ 

但是对于压力相关算子，由于需要代入壁面压力边界条件，需要对邻域中的壁面粒子的影响做特殊处理，于是对于该类算子离散需要重新计算moment矩阵，并在算子离散是使用不同的离散方法。  
基函数分为两种，对于邻域内的流体粒子，  
$\mathbf{P}_{ij} = [\frac{x}{r_e}, \frac{y}{r_e}, \frac{x^2}{r_e^{2}}, \frac{y^2}{r_e^{2}}, \frac{xy}{r_e^{2}}]$   
对于邻域内的壁面粒子，  
$\mathbf{Q}_{ij} = [n_x, n_y, 2n_x\frac{x}{r_e}, 2n_y\frac{y}{r_e}, n_x\frac{y}{r_e} + n_y\frac{x}{r_e}]$   
计算moment矩阵时，  
$\mathbf{M_i}=\sum_{j \in fluid} w_{ij} \mathbf{P}_{ij} \mathbf{P}_{ij}^{T} + \sum_{j \in wall} w_{ij} \mathbf{Q}_{ij} \mathbf{Q}_{ij}^{T}$  
梯度算子，  
$$ \nabla \phi = \frac{1}{r_e} \sum_{j \in fluid} w_{ij}(\phi_j - \phi_i)
\begin{bmatrix}
M_{i,0} \\
M_{i,1} \\
\end{bmatrix}P_{ij} 
+ \frac{1}{r_e} \sum_{j \in wall} w_{ij} r_e\rho\mathbf{g}\mathbf{n}
\begin{bmatrix}
M_{i,0} \\
M_{i,1} \\
\end{bmatrix}Q_{ij}
$$  
拉普拉斯算子，  
$$ \nabla^2 \phi = \frac{2}{r_e^2} \sum_{j \in fluid} w_{ij}(\phi_j - \phi_i)
\begin{bmatrix}
M_{i,2} + M_{i,3}
\end{bmatrix}P_{ij}
+ \frac{2}{r_e^2} \sum_{j \in wall} w_{ij}r_e\rho\mathbf{g}\mathbf{n}
\begin{bmatrix}
M_{i,2} + M_{i,3}
\end{bmatrix}Q_{ij}
$$ 




