# 速度耦合的压力梯度计算 
为了消除模拟中出现的近自由面粒子的速度和压力梯度“高-低-高-低”的棋盘格现象，需要在计算压力梯度时考虑粒子相对于周围粒子出现的速度突变情况。
## 整体流程
1、计算临时速度  
2、根据临时速度计算每个粒子的速度梯度并存储  
3、构建并求解PPE  
4、计算压力梯度，计算压力梯度时加入额外约束项，使得粒子的压力梯度能够抵消粒子间的速度突变  

## 计算临时速度
和原先的实现相同，没有特殊处理

## 计算临时速度的梯度
对于每个流体粒子，都要计算出当地的速度梯度，并保存起来，作为后续的速度插值备用  

## 构建并求解PPE
和原有实现相同

## 计算压力梯度
1、计算仅考虑流体邻域粒子和壁面边界条件的基向量以及矩矩阵
$$
p_{ij} = 
[
  \frac{x_{ij}}{r_s},
  \frac{y_{ij}}{r_s},
  \frac{x_{ij}^2}{r_s^2},
  \frac{y_{ij}^2}{r_s^2},
  \frac{x_{ij}y_{ij}}{r_s^2}
]^T  \\
q_{ij} = 
[
  n_x,
  n_y,
  \frac{2x_{ij}n_x}{r_s},
  \frac{2y_{ij}n_y}{r_s},
  \frac{x_{ij}n_y + y_{ij}n_x}{r_s}
]^T  \\
M_i^{vel} = \sum_{j \in fluid} w_{ij}p_{ij}p_{ij}^T +
\sum_{j \in wall} w_{ij}q_{ij}q_{ij}^T
$$

2、根据粒子i和邻域粒子j的速度梯度分别插值两者中心位置$\frac{r_i + r_j}{2}$的速度
假设i粒子的速度$V_i^* =[u_i^* , v_i^*]$, 位置 $r_i =[x_i , y_i]$  
假设$d_{ij} = r_j - r_i, dx = x_j - x_i$,$dy = y_j - y_i$ 
粒子间的方向向量$e_{ij}=[dx / |d_{ij}|, dy / |d_{ij}|]$  
**如果j粒子是流体粒子** 
j粒子必须在i粒子1.2倍粒子间距内的邻域粒子。  
根据i粒子的速度梯度插值算出的速度：
$$h_i(\frac{r_i + r_j}{2})= 
[
  u_i^* + \frac{\partial u_i^*}{\partial x} * \frac{dx}{2} + \frac{\partial u_i^*}{\partial y} * \frac{dy}{2},
  v_i^* + \frac{\partial v_i^*}{\partial x} * \frac{dx}{2} + \frac{\partial v_i^*}{\partial y} * \frac{dy}{2}
]
$$
根据j粒子的速度梯度插值算出的速度：
$$h_j(\frac{r_i + r_j}{2})= 
[
  u_j^* - \frac{\partial u_j^*}{\partial x} * \frac{dx}{2} - \frac{\partial u_i^*}{\partial y} * \frac{dy}{2},
  v_j^* - \frac{\partial v_j^*}{\partial x} * \frac{dx}{2} - \frac{\partial v_i^*}{\partial y} * \frac{dy}{2}
]
$$
两者的差值：
$$
\Delta h_{ij} = h_j - h_i=
[
  u_j^* - u_i^* - \frac{1}{2}(
    (\frac{\partial u_i^*}{\partial x} + \frac{\partial u_j^*}{\partial x})dx 
    +
    (\frac{\partial u_i^*}{\partial y} + \frac{\partial u_j^*}{\partial y})dy
  ),\\
  v_j^* - v_i^* - \frac{1}{2}(
    (\frac{\partial v_i^*}{\partial x} + \frac{\partial v_j^*}{\partial x})dx 
    +
    (\frac{\partial v_i^*}{\partial y} + \frac{\partial v_j^*}{\partial y})dy
  )
]
$$
得到间断的大小
$V_{jump}=e_{ij} \cdot \Delta h_{ij}$
计算基与$e_{ij}$的基函数
$$
q_{ij}^{'} = 
[
e_x,
e_y,
\frac{dx * e_x}{|d_{ij}|},
\frac{dy * e_y}{|d_{ij}|},
\frac{dx * e_x + dy * e_y}{|d_{ij}|}
]
$$
权重函数
$$
w_{ij}^{vel} = 
\begin{cases}
(1-\frac{|d_{ij}|}{1.2l_0})^2, |d_{ij}| < 1.2l_0 \\
0, |d_{ij}| >= 1.2l_0
\end{cases}
$$

**如果j粒子是壁面粒子**  
j粒子必须是距离i粒子0.6倍初始间距内的壁面粒子。  
因为壁面粒子没有计算临时速度的梯度，所以在j粒子处进行插值 
根据i粒子的速度梯度插值算出的速度：
$$h_i(r_j)= 
[
  u_i^* + \frac{\partial u_i^*}{\partial x} * dx + \frac{\partial u_i^*}{\partial y} * dy,
  v_i^* + \frac{\partial v_i^*}{\partial x} * dx + \frac{\partial v_i^*}{\partial y} * dy
]
$$
则两者差值为：
$$
\Delta h_{ij} =
[
u_j^* - (u_i^* + \frac{\partial u_i^*}{\partial x} * dx + \frac{\partial u_i^*}{\partial y} * dy),
v_j^* - (v_i^* + \frac{\partial v_i^*}{\partial x} * dx + \frac{\partial v_i^*}{\partial y} * dy)
]
$$
方向导数和之前的相同
$$
q_{ij}=
[
n_x,
n_y,
\frac{2dx*n_x}{r_s},
\frac{2dy*n_y}{r_s},
\frac{dx*n_y + dy*n_x}{r_s}
]
$$
得到间断的大小
$V_{jump}=q_{ij} \cdot \Delta h_{ij}$
权重函数
$$
w_{ij}^{vel} = 
\begin{cases}
(1-\frac{|d_{ij}|}{0.6l_0})^2, |d_{ij}| < 0.6l_0 \\
0, |d_{ij}| >= 1.2l_0
\end{cases}
$$

3、计算根据速度突变构建的矩矩阵
$$
M_i^{vel} = \sum_{j \in fluid} w_{ij}^{vel} q_{ij}^{'}{q_{ij}^{'}}^T + \sum_{j \in wall} w_{ij}^{vel} q_{ij}{q_{ij}}^T
$$


最终的lsmps矩阵为：  
$$C_i = (M_i + M_i^{vel})^{-1}$$

4、计算压力梯度
$$
\nabla p_i = \sum_{j \in fluid} w_{ij}
\begin{bmatrix}
C_1\\ 
C_2
\end{bmatrix}
p_{ij}\frac{p_j-p_i}{r_{ij}}
+ \sum_{j \in wall} w_{ij}
\begin{bmatrix}
C_1\\ 
C_2
\end{bmatrix}
q_{ij}\rho \mathbf{n}_j \mathbf{g}
+ \sum_{j \in fluid} w_{ij}^{vel}
\begin{bmatrix}
C_1\\ 
C_2
\end{bmatrix}
q_{ij}^{'}
(\frac{p_j-p_i}{r_{ij}}-\frac{\rho V_{jump}}{2\Delta t})
+ \sum_{j \in wall} w_{ij}
\begin{bmatrix}
C_1\\ 
C_2
\end{bmatrix}
q_{ij}(\rho \mathbf{n}_j \mathbf{g} - \frac{\rho V_{jump}}{2\Delta t})
$$