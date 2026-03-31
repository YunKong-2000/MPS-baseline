# 可变截断半径的PS算法
## 1、计算自由面法向向量
a、计算流体粒子的可变截断半径，大小为 $r_{ve}=\max(1.2l_0, d_{min})$，$d_{min}$ 为到壁面的最近距离。若邻域中没有壁面粒子，则$r_{ve}=r_e$。  
b、计算自由面法向向量：
$$
\mathbf{n}_i = - \frac{\sum_{j,r_{ve}} \left( \frac{\mathbf{r}_{ij}}{|\mathbf{r}_{ij}|}w_{ij} \right)}{\left| \sum_{j,r_{ve}} \left( \frac{\mathbf{r}_{ij}}{|\mathbf{r}_{ij}|}w_{ij} \right) \right|}
$$

## 2、计算位移矢量
a、计算基于可变截断半径的基本位移矢量：
$$
\delta r_{i}^{OPS} = -\frac{\lambda_{shift}dl_0^2}{n_0}\sum_{j,r_{ve}}w_{ij,PS}\frac{\mathbf{r}_{ij}}{|\mathbf{r}_{ij}|^2}
$$
其中 $\lambda_{shift}$ 为位移系数（设为0.1），$d$ 为空间维度。$r_{ve}$ 为PS算法使用的可变截断半径。PS算法专用的核函数为：
$$
w_{ij,PS}=
\begin{cases}
\frac{r_{ve}}{|\mathbf{r}_{ij}|} - 1, & |\mathbf{r}_{ij}| < r_{ve} \\
0, & |\mathbf{r}_{ij}| \ge r_{ve}
\end{cases}
$$

b、计算滑移矢量：
$$
\r_{i}^{SL}=\sum_{j,r_{ve}}\frac{1}{2}\tau_{ij}\sqrt{\max(0, r_c^2-|\mathbf{r}_{ij}|^2)}
$$
$r_c$为粒子最近距离限制，设置为$0.7*l_0$.
其中切向方向矢量 $\tau_{ij}$ 为粒子间相对速度的切向投影：
$$
\tau_{ij} = \frac{T_{ij}(\mathbf{u}_i - \mathbf{u}_j)}{|T_{ij}(\mathbf{u}_i - \mathbf{u}_j)|}
$$
$$
T_{ij} = I_d - \left(\frac{\mathbf{r}_{ij}}{|\mathbf{r}_{ij}|}\right)\left(\frac{\mathbf{r}_{ij}}{|\mathbf{r}_{ij}|}\right)^T
$$
*注意：计算 $\tau_{ij}$ 时 $|T_{ij}(\mathbf{u}_i - \mathbf{u}_j)|$ 可能为零，代码实现需增加除零保护（若分母为0则 $\tau_{ij}=\mathbf{0}$）*。

c、计算总位移矢量
$$
\delta \mathbf{r}_{i}^{PS} = 
\begin{cases}
r_{i}^{OPS} + r_{i}^{SL}, & i \in \text{Inner} \\
\frac{d-1}{d}(\mathbf{I}_d - \mathbf{n}_i\mathbf{n}_i^{T})(r_{i}^{OPS} + r_{i}^{SL}), & i \in \text{Free surface} \\
\mathbf{0}, & i \in \text{Splash}
\end{cases}
$$

d、最大位移限制 (Magnitude Limiter)：
防止单步位移矢量过大导致流场畸变，必须对位移矢量进行幅值截断：
$$
\delta \mathbf{r}_i^{PS} \leftarrow \min(0.1l_0, |\delta \mathbf{r}_i^{PS}|)\frac{\delta \mathbf{r}_i^{PS}}{|\delta \mathbf{r}_i^{PS}|}
$$

## 3、一阶泰勒展开修正及更新
1、提取在动量方程求解步骤中已经计算好的**速度梯度张量（矩阵）** $\nabla \mathbf{u}_i^*$。  
2、使用一阶泰勒展开对速度进行修正：
$$
\mathbf{u}_i = \mathbf{u}_i^* + \nabla \mathbf{u}_i^* \cdot \delta \mathbf{r}_i^{PS}
$$
其中 $\mathbf{u}_i^*$ 为压力梯度修正后的粒子临时速度。  
3、位置更新：
$$
\mathbf{r}_i = \mathbf{r}_i + \delta \mathbf{r}_i^{PS}
$$