# 不使用可变截断半径的Particle shifting算法

## 1、计算自由面粒子法向向量
1、计算自由面粒子的可变截断半径，大小为 $r_{ve}=\max(1.2l_0, d_{min})$，$d_{min}$ 为到壁面的最近距离;计算其他粒子的可变截断半径，大小为 $r_{ve}=\max(1.2l_0, d_{min})$，$d_{min}$ 为到壁面或者自由面的最近距离。若邻域中没有自由面或者壁面，则$r_{ve}=r_e$。  
2、计算自由面法向向量：
$$
\mathbf{n}_i = - \frac{\sum_{j,r_{ve}} \left( \frac{\mathbf{r}_{ij}}{|\mathbf{r}_{ij}|}w_{ij} \right)}{\left| \sum_{j,r_{ve}} \left( \frac{\mathbf{r}_{ij}}{|\mathbf{r}_{ij}|}w_{ij} \right) \right|}
$$

## 2、计算位移矢量
1、计算基本位移矢量：
$$
f_{i}^{PS} = -\frac{\lambda_{shift}dl_0^2}{n_0}\sum_{j,r_{ve}}w_{ij,PS}\frac{\mathbf{r}_{ij}}{|\mathbf{r}_{ij}|^2}
$$
其中 $\lambda_{shift}$ 为位移系数（建议设为0.01~0.1），$d$ 为空间维度。$r_{ve}$ 为PS算法使用的可变截断半径。PS算法专用的核函数为：
$$
w_{ij,PS}=
\begin{cases}
\frac{r_{ve}}{|\mathbf{r}_{ij}|} - 1, & |\mathbf{r}_{ij}| < r_{ve} \\
0, & |\mathbf{r}_{ij}| \ge r_{ve}
\end{cases}
$$
计算位移矢量时需要考虑壁面粒子。
2、计算滑移矢量：
$$
\hat{f}_{i}^{PS}=\sum_{j,r_{ve}}\frac{1}{2}\tau_{ij}\sqrt{\max(0, r_c^2-|\mathbf{r}_{ij}|^2)}
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

3、合成最终位移矢量：
根据粒子类型应用相应的约束：
$$
\delta \mathbf{r}_{i}^{PS} = 
\begin{cases}
f_{i}^{PS} + \hat{f}_{i}^{PS}, & i \in \text{Inner} \\
\frac{d-1}{d}(\mathbf{I}_d - \mathbf{n}_i\mathbf{n}_i^{T})(f_{i}^{PS} + \hat{f}_{i}^{PS}), & i \in \text{Free surface} \\
\mathbf{0}, & i \in \text{Splash}
\end{cases}
$$

4、最大位移限制 (Magnitude Limiter)：
防止单步位移矢量过大导致流场畸变，必须对位移矢量进行幅值截断：
$$
\delta \mathbf{r}_i^{PS} \leftarrow \min(0.1l_0, |\delta \mathbf{r}_i^{PS}|)\frac{\delta \mathbf{r}_i^{PS}}{|\delta \mathbf{r}_i^{PS}|}
$$

5、壁面防穿透保护 (Wall Collision Limiter)：
预测粒子新位置 $\mathbf{r}_i + \delta \mathbf{r}_i^{PS}$。遍历邻域内的壁面粒子找出离壁面最近距离，如果该距离小于安全阈值 $d_{safe}$（$0.25l_0$）或发生穿透，则必须剃除 $\delta \mathbf{r}_i^{PS}$ 中指向壁面内部的法向分量，将其强制转化为沿壁面的纯切向滑移。

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