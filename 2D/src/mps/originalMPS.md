# 原始MPS方法算子离散实现
## 目标
对于模拟中的自由面或者近自由面粒子，很容易出现因为邻域被自由面截断导致的邻域粒子缺失问题，使用lsmps方法计算矩矩阵时容易因为原矩阵病态导致求逆后的矩阵失真导致自由面粒子或者近自由面粒子出现了突变物理量。因此需要在自由面或近自由面处使用原始MPS方法保证程序的稳定
## 离散方法
对于梯度算子，基本的离散实现为
$$
\nabla \phi_i = \frac{d}{n_0}\sum_{j \neq i}\frac{\phi_j - \phi_i}{|\mathbf{r}_{ij}|^2}\mathbf{r}_{ij}w_{ij}
$$
散度算子类似
$$
\nabla \cdot \mathbf{\phi}_i = \frac{d}{n_0}\sum_{j \neq i}\frac{(\mathbf{\phi}_j - \mathbf{\phi}_i) \cdot \mathbf{r}_{ij}}{|\mathbf{r}_{ij}|^2}w_{ij}
$$
拉普拉斯算子
$$
\Delta \phi = \frac{2d}{n_0\lambda}\sum_{j \neq i}(\phi_j - \phi_i)w_{ij}
$$

## 具体实现
1、压力拉普拉斯算子
对于自由面或者近自由面粒子的压力拉普拉斯算子，
$$
\Delta p = \frac{2d}{n_0\lambda}\sum_{j \neq i}(p_j - p_i)w_{ij}
$$
2、速度散度算子
对于自由面或近自由面粒子的散度算子
$$
\nabla \cdot \mathbf{u}_i = \frac{d}{n_0}\sum_{j \neq i}\frac{(\mathbf{u}_j - \mathbf{u}_i) \cdot \mathbf{r}_{ij}}{|\mathbf{r}_{ij}|^2}w_{ij}
$$
其中，当j粒子是壁面粒子时，壁面粒子的速度也需要加上$\Delta t \mathbf{g}$，这一点和使用lsmps方法时相同。
3、压力梯度算子
计算压力梯度算子时使用稳定的梯度算法
$$
\nabla p_i = \frac{d}{n_0}\sum_{j \neq i}\frac{p_j - p_{min}}{|\mathbf{r}_{ij}|^2}\mathbf{r}_{ij}w_{ij}
$$