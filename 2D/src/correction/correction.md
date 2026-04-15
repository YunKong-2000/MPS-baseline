# 速度和位置更新
## 整体流程
1、临时速度$\mathbf{u}^{*}$的计算  
$$
\mathbf{u}^{*} = \mathbf{u}^k + \Delta t(\mathbf{g} + \nabla^2 \mathbf{u}^k)
$$
其中$\mathbf{u}^k$为当前时刻的速度。

2、计算完压力梯度后，使用压力梯度对临时速度$u^*$进行修正。
$$
\mathbf{u}^{**} = \mathbf{u}^* - \Delta t\frac{\nabla P}{\rho}
$$
同时需要计算速度梯度$\nabla \mathbf{u}^{**}$

3、更新位置  
更新位置时考虑位移过程的两个时刻，开始时刻的速度是$\mathbf{u}^k$，结束时刻的速度是$\mathbf{u}^{**}$。将整个位移过程看作匀加速度运动，所以位移$\mathbf{\Delta r}$为，
$$
\Delta r = \frac{\Delta t}{2}(\mathbf{u}^{k} + \mathbf{u}^{**})+\delta \mathbf{r}
$$
其中$\delta r$是PS算法计算出的shifting位移。

4、最终更新速度
速度更新前，添加速度平滑处理，每个粒子速度更新时考虑周围粒子的平均速度。
计算粒子邻域内邻域粒子的平均速度$\mathbf{\hat{u}}$，
$$
\mathbf{\hat{u}} = \frac{\sum_{j \neq i}\mathbf{u}_jw_{ij}}{\sum_{j \neq i}w_{ij}}
$$
计算完$\mathbf{\hat{u}}$再更新下一时间步速度
$$
\mathbf{u}^{k+1} = (1 - \lambda)\mathbf{u}^{**} + \lambda\mathbf{\hat{u}}
$$
其中$\lambda$为平滑系数，设置为0.1

## 对应程序模块的实现
1、使用当前速度$\mathbf{u}^k$求解临时速度$\mathbf{u}^*$  
2、correction中应当首先完成临时速度$\mathbf{u}^{**}$及其速度梯度的计算，由于correction之前是没有更新位置的，所以可以直接使用PPEbuild模块中的速度对应的lsmps矩矩阵进行梯度的计算。  
3、速度$\mathbf{u}^{**}$计算完后，调用PS算法，计算出shifting位移矢量$\delta \mathbf{r}$。  
4、correction模块此时再去计算整体位移$\mathbf{\Delta r}$，并完成位置的更新。  
5、最后correction模块通过计算下一时间步的速度$\mathbf{u}^{k+1}$，完成速度的更新。

