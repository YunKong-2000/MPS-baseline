# PPE方程组构建

## PPE方程

$$
\frac{1}{\rho} \langle \Delta p \rangle_i = \frac{1}{\Delta t} \nabla \cdot \mathbf{u}^k
$$

## 拉普拉斯算子离散

### 远离壁面的内部粒子

$$
\begin{align*}
    \langle \Delta p \rangle_i &= \frac{2}{r_e} \sum_{j \neq i} w_{ij} d_{ij} [\mathbf{C}_3 + \mathbf{C}_4] \mathbf{P}_{ij} \\
    &= \frac{2}{r_e} \sum_{j \neq i} w_{ij} \frac{p_j - p_i}{r_{ij}} [\mathbf{C}_3 + \mathbf{C}_4] \mathbf{P}_{ij} \\
    &= \left(-\frac{2}{r_e} \sum_{j \neq i} \frac{w_{ij}}{r_{ij}} [\mathbf{C}_3 + \mathbf{C}_4] \mathbf{P}_{ij}\right) p_i + \sum_{j \neq i} \left(\frac{2}{r_e} \frac{w_{ij}}{r_{ij}} [\mathbf{C}_3 + \mathbf{C}_4] \mathbf{P}_{ij}\right) p_j
\end{align*}
$$

### 近壁面内部粒子

$$
\begin{align*}
    \langle \Delta p \rangle_i &= \frac{2}{r_e} \sum_{j \in \mathrm{fluid}} w_{ij} d_{ij} [\mathbf{C}_3 + \mathbf{C}_4] \mathbf{P}_{ij} + \frac{2}{r_e} \sum_{j \in \mathrm{wall}} w_{ij} d_{ij} [\mathbf{C}_3 + \mathbf{C}_4] \mathbf{P}_{ij} \\
    &= \frac{2}{r_e} \sum_{j \in \mathrm{fluid}} w_{ij} \frac{p_j - p_i}{r_{ij}} [\mathbf{C}_3 + \mathbf{C}_4] \mathbf{P}_{ij} + \frac{2}{r_e} \sum_{j \in \mathrm{wall}} w_{ij} (-\rho \mathbf{n}_{\mathrm{wall}} \cdot \mathbf{g}) [\mathbf{C}_3 + \mathbf{C}_4] \mathbf{P}_{ij} \\
    &= \left(-\frac{2}{r_e} \sum_{j \in \mathrm{fluid}} \frac{w_{ij}}{r_{ij}} [\mathbf{C}_3 + \mathbf{C}_4] \mathbf{P}_{ij}\right) p_i + \sum_{j \in \mathrm{fluid}} \left(\frac{2}{r_e} \frac{w_{ij}}{r_{ij}} [\mathbf{C}_3 + \mathbf{C}_4] \mathbf{P}_{ij}\right) p_j \\
    &\quad - \frac{2}{r_e} \sum_{j \in \mathrm{wall}} w_{ij} (\rho \mathbf{n}_{\mathrm{wall}} \cdot \mathbf{g}) [\mathbf{C}_3 + \mathbf{C}_4] \mathbf{P}_{ij}
\end{align*}
$$

### 自由面粒子

$$
\begin{align*}
    \langle \Delta p \rangle_i &= \frac{4}{n_0 \lambda} \sum_{j \in \mathrm{fluid}} (p_j - p_i)w_{ij} - \frac{4}{n_0 \lambda}\left(n_{0} - n_{i}^{'}\right)p_i
\end{align*}
$$

其中：

$$
    n_{i}^{'} = \min(n_{i}^{*}, n_{0})
$$

其中$n_{i}^{*}$为当前粒子的粒子数密度，考虑壁面粒子

## 速度散度离散

### 远离壁面内部粒子

$$
\begin{align*}
    \nabla \cdot \mathbf{u} &= \sum_{j \neq i} w_{ij} d_{ij} \begin{bmatrix} \mathbf{C}_1 \\ \mathbf{C}_2 \end{bmatrix} \mathbf{P}_{ij} \\
    &= \sum_{j \neq i} w_{ij} \frac{\mathbf{u}_j - \mathbf{u}_i}{r_{ij}} \begin{bmatrix} \mathbf{C}_1 \\ \mathbf{C}_2 \end{bmatrix} \mathbf{P}_{ij}
\end{align*}
$$

### 近壁面内部粒子

$$
\begin{align*}
    \nabla \cdot \mathbf{u} &= \sum_{j \in \mathrm{fluid}} w_{ij} d_{ij} \begin{bmatrix} \mathbf{C}_1 \\ \mathbf{C}_2 \end{bmatrix} \mathbf{P}_{ij} + \sum_{j \in \mathrm{wall}} w_{ij} d_{ij} \begin{bmatrix} \mathbf{C}_1 \\ \mathbf{C}_2 \end{bmatrix} \mathbf{P}_{ij} \\
    &= \sum_{j \in \mathrm{fluid}} w_{ij} \frac{\mathbf{u}_j - \mathbf{u}_i}{r_{ij}} \begin{bmatrix} \mathbf{C}_1 \\ \mathbf{C}_2 \end{bmatrix} \mathbf{P}_{ij} + \sum_{j \in \mathrm{wall}} w_{ij} \frac{\mathbf{u}_{\mathrm{wall}} - \mathbf{u}_i}{r_{ij}} \begin{bmatrix} \mathbf{C}_1 \\ \mathbf{C}_2 \end{bmatrix} \mathbf{P}_{ij}
\end{align*}
$$

### 自由面粒子

$$
\begin{align*}
    \langle \nabla \cdot \mathbf{u} \rangle_i &= \frac{2}{n_0} \sum_{j \neq i} \left[\frac{\mathbf{u}_j - \mathbf{u}_i}{r_{ij}} \cdot \frac{\mathbf{r}_{ij}}{r_{ij}} w_{ij}\right]
\end{align*}
$$

## 代数式方程

### 远离壁面内部粒子

$$
\left(-\frac{2}{r_e \rho} \sum_{j \neq i} \frac{w_{ij}}{r_{ij}} [\mathbf{C}_3 + \mathbf{C}_4] \mathbf{P}_{ij}\right) p_i + \sum_{j \neq i} \left(\frac{2}{r_e \rho} \frac{w_{ij}}{r_{ij}} [\mathbf{C}_3 + \mathbf{C}_4] \mathbf{P}_{ij}\right) p_j = \frac{1}{\Delta t} \sum_{j \neq i} w_{ij} \frac{\mathbf{u}_j - \mathbf{u}_i}{r_{ij}} \begin{bmatrix} \mathbf{C}_1 \\ \mathbf{C}_2 \end{bmatrix} \mathbf{P}_{ij}
$$

### 近壁面内部粒子

$$
\begin{align*}
    &\left(-\frac{2}{r_e \rho} \sum_{j \in \mathrm{fluid}} \frac{w_{ij}}{r_{ij}} [\mathbf{C}_3 + \mathbf{C}_4] \mathbf{P}_{ij}\right) p_i + \sum_{j \in \mathrm{fluid}} \left(\frac{2}{r_e \rho} \frac{w_{ij}}{r_{ij}} [\mathbf{C}_3 + \mathbf{C}_4] \mathbf{P}_{ij}\right) p_j \\
    &= \frac{1}{\Delta t} \sum_{j \in \mathrm{fluid}} w_{ij} \frac{\mathbf{u}_j - \mathbf{u}_i}{r_{ij}} \begin{bmatrix} \mathbf{C}_1 \\ \mathbf{C}_2 \end{bmatrix} \mathbf{P}_{ij} \\
    &\quad + \frac{1}{\Delta t} \sum_{j \in \mathrm{wall}} w_{ij} \frac{\mathbf{u}_{\mathrm{wall}} - \mathbf{u}_i}{r_{ij}} \begin{bmatrix} \mathbf{C}_1 \\ \mathbf{C}_2 \end{bmatrix} \mathbf{P}_{ij} \\
    &\quad + \frac{2}{r_e \rho} \sum_{j \in \mathrm{wall}} w_{ij} (\rho \mathbf{n}_{\mathrm{wall}} \cdot \mathbf{g}) [\mathbf{C}_3 + \mathbf{C}_4] \mathbf{P}_{ij}
\end{align*}
$$

### 自由面粒子

$$
\begin{align*}
    &-\frac{4}{n_0 \lambda \rho}\left( \sum_{j \neq i} w_{ij} + \left(n_{i}^{'} - n_{i}^{*}\right)\right) p_i + \sum_{j \neq i} \left(\frac{4}{n_0 \lambda \rho} w_{ij}\right) p_j \\
    &\quad = \frac{2}{n_0 \Delta t} \sum_{j \neq i} \left[\frac{\mathbf{u}_j - \mathbf{u}_i}{r_{ij}} \cdot \frac{\mathbf{r}_{ij}}{r_{ij}} w_{ij}\right]
\end{align*}
$$
