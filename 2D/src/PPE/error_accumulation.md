# 误差累计原因分析  
## bias vector计算
$$
\mathbf{b}_{i} = \sum_{j \neq i}\frac{w_{ij}}{n_0} (
\begin{bmatrix}
C_1 \\
C_2 \\
\end{bmatrix}
P)
$$

## Laplacian-model coefficients 计算
$$
c_{i} = 2\sum_{j \neq i}\frac{w_{ij}}{n_0} (
\begin{bmatrix}
C_3 + C_4 \\
\end{bmatrix}
P)
$$

## Ri计算
$$
R_i = \frac{(\| \mathbf{b}_i \|)^2}{c_i}
$$