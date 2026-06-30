# H2 一阶密度矩阵响应

## 模型参数

- `L`: 4
- `N_G`: 9
- 矩阵维度: `18 x 18`
- `nk`: 181

## 响应参数

- `E_omega`: 0.03
- `omega`: 0.35
- `gamma_inter`: 0.04
- `gamma_intra`: 0.04
- `mu`: 0.0
- `beta`: 30.0

## 输出

- `h2_rho1_intraband.png`：一阶对角密度矩阵 `|rho^(1)_nn(k)|`。
- `h2_rho1_interband.png`：能带 `8` 到 `9` 的能隙、`|A_nm|` 与 `|rho^(1)_nm|`。
- `h2_rho1_matrix_at_k0.png`：`k≈0` 处的一阶密度矩阵模长热图。

## 数值检查

- max diagonal `|rho^(1)_nn|`: `3.054491e-03`
- max off-diagonal `|rho^(1)_nm|`: `2.775085e-03`

## 说明

本脚本使用频域一阶响应和简单弛豫时间近似。它是教学用最小模型，
尚未处理完整规范协变导数、二阶响应、真实材料参数或电流/极化张量积分。
