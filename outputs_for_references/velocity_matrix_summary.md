# H2 速度矩阵元

## 参数

- `L`: 4
- `N_G`: 9
- 矩阵维度: `18 x 18`
- `nk`: 181

## 输出

- `h2_intraband_velocities.png`：近零能带的对角速度矩阵元 `v_nn(k)`。
- `h2_interband_velocity.png`：能带 `8` 到 `9` 的能隙和 `|v_nm(k)|`。

## 检查

- Band-basis velocity Hermitian error: `3.553e-15`

## 说明

`v_nn(k)` 是带内速度；`v_nm(k), n != m` 是带间速度矩阵元，是后续光学响应的输入之一。
