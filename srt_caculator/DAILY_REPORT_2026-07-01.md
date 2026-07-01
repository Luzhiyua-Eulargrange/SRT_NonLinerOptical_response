# 工作日报

日期：2026年7月1日

## 一句话总结

今天在 `dev` 主线中完成了局域电流分布计算、粗细 k 网格收敛诊断和密度/电流可视化输出的集成，同时确认当前默认参数下密度收敛良好，但电流分布仍存在明显 k 网格不收敛和长度规范响应定义待确认的问题。

## 今日工作

1. 阅读并检查当前 `dev` 项目结构，确认主线已经收敛到长度规范 RDM、实空间密度和后处理诊断流程。
2. 检查 `VELOCITY_MATRIX_DIAGNOSTICS.md` 的速度矩阵推导，确认其主要公式是本征方程求导得到的恒等式，不应改写成带能量分母的微扰论形式。
3. 仔细检查 `Theory_Formula_of_Current_Space_Distribution.md` 中的局域电流推导，指出其中 `B` 层周期势符号、显式 `hbar` 约定、Fourier 归一化和长度规范代入等需要继续修正的问题。
4. 新增 `Current_Distribution.py`，从 band-basis RDM 经 plane-wave basis 重构晶胞内局域电流分布。
5. 将局域电流分布接入 `main_length_gauge.py` 主流程，默认生成 `length_current_distribution.npz`。
6. 为电流结果增加粗、细 k 网格收敛诊断，默认生成 `length_current_convergence.npz`。
7. 维护 `README.md`，移除已不存在的速度规范主线描述，补充局域电流分布、输出文件和当前已知问题。
8. 扩展 `Debug_Tools.py`，新增从 `.npz` 输出生成密度和电流静态图的工具函数。
9. 将绘图集成进主流程，目前默认输出两张全量时空热图和两组响应量图。

## 新增和修改的主要文件

1. 新增 `Current_Distribution.py`。
2. 修改 `main_length_gauge.py`，集成电流计算、电流收敛诊断和绘图输出。
3. 修改 `Debug_Tools.py`，新增 `.npz` 读取、密度/电流时空图、时间切片图和空间采样点时间曲线图。
4. 修改 `README.md`，同步当前项目结构、数据流、输出文件和已知问题。

## 当前主流程输出

运行：

```bash
python main.py
```

会生成或更新：

1. `length_band_result.npz`
2. `length_band_structure.png`
3. `length_rdm_result.npz`
4. `length_rho_convergence.npz`
5. `length_density_result.npz`
6. `length_density_convergence.npz`
7. `length_current_distribution.npz`
8. `length_current_convergence.npz`

当前默认图像输出包括：

1. `density_spacetime_total.png`
2. `current_spacetime_total.png`
3. `density_response_spacetime_total.png`
4. `density_response_time_slices.png`
5. `density_response_time_traces.png`
6. `current_response_spacetime_total.png`
7. `current_response_time_slices.png`
8. `current_response_time_traces.png`

其中响应量定义为：

```text
delta density = density(t) - density(0)
delta current = current(t) - current(0)
```

## 数值检查结果

当前默认参数：

1. `L = 2`
2. `num_k_band = 51`
3. `num_k_rdm = 21`
4. 细网格 `Nk = 42`
5. `num_time_steps = 81`
6. 实空间网格 `Nx = 201`

密度结果保持良好收敛：

1. `rho_relative_max_difference` 约为 `3.79e-5`
2. `diagonal_relative_max_difference` 约为 `1.28e-9`
3. `density_total_relative_max_difference` 约为 `4.45e-9`
4. `density_component_relative_max_difference` 约为 `7.31e-9`
5. 每晶胞电荷差为机器精度量级

电流模块通过了内部一致性检查：

1. 电流计算直接使用 plane-wave basis 中的解析矩阵元。
2. 没有使用微扰论、能量分母或本征矢的 k 导数。
3. 晶胞平均电流与解析 `dH/dk` 宏观电流在机器精度内一致，最大差约为 `6.9e-15`。

但当前默认粗、细 k 网格下电流分布不收敛：

1. `current_total_relative_max_difference` 约为 `0.999997`
2. `average_current_relative_max_difference` 约为 `1.000002`
3. 静态初态平均电流随 `Nk` 增大近似按 `1/Nk` 衰减，说明当前默认 k 网格对绝对电流分布明显不够。

## 已确认的问题

1. 局域电流分布虽然已经可以计算和保存，但默认 `Nk=21/42` 下没有收敛，结果只能作为当前公式和代码路径下的可能输出。
2. 长度规范下局域电流响应是否还需要额外极化项、规范变换项或更严格的响应定义，仍需从理论上继续确认。
3. `Theory_Formula_of_Current_Space_Distribution.md` 中仍需系统修正 `B` 层势能符号、`hbar` 约定和 Fourier 归一化。
4. 速度矩阵的 Berry 联络有限差分形式仍受近简并子空间、band tracking 和 BZ 边界 sewing matrix 影响，不应作为生产路径。
5. 当前密度和电流的时间演化响应量非常小，绝对值图中几乎不可见，因此需要优先查看响应量图。

## 后续建议

1. 系统扫描 `Nk`、`L`、`Nx` 和时间积分容差，建立密度与电流的收敛测试矩阵。
2. 优先解决电流中静态背景项和 BZ 边界贡献的来源，检查是否需要边界 sewing 或规范一致的电流定义。
3. 修订电流理论文档，使公式、代码 Hamiltonian、归一化和长度规范响应定义完全一致。
4. 在确认理论定义后，再考虑把电流结果用于物理结论或论文图。
