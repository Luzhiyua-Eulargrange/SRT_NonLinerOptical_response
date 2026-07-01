# SRT RDM Solver

本项目用于计算一维双组分连续模型在外加脉冲电场下的能带结构、约化密度矩阵
（Reduced Density Matrix, RDM）动力学，以及由 RDM 重构得到的晶胞内实空间密度
和晶胞内局域电流分布。

当前主线开发已经从速度规范电流收敛排查，转向以长度规范 RDM、实空间密度收敛
和最小耦合局域电流分布收敛为核心的验证流程。速度规范电流链条没有作为 `dev`
主线代码保留；速度矩阵的 Berry 联络有限差分形式仍只适合作为理论和数值诊断，
不应作为生产路径。

## 项目结构

- `main.py`：默认入口，转发到 `main_length_gauge.py`。
- `main_length_gauge.py`：长度规范 RDM、实空间密度重构和收敛诊断主流程。
- `config.py`：默认参数、参数归一化与校验、第一布里渊区 k 网格生成。
- `Band_Solver.py`：构造一维双组分连续模型哈密顿量，计算能带和本征矢，并进行相位平滑。
- `Geometry.py`：计算 Berry 联络和协变导数。
- `RDM_Common.py`：RDM 求解器共用的外场、初态密度矩阵、时间网格和 ODE 积分工具。
- `RDM_Length_Gauge.py`：长度规范 RDM 演化模块，对应 Phys. Rev. B 96, 035431 的式 (55)。
- `Density.py`：从 band-basis RDM 重构晶胞内实空间密度 `n_A(x,t)`、`n_B(x,t)` 和总密度。
- `Current_Distribution.py`：从 band-basis RDM 重构晶胞内局域电流 `j_A(x,t)`、`j_B(x,t)` 和总电流。
- `Debug_Tools.py`：参数打印、能带绘图、RDM 轨迹检查和 `.npz` 可视化工具。
- `Velocity_Matrix_Diagnostics.py`：历史速度矩阵诊断脚本；当前不属于主线流程。
- `outputs_for_references/`：H2 最小实验的参考图片和 Markdown 摘要，用于对照和后续整理。

## 环境依赖

建议使用 Python 3.10 或更高版本。

```bash
pip install -r requirements.txt
```

依赖列表维护在 `requirements.txt` 中：

- `numpy`：矩阵构造、线性代数、数组保存与读取。
- `scipy`：调用 `scipy.integrate.solve_ivp` 进行 RDM 时间演化；如果未安装 `scipy`，会回退到内置四阶 Runge-Kutta 方法。
- `matplotlib`：用于能带、密度和电流可视化。

## 快速运行

在项目目录下运行：

```bash
python main.py
```

当前默认流程会执行长度规范计算：

1. 使用 `num_k_band` 网格计算并保存能带。
2. 在粗、细两个嵌套 k 网格上求解长度规范 RDM。
3. 比较 RDM 轨迹的收敛性。
4. 从 band-basis RDM 重构晶胞内实空间密度。
5. 比较实空间密度和每晶胞电荷的收敛性。
6. 从 band-basis RDM 重构晶胞内局域电流分布。
7. 比较局域电流、分量电流和晶胞平均电流的收敛性。

局域电流分布也可以作为独立后处理脚本运行：

```bash
python Current_Distribution.py
```

该脚本会读取 `length_rdm_result.npz`，并重新生成 `length_current_distribution.npz`。

默认输出文件包括：

- `length_band_result.npz`：长度规范流程使用的能带数据。
- `length_band_structure.png`：能带图。
- `length_rdm_result.npz`：细 k 网格上的长度规范 RDM 结果。
- `length_rho_convergence.npz`：粗、细 k 网格 RDM 收敛诊断。
- `length_density_result.npz`：晶胞内实空间密度结果。
- `length_density_convergence.npz`：实空间密度收敛诊断。
- `length_current_distribution.npz`：晶胞内局域电流分布结果。
- `length_current_convergence.npz`：局域电流分布收敛诊断。

默认图像输出包括：

- `length_band_structure.png`：能带图。
- `density_spacetime_total.png`：全量总粒子密度 `n(x,t)` 时空热图。
- `current_spacetime_total.png`：全量总局域电流 `j(x,t)` 时空热图。
- `density_response_spacetime_total.png`：密度响应量 `delta n(x,t)=n(x,t)-n(x,0)` 时空热图。
- `density_response_time_slices.png`：密度响应量的若干时间切片 `delta n(x,t_i)`。
- `density_response_time_traces.png`：密度响应量在若干空间采样点的时间曲线 `delta n(x_i,t)`。
- `current_response_spacetime_total.png`：电流响应量 `delta j(x,t)=j(x,t)-j(x,0)` 时空热图。
- `current_response_time_slices.png`：电流响应量的若干时间切片 `delta j(x,t_i)`。
- `current_response_time_traces.png`：电流响应量在若干空间采样点的时间曲线 `delta j(x_i,t)`。

绝对值图只保留密度和电流各一张时空热图；时间切片和空间采样点时间曲线默认只对响应量输出。

## 当前默认结果

以当前 `main_length_gauge.py` 中的测试参数为例：

- `L = 2`
- `num_k_band = 51`
- `num_k_rdm = 21`
- 细网格为 `2 * num_k_rdm = 42`
- `num_time_steps = 81`
- 实空间密度网格 `Nx = 201`

已有一次运行结果显示：

- RDM 相对最大差约为 `3.79e-5`
- 对角占据相对最大差约为 `1.28e-9`
- 实空间总密度相对最大差约为 `4.45e-9`
- 实空间分量密度相对最大差约为 `7.31e-9`
- 每晶胞电荷差为机器精度量级
- 局域电流分布的晶胞平均电流与解析 `dH/dk` 宏观电流一致，最大差约为 `6.9e-15`
- 当前默认 `Nk=21/42` 下局域电流粗、细网格尚未收敛，`current_total_relative_max_difference` 约为 `1.0`

这些结果说明当前参数下实空间密度收敛很好。不过，由于 Berry 联络和能带规范处理仍在诊断中，
仍需继续检查近简并带、band tracking 和 BZ 边界 sewing matrix 的影响。局域电流分布模块
已经与当前平面波基和 k 积分归一化做过一致性检查，但默认 k 网格对绝对电流分布明显不够；
长度规范下是否需要额外极化项或规范变换项也仍需从理论上继续确认。

当前外场参数下绝对密度和绝对电流中的时间演化很小。典型量级为：

- `density_total` 约为 `5`，而 `max |density_total(t)-density_total(0)|` 约为 `9.6e-8`。
- `current_total` 约为 `-0.374`，而 `max |current_total(t)-current_total(0)|` 约为 `2.4e-7`。

因此默认图像同时保存全量热图和响应量图。判断时间演化时应优先查看响应量图。

## 主要参数

参数通过字典传入 `normalize_params`，未提供的参数会使用 `config.py` 中的默认值。

常用参数包括：

- `L`：倒格矢截断参数；总基维度 `Nb = 2 * (2L + 1)`。
- `num_k_band`：能带绘图和能带结果保存使用的 k 网格点数。
- `num_k_rdm`：RDM、Berry 联络和协变导数使用的 k 网格点数。
- `a0`：晶格常数，倒格矢大小为 `b = 2*pi/a0`。
- `m`、`hbar`：有效质量与约化普朗克常数。
- `kappa`、`v0`、`vA`、`vB`、`w1`、`w2`：模型哈密顿量参数。
- `E0`、`omega`、`pulse_duration`、`t_switch`：外加脉冲电场的强度、频率、持续时间和开关时间。
- `num_time_steps`：时间演化采样点数。
- `fermi_energy`：零温初态占据使用的费米能级。
- `rtol`、`atol`：`solve_ivp` 的相对和绝对误差容限。

示例：

```python
params = normalize_params({
    "L": 2,
    "E0": 0.02,
    "omega": 0.5,
    "pulse_duration": 20.0,
    "t_switch": 4.0,
    "num_time_steps": 81,
    "num_k_band": 51,
    "num_k_rdm": 21,
    "fermi_energy": 0.0,
    "rtol": 1e-6,
    "atol": 1e-8,
})
```

## 数据文件说明

`length_rdm_result.npz` 主要包含：

- `gauge`：规范名称，当前为 `length`。
- `k_grid`：RDM 使用的 k 网格。
- `k_weight`：k 空间积分权重。
- `time_grid`：时间网格。
- `rho_trajectory`：RDM 轨迹，形状为 `(Nt, Nk, Nb, Nb)`。
- `energies`：能带能量，形状为 `(Nk, Nb)`。
- `eigenvectors`：本征矢矩阵，形状为 `(Nk, Nb, Nb)`。

`length_density_result.npz` 主要包含：

- `x_grid`：晶胞内实空间网格。
- `time_grid`：时间网格。
- `density_components`：A/B 分量密度，形状为 `(Nt, 2, Nx)`。
- `density_total`：总粒子密度，形状为 `(Nt, Nx)`。
- `charge_density_components`：A/B 分量电荷密度。
- `charge_density_total`：总电荷密度。
- `particles_per_cell`：每个时刻的每晶胞粒子数。
- `charge_per_cell`：每个时刻的每晶胞电荷。

`length_current_distribution.npz` 主要包含：

- `x_grid`：晶胞内实空间网格。
- `time_grid`：时间网格。
- `vector_potential`：用于抗磁项的矢势采样；长度规范默认全零。
- `density_components`：A/B 分量粒子密度，形状为 `(Nt, 2, Nx)`。
- `paramagnetic_current_components`：A/B 分量顺磁电流，形状为 `(Nt, 2, Nx)`。
- `diamagnetic_current_components`：A/B 分量抗磁电流，形状为 `(Nt, 2, Nx)`。
- `current_components`：A/B 分量总局域电流，形状为 `(Nt, 2, Nx)`。
- `current_total`：总局域电流，形状为 `(Nt, Nx)`。
- `current_per_cell`：每个时刻对一个晶胞积分得到的电流。
- `average_current`：每个时刻的晶胞平均电流。

`length_current_convergence.npz` 主要包含：

- `coarse_current_total`、`fine_current_total`：粗、细 k 网格的总局域电流。
- `coarse_current_components`、`fine_current_components`：粗、细 k 网格的 A/B 分量局域电流。
- `coarse_current_per_cell`、`fine_current_per_cell`：粗、细 k 网格的晶胞积分电流。
- `coarse_average_current`、`fine_average_current`：粗、细 k 网格的晶胞平均电流。
- `current_total_relative_max_difference` 等标量收敛指标。

读取示例：

```python
import numpy as np

rdm = np.load("length_rdm_result.npz")
rho = rdm["rho_trajectory"]
time_grid = rdm["time_grid"]

density = np.load("length_density_result.npz")
density_total = density["density_total"]

current = np.load("length_current_distribution.npz")
current_total = current["current_total"]
average_current = current["average_current"]

current_conv = np.load("length_current_convergence.npz")
current_error = current_conv["current_total_relative_max_difference"]
```

## 长度规范 RDM 方程

长度规范模块使用：

```text
(i hbar d_t - epsilon_ss') rho^E_kss' = i e E(t) [D, rho^E]_kss'
```

其中一维情况下：

```text
[D, rho] = d_k rho - i [xi, rho]
```

`xi` 是 `Geometry.py` 中通过本征矢有限差分计算的 Berry 联络。

## 局域电流分布

`Current_Distribution.py` 使用 `Band_Solver.py` 和 `Density.py` 的平面波基约定。能带本征矢
`eigenvectors[k]` 的列向量是 `U_{alpha G,s}(k)`，基底顺序为先 A 分量的 `G=-L..L`，
再 B 分量的 `G=-L..L`。晶胞内 Bloch 周期部分按

```text
u_{ks,alpha}(x) = sum_G U_{alpha G,s}(k) exp(i G x)
```

重构，k 空间积分权重沿用 `k_weight = dk/(2*pi)`。

局域电流算符来自最小耦合 `j(x) = -delta H / delta A(x)`。在当前模型和归一化下，模块直接在
plane-wave basis 中计算解析矩阵元，不使用微扰论、不使用能量分母，也不使用本征矢的 k 导数。
长度规范默认取 `A(t)=0`，因此默认输出中抗磁项为零；如需速度规范或其他规范下的抗磁项，
可以通过 `vector_potential` 参数显式传入。

当前电流分布实现没有使用数值微分。如果后续改为有限差分计算空间导数、Berry 联络或其他微分量，
对应微分步长和网格必须加入收敛性验证。现有主流程已经对粗、细 k 网格上的电流分布进行收敛比较。

作为一致性检查，一个晶胞内的局域电流积分满足

```text
average_current(t) = integral_cell dx current_total(x,t) / a0
```

并与解析 `dH/dk` 得到的宏观电流在机器精度内一致。

## 诊断脚本

从已有 RDM 结果重建密度：

```bash
python Density.py
```

该脚本会从脚本所在目录读取 `length_rdm_result.npz`，并生成 `length_density_result.npz`。

从已有 RDM 结果重建局域电流：

```bash
python Current_Distribution.py
```

该脚本会从脚本所在目录读取 `length_rdm_result.npz`，并生成 `length_current_distribution.npz`。
主流程 `python main.py` 还会额外生成 `length_current_convergence.npz`。

从已有 `.npz` 结果生成可视化图像：

```python
from Debug_Tools import plot_density_standard_views_from_npz
from Debug_Tools import plot_current_standard_views_from_npz

plot_density_standard_views_from_npz(
    "length_density_result.npz",
    "density_response",
    subtract_initial=True,
)
plot_current_standard_views_from_npz(
    "length_current_distribution.npz",
    "current_response",
    subtract_initial=True,
)
```

`subtract_initial=True` 表示画响应量 `data(t)-data(0)`。如果需要全量图，可以设为 `False`，
但当前主流程只默认保存全量时空热图，不再保存全量切片和全量时间曲线。

速度矩阵诊断文档 `VELOCITY_MATRIX_DIAGNOSTICS.md` 仍可作为理论参考。当前 `dev` 代码已经清理了
`Geometry.py` 中的速度矩阵实现，因此历史脚本 `Velocity_Matrix_Diagnostics.py` 不属于可运行主线；
如需恢复，应先重新实现并验证解析 `dH/dk` 速度矩阵接口。

## 工作记录

- `DAILY_REPORT_2026-06-30.md`：6月30日工作日报。
- `DAILY_REPORT_2026-07-01.md`：7月1日工作日报，记录局域电流、收敛诊断和可视化集成。
- `ONE_LINE_SUMMARY_2026-07-01.md`：7月1日一句话总结。


## 参考输出

`outputs_for_references/` 保存了一组 H2 最小实验的参考输出：

- `h2_bands.png`
- `h2_bands_near_zero.png`
- `h2_component_weights.png`
- `h2_coupling_comparison.png`
- `h2_intraband_velocities.png`
- `h2_interband_velocity.png`
- `h2_rho1_intraband.png`
- `h2_rho1_interband.png`
- `h2_density_matrix_time_dynamics.png`
- 若干 Markdown 摘要

这些文件目前更适合作为参考结果和回归对照。后续建议将其整理为 `examples/h2_reference/`
脚本和 `reference_outputs/h2/` 数据，而不是作为独立项目维护。

## 已知问题和后续工作

1. 速度矩阵诊断中的 Berry 联络有限差分形式仍需要处理近简并子空间、band tracking 和 BZ 边界 sewing matrix。
2. 局域电流分布已经按最小耦合公式实现并接入主流程收敛检查，但长度规范响应中是否需要额外极化项或规范一致的响应定义还需要继续确认。
3. `README.md` 和历史诊断脚本中仍有速度规范旧路径的背景信息，后续如恢复速度规范，应与长度规范主线清晰分离。
4. 需要建立系统收敛测试矩阵，分别检查 `Nk`、平面波截断 `L`、时间积分容差、实空间网格 `Nx` 和局域电流分布。
