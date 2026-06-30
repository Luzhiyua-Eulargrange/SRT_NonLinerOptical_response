# SRT RDM Solver

本项目用于计算一维双组分连续模型在外加脉冲电场下的能带结构、约化密度矩阵
（Reduced Density Matrix, RDM）动力学，以及由 RDM 重构得到的晶胞内实空间密度。

当前主线开发已经从速度规范电流收敛排查，转向以长度规范 RDM 和实空间密度收敛
为核心的验证流程。速度规范电流和速度矩阵诊断仍保留在项目中，但速度矩阵的
Berry 联络有限差分形式尚未完全稳定，相关结果应作为诊断而不是最终生产路径。

## 项目结构

- `main.py`：默认入口，转发到 `main_length_gauge.py`。
- `main_length_gauge.py`：长度规范 RDM、实空间密度重构和收敛诊断主流程。
- `main_velocity_gauge.py`：速度规范 RDM 与电流收敛主流程；当前该路径对截断带空间较敏感。
- `config.py`：默认参数、参数归一化与校验、第一布里渊区 k 网格生成。
- `Band_Solver.py`：构造一维双组分连续模型哈密顿量，计算能带和本征矢，并进行相位平滑。
- `Geometry.py`：计算速度算符矩阵、能带表象速度矩阵、Berry 联络和协变导数。
- `RDM_Common.py`：RDM 求解器共用的外场、初态密度矩阵、时间网格和 ODE 积分工具。
- `RDM_Length_Gauge.py`：长度规范 RDM 演化模块，对应 Phys. Rev. B 96, 035431 的式 (55)。
- `RDM_Velocity_Gauge.py`：速度规范 RDM 演化模块，对应 Phys. Rev. B 96, 035431 的式 (56)。
- `Density.py`：从 band-basis RDM 重构晶胞内实空间密度 `n_A(x,t)`、`n_B(x,t)` 和总密度。
- `Current.py`：从 RDM 轨迹计算宏观电流；主要服务速度规范路径。
- `Debug_Tools.py`：参数打印、能带绘图、RDM 轨迹检查和电流图绘制工具。
- `Velocity_Matrix_Diagnostics.py`：比较速度矩阵直接变换公式与 Berry 联络有限差分公式。
- `outputs_for_references/`：H2 最小实验的参考图片和 Markdown 摘要，用于对照和后续整理。

## 环境依赖

建议使用 Python 3.10 或更高版本。

```bash
pip install -r requirements.txt
```

依赖列表维护在 `requirements.txt` 中：

- `numpy`：矩阵构造、线性代数、数组保存与读取。
- `scipy`：调用 `scipy.integrate.solve_ivp` 进行 RDM 时间演化；如果未安装 `scipy`，会回退到内置四阶 Runge-Kutta 方法。
- `matplotlib`：用于能带图和电流图绘制。

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

默认输出文件包括：

- `length_band_result.npz`：长度规范流程使用的能带数据。
- `length_band_structure.png`：能带图。
- `length_rdm_result.npz`：细 k 网格上的长度规范 RDM 结果。
- `length_rho_convergence.npz`：粗、细 k 网格 RDM 收敛诊断。
- `length_density_result.npz`：晶胞内实空间密度结果。
- `length_density_convergence.npz`：实空间密度收敛诊断。

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

这些结果说明当前参数下实空间密度收敛很好。不过，由于 Berry 联络和能带规范处理仍在诊断中，
仍需继续检查近简并带、band tracking 和 BZ 边界 sewing matrix 的影响。

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

读取示例：

```python
import numpy as np

rdm = np.load("length_rdm_result.npz")
rho = rdm["rho_trajectory"]
time_grid = rdm["time_grid"]

density = np.load("length_density_result.npz")
density_total = density["density_total"]
```

## 双规范 RDM 方程

长度规范模块使用：

```text
(i hbar d_t - epsilon_ss') rho^E_kss' = i e E(t) [D, rho^E]_kss'
```

其中一维情况下：

```text
[D, rho] = d_k rho - i [xi, rho]
```

`xi` 是 `Geometry.py` 中通过本征矢有限差分计算的 Berry 联络。

速度规范模块使用：

```text
(i hbar d_t - epsilon_ss') rho^A_kss' = e A(t) [v, rho^A]_kss'
```

当前速度规范生产路径中的速度矩阵使用直接解析变换：

```text
v_band(k) = U^\dagger(k) (1/hbar dH/dk) U(k)
```

这条路径比 Berry 联络有限差分形式更稳定。后者目前只建议作为诊断工具。

## 诊断脚本

速度矩阵诊断：

```bash
python Velocity_Matrix_Diagnostics.py
```

该脚本比较：

- 直接变换公式 `U^\dagger (1/hbar dH/dk) U`
- Eq. (28) 中基于 Berry 联络和能量导数的有限差分形式

已知有限差分形式对本征矢相位、近简并带的酉旋转、能带标号、BZ 边界 sewing matrix 和 k 网格误差非常敏感。

从已有 RDM 结果重建密度：

```bash
python Density.py
```

该脚本会从脚本所在目录读取 `length_rdm_result.npz`，并生成 `length_density_result.npz`。

速度规范电流流程：

```bash
python main_velocity_gauge.py
```

该流程会检查速度规范电流在粗、细 k 网格之间的收敛。如果超过容差，会显式抛出错误，而不是保存不可靠结果。

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

1. 速度矩阵诊断仍显示 Berry 联络有限差分结果与直接解析变换结果存在较大差异，需要处理近简并子空间、band tracking 和 BZ 边界 sewing matrix。
2. 长度规范下的电流定义尚未完全整理清楚，是否需要额外极化项或规范一致的响应定义还需要从理论上确认。
3. 当前入口脚本、诊断脚本、物理模块和保存逻辑仍有冗余，后续应进一步收敛项目结构。
4. 需要建立系统收敛测试矩阵，分别检查 `Nk`、平面波截断 `L`、时间积分容差和实空间网格 `Nx`。
