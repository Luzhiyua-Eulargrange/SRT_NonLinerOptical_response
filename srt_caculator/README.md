# SRT RDM Solver

本项目用于计算一维双组分连续模型在外加脉冲电场下的能带结构与宏观电流响应。程序首先在第一布里渊区的动量网格上构造并对角化截断平面波基下的哈密顿量，随后使用约化密度矩阵（Reduced Density Matrix, RDM）方法推进体系随时间的演化，并对不同初始动量的电流贡献进行积分。

## 文件结构

- `main.py`：示例运行入口，设置参数、计算能带、保存结果，并对两种 RDM 规范做轻量调试检查。
- `config.py`：默认参数、参数归一化与校验、第一布里渊区 k 网格生成。
- `Band_Solver.py`：构造一维双组分连续模型哈密顿量，计算单个或多个 k 点的本征能量和本征矢。
- `Geometry.py`：计算速度算符矩阵、能带表象速度矩阵、Berry 联络以及协变导数。
- `RDM_Common.py`：RDM 求解器共用的外场、初态密度矩阵和 ODE 积分工具。
- `RDM_Velocity_Gauge.py`：速度规范 RDM 演化模块，对应 Phys. Rev. B 96, 035431 的式 (56)。
- `RDM_Length_Gauge.py`：长度规范 RDM 演化模块，对应 Phys. Rev. B 96, 035431 的式 (55)。
- `Debug_Tools.py`：默认参数打印、能带图绘制、本征矢打印以及 RDM 轨迹检查工具。
- `Current.py`：计算单个 k 轨道的电流贡献，并在第一布里渊区积分得到总电流。

## 环境依赖

建议使用 Python 3.10 或更高版本。项目主要依赖：

```bash
pip install -r requirements.txt
```

依赖列表维护在 `requirements.txt` 中：

- `numpy`：矩阵构造、线性代数、数组保存与读取。
- `scipy`：调用 `scipy.integrate.solve_ivp` 进行 RDM 时间演化；如果未安装 `scipy`，程序会回退到内置的四阶 Runge-Kutta 方法。
- `matplotlib`：用于 `Debug_Tools.py` 中的能带图绘制。

## 运行方法

在项目目录下运行：

```bash
python main.py
```

当前主程序会先进行能带计算和调试检查，并在当前目录生成：

- `band_result.npz`：能带计算数据文件。
- `band_structure.png`：能带图。

`band_result.npz` 包含：

- `k_grid`：第一布里渊区内的 k 点网格。
- `k_weight`：k 空间积分权重。
- `energies`：每个 k 点对应的能带能量，形状为 `(Nk, Nb)`。
- `eigenvectors`：每个 k 点对应的本征矢矩阵，形状为 `(Nk, Nb, Nb)`。

随后主程序会使用一个较小的调试网格分别调用速度规范和长度规范 RDM 模块，检查密度矩阵的 Hermitian 性和迹守恒。这个步骤用于快速验证模块接口；正式计算可以直接调用对应的 `propagate_velocity_gauge_rdm` 或 `propagate_length_gauge_rdm`。

## 主要参数

参数通过字典传入 `normalize_params`，未提供的参数会使用 `config.py` 中的默认值。常用参数包括：

- `L`：倒格矢截断参数，决定平面波基大小；总基维度 `Nb = 2 * (2L + 1)`。
- `num_k_band`：能带绘图和能带结果保存使用的 k 网格点数。
- `num_k_rdm`：RDM、Berry 联络和协变导数使用的 k 网格点数，正式计算时建议做收敛测试。
- `a0`：晶格常数，倒格矢大小为 `b = 2*pi/a0`。
- `m`、`hbar`：有效质量与约化普朗克常数。
- `kappa`、`v0`、`vA`、`vB`、`w1`、`w2`：模型哈密顿量参数。
- `E0`、`omega`、`pulse_duration`、`t_switch`：外加脉冲电场的强度、频率、持续时间与开关时间。
- `num_time_steps`：时间演化采样点数。
- `fermi_energy`：零温初态占据所用的费米能级。
- `rtol`、`atol`：`solve_ivp` 的相对和绝对误差容限。

可以直接修改 `main.py` 中的参数字典，例如：

```python
params = normalize_params({
    "L": 2,
    "E0": 0.02,
    "omega": 0.5,
    "pulse_duration": 40.0,
    "num_time_steps": 201,
    "num_k_band": 51,
    "num_k_rdm": 101,
    "fermi_energy": 0.0,
})
```

## 结果读取示例

```python
import numpy as np

data = np.load("band_result.npz")
k_grid = data["k_grid"]
k_weight = float(data["k_weight"])
energies = data["energies"]
eigenvectors = data["eigenvectors"]
```

这些数据可进一步用于绘制能带图、检查本征矢或作为后续电流计算的输入。

## RDM 双规范模块

两个新 RDM 模块都在能带表象中推进密度矩阵，输入为 `k_grid`、`energies` 和 `eigenvectors`。

速度规范模块使用：

```text
(i hbar d_t - epsilon_ss') rho^A_kss' = e A(t) [v, rho^A]_kss'
```

长度规范模块使用：

```text
(i hbar d_t - epsilon_ss') rho^E_kss' = i e E(t) [D, rho^E]_kss'
```

其中一维情况下：

```text
[D, rho] = d_k rho - i [xi, rho]
```

`xi` 是 `Geometry.py` 中通过本征矢有限差分计算的 Berry 联络。
