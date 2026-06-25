# SRT RDM Solver

本项目用于计算一维双组分连续模型在外加脉冲电场下的能带结构与宏观电流响应。程序首先在第一布里渊区的动量网格上构造并对角化截断平面波基下的哈密顿量，随后使用约化密度矩阵（Reduced Density Matrix, RDM）方法推进体系随时间的演化，并对不同初始动量的电流贡献进行积分。

## 文件结构

- `main.py`：示例运行入口，设置参数、生成 k 网格、计算能带和总电流，并保存结果。
- `config.py`：默认参数、参数归一化与校验、第一布里渊区 k 网格生成。
- `Band_Solver.py`：构造一维双组分连续模型哈密顿量，计算单个或多个 k 点的本征能量和本征矢。
- `RDM_Solver.py`：定义电场、矢势、动量漂移以及 RDM 时间演化。
- `Geometry.py`：计算速度算符矩阵，支持平面波基和能带本征基。
- `Current.py`：计算单个 k 轨道的电流贡献，并在第一布里渊区积分得到总电流。

## 环境依赖

建议使用 Python 3.10 或更高版本。项目主要依赖：

```bash
pip install numpy scipy
```

其中 `scipy` 用于调用 `scipy.integrate.solve_ivp` 进行时间演化；如果未安装 `scipy`，程序会回退到内置的四阶 Runge-Kutta 方法。

## 运行方法

在项目目录下运行：

```bash
python main.py
```

运行完成后，程序会在当前目录生成 `result.npz`，其中包含：

- `k_grid`：第一布里渊区内的 k 点网格。
- `energies`：每个 k 点对应的能带能量，形状为 `(Nk, Nb)`。
- `time_grid`：时间演化网格。
- `current`：随时间变化的总电流。

## 主要参数

参数通过字典传入 `normalize_params`，未提供的参数会使用 `config.py` 中的默认值。常用参数包括：

- `L`：倒格矢截断参数，决定平面波基大小；总基维度 `Nb = 2 * (2L + 1)`。
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
    "fermi_energy": 0.0,
})
```

## 结果读取示例

```python
import numpy as np

data = np.load("result.npz")
k_grid = data["k_grid"]
energies = data["energies"]
time_grid = data["time_grid"]
current = data["current"]
```

这些数据可进一步用于绘制能带图、时间依赖电流曲线或频谱分析。
