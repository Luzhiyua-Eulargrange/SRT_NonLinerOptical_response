# 工作日报

日期：2026年6月30日

## 一句话总结

围绕一维双分量连续模型，今天完成了电流空间分布理论推导的重写和补充，并对 `dev` 分支相对 `main_0.4` 的规范模块做了收敛到长度规范与 Berry 联络主线的清理。

## 今日工作

1. 重写 `Theory_Formula_of_Current_Space_Distribution.md`，从最小耦合出发重新组织电流算符推导，并将双层 AB Hamiltonian 二次量子化到平面波/Bloch 块形式。
2. 补充动量空间产生湮灭算符 `c_{n k}` 与实空间场算符 `\Psi_A(x),\Psi_B(x)` 之间的 Fourier 关系，为后续从 Hamiltonian 变分得到局域电流打基础。
3. 新增两份更完整的电流空间分布推导稿：`Theory_Formula_of_Current_Space_Distribution_codex.md` 和 `Theory_Formula_of_Current_Space_Distribution_full_derivation.md`，分别用于整理从 `H_2` 模型到局域电流公式、以及从 `c_k` 量子化到密度矩阵表达的完整链条。
4. 对比 `main_0.4` 后，清理 `Geometry.py` 中速度矩阵相关实现，使该模块聚焦 Berry connection 与协变导数计算。
5. 对比 `main_0.4` 后，清理 `RDM_Common.py` 中速度规范、矢势漂移动量和 legacy 单 k 平面波 RDM 传播接口，使共享模块只保留时间网格、外场、初始 band-basis RDM、Hermitize 与复 ODE 求解等通用能力。
6. 清理 `Debug_Tools.py` 中电流时间演化画图函数，避免调试工具继续依赖已经拆出的电流计算路径。
7. 保留并整理一维模型说明、能带数据与图像结果，包括 `1d_model.md`、`1d_model_bands.csv` 和 `1d_model_bands.png`，作为理论推导和数值结果的背景材料。

## 与 main_0.4 的主要差异

1. `dev` 新增：`1d_model.md`、`1d_model_bands.csv`、`1d_model_bands.png`、`Theory_Formula_of_Current_Space_Distribution_codex.md`、`Theory_Formula_of_Current_Space_Distribution_full_derivation.md`。
2. `dev` 修改：`Debug_Tools.py`、`Geometry.py`、`RDM_Common.py`、`Theory_Formula_of_Current_Space_Distribution.md`。
3. `dev` 中未保留 `main_0.4` 的 `Current.py`、`main_velocity_gauge.py`、`RDM_Velocity_Gauge.py`、`DAILY_REPORT_2026-06-29.md` 和 `WORK_SUMMARY.md`。
4. 两个目录中 `Band_Solver.py`、`config.py`、`Density.py`、`main.py`、`main_length_gauge.py`、`RDM_Length_Gauge.py`、`README.md`、`VELOCITY_MATRIX_DIAGNOSTICS.md`、`Velocity_Matrix_Diagnostics.py` 以及主要长度规范输出文件哈希一致。

## 当前结论

1. 6月30日的核心进展是把电流空间分布问题从直接套用全空间平均电流公式，转向从二次量子化 Hamiltonian、Fourier 展开和最小耦合出发重新推导局域表达式。
2. 当前 `dev` 代码结构更偏向长度规范 RDM、实空间密度和 Berry 联络诊断；速度规范电流计算链条没有作为核心代码保留在 `dev` 中。
3. 理论推导中需要继续检查符号一致性，尤其是 `B` 层势能项在矩阵 Hamiltonian 中为 `-V_B(x)`，而某些动量空间写法若不带负号会导出 `+V_B(x)`。

## 后续计划

1. 将 `c_{n k}` 形式 Hamiltonian 中每一项代回 Fourier 展开的计算补全到主理论文档，并显式化简掉 `k,n,x,x'` 积分中的完备关系。
2. 基于最终实空间 Hamiltonian 对外加矢势做变分，得到 A/B 两分量的局域顺磁电流与可能的抗磁项。
3. 将局域电流公式与已有 RDM 结果连接，明确应在 plane-wave basis、band basis 还是实空间密度矩阵中计算电流分布。
4. 若后续恢复速度规范计算，应把 `Current.py`、`RDM_Velocity_Gauge.py` 与长度规范路径分离成清晰的可选模块，避免再次混入共享工具层。
