# 速度矩阵诊断说明

主程序中速度规范 RDM 的生产路径使用直接的能带表象变换：

$$
v_{\mathrm{band}}(k)
= U^\dagger(k)\,\frac{1}{\hbar}\frac{\partial H(k)}{\partial k}\,U(k)
$$

这条路径只需要同一个 \(k\) 点处的本征矢和解析的哈密顿量导数，目前相对稳定。

论文中的等价公式可以写成：

$$
v_{kss'} = \frac{1}{\hbar}\,[D_k,H(k)]_{ss'} .
$$

在能量本征表象中，它等价于：

$$
v_{kss'}
= \frac{1}{\hbar}
\left[
\delta_{ss'}\,\partial_k \epsilon_{ks}
- i\left(\epsilon_{ks'}-\epsilon_{ks}\right)\xi_{kss'}
\right],
$$

其中

$$
\xi_{kss'} = i\langle u_{ks}|\partial_k u_{ks'}\rangle
$$

是 Berry 联络。

## 理论上的等价性

从本征方程出发：

$$
H(k)|u_s(k)\rangle = \epsilon_s(k)|u_s(k)\rangle .
$$

对 \(k\) 求导：

$$
\frac{\partial H}{\partial k}|u_s\rangle
+ H|\partial_k u_s\rangle
=
(\partial_k \epsilon_s)|u_s\rangle
+ \epsilon_s|\partial_k u_s\rangle .
$$

用 \(\langle u_{s'}|\) 左乘。对角项给出：

$$
\langle u_s|\partial_k H|u_s\rangle
= \partial_k \epsilon_s .
$$

非对角项给出：

$$
\langle u_{s'}|\partial_k H|u_s\rangle
=
(\epsilon_s-\epsilon_{s'})
\langle u_{s'}|\partial_k u_s\rangle .
$$

结合 Berry 联络定义，就可以得到 Eq. (28) 中的速度矩阵表达式。因此，如果导数是精确的，而且能带规范处理完全一致，那么 Berry 联络公式和直接变换公式应该给出相同的速度矩阵。

## 数值诊断为何仍会失败

直接变换公式使用解析的 \(\partial_k H\)，只依赖同一 \(k\) 点的本征矢；而 Eq. (28) 的数值诊断需要对不同 \(k\) 点的本征矢做有限差分，因此对以下问题非常敏感：

- 数值对角化给出的任意本征矢相位；
- 能带交叉或近简并附近的能带标号变化；
- 近简并子空间内部的任意酉旋转；
- 第一布里渊区边界处的周期 sewing matrix；
- 有限 \(k\) 网格下 Berry 联络和协变导数的离散化误差。

目前已经把相位平滑机制前移到 `Band_Solver.py`，并使 `Velocity_Matrix_Diagnostics.py` 从 `Geometry.py` 调用统一的 Berry 联络实现。但最新诊断仍显示 Eq. (28) 有限差分结果和直接变换结果差异较大。这说明问题不只是逐带 U(1) 相位，还可能涉及近简并子空间处理、band tracking 和 BZ 边界 sewing。

因此，当前生产路径仍应使用直接变换：

```text
U^\dagger (1/hbar dH/dk) U
```

Eq. (28) 路径目前只作为 Berry 联络质量和能带规范质量的诊断工具，不应直接用于生产计算。

独立运行诊断：

```bash
python Velocity_Matrix_Diagnostics.py
```
