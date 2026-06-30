# 长度规范 RDM 不收敛问题阶段性总结

日期：2026-06-27

本文记录目前对长度规范（length gauge）RDM 电流计算不收敛问题的排查结论。重点文件包括：

- `main_length_gauge.py`
- `RDM_Length_Gauge.py`
- `RDM_Common.py`
- `Geometry.py`
- `Band_Solver.py`
- `Current.py`

当前结论可以先概括为一句话：长度规范结果不收敛的主要嫌疑不是时间积分器，而是 `k` 空间几何量的离散化，尤其是 Berry connection 和协变导数在第一布里渊区边界处错误地使用了周期有限差分，但本征矢规范并没有真正做到周期连续。

---

## 1. 当前长度规范计算流程

长度规范主程序在 `main_length_gauge.py` 中。主要流程是：

1. 设置参数，例如：
   - `L = 2`
   - `E0 = 0.02`
   - `omega = 0.5`
   - `pulse_duration = 20.0`
   - `num_time_steps = 81`
   - `num_k_rdm = 21`
   - `rtol = 1e-6`
   - `atol = 1e-8`

2. 用 `make_k_grid()` 生成第一布里渊区中的均匀 `k` 网格。

3. 用 `solve_bands()` 对每个 `k` 点独立对角化哈密顿量，得到：
   - `energies`
   - `eigenvectors`

4. 调用 `propagate_length_gauge_rdm()` 做长度规范 RDM 时间演化。

5. 调用 `total_current_from_band_rdm()` 从 RDM 轨迹计算总电流。

6. 通过粗网格和细网格比较电流曲线，判断是否收敛。

当前主程序中的收敛检查是：

```python
coarse_num_k = max(5, params["num_k_rdm"] // 2)
fine_num_k = params["num_k_rdm"]
convergence_error = relative_current_difference(fine["current"], coarse["current"])
```

也就是说，默认 `num_k_rdm = 21` 时，实际比较的是大约 `Nk = 10` 和 `Nk = 21` 的电流差异。

---

## 2. 已经排除的问题

### 2.1 时间积分器不是目前的首要问题

`RDM_Common.py` 中的 `solve_complex_ode()` 对复数 ODE 的求解方式是：

```python
solve_ivp(..., method="RK45", ...)
```

也就是 SciPy 的自适应 Runge-Kutta 4(5) 方法。

如果 SciPy 不可用，代码才会退回到 `_rk4_fallback()`，即手写的固定步长四阶 Runge-Kutta。

目前检查了长度规范 RDM 演化中的基本守恒量：

```text
Nk = 21:
  Hermitian error = 0.0
  trace drift     ≈ 2.13e-14

Nk = 41:
  Hermitian error = 0.0
  trace drift     ≈ 3.55e-14
```

这说明：

- 密度矩阵的厄米性没有明显破坏。
- `Tr(rho)` 基本守恒。
- 时间积分器没有表现出明显的不稳定或者漂移。

因此，目前没有证据表明“不收敛”主要来自 `RK45` 或时间步长容差。

当然，这并不代表时间积分永远没有问题。只是就当前诊断结果看，它不是最优先修复对象。

### 2.2 费米面占据跳变不是当前看到的主要问题

初始 RDM 是通过能量和费米能级判断占据：

```python
occupations = (energies < p["fermi_energy"]).astype(float)
```

如果某条能带刚好穿过费米能级，随着 `k` 变化占据数会发生跳变，`d_k rho` 会出现不光滑，这也可能导致长度规范不收敛。

但在当前参数下，检查到的占据模式是稳定的：

```text
Nk = 10, 21, 41, 81:
  occ_counts = 一致
  max_occ_jump = 0
```

也就是说，当前测试参数下没有看到明显的占据跳变。因此，费米面不连续不是这次最直接的原因。

---

## 3. 已经发现的核心问题

### 3.1 `Geometry.py` 中使用了周期中心差分

`Geometry.py` 中 Berry connection 的计算使用了：

```python
derivatives = (np.roll(vectors, -1, axis=0) - np.roll(vectors, 1, axis=0)) / (2.0 * dk)
```

协变导数也使用了类似形式：

```python
derivative = (np.roll(matrices, -1, axis=0) - np.roll(matrices, 1, axis=0)) / (2.0 * dk)
```

`np.roll` 的含义是把数组首尾相接。因此这段代码默认认为：

```text
k 网格最左端和最右端可以直接相邻相减。
```

也就是说，代码假设第一布里渊区边界处是周期连续的。

这个假设在物理上不是完全错误，但在数值实现中必须非常小心。因为能带本征矢 `eigenvectors` 本身带有任意相位，而且在多能带问题中还可能有能带重排、近简并子空间旋转等问题。只有本征矢规范已经处理成真正的周期连续规范时，才可以安全地用 `np.roll` 跨边界做差分。

当前代码没有做到这一点。

### 3.2 当前的相位平滑不保证 BZ 边界周期连续

`Geometry.py` 中有一个函数：

```python
def smooth_eigenvector_gauge(eigenvectors):
```

它的作用是沿着 `k` 网格逐点修正相邻本征矢的相位，使相邻点看起来更平滑。

但它只做了从第 0 个点到最后一个点的顺序平滑，没有强制最后一个点和第 0 个点在第一布里渊区边界处也连续。

更重要的是，它只逐条能带修相位，不会解决以下问题：

- 近简并能带之间的子空间旋转；
- 能带交叉或避免交叉附近的带标号交换；
- 第一布里渊区边界处的周期 sewing matrix；
- 截断平面波基下 `k = -b/2` 和 `k = +b/2` 的基底重标号问题。

因此，即使调用了 `smooth_eigenvector_gauge()`，也不能保证 `np.roll` 跨边界差分是合法的。

### 3.3 Berry connection 随 `Nk` 增大反而变大

做过一个轻量诊断，结果如下：

```text
Nk = 10: max|xi| ≈ 1.12
Nk = 21: max|xi| ≈ 2.00
Nk = 41: max|xi| ≈ 4.20
Nk = 81: max|xi| ≈ 7.48
```

这里的 `xi` 是 Berry connection。

如果离散化是正常收敛的，随着 `Nk` 增大、`dk` 变小，`max|xi|` 不应该这样持续变大。当前结果更像是：

```text
某个不连续跳变 / dk
```

由于 `dk` 越来越小，这个假尖峰就越来越大。

这非常符合“在 BZ 边界处把不连续的本征矢拿去做周期差分”的症状。

### 3.4 BZ 边界处本征矢几乎不连续

进一步检查了平滑后的相邻同带重叠，发现 BZ 边界处的最小重叠约为：

```text
boundary min overlap ≈ 1.7e-11
```

这个数非常接近 0，说明代码认为应该周期相邻的两个点，在当前本征矢表示下几乎正交。

这不是一个小误差，而是结构性问题。

如果这两个点被 `np.roll` 当成相邻点做差分，就会得到一个很大的、非物理的导数项。

### 3.5 电流确实没有随 `Nk` 收敛

当前长度规范电流的相对差异诊断结果是：

```text
Nk 10 vs 21: relative difference ≈ 1.10
Nk 21 vs 41: relative difference ≈ 0.95
Nk 41 vs 81: relative difference ≈ 0.98
```

这说明增大 `Nk` 并没有让电流曲线靠近某个稳定结果。

结合 Berry connection 的发散趋势，当前最合理的解释是：

```text
k 空间几何量离散化错误污染了长度规范演化方程，导致电流不收敛。
```

---

## 4. 为什么这个问题主要影响长度规范

长度规范 RDM 方程中有协变导数项。代码在 `RDM_Length_Gauge.py` 中写成：

```python
covariant_derivative = covariant_derivative_matrix(
    rho_grid,
    berry_connection_grid,
    k_grid,
)

field_term = p["e_charge"] * electric_field(t, p) * covariant_derivative / p["hbar"]
```

也就是说，长度规范直接依赖：

- `d_k rho`
- Berry connection `xi`
- `-i[xi, rho]`

这些都是 `k` 空间几何量，对本征矢规范、能带连续性、边界处理非常敏感。

相比之下，速度规范主要用速度矩阵和对易子，不需要直接对 `rho(k)` 做 `k` 空间导数。因此同样的本征矢规范问题在速度规范里通常没有这么直接、这么严重。

这就是为什么现在应该专注长度规范的 `Geometry.py` 和相关离散化，而不是先去改 ODE 求解器。

---

## 5. 可能还存在的问题

下面这些问题目前还没有完全确认，但很可能会影响后续收敛性。

### 5.1 能带排序只按本征值排序，可能存在带标号不连续

`Band_Solver.py` 中每个 `k` 点都是独立调用：

```python
np.linalg.eigh(H)
```

`np.linalg.eigh` 会按本征值排序。但在多能带系统中，仅按能量排序不一定能保证“第 n 条带”随 `k` 连续。

当存在近简并、避免交叉或真实交叉时，两个相邻 `k` 点上的能带索引可能发生交换。

如果带标号不连续，那么后面计算：

```python
<u_n(k) | d_k u_m(k)>
```

就会出现非物理尖峰。

解决这类问题通常需要基于相邻 `k` 点本征矢 overlap 的 band matching，而不是只依赖本征值排序。

### 5.2 近简并子空间需要整体处理

如果两条或多条能带非常接近，逐条能带修相位是不够的。

原因是：在近简并子空间内部，数值对角化返回的本征矢可以发生任意酉旋转。这个旋转并不一定有物理意义，但会严重影响 Berry connection 的有限差分。

更稳的做法是对子空间做 parallel transport，或者使用 gauge-covariant 的离散公式。

### 5.3 BZ 边界可能需要 sewing matrix，而不是简单相等

对周期系统来说，`k` 和 `k + G` 的物理态等价，但在平面波截断基底中，矩阵表示未必是简单相等。

特别是这里使用的是有限数量的 reciprocal lattice components。BZ 边界两侧的基底可能需要重标号才能对应起来。

因此，即使物理上是周期的，代码中的 `eigenvectors[0]` 和 `eigenvectors[-1]` 也不一定可以直接相减。

如果要做严格的周期差分，可能需要构造边界处的 sewing matrix。

### 5.4 当前收敛判据比较粗糙

现在的收敛判据是比较：

```text
coarse Nk = num_k_rdm // 2
fine Nk = num_k_rdm
```

并使用电流最大范数的相对误差。

这可以作为初步检查，但还不够系统。更好的做法是比较一组网格：

```text
Nk = 21, 41, 81, 161, ...
```

同时观察：

- `max|current|`
- `max|J_Nk - J_2Nk|`
- Berry connection 的最大值和分布
- RDM trace drift
- Hermitian error
- 电流曲线形状是否稳定

### 5.5 平面波截断 `L` 也可能影响结果

当前长度规范主程序里使用了：

```python
L = 2
```

这意味着平面波基底比较小。即使修复了 `k` 空间规范问题，也应该检查 `L` 收敛：

```text
L = 2, 3, 4, ...
```

否则可能出现 `k` 网格收敛了，但基底截断还没有收敛的情况。

---

## 6. 建议的解决方案

这里按优先级从低风险到高质量修复排列。

### 6.1 第一阶段：增加诊断，不急着改物理公式

建议先在代码中加入诊断函数，输出以下量：

1. 相邻 `k` 点本征矢 overlap：

   ```text
   min |<u_n(k_i) | u_n(k_{i+1})>|
   ```

2. BZ 边界 overlap：

   ```text
   min |<u_n(k_last) | u_n(k_0)>|
   ```

3. Berry connection 的尺度：

   ```text
   max|xi|
   mean|xi|
   max|xi| 出现在哪个 k 点、哪两个带之间
   ```

4. 协变导数的尺度：

   ```text
   max|D_k rho|
   ```

5. RDM 守恒量：

   ```text
   Hermitian error
   trace drift
   ```

这样可以避免后续修改时只看最终电流，导致不知道问题到底变好了还是被掩盖了。

### 6.2 第二阶段：暂时避免跨 BZ 边界的错误周期差分

作为短期诊断，可以先把 `np.roll` 周期差分替换成非周期差分：

- 内点使用中心差分；
- 左边界使用前向差分；
- 右边界使用后向差分。

这不是最终最优解，因为物理上 BZ 是周期的；但它可以帮助判断：

```text
不收敛是否主要来自首尾边界差分。
```

如果去掉跨边界差分后 `max|xi|` 不再随 `1/dk` 发散，电流也明显改善，就可以进一步确认问题来源。

注意：这个方案只是诊断和临时稳定用，不是严格的周期规范实现。

### 6.3 第三阶段：实现基于 overlap 的能带匹配

现在每个 `k` 点的能带只按能量排序。建议增加一个 band tracking 步骤：

1. 从第一个 `k` 点开始。
2. 对相邻两个 `k` 点计算 overlap 矩阵：

   ```python
   O_mn = <u_m(k_i) | u_n(k_{i+1})>
   ```

3. 根据 overlap 最大原则重新排列 `k_{i+1}` 点的能带顺序。
4. 再做相位平滑。

如果有近简并，单纯最大 overlap 还不一定够，但它比只按本征值排序要稳很多。

### 6.4 第四阶段：处理周期规范或 sewing matrix

如果要继续使用周期差分，就必须让 BZ 首尾真正连续。

可以考虑两种路线：

#### 路线 A：构造周期平滑规范

目标是让：

```text
u_n(k_last) 和 u_n(k_0)
```

在相位和带标号上尽可能连续。

这通常需要：

- 沿整个 BZ 做 parallel transport；
- 计算首尾 mismatch；
- 把 mismatch 分散回整个 `k` 网格，使规范变成周期的。

#### 路线 B：在边界使用 sewing matrix

不强行要求本征矢数组首尾相同，而是在跨边界差分时显式加入一个变换矩阵，把 `k_last` 的表示转换到 `k_0` 的规范下。

这种方式更物理，但实现难度更高。

### 6.5 第五阶段：使用 gauge-covariant 离散导数

更推荐的长期方案是避免直接差分本征矢，改用 link variable / Wilson line 形式。

核心思想是：

- 不直接比较不同 `k` 点上任意规范下的矩阵；
- 先用 overlap 矩阵建立相邻 `k` 点之间的平行移动；
- 再定义协变差分。

这样可以显著降低结果对本征矢相位选择的敏感性。

这类方法比直接有限差分 Berry connection 更适合多能带数值计算。

### 6.6 第六阶段：重新做系统收敛测试

修复几何离散化之后，需要重新做收敛测试。建议至少检查：

```text
Nk = 21, 41, 81, 161
L  = 2, 3, 4
rtol/atol = 1e-6/1e-8 和 1e-8/1e-10
```

判断顺序建议是：

1. 先固定 `L` 和 ODE 容差，看 `Nk` 是否收敛。
2. 再固定足够大的 `Nk`，看 `L` 是否收敛。
3. 最后检查 ODE 容差是否影响结果。

不要一开始同时改变所有参数，否则很难判断是哪一类误差在主导。

---

## 7. 建议的近期修改顺序

建议下一步按下面顺序做：

1. 在 `Geometry.py` 或单独的 debug 脚本里加入 overlap、Berry connection、协变导数诊断。
2. 临时实现非周期边界差分版本，用来验证 BZ 边界是否是主要污染源。
3. 增加基于 overlap 的 band tracking，减少能带标号不连续。
4. 重新检查 `max|xi|` 是否还随 `Nk` 增大而发散。
5. 如果诊断改善，再重新跑长度规范电流收敛。
6. 最后再实现更严格的周期规范或 gauge-covariant link 方案。

---

## 8. 当前判断

目前最可信的判断是：

```text
长度规范不收敛主要来自 k 空间有限差分和本征矢规范处理不一致。
```

具体表现为：

- `Geometry.py` 使用 `np.roll` 做周期差分；
- 但 `eigenvectors` 没有被处理成 BZ 周期连续规范；
- BZ 边界处本征矢重叠几乎为 0；
- Berry connection 最大值随 `Nk` 增大而增大；
- 电流随 `Nk` 没有收敛趋势；
- RDM 的 Hermiticity 和 trace 基本正常，因此 ODE 时间积分不是当前首要矛盾。

因此，不建议现在优先调 `rtol`、`atol` 或单纯增大 `num_k_rdm`。这些操作很可能只会让计算更慢，甚至让边界假尖峰更明显。

更应该先修复或替换长度规范中的 `k` 空间几何离散化。
