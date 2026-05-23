# 一维连续模型与数值求解入门

本文档先介绍一个一维连续模型的形式，再讨论如何用平面波基组对该模型进行数值求解。

## 目录

- [一维连续模型与数值求解入门](#一维连续模型与数值求解入门)
  - [目录](#目录)
  - [1. 研究目标](#1-研究目标)
  - [2. 1D 连续模型的形式](#2-1d-连续模型的形式)
  - [3. 连续模型的求解方法](#3-连续模型的求解方法)
    - [3.1 平面波基组与周期性](#31-平面波基组与周期性)
    - [3.2 不变子空间与布里渊区](#32-不变子空间与布里渊区)
    - [3.3 两分量波函数与有限矩阵](#33-两分量波函数与有限矩阵)
    - [3.4 矩阵元](#34-矩阵元)
- [H\_{\\alpha\\beta}(\\ell,\\ell')](#h_alphabetaellell)
- [H\_{AA}(\\ell,\\ell')](#h_aaellell)
- [H\_{BB}(\\ell,\\ell')](#h_bbellell)
- [H\_{AB}(\\ell,\\ell')](#h_abellell)
- [H\_{BA}(\\ell,\\ell')](#h_baellell)
- [\\langle k+G\_\\ell|e^{ibx}|k+G\_{\\ell'}\\rangle](#langle-kg_elleibxkg_ellrangle)
    - [3.5 平面波截断](#35-平面波截断)
    - [3.6 数值流程](#36-数值流程)
    - [3.7 参考结果](#37-参考结果)

## 1. 研究目标

本文档的目标是理解一个一维连续模型的基本形式，并掌握数值求解这类模型的基本流程。我们将从单粒子哈密顿量出发，说明晶格周期性如何简化问题，并以 $H_2$ 为例（定义见下文）写出平面波基组下的有限维矩阵。

## 2. 1D 连续模型的形式

考虑电子在一维晶格中运动，晶格常数为 $a_0$。这里忽略电子之间的相互作用。

![电子一维运动示意图](figures/1d_motion.png)

由于已经忽略了电子之间的相互作用，我们可以先求解单电子本征态，再按照泡利不相容原理填充这些态，从而得到多电子体系的基态或激发态。

最简单的单电子哈密顿量形式如下：

$$
\tag{1}
H_1 = \frac{p^2}{2 m} + V(x),
$$

其中 $V(x+a_0) = V(x)$ 是电子受到的周期性晶格势。例如，$V(x)=V_0 \cos(2\pi x/a_0)+ V_1 \sin(4\pi x/a_0)$ 就是一个满足周期性的势能。

本文档重点讨论的是下面这个稍复杂的单电子哈密顿量：

$$
\tag{2}
H_2 = 
  \begin{pmatrix}
  \frac{(p-\hbar\kappa)^2}{2m}+V_A(x)+v_0 & W(x) \\
  W^*(x) & -\frac{p^2}{2m}-V_B(x)-v_0
  \end{pmatrix},
$$
其中 $V_A(x)=v_A e^{ibx}+v^*_A e^{-ibx}$，$V_B(x)=v_B e^{ibx}+v^*_B e^{-ibx}$，$W(x)=w_1 e^{ibx}+w_2 e^{-ibx}$，$b=2\pi/a_0$。

讨论 $H_2$ 而不是 $H_1$ 的原因如下。

1. 求解 $H_1$ 和求解 $H_2$ 的基本思路是相通的。学会了其中一个模型的处理方法，就能理解另一个模型的大部分求解步骤。这里选择讨论更复杂、物理内容更丰富的 $H_2$。如果遇到概念上想不清楚的地方，也可以先退回 $H_1$，思考在更简单的模型中应该如何处理。

2. $H_2$ 是一个 $2\times 2$ 矩阵哈密顿量，对应的波函数有两个分量：

$$
\psi(x)=
\begin{pmatrix}
\psi_A(x) \\
\psi_B(x)
\end{pmatrix}.
$$

可以把这两个分量理解为上下两层、两个子晶格，或者两种自旋态。对角元描述两个分量各自的运动，非对角元 $W(x)$ 描述它们之间的耦合。因此，$H_2$ 比 $H_1$ 更接近莫尔转角材料等实际体系中的有效模型。

3. 更重要的是，我们的最终目标是在加入光场后，模拟电荷和电流随时间的演化。这里重点关注的注入电流（injection current）在具有时间反演对称性的体系中通常被禁止，例如最简单的 $H_1$ 模型。在 $H_2$ 中，非零的 $\kappa$ 可以破坏这种时间反演对称性，从而允许注入电流出现。

关于 $H_2$ 形式的几点说明：

1. $H_2$ 不是单分量的标量算符，而是一个作用在两分量波函数上的算符矩阵。也就是说，在每一个位置 $x$ 上，波函数都有两个分量 $\psi_A(x)$ 和 $\psi_B(x)$。哈密顿量作用在波函数上时，既会分别改变这两个分量，也会通过非对角项把一个分量耦合到另一个分量。

2. 矩阵的两个对角元描述两个分量各自的运动。左上角

$$
\frac{(p-\hbar\kappa)^2}{2m}+V_A(x)+v_0
$$

可以理解为 $A$ 分量的动能、周期势能和整体能量偏移；右下角

$$
-\frac{p^2}{2m}-V_B(x)-v_0
$$

可以理解为 $B$ 分量的对应项。两者前面的符号不同，表示这两个分量对应的能带弯曲方向不同；粗略地说，一个对应“向上弯”的能带，另一个对应“向下弯”的能带。

3. 动量算符为 $p=-i\hbar\partial_x$。左上角的 $p-\hbar\kappa$ 表示 $A$ 分量的能带在动量空间中发生了平移，平移量由 $\kappa$ 控制。

4. $V_A(x)$ 和 $V_B(x)$ 是周期势，周期都是 $a_0$。写成 $e^{ibx}$ 和 $e^{-ibx}$ 的形式，是因为周期函数可以用傅里叶分量展开。这里为了得到一个尽量简单的模型，只保留了最基本的两个傅里叶分量，其中 $b=2\pi/a_0$ 是倒格矢。

5. 这个模型中一共有 8 个可调参数：

$$
a_0,\quad m,\quad \kappa,\quad v_0,\quad v_A,\quad v_B,\quad w_1,\quad w_2.
$$

其中 $a_0$、$m$、$\kappa$、$v_0$ 是实数参数。$a_0$ 是晶格常数，通常取正数；$m$ 是有效质量，也通常取正数；$\kappa$ 控制两个分量在动量空间中的相对平移；$v_0$ 控制两个分量之间的整体能量偏移。$v_A$、$v_B$、$w_1$、$w_2$ 是复数参数。由于 $V_A(x)$ 和 $V_B(x)$ 中同时出现了 $v$ 和 $v^*$，所以 $V_A(x)$、$V_B(x)$ 本身是实函数；而 $W(x)$ 一般是复函数，用来描述两个分量之间的复耦合。

6. 非对角元 $W(x)$ 和 $W^*(x)$ 描述 $A$、$B$ 两个分量之间的耦合。如果 $W(x)=0$，两个分量彼此独立，模型可以分成两个互不影响的问题；如果 $W(x)\neq 0$，两个分量会发生混合。

7. 左下角写成 $W^*(x)$ 是为了保证 $H_2$ 是厄米的。厄米性保证本征能量是实数，这是哈密顿量必须满足的基本要求。

8. 这里的 $H_2$ 是一个便于教学和后续光学响应计算的模型，并不是唯一可能的模型。理解求解原理之后，可以根据具体研究目标对它进行修改、简化，或者构造新的模型。

## 3. 连续模型的求解方法

连续模型通常难以直接得到解析解，因此需要使用数值方法。这里以 $H_2$ 为例，介绍一种常用的平面波基组方法。

### 3.1 平面波基组与周期性

数值求解的核心任务，是把连续哈密顿量写成某个基组下的矩阵，然后对角化这个矩阵。这里选择平面波基组。由于 $H_2$ 作用在两分量波函数上，基矢不仅要标记波矢，还要标记分量（$A$ 或 $B$）：

$$
|A,q\rangle =
\begin{pmatrix}
e^{iqx} \\
0
\end{pmatrix},
\quad
|B,q\rangle =
\begin{pmatrix}
0 \\
e^{iqx}
\end{pmatrix}.
$$

其中 $q$ 是平面波波矢（与动量相差系数 $\hbar$，下文不再区分两者）。动量算符在这个基组中是对角的，

$$
p\,e^{iqx}=\hbar q\,e^{iqx}.
$$

重要的是周期势如何作用在平面波上。$H_2$ 中的周期势和耦合项满足

$$
V_A(x+a_0)=V_A(x),\quad
V_B(x+a_0)=V_B(x),\quad
W(x+a_0)=W(x),
$$

并且只包含 $e^{ibx}$ 和 $e^{-ibx}$ 这两个傅里叶分量，其中

$$
b=\frac{2\pi}{a_0}.
$$

因此，周期势不会把一个平面波耦合到任意动量，而只会把动量 $q$ 耦合到 $q+b$ 或 $q-b$。例如，

$$
e^{ibx}e^{iqx}=e^{i(q+b)x},
\quad
e^{-ibx}e^{iqx}=e^{i(q-b)x}.
$$

这说明：如果从某个动量 $k$ 出发，哈密顿量只会把它耦合到

$$
k,\quad k\pm b,\quad k\pm 2b,\quad \cdots
$$

这些相差倒格矢的平面波，而不会耦合到其他无关的动量。

### 3.2 不变子空间与布里渊区

令

$$
G_\ell=\ell b,\quad \ell\in\mathbb{Z}.
$$

对于任意给定的 $k$，考虑下面这个由相差倒格矢的平面波张成的子空间：

$$
\mathcal{S}_k=
\mathrm{span}\left\{
|A,k+G_\ell\rangle,\ |B,k+G_\ell\rangle
\; \middle|\; \ell\in\mathbb{Z}
\right\}.
$$

由于 $H_2$ 中的动能项不改变平面波动量，而周期势只把动量改变整数倍的 $b$，所以

$$
H_2\mathcal{S}_k\subseteq \mathcal{S}_k.
$$

也就是说，$\mathcal{S}_k$ 是哈密顿量的不变子空间。这个结论是数值求解的关键：我们不需要一次性处理所有平面波，而可以对每一个 $k$，只在对应的 $\mathcal{S}_k$ 中写出哈密顿量矩阵并对角化。

此外，$k$ 和 $k+b$ 给出的子空间其实相同：

$$
\mathcal{S}_{k+b}=\mathcal{S}_k.
$$

因此 $k$ 只需要取一个长度为 $b$ 的区间。这个区间称为第一布里渊区。在一维情况下，一个常用选择是

$$
k\in \left[-\frac{b}{2},\frac{b}{2}\right).
$$

实际数值计算时，我们在第一布里渊区中取一组离散的 $k$ 点。对每一个 $k$，构造 $H_2$ 在 $\mathcal{S}_k$ 中的矩阵并求本征值。把不同 $k$ 点得到的本征能量连起来，就得到能带结构 $E_n(k)$。

从这个角度看，Bloch 定理不是额外假设，而是上面分块结构的直接结果。$\mathcal{S}_k$ 中的一般两分量态可以写成

$$
\psi_{n k}(x)=
\sum_\ell
\begin{pmatrix}
c_{A,\ell}^{(n k)} \\
c_{B,\ell}^{(n k)}
\end{pmatrix}
e^{i(k+G_\ell)x}.
$$

把公共因子 $e^{ikx}$ 提出来，就得到

$$
\psi_{n k}(x)
=e^{ikx}
\begin{pmatrix}
u_{A,n k}(x) \\
u_{B,n k}(x)
\end{pmatrix},
$$

其中

$$
u_{A,n k}(x)=\sum_\ell c_{A,\ell}^{(n k)}e^{iG_\ell x},
\quad
u_{B,n k}(x)=\sum_\ell c_{B,\ell}^{(n k)}e^{iG_\ell x}.
$$

因为 $G_\ell=\ell b$，所以 $u_{A,n k}(x)$ 和 $u_{B,n k}(x)$ 都满足周期性

$$
u_{A,n k}(x+a_0)=u_{A,n k}(x),\quad
u_{B,n k}(x+a_0)=u_{B,n k}(x).
$$

这就是两分量形式的 Bloch 波函数。

### 3.3 两分量波函数与有限矩阵

在固定的 $k$ 下，$\mathcal{S}_k$ 中的自然基组是

$$
|A,k+G_\ell\rangle,\quad |B,k+G_\ell\rangle,
\quad \ell\in\mathbb{Z}.
$$

如果保留所有整数 $\ell$，矩阵是无限维的。数值计算中需要先截断，只保留有限个 $\ell$。不过在写矩阵元之前，先记住一个重要事实：对于每一个保留下来的 $\ell$，都有 $A$ 和 $B$ 两个基函数。因此如果保留 $N_G$ 个不同的 $G_\ell$，矩阵维度是 $2N_G$，而不是 $N_G$。

### 3.4 矩阵元

为了写出数值计算中的矩阵元，我们令

$$
G_\ell=\ell b.
$$

在固定的 $k$ 下，动量算符作用在平面波上给出

$$
p\,e^{i(k+G_\ell)x}
=\hbar(k+G_\ell)e^{i(k+G_\ell)x}.
$$

因此，动能项在平面波基组中是对角的。下面把矩阵元记为（把 $\alpha,\ell$ 理解为矩阵行指标，把$\beta,\ell'$ 理解为矩阵列指标）

$$
H_{\alpha\beta}(\ell,\ell')
=
\langle \alpha,k+G_\ell|H_2|\beta,k+G_{\ell'}\rangle,
\quad
\alpha,\beta=A,B.
$$

矩阵的两个对角块为

$$
H_{AA}(\ell,\ell')
=
\left[
\frac{\hbar^2(k+G_\ell-\kappa)^2}{2m}
+v_0
\right]\delta_{\ell,\ell'}
+v_A\delta_{\ell,\ell'+1}
+v_A^*\delta_{\ell,\ell'-1},
$$

$$
H_{BB}(\ell,\ell')
=
\left[
-\frac{\hbar^2(k+G_\ell)^2}{2m}
-v_0
\right]\delta_{\ell,\ell'}
-v_B\delta_{\ell,\ell'+1}
-v_B^*\delta_{\ell,\ell'-1}.
$$

两个非对角块来自 $W(x)$ 和 $W^*(x)$：

$$
H_{AB}(\ell,\ell')
=
w_1\delta_{\ell,\ell'+1}
+w_2\delta_{\ell,\ell'-1},
$$

$$
H_{BA}(\ell,\ell')
=
w_2^*\delta_{\ell,\ell'+1}
+w_1^*\delta_{\ell,\ell'-1}.
$$

这里的 Kronecker delta $\delta_{\ell,\ell'}$ 表示

$$
\delta_{\ell,\ell'}=
\begin{cases}
1, & \ell=\ell',\\
0, & \ell\neq \ell'.
\end{cases}
$$

这些 delta 函数来自平面波之间的正交性。考虑周期势中的一个傅里叶分量 $e^{ibx}$。在固定的 $\mathcal{S}_k$ 中，它的矩阵元为

$$
\langle k+G_\ell|e^{ibx}|k+G_{\ell'}\rangle
=\frac{1}{a_0}\int_0^{a_0}
e^{-i(k+G_\ell)x}e^{ibx}e^{i(k+G_{\ell'})x}\,dx.
$$

其中 $k$ 在指数中相互抵消。代入 $G_\ell=\ell b$ 和 $G_{\ell'}=\ell' b$，得到

$$
\langle k+G_\ell|e^{ibx}|k+G_{\ell'}\rangle
=
\frac{1}{a_0}\int_0^{a_0}
e^{i(\ell'+1-\ell)bx}\,dx
=\delta_{\ell,\ell'+1}.
$$

类似地，

$$
\langle k+G_\ell|e^{-ibx}|k+G_{\ell'}\rangle
=\delta_{\ell,\ell'-1}.
$$

因此，$e^{ibx}$ 只把 $k+G_{\ell'}$ 耦合到 $k+G_{\ell'+1}$，而 $e^{-ibx}$ 只把 $k+G_{\ell'}$ 耦合到 $k+G_{\ell'-1}$。这正是周期势在平面波基组中的作用：它不会任意耦合两个动量，而是只耦合相差一个倒格矢的平面波。对于包含更高傅里叶分量的周期势，同样的逻辑会给出相差多个倒格矢的耦合。

### 3.5 平面波截断

上面的展开原则上包含无穷多个倒格矢 $G=\ell b$，因此矩阵维度是无穷大的。数值计算中必须进行截断。常见做法是只保留

$$
\ell=-L,-L+1,\cdots,L-1,L.
$$

这样一共有

$$
N_G=2L+1
$$

个平面波。由于波函数有 $A$、$B$ 两个分量，所以最终需要对角化的矩阵大小是

$$
2N_G\times 2N_G=2(2L+1)\times 2(2L+1).
$$

$L$ 越大，保留的平面波越多，结果通常越精确，但计算量也越大。实际计算时应该检查收敛性：逐渐增大 $L$，如果关心的低能能带几乎不再变化，就说明截断已经足够。截断边界之外的平面波被舍去，因此靠近截断边界的高能本征态通常不如低能本征态可靠。

### 3.6 数值流程

完整的数值求解流程可以概括如下。伪代码中的 `A,ell` 和 `B,ell` 表示矩阵中的两个分量索引；实际写程序时，可以把它们映射成普通的整数矩阵下标。

```text
输入模型参数：
    a0, m, kappa, v0, vA, vB, w1, w2

计算倒格矢：
    b = 2*pi/a0

选择平面波截断：
    ell_list = [-L, -L+1, ..., L]
    NG = 2*L + 1
    matrix_size = 2*NG

选择第一布里渊区中的 k 点：
    k_list = linspace(-b/2, b/2, Nk, endpoint=False)

对每一个 k:
    初始化 H 为 matrix_size x matrix_size 的零矩阵

    对每一对 ell, ell':
        G  = ell  * b

        填入 AA 块：
            kinetic_A = hbar^2 * (k + G - kappa)^2 / (2*m)
            H[A,ell; A,ell'] += (kinetic_A + v0) * delta(ell, ell')
            H[A,ell; A,ell'] += vA      * delta(ell, ell' + 1)
            H[A,ell; A,ell'] += conj(vA)* delta(ell, ell' - 1)

        填入 BB 块：
            kinetic_B = -hbar^2 * (k + G)^2 / (2*m)
            H[B,ell; B,ell'] += (kinetic_B - v0) * delta(ell, ell')
            H[B,ell; B,ell'] += -vB       * delta(ell, ell' + 1)
            H[B,ell; B,ell'] += -conj(vB) * delta(ell, ell' - 1)

        填入 AB 块：
            H[A,ell; B,ell'] += w1 * delta(ell, ell' + 1)
            H[A,ell; B,ell'] += w2 * delta(ell, ell' - 1)

        填入 BA 块：
            H[B,ell; A,ell'] += conj(w2) * delta(ell, ell' + 1)
            H[B,ell; A,ell'] += conj(w1) * delta(ell, ell' - 1)

    对角化 H:
        eigenvalues, eigenvectors = diagonalize(H)

    保存 eigenvalues 作为该 k 点的能带能量
    保存 eigenvectors 作为该 k 点的本征波函数系数

输出：
    E_n(k) 以及对应的两分量本征态
```

这里的本征值给出能带能量，本征矢给出平面波展开系数

$$
\left\{c_{A,G}^{(n k)}, c_{B,G}^{(n k)}\right\}.
$$

有了这些系数，就可以进一步计算波函数、速度矩阵元、贝里联络等物理量。

### 3.7 参考结果

![能带参考结果](figures/band.png)

说明：

1. 图中使用的参数单位如下：$a_0$ 的单位是 $\text{\AA}$，$m$ 以电子质量为单位，$\kappa$ 以 $b$ 为单位，其余能量参数的单位均为 meV。

2. 参数的选取具有一定任意性，可以根据后续研究目标进行调整。这里给出的结果主要用于检查程序流程和能带形状是否合理。

3. 图中的能带是能量最靠近 $0$ meV 的部分能带。当 $N_G$ 无穷大时，图中未画出的更高能量、更低能量能带都有无限多条。
