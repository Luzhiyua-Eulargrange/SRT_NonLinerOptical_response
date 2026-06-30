# 从 H2 模型推导空间电流分布

本文从 `1d_model.md` 中的一维两分量连续模型出发，推导适用于该模型的局域电流算符和空间电流分布公式。

模型哈密顿量为

$$
H_2=
\begin{pmatrix}
\dfrac{(p-\hbar\kappa)^2}{2m}+V_A(x)+v_0 & W(x)\\
W^*(x) & -\dfrac{p^2}{2m}-V_B(x)-v_0
\end{pmatrix},
\qquad p=-i\hbar\partial_x .
$$

其中

$$
V_A(x)=v_Ae^{ibx}+v_A^*e^{-ibx},\quad
V_B(x)=v_Be^{ibx}+v_B^*e^{-ibx},
$$

$$
W(x)=w_1e^{ibx}+w_2e^{-ibx},\qquad b=\frac{2\pi}{a_0}.
$$

下面约定 $e>0$，电子电荷为 $-e$。自旋自由度已经忽略，因此不额外乘自旋简并因子 2。

## 1. 为什么不能直接使用普通电流算符

普通单分量模型

$$
H=\frac{p^2}{2m}+V(x)
$$

对应的电子电荷流为

$$
j(x)=
i\frac{e\hbar}{2m}
\left[
\psi^*(x)\partial_x\psi(x)
-
\left(\partial_x\psi^*(x)\right)\psi(x)
\right].
$$

这个公式只适用于正号动能 $p^2/(2m)$。当前 $H_2$ 模型中：

1. A 分量的动量是 $p-\hbar\kappa$，不是 $p$。
2. B 分量的动能是 $-p^2/(2m)$，符号与普通电子相反。
3. $V_A,V_B,W$ 都是局域位置算符，不含动量导数，因此它们不直接产生输运电流。

所以必须从 $H_2$ 的动量依赖项重新推导电流算符。

## 2. 用最小耦合定义局域电流

引入外加矢势 $A_{em}(x,t)$。电子电荷为 $-e$，所以最小耦合为

$$
p\rightarrow p+eA_{em}(x,t).
$$

因此

$$
H_A[A_{em}]
=
\frac{(p+eA_{em}-\hbar\kappa)^2}{2m},
$$

$$
H_B[A_{em}]
=
-\frac{(p+eA_{em})^2}{2m}.
$$

电荷流密度算符定义为

$$
\hat j(x)=-\left.\frac{\delta \hat H[A_{em}]}{\delta A_{em}(x)}\right|_{A_{em}=0}.
$$

写成场算符 $\Psi(x)=(\Psi_A(x),\Psi_B(x))^T$ 后，总局域电流为

$$
\hat j(x)=\hat j_A(x)+\hat j_B(x),
$$

其中

$$
\boxed{
\hat j_A(x)=
i\frac{e\hbar}{2m}
\left[
\Psi_A^\dagger\partial_x\Psi_A
-
\left(\partial_x\Psi_A^\dagger\right)\Psi_A
\right]
+
e\frac{\hbar\kappa}{m}\Psi_A^\dagger\Psi_A
}
$$

以及

$$
\boxed{
\hat j_B(x)=
-i\frac{e\hbar}{2m}
\left[
\Psi_B^\dagger\partial_x\Psi_B
-
\left(\partial_x\Psi_B^\dagger\right)\Psi_B
\right].
}
$$

等价地，引入

$$
\eta_A=+1,\quad \eta_B=-1,
\qquad
\kappa_A=\kappa,\quad \kappa_B=0,
$$

则两个分量可以合写为

$$
\boxed{
\hat j_\alpha(x)=
\eta_\alpha i\frac{e\hbar}{2m}
\left[
\Psi_\alpha^\dagger\partial_x\Psi_\alpha
-
\left(\partial_x\Psi_\alpha^\dagger\right)\Psi_\alpha
\right]
+
\eta_\alpha e\frac{\hbar\kappa_\alpha}{m}
\Psi_\alpha^\dagger\Psi_\alpha ,
\quad \alpha=A,B.
}
$$

这个公式是当前模型的基本局域电流算符。

## 3. 连续性方程检查

定义粒子数密度

$$
\hat n(x)=\Psi_A^\dagger\Psi_A+\Psi_B^\dagger\Psi_B,
$$

电荷密度为

$$
\hat\rho_q(x)=-e\hat n(x).
$$

由 Heisenberg 方程可得

$$
\partial_t\hat\rho_q(x)+\partial_x\hat j(x)=0.
$$

局域势 $V_A,V_B$ 与密度对易，不产生电流。非对角耦合 $W(x)$ 会在 A、B 分量各自的连续性方程中产生相反的局域转移项，但在总电荷连续性方程中互相抵消，因此它也不作为额外的空间输运电流出现。

## 4. Bloch band basis 中的空间电流

Bloch 本征态写为

$$
\psi_{kn}(x)=e^{ikx}
\begin{pmatrix}
u_{A,kn}(x)\\
u_{B,kn}(x)
\end{pmatrix},
\qquad
u_{\alpha,kn}(x+a_0)=u_{\alpha,kn}(x).
$$

场算符展开为

$$
\Psi_\alpha(x)=
\sum_n\int_{\mathrm{FBZ}}\frac{dk}{2\pi}\,
e^{ikx}u_{\alpha,kn}(x)c_{kn}.
$$

采用如下 band-basis 单体密度矩阵约定：

$$
\Gamma_{nm}(k,t)=\langle c^\dagger_{kn}c_{km}\rangle .
$$

对于空间均匀外场和平移不变初态，不同 $k$ 块之间不相干：

$$
\langle c^\dagger_{kn}c_{k'm}\rangle
=2\pi\delta(k-k')\Gamma_{nm}(k,t).
$$

把 Bloch 展开代入局域电流算符，得到

$$
\boxed{
j(x,t)=
\int_{\mathrm{FBZ}}\frac{dk}{2\pi}
\sum_{n,m}\Gamma_{nm}(k,t)
\sum_{\alpha=A,B}
\mathcal J^{(\alpha)}_{nm}(k,x)
}
$$

其中局域电流核为

$$
\boxed{
\mathcal J^{(\alpha)}_{nm}(k,x)=
\eta_\alpha i\frac{e\hbar}{2m}
\left[
u^*_{\alpha,kn}\partial_xu_{\alpha,km}
-
\left(\partial_xu^*_{\alpha,kn}\right)u_{\alpha,km}
\right]
-
\eta_\alpha e\frac{\hbar}{m}(k-\kappa_\alpha)
u^*_{\alpha,kn}u_{\alpha,km}.
}
$$

上式中的 $-e\hbar(k-\kappa_\alpha)/m$ 项来自 Bloch 因子 $e^{ikx}$ 和 A 分量的动量平移。它不能省略。

## 5. 平面波基中的可计算公式

展开周期部分：

$$
u_{\alpha,kn}(x)=
\sum_\ell U_{\alpha\ell,n}(k)e^{iG_\ell x},
\qquad
G_\ell=\ell b.
$$

定义平面波基密度矩阵

$$
\rho^{\mathrm{pw}}_{\alpha\ell,\beta\ell'}(k,t)
=
\sum_{n,m}
U_{\alpha\ell,m}(k)\,
\Gamma_{nm}(k,t)\,
U^*_{\beta\ell',n}(k).
$$

空间电流只需要同一分量内的平面波密度矩阵块：

$$
\rho^{\mathrm{pw}}_{\alpha\ell,\alpha\ell'}(k,t).
$$

代入后得到适合数值计算的晶胞内电流分布：

$$
\boxed{
j(x,t)=
-e\frac{\hbar}{2m}
\int_{\mathrm{FBZ}}\frac{dk}{2\pi}
\sum_{\alpha=A,B}
\sum_{\ell,\ell'}
\eta_\alpha
\left(
2k+G_\ell+G_{\ell'}-2\kappa_\alpha
\right)
e^{i(G_\ell-G_{\ell'})x}
\rho^{\mathrm{pw}}_{\alpha\ell,\alpha\ell'}(k,t).
}
$$

数值上 $j(x,t)$ 应取实部：

$$
j_{\mathrm{num}}(x,t)=\mathrm{Re}\,j(x,t).
$$

若平面波基函数采用 $e^{i(k+G)x}/\sqrt{a_0}$ 的归一化，上式还需要整体乘 $1/a_0$。如果使用未显式归一化的胞内周期函数，需保证密度和电流使用同一归一化约定。

## 6. 晶胞平均检查

对一个晶胞平均：

$$
\frac{1}{a_0}\int_0^{a_0}dx\,
e^{i(G_\ell-G_{\ell'})x}
=\delta_{\ell,\ell'}.
$$

因此

$$
\boxed{
\bar j(t)=
-e
\int_{\mathrm{FBZ}}\frac{dk}{2\pi}
\sum_{\alpha=A,B}
\sum_\ell
\eta_\alpha
\frac{\hbar}{m}(k+G_\ell-\kappa_\alpha)
\rho^{\mathrm{pw}}_{\alpha\ell,\alpha\ell}(k,t).
}
$$

这与从 $H_2(k)$ 的动量依赖项取

$$
v(k)=\frac{1}{\hbar}\frac{\partial H_2(k)}{\partial k}
$$

得到的总电流一致。具体地：

$$
v_{A\ell}=\frac{\hbar}{m}(k+G_\ell-\kappa),
\qquad
v_{B\ell}=-\frac{\hbar}{m}(k+G_\ell).
$$

乘以电子电荷 $-e$ 后，正好给出上面的 $\bar j(t)$。

## 7. 实现时的输入输出

若已有 band basis RDM $\Gamma(k,t)$ 和本征矢 $U(k)$，计算空间电流分布的步骤是：

1. 对每个 $k,t$，先变换到平面波基：

   $$
   \rho^{\mathrm{pw}}(k,t)=U(k)\Gamma(k,t)U^\dagger(k).
   $$

2. 分别取 A/A 和 B/B 块。

3. 在给定 $x$ 网格上累加

   $$
   \eta_\alpha
   \left(
   2k+G_\ell+G_{\ell'}-2\kappa_\alpha
   \right)
   e^{i(G_\ell-G_{\ell'})x}
   \rho^{\mathrm{pw}}_{\alpha\ell,\alpha\ell'}(k,t).
   $$

4. 乘整体因子 $-e\hbar/(2m)$、k 空间积分权重，并取实部。

这个公式才是当前 $H_2$ 模型下的空间电流分布公式；普通单分量 $p^2/2m$ 的电流算符只是在 $\eta_A=1$、无 B 分量、$\kappa=0$ 时的特殊情形。
