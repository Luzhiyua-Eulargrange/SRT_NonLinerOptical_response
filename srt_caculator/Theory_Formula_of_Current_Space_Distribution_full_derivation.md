# 电流空间分布公式：从 $c_k$ 量子化到密度矩阵

本文保留原始单粒子 Hamiltonian

$$
H_2=
\begin{pmatrix}
\dfrac{(p-\hbar\kappa)^2}{2m}+V_A(x)+v_0 & W(x)\\
W^*(x) & -\dfrac{p^2}{2m}-V_B(x)-v_0
\end{pmatrix},
\qquad p=-i\hbar\partial_x ,
$$

其中

$$
V_A(x)=v_Ae^{ibx}+v_A^*e^{-ibx},\quad
V_B(x)=v_Be^{ibx}+v_B^*e^{-ibx},
$$

$$
W(x)=w_1e^{ibx}+w_2e^{-ibx},\qquad b=\frac{2\pi}{a_0}.
$$

约定 $e>0$，电子电荷为 $-e$。自旋已经忽略。下面先用平面波产生湮灭算符 $c_k$ 做二次量子化，再变换到能带基，最后与单体密度矩阵合并得到空间电流分布。

## 1. 用 $c_k$ 量子化 Hamiltonian

引入两分量场算符

$$
\hat\Psi(x)=
\begin{pmatrix}
\hat\psi_A(x)\\
\hat\psi_B(x)
\end{pmatrix},
\qquad
\hat\psi_\alpha(x)=\frac{1}{\sqrt{L_x}}\sum_q e^{iqx}\hat c_{\alpha q},
\quad \alpha=A,B ,
$$

满足

$$
\{\hat c_{\alpha q},\hat c^\dagger_{\beta q'}\}
=\delta_{\alpha\beta}\delta_{q,q'},\qquad
\{\hat c_{\alpha q},\hat c_{\beta q'}\}=0 .
$$

二次量子化 Hamiltonian 为

$$
\hat H_2=\int dx\,\hat\Psi^\dagger(x)H_2\hat\Psi(x).
$$

先看对角动能和常数项。因为 $p e^{iqx}=\hbar q e^{iqx}$，有

$$
\hat H_{\rm diag}
=
\sum_q
\left[
\frac{\hbar^2(q-\kappa)^2}{2m}+v_0
\right]
\hat c^\dagger_{Aq}\hat c_{Aq}
+
\sum_q
\left[
-\frac{\hbar^2q^2}{2m}-v_0
\right]
\hat c^\dagger_{Bq}\hat c_{Bq}.
$$

周期势只改变动量一个倒格矢 $b$：

$$
\begin{aligned}
\hat H_{V_A}
&=
\sum_q
\left(
v_A\hat c^\dagger_{A,q+b}\hat c_{Aq}
+v_A^*\hat c^\dagger_{A,q-b}\hat c_{Aq}
\right),\\
\hat H_{V_B}
&=
-\sum_q
\left(
v_B\hat c^\dagger_{B,q+b}\hat c_{Bq}
+v_B^*\hat c^\dagger_{B,q-b}\hat c_{Bq}
\right).
\end{aligned}
$$

层间耦合为

$$
\begin{aligned}
\hat H_W
=\sum_q\big(
&w_1\hat c^\dagger_{A,q+b}\hat c_{Bq}
+w_2\hat c^\dagger_{A,q-b}\hat c_{Bq}\\
&+w_1^*\hat c^\dagger_{B,q-b}\hat c_{Aq}
+w_2^*\hat c^\dagger_{B,q+b}\hat c_{Aq}
\big).
\end{aligned}
$$

因此

$$
\hat H_2=\hat H_{\rm diag}+\hat H_{V_A}+\hat H_{V_B}+\hat H_W .
$$

为了得到 Bloch 块形式，把动量写成

$$
q=k+G_n,\qquad G_n=nb,\qquad k\in{\rm BZ}.
$$

定义

$$
\hat c_{A,nk}\equiv \hat c_{A,k+G_n},\qquad
\hat c_{B,nk}\equiv \hat c_{B,k+G_n}.
$$

对每个 $k$，取平面波基列向量

$$
\hat C_k=
\begin{pmatrix}
\cdots\\
\hat c_{A,nk}\\
\cdots\\
\hat c_{B,nk}\\
\cdots
\end{pmatrix}.
$$

这里采用代码中的排序：先所有 $A$ 分量的 $n=-L,\ldots,L$，再所有 $B$ 分量的 $n=-L,\ldots,L$。于是

$$
\hat H_2=\sum_{k\in{\rm BZ}}\hat C_k^\dagger\,\mathcal H(k)\,\hat C_k .
$$

非零矩阵元为

$$
\begin{aligned}
\mathcal H_{A n,A n}(k)
&=\frac{\hbar^2(k+G_n-\kappa)^2}{2m}+v_0,\\
\mathcal H_{B n,B n}(k)
&=-\frac{\hbar^2(k+G_n)^2}{2m}-v_0,\\
\mathcal H_{A,n+1;A,n}(k)&=v_A,\qquad
\mathcal H_{A,n-1;A,n}(k)=v_A^*,\\
\mathcal H_{B,n+1;B,n}(k)&=-v_B,\qquad
\mathcal H_{B,n-1;B,n}(k)=-v_B^*,\\
\mathcal H_{A,n+1;B,n}(k)&=w_1,\qquad
\mathcal H_{A,n-1;B,n}(k)=w_2,\\
\mathcal H_{B,n-1;A,n}(k)&=w_1^*,\qquad
\mathcal H_{B,n+1;A,n}(k)=w_2^* .
\end{aligned}
$$

这就是用 $c_k$ 算符量子化后的二次型费米 Hamiltonian。

## 2. 从最小耦合得到平面波基电流矩阵

电子电荷为 $-e$，因此最小耦合取

$$
p\rightarrow p+eA(x,t).
$$

电流密度算符定义为

$$
\hat j(x)=-\left.\frac{\delta \hat H[A]}{\delta A(x)}\right|_{A=0}.
$$

只有含 $p$ 的动能项直接产生输运电流；$V_A,V_B,W$ 都是不含导数的局域势，不直接贡献电流算符。令

$$
\eta_A=+1,\quad \eta_B=-1,\qquad
\kappa_A=\kappa,\quad \kappa_B=0 .
$$

则局域电流算符可以写成

$$
\boxed{
\hat j_\alpha(x)=
\eta_\alpha i\frac{e\hbar}{2m}
\left[
\hat\psi_\alpha^\dagger\partial_x\hat\psi_\alpha
-\left(\partial_x\hat\psi_\alpha^\dagger\right)\hat\psi_\alpha
\right]
+
\eta_\alpha e\frac{\hbar\kappa_\alpha}{m}
\hat\psi_\alpha^\dagger\hat\psi_\alpha ,
}
$$

并且

$$
\hat j(x)=\hat j_A(x)+\hat j_B(x).
$$

代入平面波展开，在固定 $k$ 的 Bloch 块中可写成单粒子矩阵形式。有限长度体系中局域算符带有场算符归一化给出的整体因子 $1/L_x$；在热力学极限中该因子与 $\sum_k$ 合并为 $\int dk/(2\pi)$。下面的 $\mathcal J^{\rm pw}$ 表示不含这个整体 $k$ 空间权重的电流核：

$$
\hat j(x)=\frac{1}{L_x}\sum_{k\in{\rm BZ}}\hat C_k^\dagger\,\mathcal J^{\rm pw}(k,x)\,\hat C_k .
$$

平面波基电流矩阵的非零元只在同一分量内：

$$
\boxed{
\mathcal J^{\rm pw}_{\alpha n,\beta n'}(k,x)
=
-\delta_{\alpha\beta}\,
\eta_\alpha\frac{e\hbar}{2m}
\left(
2k+G_n+G_{n'}-2\kappa_\alpha
\right)
e^{i(G_{n'}-G_n)x}.
}
$$

这个矩阵是局域电流矩阵。对一个晶胞取平均时，

$$
\frac{1}{a_0}\int_0^{a_0}dx\,e^{i(G_{n'}-G_n)x}
=\delta_{nn'},
$$

于是得到平均电流矩阵

$$
\bar{\mathcal J}^{\rm pw}_{\alpha n,\beta n'}(k)
=
-e\,\delta_{\alpha\beta}\delta_{nn'}\,v_{\alpha n}(k),
$$

其中

$$
v_{A n}(k)=\frac{\hbar}{m}(k+G_n-\kappa),
\qquad
v_{B n}(k)=-\frac{\hbar}{m}(k+G_n).
$$

也就是说

$$
\bar{\mathcal J}^{\rm pw}(k)
=-e\,v^{\rm pw}(k),
\qquad
v^{\rm pw}(k)=\frac{1}{\hbar}\frac{\partial \mathcal H(k)}{\partial k}.
$$

这与代码中 `Geometry.velocity_matrix` 使用的速度矩阵完全一致。

## 3. 带入能带基得到电流矩阵

对每个 $k$ 对角化平面波 Hamiltonian：

$$
\mathcal H(k)U(k)=U(k)\varepsilon(k),
\qquad
U^\dagger(k)U(k)=1.
$$

其中 $U_{\mu s}(k)$ 是第 $s$ 条能带在平面波基 $\mu=(\alpha,n)$ 上的本征矢分量，$\varepsilon(k)$ 为对角矩阵。定义能带基湮灭算符

$$
\hat d_{sk}=\sum_\mu U^*_{\mu s}(k)\hat C_{\mu k},
\qquad
\hat C_{\mu k}=\sum_s U_{\mu s}(k)\hat d_{sk}.
$$

则 Hamiltonian 变为

$$
\hat H_2=\sum_{k,s}\varepsilon_s(k)\hat d^\dagger_{sk}\hat d_{sk}.
$$

局域电流矩阵在能带基中为

$$
\boxed{
\mathcal J^{\rm band}(k,x)
=U^\dagger(k)\mathcal J^{\rm pw}(k,x)U(k).
}
$$

写成指标形式：

$$
\boxed{
\mathcal J^{\rm band}_{ss'}(k,x)
=
\sum_{\mu,\nu}
U^*_{\mu s}(k)
\mathcal J^{\rm pw}_{\mu\nu}(k,x)
U_{\nu s'}(k).
}
$$

由于 $\mathcal J^{\rm pw}$ 只含同一分量内的矩阵元，也可以显式写成

$$
\boxed{
\begin{aligned}
\mathcal J^{\rm band}_{ss'}(k,x)
=
-\frac{e\hbar}{2m}
\sum_{\alpha=A,B}\eta_\alpha
\sum_{n,n'}
&U^*_{\alpha n,s}(k)
\left(
2k+G_n+G_{n'}-2\kappa_\alpha
\right)\\
&\times
e^{i(G_{n'}-G_n)x}
U_{\alpha n',s'}(k).
\end{aligned}
}
$$

晶胞平均后的能带基电流矩阵为

$$
\boxed{
\bar{\mathcal J}^{\rm band}(k)
=
U^\dagger(k)\bar{\mathcal J}^{\rm pw}(k)U(k)
=
-e\,U^\dagger(k)v^{\rm pw}(k)U(k).
}
$$

等价地，在非简并且规范光滑的情形下，对角元满足

$$
\bar{\mathcal J}^{\rm band}_{ss}(k)
=
-\frac{e}{\hbar}\partial_k\varepsilon_s(k),
$$

而非对角元由 $U^\dagger v^{\rm pw}U$ 给出，不能只用能量导数表示。

## 4. 与密度矩阵合并

采用代码中常用的单体密度矩阵约定：

$$
\rho^{\rm band}_{ss'}(k,t)
=
\langle \hat d^\dagger_{s'k}(t)\hat d_{sk}(t)\rangle .
$$

在这个约定下，任意单体算符

$$
\hat O=\sum_{k,s,s'}
\hat d^\dagger_{sk}O^{\rm band}_{ss'}(k)\hat d_{s'k}
$$

的期望值为

$$
\langle \hat O\rangle
=
\sum_k {\rm Tr}\left[
O^{\rm band}(k)\rho^{\rm band}(k,t)
\right].
$$

因此空间电流分布为

$$
\boxed{
j(x,t)
=
\int_{\rm BZ}\frac{dk}{2\pi}\,
{\rm Re}\,
{\rm Tr}\left[
\mathcal J^{\rm band}(k,x)\rho^{\rm band}(k,t)
\right].
}
$$

数值离散时，若 $k_{\rm weight}=\Delta k/(2\pi)$，则

$$
\boxed{
j(x,t)
=
\sum_{k_i}k_{\rm weight}\,
{\rm Re}\,
{\rm Tr}\left[
\mathcal J^{\rm band}(k_i,x)\rho^{\rm band}(k_i,t)
\right].
}
$$

也可以先把密度矩阵从能带基变回平面波基：

$$
\rho^{\rm pw}(k,t)
=U(k)\rho^{\rm band}(k,t)U^\dagger(k).
$$

这里

$$
\rho^{\rm pw}_{\mu\nu}(k,t)
=\langle \hat C^\dagger_{\nu k}(t)\hat C_{\mu k}(t)\rangle ,
$$

所以把 $\mathcal J^{\rm pw}_{\mu\nu}$ 与 $\rho^{\rm pw}_{\nu\mu}$ 收缩后，可将求和指标互换，写成和密度分布代码相同的相位约定 $e^{i(G_n-G_{n'})x}\rho^{\rm pw}_{\alpha n,\alpha n'}$。

此时最终公式可以直接写成

$$
\boxed{
\begin{aligned}
j(x,t)
=
-\frac{e\hbar}{2m}
\int_{\rm BZ}\frac{dk}{2\pi}\,
{\rm Re}
\sum_{\alpha=A,B}\eta_\alpha
\sum_{n,n'}
&
\left(
2k+G_n+G_{n'}-2\kappa_\alpha
\right)\\
&\times
e^{i(G_n-G_{n'})x}
\rho^{\rm pw}_{\alpha n,\alpha n'}(k,t).
\end{aligned}
}
$$

其中

$$
\rho^{\rm pw}_{\alpha n,\alpha n'}(k,t)
=
\left[
U(k)\rho^{\rm band}(k,t)U^\dagger(k)
\right]_{\alpha n,\alpha n'}.
$$

离散形式为

$$
\boxed{
\begin{aligned}
j(x,t)
=
-\frac{e\hbar}{2m}
\sum_{k_i}k_{\rm weight}\,
{\rm Re}
\sum_{\alpha=A,B}\eta_\alpha
\sum_{n,n'}
&
\left(
2k_i+G_n+G_{n'}-2\kappa_\alpha
\right)\\
&\times
e^{i(G_n-G_{n'})x}
\left[
U(k_i)\rho^{\rm band}(k_i,t)U^\dagger(k_i)
\right]_{\alpha n,\alpha n'} .
\end{aligned}
}
$$

这就是所需的电流空间分布理论公式。它与密度分布公式

$$
n(x,t)=
\int_{\rm BZ}\frac{dk}{2\pi}\,
{\rm Re}
\sum_{\alpha=A,B}\sum_{n,n'}
e^{i(G_n-G_{n'})x}
\rho^{\rm pw}_{\alpha n,\alpha n'}(k,t)
$$

只差一个由最小耦合给出的电流权重

$$
-\eta_\alpha\frac{e\hbar}{2m}
\left(
2k+G_n+G_{n'}-2\kappa_\alpha
\right).
$$

对一个晶胞平均后，空间相位只保留 $n=n'$ 项，得到总电流

$$
\boxed{
\bar j(t)
=
-e
\int_{\rm BZ}\frac{dk}{2\pi}\,
{\rm Re}\,
{\rm Tr}
\left[
v^{\rm band}(k)\rho^{\rm band}(k,t)
\right],
}
$$

其中

$$
v^{\rm band}(k)=U^\dagger(k)v^{\rm pw}(k)U(k).
$$

这就是 `Current.py` 中总电流公式的理论来源；上面的 $j(x,t)$ 则是其未做晶胞平均之前的空间分布版本。
