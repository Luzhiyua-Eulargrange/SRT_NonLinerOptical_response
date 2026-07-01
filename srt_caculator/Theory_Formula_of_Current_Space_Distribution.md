# 电流分布理论公式的推导

## 最小耦合法

首先引入推导电流的方法。由于规范矢量场的性质，场量必然耦合到一个守恒流上，于是有
$$
    j_\mu = \frac{\delta \mathcal{L}}{\delta A^\mu}
$$
从而空间分量有
$$
    j = -\frac{\delta H}{\delta A}
$$
譬如对于最简单的非相对论多体模型

$$
    H = \frac{k^2}{2m} c^\dagger_{k} c_{k} + \hat{V}(x) = \frac{1}{2m}(\nabla \Psi(x))^\dagger \nabla\Psi(x) + V(x) \Psi(x)^\dagger\Psi(x)
$$

其中对应的产生湮灭算符
$$
    \Psi(x) = \int dk e^{ikx} u_k(x) c_k \, c_k = \int dx e^{-ikx} u_k^*(x) \Psi(x)
$$
其中$e^{ikx}u_k(x)$是满足正交关系的单粒子波函数。
引入最小耦合得到
$$
    \nabla \rightarrow \nabla + ieA
$$
带入Hamiltonian有
$$
    H = \frac{1}{2m}(\nabla -ieA)\Psi(x)^\dagger (\nabla+ieA)\Psi(x) + V(x) \Psi(x)^\dagger\Psi(x)
    = H_0 -i\frac{e}{2m}[\Psi(x)^\dagger\nabla\Psi(x)-(\nabla\Psi(x)^\dagger)\Psi(x)]A + \frac{e^2}{2m}\Psi(x)^\dagger\Psi(x) A^2
$$
变分，即得到与量子力学相同的顺磁电流和与高阶电磁场相关的抗磁电流
$$
    j(x) = i\frac{e}{2m}[\Psi(x)^\dagger\nabla\Psi(x)-(\nabla\Psi(x)^\dagger)\Psi(x)] - \frac{e^2A}{m}\Psi(x)^\dagger\Psi(x)
$$

## 双层(AB)模型的电流密度算符

现有的双层Hamiltonian模型如下
$$
H_{qm}=
\begin{pmatrix}
\dfrac{(p - \kappa)^2}{2m}+V_A(x)+v_0 & W(x)\\
W^*(x) & -\dfrac{p^2}{2m}-V_B(x)-v_0
\end{pmatrix},
$$

其中

$$
V_A(x)=v_Ae^{ibx}+v_A^*e^{-ibx},\quad
V_B(x)=v_Be^{ibx}+v_B^*e^{-ibx},
$$

$$
W(x)=w_1e^{ibx}+w_2e^{-ibx},\qquad b=\frac{2\pi}{a_0}.
$$

下面约定 $e>0$，电子电荷为 $-e$。电子的自旋已经被忽略，用标量场描述。为了应用场论方法，首先做二次量子化。注意，这里为了使动量的定义明确，便于引入最小耦合方法，采取平面波基下的量子化。观察Hamiltonian不难发现

$$
\begin{aligned}
  H &= \left(\frac{(p-\kappa)^2}{2m}+v_0\right) c_{Ap}^\dagger c_{Ap} -(\dfrac{p^2}{2m} + v_0) c_{Bp}^\dagger c_{Bp}\\
    &+\left( v_A\hat c^\dagger_{A,q+b}\hat c_{Aq} + v_A^*\hat c^\dagger_{A,q-b}\hat c_{Aq} \right) + \left(v_B\hat c^\dagger_{B,q+b}\hat c_{Bq} + v_B^*\hat c^\dagger_{B,q-b}\hat c_{Bq} \right)\\
    &+\left(w_1\hat c^\dagger_{A,q+b}\hat c_{Bq}+w_2\hat c^\dagger_{A,q-b}\hat c_{Bq}+w_1^*\hat c^\dagger_{B,q-b}\hat c_{Aq}+w_2^*\hat c^\dagger_{B,q+b}\hat c_{Aq}\right)
\end{aligned}
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
Hamiltonian就写作：

$$
\begin{aligned}
H &= \sum_{n} \int dk \left[
\epsilon_A(n,k)\, \hat c^\dagger_{A,nk}\hat c_{A,nk}
+ \epsilon_B(n,k)\, \hat c^\dagger_{B,nk}\hat c_{B,nk}
\right] \\
&+ \sum_{n} \int dk \Big[
v_A\, \hat c^\dagger_{A,n+1,k}\hat c_{A,nk}
+ v_A^*\, \hat c^\dagger_{A,n-1,k}\hat c_{A,nk}
+ v_B\, \hat c^\dagger_{B,n+1,k}\hat c_{B,nk}
+ v_B^*\, \hat c^\dagger_{B,n-1,k}\hat c_{B,nk}
\Big] \\
&+ \sum_{n} \int dk \Big[
w_1\, \hat c^\dagger_{A,n+1,k}\hat c_{B,nk}
+ w_2\, \hat c^\dagger_{A,n-1,k}\hat c_{B,nk}
+ w_1^*\, \hat c^\dagger_{B,n-1,k}\hat c_{A,nk}
+ w_2^*\, \hat c^\dagger_{B,n+1,k}\hat c_{A,nk}
\Big],
\end{aligned}
$$

其中单粒子能量为：
$$
\epsilon_A(n,k) = \frac{(k + G_n - \kappa)^2}{2m} + v_0,
$$
$$
\epsilon_B(n,k) = -\frac{(k + G_n)^2}{2m} - v_0.
$$

Fourier变换为
$$
\Psi_A(x) = \sum_{n} \int_{k \in FBZ} e^{i(k+G_n)x}\, u_{A,nk}\, \hat c_{A,nk},\\
\Psi_B(x) = \sum_{n} \int_{k \in FBZ} e^{i(k+G_n)x}\, u_{B,nk}\, \hat c_{B,nk},
$$

$$
\hat c_{A,nk} = \int_0^{a_0} dx\, e^{-i(k+G_n)x} u_{A,nk}^*\, \Psi_A(x),\\
\hat c_{B,nk} = \int_0^{a_0} dx\, e^{-i(k+G_n)x} u_{B,nk}^*\, \Psi_B(x).
$$
由于是标量场平面波，单粒子波函数$u_{AB,nk} = u$只是一个常数，有
$$
-i(k+G_n) \hat c_{A,nk} = \int_0^{a_0} dx\, \nabla e^{-i(k+G_n)x} u^*\, \Psi_A(x) = -\int_0^{a_0} dx\,  e^{-i(k+G_n)x} u^*\, \nabla \Psi_A(x)
$$
$$
-i(k+G_n) \hat c_{B,nk} = \int_0^{a_0} dx\, \nabla e^{-i(k+G_n)x} u^*\, \Psi_B(x) = -\int_0^{a_0} dx\,  e^{-i(k+G_n)x} u^*\, \nabla \Psi_B(x)
$$

带入Hamiltonian，得到

$$
    H_{diag} = \frac{1}{2m} \int dx [(-i\nabla -\kappa)\Psi_A(x)]^\dagger[(-i\nabla -\kappa)\Psi_A(x)] + \int dx \, v_0 \Psi_A(x)^\dagger \Psi_A(x) - \frac{1}{2m} \int dx (\nabla\Psi_B(x))^\dagger (\nabla\Psi_B(x)) - \int dx \, v_0 \Psi_B(x)^\dagger \Psi_B(x)
$$

$$
    H_V = v_A \int dx \, e^{ibx} \Psi(x)_A^\dagger \Psi_A(x) + \cdots;\ H_W = w_1 \int dx \, e^{ibx} \Psi_A(x)^\dagger \Psi_B(x) + \cdots
$$

现在动量只集中在对角块上，可以方便地引入耦合 $\nabla \rightarrow \nabla + ieA$, 从而可以得到Hamiltonian density

$$
    \mathcal{H}(A) = \mathcal{H}_0 + \frac{eA}{2m} [-i(\Psi_A^\dagger \nabla \Psi_A - (\nabla\Psi_A^\dagger)  \Psi_A) - 2\kappa\Psi_A^\dagger\Psi_A] + \frac{e^2 A^2}{2m} \Psi_A^\dagger \Psi_A \qquad\\
    -\frac{eA}{2m} [-i(\Psi_B^\dagger \nabla \Psi_B - (\nabla\Psi_B^\dagger)  \Psi_B)] - \frac{e^2 A^2}{2m} \Psi_B^\dagger \Psi_B
$$

变分即得到电流密度算符

$$
    j(x) = \frac{e}{2m} [i(\Psi_A^\dagger \nabla \Psi_A - (\nabla\Psi_A^\dagger)  \Psi_A) + 2\kappa\Psi_A^\dagger\Psi_A] - \frac{e^2 A}{m} \Psi_A^\dagger \Psi_A - \frac{e}{2m} i(\Psi_B^\dagger \nabla \Psi_B - (\nabla\Psi_B^\dagger)  \Psi_B) + \frac{e^2 A}{m} \Psi_B^\dagger \Psi_B
$$

## 能带基展开

考虑傅里叶变换
$$
    c^\dagger_{ks} = \int dx \psi(x) \Psi^\dagger(x) = \int dx e^{ikx} (u_{ks.A}(x) \Psi_A^\dagger(x) + u_{ks.B}(x) \Psi_B^\dagger(x))
$$
相应的逆变换
$$
\Psi_A(x) = \sum_{s} \int_{\mathrm{FBZ}} dk \, e^{ikx} u_{\mathbf{k}s,A}(x) \, c_{\mathbf{k}s},
$$

$$
\Psi_B(x) = \sum_{s} \int_{\mathrm{FBZ}} dk \, e^{ikx} u_{\mathbf{k}s,B}(x) \, c_{\mathbf{k}s}.
$$

对于能带基展开，得到电流算符
$$
    j(x) = \sum\limits_{ss'} \int dk dk' e^{i(k'-k)x} c^\dagger_{ks}c_{k's'} \big[\qquad\\
    \frac{e}{2m} [-(k+k')(u_{ks,A}^*u_{k's',A} - u_{ks,B}^*u_{k's',B}) + i u_{ks,A}^* \overset{\leftrightarrow}{\nabla} u_{k's',A}- i u_{ks,B}^* \overset{\leftrightarrow}{\nabla} u_{k's',B}]\\
    +\frac{e\kappa}{m} (u_{ks,A}^*u_{k's',A})
    - \frac{e^2 A}{m} (u_{ks,A}^*u_{k's',A} - u_{ks,B}^*u_{k's',B})
    \big]
$$
密度矩阵
$$
    \rho_{kk'ss'} = \langle c_{ks}^\dagger c_{k's'} \rangle = \delta(k-k') \langle c_{ks}^\dagger c_{ks'} \rangle = \delta(k-k')\rho_{kss'}
$$
从而电流密度的期望值
$$
    \langle j(x, t) \rangle = \sum\limits_{ss'} \int dk \rho_{kss'}(t) \big[\qquad\\
    \frac{e}{2m} [-2k(u_{ks,A}^*u_{ks',A} - u_{ks,B}^*u_{ks',B}) + i u_{ks,A}^* \overset{\leftrightarrow}{\nabla} u_{ks',A}- i u_{ks,B}^* \overset{\leftrightarrow}{\nabla} u_{ks',B}]\\
    +\frac{e\kappa}{m} (u_{ks,A}^*u_{ks',A})
    - \frac{e^2 A}{m} (u_{ks,A}^*u_{ks',A} - u_{ks,B}^*u_{ks',B})
    \big]
$$
从这个式子就可以推导电流的分布了。对于主程序使用的长度规范，只需把A取为0，同时修改对应产生湮灭算符的定义即可。
之前的文献中给出的是对一个晶格内电流密度的积分
$$
    \langle J^E(t) \rangle = -e \int \frac{dk}{2\pi} \text{Tr}[V^E(k)\rho^E(k)]
$$
不够精细。故做如上的推导。
