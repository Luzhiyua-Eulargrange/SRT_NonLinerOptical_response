# Velocity Matrix Diagnostics

The main code uses the direct band-basis transform

$$
v_{\mathrm{band}}(k)
= U^\dagger(k)\,\frac{1}{\hbar}\frac{\partial H(k)}{\partial k}\,U(k)
$$

for velocity-gauge time evolution. This is the stable production path.

The paper formula is

$$
v_{kss'} = \frac{1}{\hbar}\,[D_k,H(k)]_{ss'} .
$$

In the energy eigenbasis it can be written as

$$
v_{kss'}
= \frac{1}{\hbar}
\left[
\delta_{ss'}\,\partial_k \epsilon_{ks}
- i\left(\epsilon_{ks'}-\epsilon_{ks}\right)\xi_{kss'}
\right],
$$

where

$$
\xi_{kss'} = i\langle u_{ks}|\partial_k u_{ks'}\rangle
$$

is the Berry connection.

## Why The Two Formulas Are Equivalent

Start from the eigenvalue equation

$$
H(k)|u_s(k)\rangle = \epsilon_s(k)|u_s(k)\rangle .
$$

Differentiate with respect to \(k\):

$$
\frac{\partial H}{\partial k}|u_s\rangle
+ H|\partial_k u_s\rangle
=
(\partial_k \epsilon_s)|u_s\rangle
+ \epsilon_s|\partial_k u_s\rangle .
$$

Project with \(\langle u_{s'}|\). For \(s'=s\),

$$
\langle u_s|\partial_k H|u_s\rangle
= \partial_k \epsilon_s .
$$

For \(s'\ne s\),

$$
\langle u_{s'}|\partial_k H|u_s\rangle
=
(\epsilon_s-\epsilon_{s'})
\langle u_{s'}|\partial_k u_s\rangle .
$$

Using

$$
\xi_{ks's}=i\langle u_{ks'}|\partial_k u_{ks}\rangle ,
$$

one obtains the off-diagonal part of Eq. (28). Combining the diagonal and
off-diagonal cases gives the same object as

$$
U^\dagger(k)\,\frac{1}{\hbar}\frac{\partial H(k)}{\partial k}\,U(k).
$$

Thus, with exact derivatives and a consistent smooth gauge, the two formulas
are theoretically equivalent.

## Why The Numerical Eq. (28) Diagnostic Can Fail

The direct transform uses the analytic derivative of the Hamiltonian and only
requires eigenvectors at the same \(k\).

The Eq. (28) diagnostic needs finite differences of eigenvectors across
neighboring \(k\) points. This makes it sensitive to:

- arbitrary eigenvector phases from numerical diagonalization;
- band sorting changes near crossings or near degeneracies;
- insufficient \(k\)-grid density;
- imperfect periodic matching at the Brillouin-zone boundary;
- using single-band phase smoothing where a degenerate subspace treatment is needed.

For that reason, the Eq. (28) implementation is kept in
`Velocity_Matrix_Diagnostics.py` only. It is useful for studying Berry
connection quality, but it is not used by the production velocity-gauge RDM
evolution.

Run the diagnostic independently with:

```bash
python Velocity_Matrix_Diagnostics.py
```
