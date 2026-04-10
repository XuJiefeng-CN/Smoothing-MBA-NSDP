# SmoothingMBA-NSDP-ℓ₁
This is a MATLAB package that implements the smoothing moving ball approximation (*s*MBA) method in [[1]] for solving the following convex nonlinear semidefinite programs with $\ell_1$ regularization (NSDP-ℓ₁):

$$\eqalign{
	\min\limits_{x\in ℝ^{n}} & \sum_{i=1}^{n}\left( \frac{1}{4} d_i x_i^4 + \frac{1}{3} c_i |x_i|^3 \right)+ x^\top Q x + b^{\top}x + \rho ‖x‖_1 \\
	\text{s.t.} & A_0 + x_1 A_1 + \cdots + x_n A_0  \succeq 0,
}$$

where $\rho$ is a nonnegative number, $b$ is a real-valued vector of dimension $n$, $c$ and $d$ are real-valued nonnegative vectors of dimension $n$, $Q$ is a symmetric positive semidefinite matrix of dimension $n$ by $n$, and $A_i$ are symmetric positive semidefinite matrices of dimension $m$ by $m$ for $i=0,...,n$.
The performance of *s*MBA is compared with CVX using the SDPT3 solver.
For further details, please refer to our paper in [[1]].

# Matlab source codes
- **demo_NSDP.m**\
A demo of the numerical experiments in [[1]].

- **NSDP.m**\
A function that generates the data of a convex QCQP.

- **sMBA_NSDP.m**\
The implementation of sMBA for solving QCQP-ℓ₁.

- **CaseSg.m** and **SubP_alpha.m**\
The implementations of the root-finding scheme described in [[2], Appendix A], which we use to solve the subproblem of *s*MBA (i.e., [[1], (3.1)]). The codes are hosted by Pong T. K. and are available at [https://www.polyu.edu.hk/ama/profile/pong/MBA_l1vl2/.](https://www.polyu.edu.hk/ama/profile/pong/MBA_l1vl2/)

# References
[1]: https://arxiv.org/pdf/2505.12314 "J. Xu, T. K. Pong and N. S. Sze. A smoothing moving balls approximation method for a class of conic-constrained difference-of-convex optimization problems. Preprint (2025)."
\[1\] [J. Xu, T. K. Pong and N. S. Sze. A smoothing moving balls approximation method for a class of conic-constrained difference-of-convex optimization problems. Preprint (2025).](https://arxiv.org/pdf/2505.12314)

[2]: https://epubs.siam.org/doi/abs/10.1137/20M1314057 "P. Yu, T. K. Pong and Z. Lu. Convergence rate analysis of a sequential convex programming 
method with line search for a class of constrained difference-of-convex optimization problems. *SIAM J. Optim.* 31, pp. 2024--2054 (2021)."
\[2\] [P. Yu, T. K. Pong and Z. Lu. Convergence rate analysis of a sequential convex programming method with line search for a class of constrained difference-of-convex optimization problems. *SIAM J. Optim.* 31, pp. 2024--2054 (2021).](https://epubs.siam.org/doi/abs/10.1137/20M1314057)



