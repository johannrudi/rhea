# Rhea

[![Documentation Status](https://readthedocs.org/projects/rhea/badge/?version=latest&style=flat)](https://rhea.readthedocs.io/en/latest/)

<!-- start include-docs-index -->

Extreme-scale Simulation and Inference for Earth's Mantle Convection.

The adaptive nonlinear Stokes solver for mantle convection and plate tectonics.
Forward models of Earth's mantle are solved at a global scale while resolving
plate-boundaries. Inverse problems of finding rheological parameters of the
mantle are computed with adjoint-based gradients and Hessians.

The discretization of Earth's mantle is carried out by finite elements on
aggressively adaptively refined hexahedral meshes. We use a quadratic finite
element approximation for the velocity and linear elements for the pressure. To
distribute the discrete problem onto distributed-memory parallel computing
clusters, parallel forest-of-octrees algorithms are used for efficient, scalable
mesh refinement/coarsening, mesh balancing and repartitioning. The large
implicit nonlinear systems to be solved are poorly conditioned, and specifically
tailored iterative numerical methods including advanced preconditioning
techniques are required. Our hybrid spectral-geometric-algebraic multigrid (HMG)
method constitutes a core preconditioning component of the solver. HMG is
essential for preconditioning efficacy (reducing iteration counts) and
algorithmic and parallel scalability to extreme scales of 106 processor cores.
In addition, we used a preconditioner for the Schur complement that is robust in
the presence of the highly heterogeneous viscosities present in the models.
Combining this Schur complement preconditioner with HMG results in a scalable
and highly robust implicit linear solver. The implicit linear solver was
embedded into an inexact Newton-Krylov method to deal with severely nonlinear
rheologies. We developed a nonlinear preconditioner to avoid prohibitively
stagnating convergence of Newton's method due to highly nonlinear physics at
plate fault zones.

## Contributors

Johann Rudi, Max Heldman, Leonid Pereiaslov, Jiaqi Fang, Jiashun Hu

<!-- end include-docs-index -->
