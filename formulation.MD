The electric ($\mathbf{E}$) and magnetic ($\mathbf{H}$) field components of the electromagnetic wave propagating inside the dielectric ring can be represented as sums of three orthogonal vectors as follows:

$$
\begin{align}
\mathbf{E}(\rho, \phi, z) &  = \left\{ \hat{\rho} E_\rho(\rho, z) +\hat{\phi} E_\phi(\rho, z) + \hat{z}E_z(\rho, z) \right\} e^{-j \nu \phi}, 
\label{Efield} \\
\mathbf{H}(\rho, \phi, z) & = \left\{ \hat{\rho} H_\rho(\rho, z) +\hat{\phi} H_\phi(\rho, z) + \hat{z}H_z(\rho, z) \right\} e^{-j \nu \phi},    
\label{Hfield}
\end{align}
$$

assuming $e^{-j \nu \phi}$ dependence, where $\nu$ is the angular propagation constant, $\rho$ is the radial distance from the origin to the point projected onto the $xy$ plane, $\phi$ is the azimuthal angle, and $z$ is the height or vertical distance from the $xy$ plane.

## Modal Analysis

To carry out modal analysis, we start with Maxwell's equations to derive the Helmholtz wave equations. The $\hat{\rho}$ direction:

$$
\left( \mathcal{L} - \frac{1}{R_c^2} + \frac{\rho^2}{R_c^2}\frac{\partial}{\partial \rho}  \frac{1}{\varepsilon} \frac{\partial \varepsilon}{\partial \rho}\right) E_\rho  
    + \frac{2j \beta}{R_c}E_\phi + \frac{\rho^2}{R_c^2} \frac{\partial}{\partial \rho}\frac{1}{\varepsilon} \frac{\partial \varepsilon}{\partial z} E_z  = \beta^2 E_\rho
$$

where

$$
\mathcal{L} = \frac{\rho^2}{R_c^2} \left( \frac{\partial^2}{\partial \rho^2}  + \frac{1}{\rho}\frac{\partial}{\partial \rho} + \frac{\partial^2}{\partial z^2} +k_0^2 \varepsilon_r \right)
$$

## Matrix Formulation

The equations can be cast into a matrix equation:

$$
\overline{\overline{M}}\mathbf{A} + \beta \overline{\overline{L}} \mathbf{A} - \beta^2 \overline{\overline{I}} \mathbf{A} = 0
$$

where $\mathbf{A} = [E_{\rho} \ E_{z} \  H_{\rho} \ H_{z}]^T$.
