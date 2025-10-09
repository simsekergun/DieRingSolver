The electric ($\mathbf{E}$) and magnetic ($\mathbf{H}$) field components of the electromagnetic wave propagating inside the dielectric ring can be represented as sums of three orthogonal vectors as follows
\begin{align}
\mathbf{E}(\rho, \phi, z) &  = \left\{ \hat{\rho} E_\rho(\rho, z) +\hat{\phi} E_\phi(\rho, z) + \hat{z}E_z(\rho, z) \right\} e^{-j \nu \phi}, 
\label{Efield} \\
\mathbf{H}(\rho, \phi, z) & = \left\{ \hat{\rho} H_\rho(\rho, z) +\hat{\phi} H_\phi(\rho, z) + \hat{z}H_z(\rho, z) \right\} e^{-j \nu \phi},    
\label{Hfield}
\end{align}    
assuming $e^{-j \nu \phi}$ dependence, $\nu$ is the angular propagation constant, $\rho$ is the radial distance from the origin to the point projected onto the $xy$ plane, $\phi$ is the azimuthal angle, and $z$ is the height or vertical distance from the $xy$ plane. Our goal is to determine the $\nu$ and to compute the propagating modes ($\mathbf{E}$ and $\mathbf{H}$) along the $\phi$ direction inside the dielectric ring for a given frequency ($f$), corresponding to the wavelength $\lambda = c_0/f$, where $c_0$ is the velocity of electromagnetic waves in vacuum.


To carry out modal analysis, we start with Maxwell's equations to derive the Helmholtz wave equations for the electric and magnetic fields in a source-free, lossless, and non-magnetic medium. After expanding each term of the electric field equation, we obtain the following three sets of scalar equations by enforcing each side of the equation to be the same in the $\hat{\rho}$, $\hat{\phi}$, and $\hat{z}$ directions: \\
\noindent $\hat{\rho}$ direction:
\begin{equation}
\left( \mathcal{L} - \frac{1}{R_c^2} + \frac{\rho^2}{R_c^2}\frac{\partial}{\partial \rho}  \frac{1}{\varepsilon} \frac{\partial \varepsilon}{\partial \rho}\right) E_\rho  
    + \frac{2j \beta}{R_c}E_\phi + \frac{\rho^2}{R_c^2} \frac{\partial}{\partial \rho}\frac{1}{\varepsilon} \frac{\partial \varepsilon}{\partial z} E_z  = \beta^2 E_\rho,
\label{EformRho}
\end{equation}%
%%%%%%%%
$\hat{\phi}$ direction:
\begin{equation}
\left( -\frac{2j \beta}{R_c} -\frac{j \beta \rho}{R_c \varepsilon} \frac{\partial \varepsilon}{\partial \rho} \right) E_\rho  
    + \left( \mathcal{L} - \frac{1}{R_c^2}  \right) E_\phi
    -\frac{j \beta\rho}{R_c \varepsilon} \frac{\partial \varepsilon}{\partial z} E_z = \beta^2 E_\phi,
\end{equation}
%%%%%%%%
%
%%%%%%%%
$\hat{z}$ direction:
\begin{equation}
\frac{\rho^2}{R_c^2}\frac{\partial}{\partial z} \frac{1}{\varepsilon} \frac{\partial \varepsilon}{\partial \rho} E_\rho + \mathcal{L} E_z  +  \frac{\rho^2}{R_c^2}\frac{\partial}{\partial z}  \frac{1}{\varepsilon} \frac{\partial \varepsilon}{\partial z} E_z = \beta^2 E_z, 
\label{EformZ}
\end{equation}
where 
%
\begin{equation}
\mathcal{L} = \frac{\rho^2}{R_c^2} \left( \frac{\partial^2}{\partial \rho^2}  + \frac{1}{\rho}\frac{\partial}{\partial \rho} + \frac{\partial^2}{\partial z^2} +k_0^2 \varepsilon_r \right). 
\label{Loperator}
\end{equation}
Following a similar procedure on the Helmholtz wave equation for the magnetic field, we obtain three other equations. However, since we will only need the one along the $\hat{z}$ direction, we provide it below for the sake of completeness. 
\begin{equation}  
\frac{\rho^2}{R_c^2} \frac{1}{\varepsilon} \frac{\partial \varepsilon}{\partial \rho} \frac{\partial H_\rho}{\partial z}
+ \left( \mathcal{L} - \frac{\rho^2}{R_c^2}  \frac{1}{\varepsilon} \frac{\partial \varepsilon}{\partial \rho} \frac{\partial }{\partial \rho} \right) H_z = \beta^2 H_z.
\label{HformZ}
\end{equation}

It is known that mixed electric field–magnetic field formulations allow for direct enforcement of boundary conditions on both electric and magnetic fields at interfaces between materials \cite{bane}. In single-field formulations, deriving the secondary field (electric from magnetic or vice versa) often complicates boundary condition enforcement, introducing errors at material interfaces\cite{Simsek:25}. 
From Maxwell's equations, one can derive the following expressions:\\
\begin{minipage}{0.42\linewidth}
\begin{equation}
E_\phi = \frac{1}{j \omega \varepsilon} \left(\frac{\partial H_\rho}{\partial z} - \frac{\partial H_z}{\partial \rho} \right), 
\label{ephi} 
\end{equation}
\end{minipage}
\hspace{0.15\linewidth}
\begin{minipage}{0.42\linewidth}
\begin{equation}
H_\phi = \frac{1}{j\omega\mu_0}  \left( \frac{\partial E_z}{\partial \rho} - \frac{\partial E_\rho}{\partial z} \right).
\label{hphi} 
\end{equation}
\end{minipage}

\vspace{2mm}
\noindent By using Eq. \eqref{ephi}, we rewrite Eq. \eqref{EformRho} in the following form 
\begin{equation}  
 \left(\mathcal{L} - \frac{1}{R_c^2} + \frac{\rho^2}{R_c^2}\frac{\partial}{\partial \rho}  \frac{1}{\varepsilon} \frac{\partial \varepsilon}{\partial \rho} \right) E_\rho 
+ \left( \frac{\rho^2}{R_c^2} \frac{\partial}{\partial \rho}\frac{1}{\varepsilon} \frac{\partial \varepsilon}{\partial z} \right) E_z 
 + \frac{2\beta}{\omega \varepsilon R_c} \left(\frac{\partial H_\rho}{\partial z} - \frac{\partial H_z}{\partial \rho} \right)   = \beta^2 E_\rho. 
\label{mixedErho}
\end{equation}  
Similarly,  $\hat{\rho}$ component of the $H$-formulation becomes
\begin{equation}  
\left(\mathcal{L} -\frac{\rho^2}{R_c^2}\frac{1}{\varepsilon} \frac{\partial \varepsilon}{\partial z}\frac{\partial }{\partial z} -\frac{1}{R_c^2}\right) H_\rho
+ \frac{\rho^2}{R_c^2}\frac{1}{\varepsilon} \frac{\partial \varepsilon}{\partial z}\frac{\partial }{\partial z} H_z + \frac{2\beta}{\omega\mu_0 R_c}  \left( \frac{\partial E_z}{\partial \rho} - \frac{\partial E_\rho}{\partial z} \right) = \beta^2 H_\rho.
\label{mixedHrho}
\end{equation}  
Eqs. \eqref{mixedErho}, \eqref{mixedHrho}, \eqref{EformZ}, and \eqref{HformZ} can be cast into a matrix equation such as
\begin{equation}
\overline{\overline{{M}}}\mathbf{A} + \beta \overline{\overline{{L}}} \mathbf{A} - \beta^2 \overline{\overline{{I}}} \mathbf{A} = 0,
\end{equation}
where $\mathbf{A} = [E_{\rho} \ E_{z} \  H_{\rho} \ H_{z}]^T$, $\overline{\overline{{M}}}$ is the matrix independent of $\beta$, $\overline{\overline{{L}}}$ is the linear term in $\beta$, and the last term is the quadratic term in $\beta$. 

\begin{minipage}{0.4\linewidth}
\begin{equation}
\overline{\overline{{M}}}  = 
\begin{bmatrix}
M_1 & M_2 & 0 & 0 \\
M_3 & M_4 & 0 & 0  \\
0 & 0  & M_{5} & M_{6} \\
0 & 0  & M_{7} & M_{7} \\
\end{bmatrix}
\end{equation}    
\end{minipage}
\hspace{0.15\linewidth}
\begin{minipage}{0.4\linewidth}
\begin{equation}
\overline{\overline{{L}}}  = 
\begin{bmatrix}
0 & 0 & L_1 & L_2\\
0 &  0 & 0 & 0 \\
L_3 & L_{4} & 0 & 0 \\
0 &  0 & 0 & 0 \\
\end{bmatrix}
\end{equation}
\end{minipage}

\begin{minipage}{0.4\linewidth}
\begin{align}
    M_1 = & \mathcal{L} - \frac{1}{R_c^2} + \frac{\rho^2}{R_c^2}\frac{\partial}{\partial \rho}  \frac{1}{\varepsilon} \frac{\partial \varepsilon}{\partial \rho} \\
    M_2 = & \frac{\rho^2}{R_c^2} \frac{\partial}{\partial \rho}\frac{1}{\varepsilon} \frac{\partial \varepsilon}{\partial z}  \\
    M_3 = & \frac{\rho^2}{R_c^2}\frac{\partial}{\partial z} \frac{1}{\varepsilon} \frac{\partial \varepsilon}{\partial \rho}   \\
    M_4 = & \mathcal{L} +  \frac{\rho^2}{R_c^2}\frac{\partial}{\partial z}  \frac{1}{\varepsilon} \frac{\partial \varepsilon}{\partial z}
\end{align}    
\end{minipage} \hspace{0.12\linewidth}
\begin{minipage}{0.4\linewidth}
\begin{align}    
    M_{5} = & \mathcal{L} -\frac{\rho^2}{R_c^2}\frac{1}{\varepsilon} \frac{\partial \varepsilon}{\partial z}\frac{\partial }{\partial z} -\frac{1}{R_c^2}\\
    M_{6} = & \frac{\rho^2}{R_c^2}\frac{1}{\varepsilon} \frac{\partial \varepsilon}{\partial z}\frac{\partial }{\partial z} \\
    M_{7} = & \frac{\rho^2}{R_c^2} \frac{1}{\varepsilon} \frac{\partial \varepsilon}{\partial \rho} \frac{\partial }{\partial z} \\
    M_{8} = &  \mathcal{L} - \frac{\rho^2}{R_c^2}  \frac{1}{\varepsilon} \frac{\partial \varepsilon}{\partial \rho} \frac{\partial }{\partial \rho}
\end{align}
\end{minipage}

\begin{minipage}{0.4\linewidth}
\begin{align}
    L_1 = & \frac{2\beta}{\omega \varepsilon R_c} \frac{\partial H_\rho}{\partial z} \\
    L_2 = & -\frac{2\beta}{\omega \varepsilon R_c}  \frac{\partial H_z}{\partial \rho}
\end{align}    
\end{minipage} \hspace{0.12\linewidth}
\begin{minipage}{0.4\linewidth}
\begin{align}        
    L_3 = &  - \frac{2\beta}{\omega\mu_0 R_c} \frac{\partial}{\partial z}  \\
    L_4 = &  \frac{2\beta}{\omega\mu_0 R_c} \frac{\partial}{\partial \rho} 
\end{align}
\end{minipage}

\vspace{2mm}
This can be rewritten as a linear generalized eigenvalue problem by introducing an auxiliary variable $\mathbf{F} = \beta \mathbf{A}$, leading to $\overline{\overline{{M}}} \mathbf{A} + \overline{\overline{{L}}} \mathbf{F} - \mathbf{F} \beta = 0$, which gives the augmented system:
\begin{equation}
    \begin{pmatrix}
        \overline{\overline{{M}}} & \overline{\overline{{L}}} \\
        0 & \overline{\overline{{I}}}
    \end{pmatrix}
    \begin{pmatrix}
        \mathbf{A} \\
        \mathbf{F}
    \end{pmatrix}
    = \beta
    \begin{pmatrix}
        0 & \overline{\overline{{I}}} \\
        \overline{\overline{{I}}} & 0
    \end{pmatrix}
    \begin{pmatrix}
        \mathbf{A} \\
        \mathbf{F}
    \end{pmatrix}.
\end{equation}
We solve this linear generalized eigenvalue problem using MATLAB's \textsf{eigs} function, which is optimized for large sparse matrices. After solving, we retain only the physically meaningful modes (those with $n_{\rm eff}^{\rm ring}<n_{\rm{ring}}$, where $n_{\rm eff}^{\rm ring} = \beta/k_0$ is the effective index of the ring and $n_{\rm{ring}}$ is the refractive index of the material which ring is made from.)
