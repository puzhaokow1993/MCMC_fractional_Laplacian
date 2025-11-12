> [!NOTE] 
> For proper equation rendering, please view this documentation in day mode instead of night mode.

This repository provides the MATLAB implementation of the MCMC algorithm presented in the following paper: 

[Pu-Zhao Kow](https://puzhaokow1993.github.io/homepage/), [Janne Nurminen](https://users.jyu.fi/~jasanurm/) and [Jesse Railo](https://sites.google.com/view/jesserailo), *Bayesian inference for the fractional Calderón problem with a single measurement*, manuscript 

We will not going to explain all notations here, please refer to our paper for more details. 

Let ![Omega](https://latex.codecogs.com/png.image?\dpi{110}\Omega) be a bounded smooth domain in ![Rd](https://latex.codecogs.com/png.image?\dpi{110}\mathbb{R}^{d}), and let ![Laplacian](https://latex.codecogs.com/png.image?\dpi{110}(-\Delta)^{s}) denote the fractional Laplacian of order ![s\in(0,1)](https://latex.codecogs.com/png.image?\dpi{110}s\in(0,1)). 
We fix any ![phi](https://latex.codecogs.com/png.image?\dpi{110}0\not\equiv0\in%20C_{c}^{\infty}(\Omega_{e})) with ![Omega_e](https://latex.codecogs.com/png.image?\dpi{110}\Omega_e:=\mathbb{R}^{d}\setminus\overline{\Omega}), and consider the following Dirichlet problem: 
<div align="center">
  
![((-\Delta)^s+f)u=0](https://latex.codecogs.com/png.image?\dpi{110}((-\Delta)^s+f)u=0) in ![Omega](https://latex.codecogs.com/png.image?\dpi{110}\Omega) with ![u|_{\Omega_e}=\phi](https://latex.codecogs.com/png.image?\dpi{110}u|_{\Omega_e}=\phi), 
</div>

It is well-known that the potential ![f](https://latex.codecogs.com/png.image?\dpi{110}f) from the Dirichlet-to-Neumann map (DN map) 
<div align="center">
  
![\Lambda_{f}\phi \equiv (-\Delta)^{s}u_{f}|_{\mathcal{D}}](https://latex.codecogs.com/png.image?\dpi{110}\Lambda_{f}\phi\equiv(-\Delta)^{s}u_{f}|_{\mathcal{D}}), 
</div>

where ![\mathcal{D}\subset\Omega_e](https://latex.codecogs.com/png.image?\dpi{110}\mathcal{D}\subset\Omega_e) is a fixed open set. We adopt a Bayesian approach to this problem for the case when ![n=1](https://latex.codecogs.com/png.image?\dpi{110}n=1), providing not only rigorous theoretical justifications but also supporting numerical simulations. 

# Algorithm # 

**Require:** an initial guess ![f^{(0)}](https://latex.codecogs.com/png.image?\dpi{110}f^{(0)}) and a resolution parameter ![J_{(0)}](https://latex.codecogs.com/png.image?\dpi{110}J_{(0)})

1. Set ![{\rm%20status}(0)={\rm%20accept}](https://latex.codecogs.com/png.image?\dpi{110}{\rm%20status}(0)={\rm%20accept})
2. Creating an empty sequence ![({\rm%20status}(\tau))_{\tau=1}^{\infty}](https://latex.codecogs.com/png.image?\dpi{110}({\rm%20status}(\tau))_{\tau=1}^{\infty})
3. **for ![\tau=0,1,2,\cdots](https://latex.codecogs.com/png.image?\dpi{110}\tau=0,1,2,\cdots) do**
4. $~~~~$ Generate ![f(x)=\sum_{r=-2^{J_0}}^{2^{J_0-1}}f_{r}\chi_{(0,1)}(2(2^{J_0}x-r))](https://latex.codecogs.com/png.image?\dpi{110}f(x)=\sum_{r=-2^{J_0}}^{2^{J_0-1}}f_{r}\chi_{(0,1)}(2(2^{J_0}x-r))) with randomly chosen ![f_{r}\sim\mathcal{N}(0,1)](https://latex.codecogs.com/png.image?\dpi{110}f_{r}\sim\mathcal{N}(0,1))
5. $~~~~$ Propose ![f^{(\tau+1)}=1+\sqrt{(1-\beta^2)}(f^{(\tau)}-1)+\beta%20f](https://latex.codecogs.com/png.image?\dpi{110}f^{(\tau+1)}=1+\sqrt{(1-\beta^2)}(f^{(\tau)}-1)+\beta%20f)
6. $~~~~$ **if ![{\rm%20status}(\tau)={\rm%20accept}](https://latex.codecogs.com/png.image?\dpi{110}{\rm%20status}(\tau)={\rm%20accept}) then**
7. $~~~~~~~~$ ![\ell_{\rm%20current}=\tilde\ell^{(N)}(F^{(\tau)})](https://latex.codecogs.com/png.image?\dpi{110}\ell_{\rm%20current}=\tilde\ell^{(N)}(F^{(\tau)})), where ![\tilde\ell^{(N)}(f)=-\frac{1}{2\sigma^2}\sum_{i=1}^{N}(Y_{i}-G(f)(X_{i}))^{2}](https://latex.codecogs.com/png.image?\dpi{110}\tilde\ell^{(N)}(f)=-\frac{1}{2\sigma^2}\sum_{i=1}^{N}(Y_{i}-G(f)(X_{i}))^{2}) is the log-likelihood function 
8. $~~~~$ **end if**
9. $~~~~$ **if ![\tilde\ell^{(N)}(f^{(\tau+1)})>\ell_{\rm%20current}](https://latex.codecogs.com/png.image?\dpi{110}\tilde\ell^{(N)}(f^{(\tau+1)})>\ell_{\rm%20current}) then**
10. $~~~~~~~~$ Set ![{\rm%20status}(\tau+1)={\rm%20accept}](https://latex.codecogs.com/png.image?\dpi{110}{\rm%20status}(\tau+1)={\rm%20accept})
11. $~~~~$ **else**
12. $~~~~~~~~$ Set ![f^{(\tau+1)}=f^{(\tau)}](https://latex.codecogs.com/png.image?\dpi{110}f^{(\tau+1)}=f^{(\tau)}) and ![{\rm%20status}(\tau+1)={\rm%20reject}](https://latex.codecogs.com/png.image?\dpi{110}{\rm%20status}(\tau+1)={\rm%20reject})
13. $~~~~$ **end if**
14. **end for**

**Return:** ![(f^{(\tau)})_{\tau=1}^{\infty}](https://latex.codecogs.com/png.image?\dpi{110}(f^{(\tau)})_{\tau=1}^{\infty}) by removing all entries ![F^{(\tau)}](https://latex.codecogs.com/png.image?\dpi{110}F^{(\tau)}) each corresponds to ![{\rm%20status}(\tau)={\rm%20reject}](https://latex.codecogs.com/png.image?\dpi{110}{\rm%20status}(\tau)={\rm%20reject})
