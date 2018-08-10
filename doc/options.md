# Options

Here we describe the available options that can be changed in BOUT.inp.
Storm-specific options are in the [storm] section. We mention a few useful
BOUT++ options as well.

**Contents**
* [Storm](#storm)
* [Storm2D](#storm2d)

Storm
-----

| **name** | **description** | **default**
|---|---|---|
| timestep | time between outputs, in cyclotron times | 1
| nout | number of output timesteps | 1
| mz | number of z-points | 64
| zmin, zmax | length of z-grid is `2\*pi\*(zmax-zmin)` | 0, 1
| **[Storm]** |  |
| B\_0 | normalising magnetic field (T) | 0.5
| T\_e0 | normalizing electron temperature (eV) | 40
| T\_i0 | ion temperature (eV) used to calculate dissipation parameters | 40
| m\_i | ion mass (amu) | 2
| q | 'safety factor' used to enhance dissipation parameters | 7
| R\_c | radius of curvature of the magnetic field (m), used to calculate `g0` | 1.5
| n\_0 | normalizing density (m\^-3) | 8e18
| loglambda | Coulomb logarithm, calculated from other parameters if negative value set | -1
| phi\_wall | electric potential of the wall (V/`T_e0`) | 0
| mu\_n0 | density diffusion coefficient, calculated from normalization parameters if negative | -1
| mu\_vort0 | vorticity diffusion coefficient, calculated from normalization parameters if negative | -1
| mu | ion-electron mass ration m\_i/m\_e, calculated if negative | -1
| nu\_parallel0 | normalized resistivity, calculated from normalization parameters if negative | -1
| kappa0 | constant prefactor of parallel heat conduction coefficient, calculated from normalization parameters if negative; not affected by uniform\_diss\_paras | -1
| kappa0\_perp | heat diffusion coefficient, calculated from normalization parameters if negative | -1
| g0 | curvature drive term, calculated from `R_c` if negative | -1
| isothermal | switch for evolution of electron temperature | false
| uniform\_diss\_paras | if true, dissipation coefficients are just calculated from normalization parameters, variation with `n` and `T` is neglected | false
| verbose | print current simulation time in stdout, updated each internal timestep | false
| bracket | choose which BOUT++ option to use for the Poisson bracket operators: 0 - default finite differencing; 2 - Arakawa scheme; other options should not be used | 2
| add\_blob | if true, add a blob perturbation using parameters in `[blob]` | true
| symmetry\_plane | use reflection-symmetry boundary conditions at the lower y-boundary | true
| run\_1d | if true, run in 1d (y-direction only) mode for setting up a steady background; for use by `create_bg_1d` script | false
| run\_1d\_T\_slowdown | reduce time derivative of electron temperature `T` by this factor to allow longer time steps and speed up finding steady profiles | 20
| equilibrium\_source | where to get the initial 'background' profiles from: '1d\_profiles' - binary files created by `create_bg_1d`; 'input\_file' - standard BOUT++ variable initialization from the input file; 'profiles\_file' - read arrays from a binary 'grid' file (netCDF or HDF5) given by `equilibrium_data_file` | "1d\_profiles"
| equilibrium\_file\_path | relative or absolute path to directory where files for `equilibrium_source=1d_profiles` or `equilibrium_source=profiles_file` options are located - if not given, the subdirectory given by `equilibrium_directory` of the directory specified for input/output | ""
| equilibrium\_directory | name of the subdirectory to read background profiles from | "equilibrium"
| equilibrium\_data\_file | see `equilibrium_source` | ""
| **[blob]** | |
| delta\_perp | width of density perturbation in units of rho\_s | 10
| elongation | ellipticity of density perturbation (ratio of long to short axis) | 1
| angle\_blob | tilt of axis of density perturbation relative to x-direction (radians) | 0
| A | amplitude of density perturbation | 2
| A\_relative\_to\_bg | if true `A` is a multiple of the background density and `A_T` a multiple of the background temperature, if false `A` is a density normalized to `n_0` and `A_T` a temperature normalized to `T_e0` | true
| L\_b | parallel extent of density perturbation as a fraction of the domain size | 0.5
| xoffset | initial x-position of the centre of the density perturbation as a fraction of the domain size | 0.25
| zoffset | initial z-position of the centre of the density perturbation as a fraction of the domain size | 0.5
| boltzmann | vorticity is initialized so that phi=log(n) | false
| delta\_front | scale length of the parallel gradient at the ends of the density perturbation, as a fraction of the domain size | 0.3
| conserve\_momentum | `U` and `V` are reduced so that the momentum density after adding the filament is the same as the momentum density of the background | false
| delta\_perp\_T | width of temperature perturbation in units of rho\_s | 10
| elongation\_T | ellipticity of temperature perturbation | 1
| angle\_blob\_T | tilt of axis of temperature perturbation relative to x-direction (radians) | 0
| A\_T | amplitude of temperature perturbation | 2
| L\_b\_T | parallel extent of temperature perturbation as a fraction of the domain size | 0.5
| xoffset\_T | initial x-position of the centre of the temperature perturbation as a fraction of the domain size | 0.25
| zoffset\_T | initial z-position of the centre of the temperature perturbation as a fraction of the domain size | 0.5
| delta\_front\_T | scale length of the parallel gradient at the ends of the temperature perturbation, as a fraction of the domain size | 0.3
| **[mesh]** | |
| staggergrids | switch for staggered grids, must be set to true for STORM | false
| Lx | length of x-grid in units of rho\_s (just used to calculate dx) | N/A
| Ly | length of y-grid in units of rho\_s (just used to calculate dy) | N/A
| ixseps1, ixseps2 | x-index of separatrices, should both be negative or zero for SOL simulations | nx
| **[S]** | |
| function | expression in terms of x, y, z giving particle source density multiplied by Ly (multiplying by Ly means if this number is kept fixed, the total particle source in the domain is constant as Ly changes, which keeps a similar density profile since the sink at the sheath also does not change with Ly) [0\<x\<1, 0\<y\<2*pi, 0\<z\<2*pi] |
| **[S\_E]** | |
| function | expression in terms of x, y, z giving energy source density multiplied by Ly (multiplying by Ly means if this number is kept fixed, the total particle source in the domain is constant as Ly changes, which keeps a similar density profile since the sink at the sheath also does not change with Ly) [0\<x\<1, 0\<y\<2*pi, 0\<z\<2*pi] |
| **[n]** | |
| function | expression giving initial profile of `n` if `equilibrium_source=input_file` |
| bndry\_yup | boundary condition at the upper target, should be `free_o3` |
| bndry\_ydown | if `symmetry_plane=true` boundary condition at the symmetry plane, should be `neumann`; otherwise boundary coundition at lower target, should be `free_o3` |
| bndry\_xin | boundary condition at the inner radial boundary, should be `neumann` |
| bndry\_xout | boundary condition at the outer radial boundary, should be `neumann` |
| **[T]** | |
| function | expression giving initial profile of `T` if `equilibrium_source=input_file` |
| bndry\_yup | boundary condition at the upper target, should be `free_o3` |
| bndry\_ydown | if `symmetry_plane=true` boundary condition at the symmetry plane, should be `neumann`; otherwise boundary coundition at lower target, should be `free_o3` |
| bndry\_xin | boundary condition at the inner radial boundary, should be `neumann` |
| bndry\_xout | boundary condition at the outer radial boundary, should be `neumann` |
| **[vort]** | |
| function | expression giving initial profile of `vort` if `equilibrium_source=input_file` |
| bndry\_yup | boundary condition at the upper target, should be `free_o3` |
| bndry\_ydown | if `symmetry_plane=true` boundary condition at the symmetry plane, should be `neumann`; otherwise boundary coundition at lower target, should be `free_o3` |
| bndry\_xin | boundary condition at the inner radial boundary, should be `neumann` |
| bndry\_xout | boundary condition at the outer radial boundary, should be `neumann` |
| **[U]** | |
| function | expression giving initial profile of `U` if `equilibrium_source=input_file` |
| bndry\_yup | boundary condition at the upper target, should be `free_o3` |
| bndry\_ydown | if `symmetry_plane=true` boundary condition at the symmetry plane, should be `dirichlet`; otherwise boundary coundition at lower target, should be `none` |
| bndry\_xin | boundary condition at the inner radial boundary, should be `neumann` |
| bndry\_xout | boundary condition at the outer radial boundary, should be `neumann` |
| **[V]** | |
| function | expression giving initial profile of `V` if `equilibrium_source=input_file` |
| bndry\_yup | boundary condition at the upper target, should be `free_o3` |
| bndry\_ydown | if `symmetry_plane=true` boundary condition at the symmetry plane, should be `dirichlet`; otherwise boundary coundition at lower target, should be `none` |
| bndry\_xin | boundary condition at the inner radial boundary, should be `neumann` |
| bndry\_xout | boundary condition at the outer radial boundary, should be `neumann` |
| **[phi]** | |
| function | expression giving initial profile of `phi` if `equilibrium_source=input_file` |
| bndry\_yup | boundary condition at the upper target, should be `free_o3` |
| bndry\_ydown | if `symmetry_plane=true` boundary condition at the symmetry plane, should be `neumann`; otherwise boundary coundition at lower target, should be `free_o3` |
| bndry\_xin | boundary condition at the inner radial boundary, should be `none` |
| bndry\_xout | boundary condition at the outer radial boundary, should be `none` |
| **[phi_stag]** | |
| bndry\_yup | boundary condition at the upper target, should be `none` |
| bndry\_ydown | boundary condition at the symmetry plane or lower target, should be `none` |
| bndry\_xin | boundary condition at the inner radial boundary, should be `free_o3` |
| bndry\_xout | boundary condition at the outer radial boundary, should be `free_o3` |
| **[qpar]** | |
| bndry\_yup | boundary condition at the upper target, should be `none` |
| bndry\_ydown | if `symmetry_plane=true` boundary condition at the symmetry plane, should be `dirichlet`; otherwise boundary coundition at lower target, should be `none` |
| bndry\_xin | boundary condition at the inner radial boundary, should be `none` |
| bndry\_xout | boundary condition at the outer radial boundary, should be `none` |
| **[laplace]** | |
| inner\_boundary\_flags | needs to be set to `16` to use dirichlet boundary conditions set to the equilibrium value of phi at the radial boundaries
| outer\_boundary\_flags | needs to be set to `16` to use dirichlet boundary conditions set to the equilibrium value of phi at the radial boundaries
| **[solver]** | |
| type | `pvode` - adaptive timestep, adaptive order implicit scheme; `cvode` - newer version of pvode, if you have SUNDIALS library installed; `rk4` 4th order Runke-Kutta scheme (fixed timestep by default) [for more options see BOUT++ documentation] | pvode
| mxstep | maximum number of internal timesteps between output timesteps, usually good to set this to a very large number, like 100000000 | 500
| timestep | internal timestep for fixed-timestep solvers like `rk4` | output timestep
| rtol | relative tolerance for adaptive timestep | 1e-5
| atol | absolute tolerance for adaptive timestep | 1e-12
| **[mesh:ddx]** | |
| first | differencing scheme for first derivatives in the x-direction: C2 - 2nd order central difference; C4 - 4th order central difference | C2
| second | differencing scheme for second derivatives in the x-direction: C2 - 2nd order central difference; C4 - 4th order central difference | C2
| upwind | differencing scheme for upwind derivatives in the x-direction: C2 - 2nd order central difference; C4 - 4th order central difference; U1 - 1st order upwind; U2 - 2nd order upwind; W3 - WENO scheme | U1
| **[mesh:ddy]** | |
| first | differencing scheme for first derivatives in the y-direction: C2 - 2nd order central difference; C4 - 4th order central difference | C2
| firststag | differencing scheme for first derivatives output on a staggered grid in the y-direction: C2 - 2nd order central difference; C4 - 4th order central difference | same as first
| second | differencing scheme for second derivatives in the y-direction: C2 - 2nd order central difference; C4 - 4th order central difference | C2
| secondstag | differencing scheme for second derivatives output on a staggered grid in the y-direction: C2 - 2nd order central difference; C4 - 4th order central difference | same as second
| upwind | differencing scheme for upwind derivatives in the y-direction: C2 - 2nd order central difference; C4 - 4th order central difference; U1 - 1st order upwind; U2 - 2nd order upwind; W3 - WENO scheme | U1
| upwindstag | differencing scheme for upwind derivatives involving staggered fields in the y-direction: C2 - 2nd order central difference; C4 - 4th order central difference; U1 - 1st order upwind; U2 - 2nd order upwind; W3 - WENO scheme | same as upwind
| **[mesh:ddz]** | |
| first | differencing scheme for first derivatives in the z-direction: C2 - 2nd order central difference; C4 - 4th order central difference; FFT - fast-Fourier-transform scheme | C2
| firststag | not used in STORM, but may need to be set explicitly if using FFT derivatives since FFT is not available for staggered derivatives: C2 - 2nd order central difference; C4 - 4th order central difference | same as first
| second | differencing scheme for second derivatives in the z-direction: C2 - 2nd order central difference; C4 - 4th order central difference; FFT - fast-Fourier-transform scheme  | C2
| secondstag | not used in STORM, but may need to be set explicitly if using FFT derivatives since FFT is not available for staggered derivatives: C2 - 2nd order central difference; C4 - 4th order central difference | same as second
| upwind | differencing scheme for upwind derivatives in the z-direction: C2 - 2nd order central difference; C4 - 4th order central difference; U1 - 1st order upwind; U2 - 2nd order upwind; W3 - WENO scheme | U1
| **[output]** | |
| floats | if true output single precision floats instead of double precision | false

Storm2D
-------

| **name** | **description** | **default**
|---|---|---|
| timestep | time between outputs, in cyclotron times | 1
| nout | number of output timesteps | 1
| mz | number of z-points | 64
| zmin, zmax | length of z-grid is `2\*pi\*(zmax-zmin)` | 0, 1
| **[Storm]** |  |
| B\_0 | normalising magnetic field (T) | 0.5
| T\_e0 | normalizing electron temperature (eV) | 40
| T\_i0 | ion temperature (eV) used to calculate dissipation parameters | 40
| m\_i | ion mass (amu) | 2
| q | 'safety factor' used to enhance dissipation parameters | 7
| R\_c | radius of curvature of the magnetic field (m), used to calculate `g0` | 1.5
| n\_0 | normalizing density (m\^-3) | 8e18
| loglambda | Coulomb logarithm, calculated from other parameters if negative value set | -1
| mu\_n0 | density diffusion coefficient, calculated from normalization parameters if negative | -1
| mu\_vort0 | vorticity diffusion coefficient, calculated from normalization parameters if negative | -1
| kappa0 | constant prefactor of parallel heat conduction coefficient, calculated from normalization parameters if negative; used by `SOL_closure=vort_adv` | -1
| kappa0\_perp | heat diffusion coefficient, calculated from normalization parameters if negative | -1
| g0 | curvature drive term, calculated from `R_c` if negative | -1
| isothermal | switch for evolution of electron temperature | false
| SOL_closure | choice of SOL closure model: "sheath_diss", "vort_adv" or "ESEL_like" | "sheath_diss"
| sheath_linear | in an isothermal simulation with `SOL_closure="sheath_diss"`, calculate the sheath potential as (1-phi) instead of exp(-phi) | false
| n_bg | set a constant background density, a density source will be set to restore the density to this value | 1
| T_bg | set a constant background temperature, a heat source will be set to restore the temperature to this value | 1
| bracket | choose which BOUT++ option to use for the Poisson bracket operators: 0 - default finite differencing; 2 - Arakawa scheme; other options should not be used | 2
| **[mesh]** | |
| Lx | length of x-grid in units of rho\_s (just used to calculate dx) | N/A
| Ly | parallel connection length in units of rho\_s (used by `SOL_closure="sheath_diss"`) | N/A
| **[All]** | |
| bndry_all | set boundary conditions for all variables, should be `neumann` |
| **[n]** | |
| function | expression in terms of x, y, z giving initial profile of `n` [0\<x\<1, 0\<y\<2*pi, 0\<z\<2*pi] |
| **[T]** | |
| function | expression in terms of x, y, z giving initial profile of `T` [0\<x\<1, 0\<y\<2*pi, 0\<z\<2*pi] |
| **[vort]** | |
| function | expression in terms of x, y, z giving initial profile of `vort` [0\<x\<1, 0\<y\<2*pi, 0\<z\<2*pi] |
| **[sigma_n]** | |
| function | expression in terms of x, y, z giving profile of density sink if `SOL_closure="ESEL_like"` [0\<x\<1, 0\<y\<2*pi, 0\<z\<2*pi] |
| **[sigma_T]** | |
| function | expression in terms of x, y, z giving profile of temperature sink if `SOL_closure="ESEL_like"` [0\<x\<1, 0\<y\<2*pi, 0\<z\<2*pi] |
| **[sigma_vort]** | |
| function | expression in terms of x, y, z giving profile of vorticity sink if `SOL_closure="ESEL_like"` [0\<x\<1, 0\<y\<2*pi, 0\<z\<2*pi] |
| **[laplace]** | |
| inner\_boundary\_flags | needs to be set to `16` to use dirichlet boundary conditions with non-zero value if the background phi is floating potential instead of zero
| outer\_boundary\_flags | needs to be set to `16` to use dirichlet boundary conditions with non-zero value if the background phi is floating potential instead of zero
| **[solver]** | |
| type | `pvode` - adaptive timestep, adaptive order implicit scheme; `cvode` - newer version of pvode, if you have SUNDIALS library installed; `rk4` 4th order Runke-Kutta scheme (fixed timestep by default) [for more options see BOUT++ documentation] | pvode
| mxstep | maximum number of internal timesteps between output timesteps, usually good to set this to a very large number, like 100000000 | 500
| timestep | internal timestep for fixed-timestep solvers like `rk4` | output timestep
| rtol | relative tolerance for adaptive timestep | 1e-5
| atol | absolute tolerance for adaptive timestep | 1e-12
| **[mesh:ddx]** | |
| first | differencing scheme for first derivatives in the x-direction: C2 - 2nd order central difference; C4 - 4th order central difference | C2
| second | differencing scheme for second derivatives in the x-direction: C2 - 2nd order central difference; C4 - 4th order central difference | C2
| upwind | differencing scheme for upwind derivatives in the x-direction: C2 - 2nd order central difference; C4 - 4th order central difference; U1 - 1st order upwind; U2 - 2nd order upwind; W3 - WENO scheme | U1
| **[mesh:ddz]** | |
| first | differencing scheme for first derivatives in the z-direction: C2 - 2nd order central difference; C4 - 4th order central difference; FFT - fast-Fourier-transform scheme | C2
| firststag | not used in STORM, but may need to be set explicitly if using FFT derivatives since FFT is not available for staggered derivatives: C2 - 2nd order central difference; C4 - 4th order central difference | same as first
| second | differencing scheme for second derivatives in the z-direction: C2 - 2nd order central difference; C4 - 4th order central difference; FFT - fast-Fourier-transform scheme  | C2
| secondstag | not used in STORM, but may need to be set explicitly if using FFT derivatives since FFT is not available for staggered derivatives: C2 - 2nd order central difference; C4 - 4th order central difference | same as second
| upwind | differencing scheme for upwind derivatives in the z-direction: C2 - 2nd order central difference; C4 - 4th order central difference; U1 - 1st order upwind; U2 - 2nd order upwind; W3 - WENO scheme | U1
| **[output]** | |
| floats | if true output single precision floats instead of double precision | false


