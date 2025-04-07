Code to reconstruct the lm-modes of the metric perturbation in Lorenz gauge, in the frequency domain, for a circular orbit on Kerr spacetime, using the method detailed in the two papers: 
 - “Metric perturbations of Kerr spacetime in Lorenz gauge: Circular equatorial orbits”, Sam R. Dolan, Leanne Durkan, Chris Kavanagh, Barry Wardell, Class. Quantum Grav. 41, 155011 (2024) [arXiv:2306.16459].
 - "Gravitational Perturbations of Rotating Black Holes in Lorenz Gauge", Sam R. Dolan, Chris Kavanagh, Barry Wardell, Phys. Rev. Lett. 128 (2022) 15, 151101 [arXiv:2108.06344].
  
If you make use of this code for your project, please cite the papers above and add a line in the acknowledgements: “This paper makes use of the Kerr-Lorenz-Circ Mathematica code written by Sam Dolan.” (or similar). This code requires the Black Hole Perturbation Toolkit.

With thanks to Barry Wardell, Patrick Bourg and Conor Dyson for code-testing and feedback.

There are three main Mathematica notebooks here:
 - metric_reconstruction_calc_radiative.nb  does the calculation of the metric perturbation for m != 0 modes, and exports the data.
 - metric_reconstruction_calc_static.nb  does the calculation of the metric perturbation for the m = 0 modes, and exports the data.
 - metric_reconstruction_visualisation.nb  reads the data in and shows how to manipulate and plot it.
   
I recommend to run the notebooks on different kernels, since there are likely to be clashes between variable names.

There is also a .m Mathematica script here (metric_reconstruction_calc_radiative.m) which can be run from the command line using:

 wolframscript -file metric_reconstruction_calc_radiative.m -args config_string

where config_string is an appropriate string labelling the config file in config/ (see below). 
The script is auto-generated from the Initialisation Cells of the corresponding notebook, each time the latter is saved. It does not use section B (see below).

There are some additional notebooks that pre-generate input files for the static (zero-frequency) m=0 calculation. These are:
 - CalcProjectionCoeffs.nb: Calculates coefficients such as <l | cos^2 theta | l’>
 - CalcStaticRadialFunctions.nb: Calculates the various radial functions in closed form.
 - CalcStaticTensorMode.nb: Calculates the components of the homogeneous metric perturbation for m=0.
The output is stored in the folders static/ and completion/.

*Metric components and spin weight*
The 10 components of the metric perturbation (MP) that I use are h_{ab} projected onto the vector legs of an unnormalised version of the Kinnerley (or Carter) tetrad defined in Eq. (7)-(8) of arXiv:2306.16459. These are numbered 1 to 10 and always in the following order:

[ h_{l+l+}, h_{l-l-}, h_{m+ m+}, h_{m- m-}, \[Rho] h_{l+ m+}, \[Rho]* h_{l+ m-}, \[Rho]* h_{l- m+}, \[Rho] h_{l- m-}, \[CapitalSigma] \[CapitalDelta] h_{l+ l-}, h (the trace) ].
(see Eq. (128) in DDKW).

The spin-weight of each MP projected component follows from counting the number of m+ projections and subtracting the number of m- projections. The radial functions labelled by l,m are the MP components projected onto a spin-weighted spherical harmonic of the appropriate spin-weight.

*Configuration parameters*
The input parameters are read from files: config/config<iConfig>.txt. Here <iConfig> is a string. Note that you will need one config file for each value of m (the azimuthal mode number). 

Some of the parameters in the config files are irrelevant; they are used only in the 2+1D time-domain calculation. Below is a list of the parameters that are relevant in metric reconstruction:

 a: The spin parameter of the black hole.
 r0: The orbital radius (in Boyer-Lindquist coords).
 m: The azimuthal number of the mode.
 n: The grid spacing in the tortoise  coordinate r* , with Delta r*  =  M / n.
 angres: Angular resolution. The number of grid points in the theta direction is  n*angres.
 lmax: The upper bound on the number of l modes of radial/angular functions to use.
 lplot: After projection onto spherical harmonics, the number of l-modes reliable enough for further use (e.g. plotting).
 nterms: The number of terms to include in the expansion of spheroidal harmonics in the spherical basis.
 inford: The order of the series expansion of the UP solutions at infinity for the kappa(r) functions.
 rgrid: 0 for a linearly-spaced grid in r* (the tortoise coordinate), or 1 for a linearly-spaced grid in r.
 horord: The order of the series expansion of the IN solutions at the horizon for the kappa(r) functions.
 rinf: The value of r at which the UP solution for kappa is evaluated; this is the initial condition for integration inwards.
 xhor: The IN solutions for kappa start at  r = r_+ + xhor.
 rmax: The maximum value of r for the Interpolation Functions for kappa.
 accgoal: Accuracy Goal and Precision Goal for the numerical integrator.
 kapord: The kappa_l function is expanded in a basis of spheroidal harmonics in the range {l-kapord, l+kapord}. Should be at least as large as nterms.
 prec: The number of digits in Working Precision. Should be at least 32.
 rstmin: The minimum value of the tortoise coordinate r* for the grid. 
 rstmax: The maximum value of the tortoise coordinate r* in the grid.
 dformat: The format for exported data. Either Real32 or Real64. The latter is recommended.

*Naming convention*
There is a naming convention for <iConfig>: it should be in the form XX-YYY-ZZZ. 
XX are the first two digits of "a" after the decimal point (i.e. 00 for Schwarzschild, 60 for a=0.6M, 99 for a=0.99M).
YYY are the digits of r0 (with 060 implying r0 = 6.0).
ZZZ are the run numbers, with the first digit typically labelling the internal parameter choices (such as grid size), and the latter two digits indicating the value of 'm'. The numbers 0ZZ are reserved for testing, whereas the numbers zZZ with z >= 1 are for production runs.

*Example config files*
Example config files are provided for r0 = 6M, and black hole spins a = 0, 0.6 and 0.99, with two choices of grid parameters:
 - on a grid spacing of M/4, in the r* domain (-20, 100), up to l_max = 20. 
 - on a grid of spacing M/8, in the range r* \in (-30, 250) up to l_max = 30.
The first set is quicker to work with. The first set starts at ZZZ=001, and the second set at ZZZ=021, like this:

Set 1:  M/4 spacing, (-20,100) domain, l_max = 20.

00-060-001: a=0, r0=6, m=1
60-060-021: a=0.6, m=1
99-060-021: a=0.99, m=1

00-060-002: a=0, r0=6, m=2
60-060-002: a=0.6, m=2
99-060-002: a=0.99, m=2

Set 2:  M/8 spacing, (-30,250) domain, l_max = 30.

00-060-021: a=0, r0=6, m=1
60-060-021: a=0.6, m=1
99-060-021: a=0.99, m=1

00-060-022: a=0, r0=6, m=2
60-060-022: a=0.6, m=2
99-060-022: a=0.99, m=2

*Calculation for m != 0 (radiative)*

The notebook for calculating m != 0 modes is:
- metric_reconstruction_calc_radiative.nb
Sections 1 to 4 of the calculation are necessary setup. These sections should run in circa 10 minutes, depending on the parameters in "config".

Sections A and B calculate the MP components on a grid. Sections A and B typically take longer to run, depending on the grid parameters. You may choose to run either A or B, instead of both. Now that it is working efficiently, I think that A is the better approach to use. The command-line script only runs part A. 

The output of Section A (1D grid) are functions of r only, for the 10 components projected onto spin-weighted spherical harmonics. The output of Section B (2D grid) are functions of (r,\[Theta]) for the 10 components. 

*Calculation for m = 0 (static)*
For the static sector, use metric_reconstruction_calc_static.nb, after first running the three "Calc" notebooks for initial setup.
Sections 1 to 6 are necessary to produce output. The later sections are legacy code, used at some point for testing (and I can’t guarantee they will run).
Note that it is necessary to use a high number of digits ('prec') to solve accurately at high l. For example, choose parameters prec=100 for lmax=20.

Note that the metric perturbation calculated by this notebook will have non-zero mass and non-zero angular momentum in the interior region (r < r0). 
The mass and a.m. at infinity is determined (numerically) from the coE and coAM coefficients (in the notebook, see section 4. Jumps -> Format). 
This is the (Kerr analogue of) the Berndtson solution: a Lorenz-gauge MP that is regular at both the horizon and infinity, and with the correct “jump” in mass and AM across the particle radius. 
One can shift to a MP with the “right” mass and AM at infinity by adding a homogeneous Lorenz-gauge perturbation, but not while maintaining regularity.

In the Schwarzschild limit, this code produces an l=1, m=0 mode that, in Barack-Lousto-Sago coordinates, has h9 = 0 everywhere and h8 = 0 at the horizon. 
By contrast, standard data that is regular on the future horizon has h9 non-zero. To convert between the two, add a homogeneous Lorenz-gauge perturbation. 
The following note from Adam Pound (Aug 2024) explains this more clearly:
  For l=1, m=0, the unique Lorenz-gauge solution that is regular (in a regular coordinate system) and has zero angular momentum on the black hole is given in Eq. (D1) of https://arxiv.org/pdf/2006.11263. 
  This is the solution in Niels' data; it isn't the traditional Zerilli solution, which (as you say) is singular at the future horizon. 
  In Schwarzschild in the Lorenz gauge you can freely modify this solution, while maintaining regularity at the horizon, by adding a homogeneous solution corresponding to a spin on the black hole; that solution is given in (D3). 

*Visualisation*

The notebook for plotting the results is:
 - metric_reconstruction_visualisation.nb
The section on 1D plotting shows the radial functions; their values at r=r0; and the metric perturbation in Barack & Lousto variables. 

An important plot to check is that showing the “regularization” at r=r0 (example below). The individual spin parts (s=0,s=1,s=2) grow with (l+½) in power law fashion, whereas the combination falls off with a power law. 
This means there is a delicate cancellation at r=r0. On the other hand, away from r=r0, the convergence with l is exponential.

I recommend to check for evidence of deviation away from correct power-law behaviour at large l, and then adjust the parameter lplot downwards accordingly. 
In due course, we will apply a puncture scheme in the region (r0 - X, r0 + X), where X is some appropriate radial width.
