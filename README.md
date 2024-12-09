Four Circles for Determining Roundness (ISO 4291:1985 (current version as of 2020)):

least squares circle (LSC)
minimum zone circle (MZC)
minimum circumscribed circle (MCC)
maximum inscribed circle (MIC)

These can all be solved as non-linear minimization problems.

# Least Squares Circle (LSC)

https://github.com/mdoube/BoneJ/blob/17ee483603afa8a7efb745512be60a29e093c94e/src/org/doube/geometry/FitCircle.java#L43

Minimize:

Σ ( (x_i - a)^2 + (y_i - b)^2 - r^2 )

where:

    Σ: summation over all data points (i = 1 to N)
    (a, b): center coordinates of the circle
    (x_i, y_i): coordinates of the i-th data point
    N: total number of data points
    r: radius of the circle

This objective function represents the sum of squared distances between each data point (x_i, y_i) and the fitted circle with center (a, b) and radius r. The minimization process aims to find the circle parameters that minimize this total squared distance, essentially providing the best fit based on minimizing these individual point-to-circle distances.

Will use Levenberg–Marquardt algorithm to solve this.

# Minimum Zone Circle (MZC)

Minimize: r

Subject to:

(x_i - a)^2 + (y_i - b)^2 <= r^2  for all data points (i = 1 to N)

where:

    r: radius of the circle
    (a, b): center coordinates of the circle
    (x_i, y_i): coordinates of the i-th data point
    N: total number of data points

This formulation minimizes the radius (r) while ensuring all data points are inside or on the circle through the squared distance constraint.

# Minimum Circumscribed Circle (MCC)

Also referred to as the "Smallest-circle problem".

See https://github.com/hbf/miniball for an example implementation.

Minimize:

 r

Subject to:

min( (x_i - a)^2 + (y_i - b)^2 ) = r^2  for at least one data point (i)

where the minimization is taken over all data points.

Here, the objective is still to minimize the radius (r). However, the constraint ensures that at least one data point lies exactly on the circle (i.e., the minimum squared distance from a data point to the center is equal to the square of the radius).

# Maximum Inscribed Circle (MIC)

Also referred to as the "Maximum empty circle problem".

Maximize:

 r

Subject to:

(x_i - a)^2 + (y_i - b)^2 >= r^2  for all data points (i = 1 to N)

This formulation maximizes the radius (r) while guaranteeing that all data points are outside or on the circle through the squared distance constraint (reversed inequality compared to MZC).

These objective functions, along with the associated constraints, are used within non-linear optimization algorithms to find the circle parameters (a, b, r) that satisfy the specific geometric conditions for each case (MZC, MCC, MIC).

## Surface Metrology

### Surface Texture Parameters
    Ft Pa Pdsm Phsc Pku Pm0 Pm2 Pm4
     Pp Ppc Pq Ps Psk Psm Pt Pv PVc
    Ra Ra1 Ra1l Ra7 Ra7l Rc Rcl Rdmd Rdmn Rdq
    Rdsk Rhsc Rk Rk+vk Rku Rm0 Rm2 Rm4 Rmax
    Rmq Rmr1 Rmr2 Rp Rpc Rpk Rpk* Rpk/k Rpk+k Rpm
    Rpm/3z Rpm7 Rpm7l Rpq Rq Rs Rsk Rsm
    Rt Rtwi Rv RVc Rvk Rvk* Rvk/k Rvm Rvo Rvq
    Ry Rz1max RzDIN RzJIS R3z
    Wa Wc Wcvx Wcvxl Wcvxm Wdq
    Weslp Weslpl Wlslp Wlslpl Wp Wpc Wpl Wpr Wq
    Ws Wseg Wsegl Wsm Wt Wtc Wv Wvda Wvdas
    Wvdc Wvdd Wvddl Wvdm Wvdmp Wvoid
    AR AW R Rx W Wte Wx
    Bearing Ratios (tpa/tpi, Pmr/Rmr): 10 Primary & 10 Roughness
    Htp Values (PHtp/RHtp): 10 Primary & 10 Roughness

### Filters

    Gaussian
    Spline Based Gaussian (with adjustable tension)
    Valley Suppression (ISO 13565-1 – 1996)
    Robust Spline-Based Gaussian (based on robust regression)
    Morphological Closing Filter (circular element applied to waviness profile)
    Morphological Opening Filter (circular element applied to waviness profile)

### Form removal

    Instrument reference (mean suppression)
    Least Squares Line
    Least Squares Arc
    Fixed Radius
    Least Squares Polynomial (user specified order)
    Spline Filter (for bandpass waviness with a user specified cutoff)
    Asphere (user defined coefficients with optional radius optimization)
    Free Form (user defined coordinates with optional pre-filtering)

### Surface Measurement Data Types

    SigmaSurf(*.sig)
    Jenoptik/Hommelwerke (*.pro; *.pip; *.asc; *.smd; *.hwp; *.hfm) 
    Mitutoyo (*.csv; *.dat; *.mes) 
    Form Talysurf (*.prf; *.fts; *.ruf, *.mod) 
    Talysurf 10 (*.ten) 
    Talysurf 6 (*.six) 
    Talyrond (*.str) 
    Surtronic 3+ (*.stp) 
    Renishaw (*.xls; *.zpx; *.rsf; *.out; *.txt) 
    Metrex (*.prf) 
    Mahr (*.pcd; *.prf; *.txt; *.pra; *.s2p) 
    TSK/Zeiss (*.rst; *.txt; *.tx1; *.tx2; *.nc; *.ynz) 
    Federal (*.dir) 
    Somicronic (*.pro; *.smd) 
    2Pros (*.hdr) 
    BTI (*.csv) 
    BM ASCII (*.pr) 
    FRT X-axis Text (*.txt) 
    Precision Devices (*.pdi; *.dat; *.nds) 
    Standard Data (*.sdf) 
    Hexagon Coordinate (*.yz) 
    CellView Text (*.txt) 
    Jenoptik Evovis XML (*.xml) 
    Zygo Profile (*.txt) 
    CMM Data (*.out; *.txt) 
    OmniSurf Settings (*.ini) 
    2-Column Radius (*.sip) 
    Insitutec (*.ist) 
    Multi-column SIG (*.sigh)