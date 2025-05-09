#  Four Circles for Determining Roundness

These circles quantify roundness using different geometric criteria. All can be solved as **non-linear minimization problems**.

---

##  Least Squares Circle (LSC)

> A best-fit circle that minimizes the sum of squared deviations from all points.

### 🔹 Objective

Minimize:
```
Σ [ (xᵢ - a)² + (yᵢ - b)² - r² ]
```

### 🔹 Where:
- `Σ`: summation over all data points (i = 1 to N)
- `(a, b)`: center coordinates of the circle
- `(xᵢ, yᵢ)`: coordinates of the i-th data point
- `r`: radius of the circle

This minimizes the squared difference between the actual radial distance of each point and the estimated radius.
**Solution method**: Levenberg–Marquardt algorithm.
🔗 [Java Implementation (BoneJ)](https://github.com/mdoube/BoneJ/blob/17ee483603afa8a7efb745512be60a29e093c94e/src/org/doube/geometry/FitCircle.java#L43)

---

##  Minimum Zone Circle (MZC)

> Finds the thinnest circular band (between two concentric circles) that fully encloses the profile.

### 🔹 Objective

Minimize:
```
r
```

### 🔹 Subject to:
```
(xᵢ - a)² + (yᵢ - b)² ≤ r²  ∀ i ∈ [1, N]
```

This ensures all data points lie **inside or on** the circle. The objective minimizes the radius to form the tightest possible enclosing zone.

---

##  Minimum Circumscribed Circle (MCC)

> Also known as the "Smallest-circle problem".

### 🔹 Objective

Minimize:
```
r
```

### 🔹 Subject to:
```
min( (xᵢ - a)² + (yᵢ - b)² ) = r²  for at least one i
```

Ensures **at least one point lies exactly** on the circle boundary while the rest are inside.
**Implementation example**:
🔗 [miniball library](https://github.com/hbf/miniball)

---

## 🟠 Maximum Inscribed Circle (MIC)

> Also known as the "Maximum empty circle problem".

### 🔹 Objective

Maximize:
```
r
```

### 🔹 Subject to:
```
(xᵢ - a)² + (yᵢ - b)² ≥ r²  ∀ i ∈ [1, N]
```

All data points must lie **outside or on** the circle. The optimization maximizes the radius of the largest possible empty circle that fits within the profile.

---

##  Summary

| Circle Type | Objective       | Constraint Type                            | Description                                |
|-------------|------------------|--------------------------------------------|--------------------------------------------|
| **LSC**     | Minimize squared radial error | No constraints | Best-fit circle using all points            |
| **MZC**     | Minimize `r`     | All points inside or on                    | Smallest enclosing zone width              |
| **MCC**     | Minimize `r`     | At least one point on circle               | Smallest enclosing circle                  |
| **MIC**     | Maximize `r`     | All points outside or on                   | Largest circle inside the profile          |


##  Deviation Types

### **`RONt` – Total Roundness (Total Runout)**
- **Definition**: Maximum radial distance between the highest peak and lowest valley from the reference circle.
- **Formula**: `RONt = max deviation - min deviation`
- **Also known as**: Total roundness deviation, roundness error.

---

### **`RONp` – Peak Roundness**
- **Definition**: The largest single peak (positive or negative) deviation from the reference circle.
- **Formula**: `RONp = max(|positive deviation|, |negative deviation|)`
- **Use case**: Indicates the worst-case point of form error.

---

### **`RONp Pos` – Peak Position**
- **Definition**: The position (in degrees) of the largest positive peak.

---

### **`RONv` – Valley Roundness**
- **Definition**: The largest single valley (negative) deviation from the reference circle.
- **Formula**: `RONv = min deviation` (in the negative direction)

---

### **`RONv Pos` – Valley Position**
- **Definition**: The position (in degrees) of the largest negative peak.

---

## ️ Reference Circle Types

### ** MIC – Maximum Inscribed Circle**
> The largest circle that fits entirely inside the profile.
> Acts like a *plug gauge*.

| **Parameter**  | **Meaning (from MIC)**                                                                       |
|----------------|-----------------------------------------------------------------------------------------------|
| `RONt`         | Max distance **outside** the MIC (how far the profile deviates from the inner circle)        |
| `RONp`         | Largest single point **outside** the MIC                                                     |
| `RONv`         | Typically `0` or undefined (no valleys **inside** the MIC)                                   |

---

### ** MCC – Minimum Circumscribed Circle**
> The smallest circle that fully encloses the profile.
> Acts like a *ring gauge*.

| **Parameter**  | **Meaning (from MCC)**                                                    |
|----------------|---------------------------------------------------------------------------|
| `RONt`         | Max distance **inside** the MCC (how far the profile retreats inward)     |
| `RONp`         | Largest positive deviation (farthest point from center of circle)         |
| `RONv`         | Largest negative deviation (most retracted point inside the circle)       |

---

### ** MZC – Minimum Zone Circle**
> A pair of concentric circles that bound the profile with the smallest radial separation.
> Gives **true form error**, minimizing both inward and outward deviations.

| **Parameter**  | **Meaning (from MZC band)**                                               |
|----------------|---------------------------------------------------------------------------|
| `RONt`         | Radial distance between the outer and inner circle (zone width)           |
| `RONp`         | Largest outward deviation from center (positive)                          |
| `RONv`         | Largest inward deviation from center (negative)                           |


#### ** LSC - Least-Squares circle

> The **Least Squares Circle (LSC)** is the best-fit circle that minimizes the sum of squared radial deviations of profile points.

| **Parameter**   | **Description** |
|------------------|-----------------|
| **RONt**         | **Total roundness deviation**: The difference between the **maximum** and **minimum** radial deviations from the LSC.<br> `RONt = max(deviation) - min(deviation)` |
| **RONp**         | **Peak deviation**: The **most positive** radial deviation (point farthest **outside** the LSC). |
| **RONv**         | **Valley deviation**: The **most negative** radial deviation (point farthest **inside** the LSC). |

---

### Surface Texture Parameters

#### P: Primary profile (unfiltered)

* **Pa**
    * **Name:** **Arithmetic Mean Deviation of the Primary Profile**
    * **Description:** Pa is the average of the absolute vertical deviations of the primary profile from its mean line over the evaluation length. It gives a general measure of the magnitude of texture variations on the unfiltered profile.
    * *Formula concept: Pa = (1/l) * integral from 0 to l of |Z(x)| dx (where l is the evaluation length and Z(x) is the profile height relative to the mean line).*

* **Pdsm**
    * **Name:** **Mean Width of the Primary Profile Elements** (*Note: `Psm` is the more standard term, `Pdsm` might imply specific element definition or be less common*)
    * **Description:** This parameter represents the average width (horizontal distance) of the profile elements (peak-valley pairs) found along the primary profile within the evaluation length. It gives an indication of the typical spacing of texture features on the unfiltered profile.

* **Phsc**
    * **Name:** **High Spot Count (Primary Profile)**
    * **Description:** Phsc measures the number of primary profile peaks projecting above a predefined height level (or threshold) within the evaluation length. It's used to quantify the density of significant peaks on the unfiltered surface.

* **Pku**
    * **Name:** **Kurtosis of the Primary Profile**
    * **Description:** Pku measures the "peakedness" or sharpness of the amplitude distribution curve (histogram of profile heights) for the primary profile.
        * `Pku > 3`: Indicates a spiky profile with sharp peaks/valleys.
        * `Pku < 3`: Indicates a bumpy profile with relatively flat peaks/valleys.
        * `Pku = 3`: Indicates a Gaussian distribution of heights, typical of random processes.

* **Pm0, Pm2, Pm4**
    * **Name:** **Profile Spectral Moments (Primary Profile)**
    * **Description:** These are statistical moments calculated from the Power Spectral Density (PSD) of the primary profile. They characterize the shape and frequency content of the profile.
        * **Pm0:** Zeroth moment. It is equal to the variance of the profile heights (`Pq^2`).
        * **Pm2:** Second moment. Related to the variance of the profile *slope*. Higher Pm2 suggests steeper features.
        * **Pm4:** Fourth moment. Related to the variance of the profile *curvature*. Higher Pm4 suggests sharper peaks and valleys.

* **Pp**
    * **Name:** **Maximum Peak Height of the Primary Profile**
    * **Description:** Pp is the height of the highest peak above the mean line within the evaluation length of the primary profile.

* **Ppc**
    * **Name:** **Peak Count (Primary Profile)**
    * **Description:** Ppc measures the number of peak/valley pairs per unit length on the primary profile, defined by crossing specific amplitude thresholds (hysteresis). It quantifies the density of local peaks and valleys exceeding certain criteria.

* **Pq**
    * **Name:** **Root Mean Square Deviation of the Primary Profile**
    * **Description:** Pq is the root mean square (RMS) average of the vertical deviations of the primary profile from its mean line over the evaluation length. Similar to Pa, but more sensitive to large peaks and valleys because it squares the deviations before averaging.
    * *Formula concept: Pq = sqrt[(1/l) * integral from 0 to l of Z^2(x) dx]*

* **Ps**
    * **Name:** **Mean Spacing of Local Peaks**
    * **Description:** Ps is the mean spacing between local peaks in the roughness profile that exceed a specified threshold relative to the mean line. Unlike Rsm (which includes both peaks and valleys), Ps focuses specifically on the frequency of peak occurrences. It is useful in assessing how often peaks occur along a surface, which is important in applications like friction or lubrication analysis.
    * *Common computation: Ps = (Evaluation Length) / (Number of peaks above threshold)*

* **Psk**
    * **Name:** **Skewness of the Primary Profile**
    * **Description:** Psk measures the asymmetry of the amplitude distribution curve (histogram of profile heights) for the primary profile.
        * `Psk > 0`: Indicates a profile dominated by peaks (e.g., heights skewed above the mean line).
        * `Psk < 0`: Indicates a profile dominated by valleys (e.g., heights skewed below the mean line).
        * `Psk = 0`: Indicates a symmetrical distribution of heights around the mean line.

* **Psm**
    * **Name:** **Mean Width of the Primary Profile Elements**
    * **Description:** Psm is the average horizontal distance (width) of the profile elements (a consecutive peak and the adjacent valley) along the mean line within the evaluation length of the primary profile. It characterizes the mean spacing of the texture features on the unfiltered profile.

* **Pt**
    * **Name:** **Total Height of the Primary Profile**
    * **Description:** Pt is the absolute vertical distance between the highest peak (Pp) and the deepest valley (Pv) within the evaluation length of the primary profile. It represents the overall peak-to-valley range of the unfiltered profile.
    * *Formula: Pt = Pp + Pv*

* **Pv**
    * **Name:** **Maximum Valley Depth of the Primary Profile**
    * **Description:** Pv is the absolute depth of the deepest valley below the mean line within the evaluation length of the primary profile.

* **PVc**
    * **Name:** **Profile Valley Core Fluid Retention Volume** (*Often related to the Rk family parameters, but applied to the Primary Profile*)
    * **Description:** This parameter typically relates to the Abbott-Firestone curve (material ratio curve) of the primary profile. It quantifies the void volume within the "core" region of the profile's valleys (often defined between material ratios like 10% and 80%), representing the profile's capacity to retain fluid (like lubricants) in the main working zone of the surface texture.

#### R: Roughness profile (after filtering to remove waviness)

* **Ra**
    * **Name:** **Arithmetic Mean Deviation of the Roughness Profile**
    * **Description:** Ra is the average of the absolute vertical deviations of the roughness profile from its mean line over the evaluation length. It is the most widely used roughness parameter globally, providing a general measure of texture magnitude.

* **Ra1**
    * **Name:** **Mean Height of the Highest Peak**
    * **Description:** Ra1 calculates the average height of the single highest peak within the evaluation length of the surface profile. It provides a focused measure of the most prominent surface asperity.

* **Ra1l**
    * **Name:** **Mean Depth of the Deepest Valley**
    * **Description:** Ra1l determines the average depth of the single deepest valley within the evaluation length. It offers insight into the most significant surface depression.

* **Ra7**
    * **Name:** **Mean Height of the Seven Highest Peaks**
    * **Description:** Ra7 computes the average height of the seven tallest peaks on the surface profile. This parameter provides a broader understanding of the surface's prominent features beyond just the single highest peak.

* **Ra7l**
    * **Name:** **Mean Depth of the Seven Deepest Valleys**
    * **Description:** Ra7l calculates the average depth of the seven deepest valleys on the surface profile, offering a comprehensive view of the most significant surface depressions.

* **Rc**
    * **Name:** **Mean Height of Profile Elements (Roughness)**
    * **Description:** Rc is the average height of the roughness profile elements (defined by consecutive peaks and valleys crossing the mean line) within the evaluation length.

* **Rcl**
    * **Name:** **Correlation Length**
    * **Description:** Rcl measures the horizontal distance over which the surface profile's autocorrelation function decays to a specified value (often 0.1 or 0.2). It indicates the scale of the surface texture features, with longer correlation lengths signifying more extended surface features.

* **Rdmd**
    * **Name:** **Mean Depth of Roughness Motifs**
    * **Description:** Rdmd represents the average depth of the individual roughness motifs found along a surface profile. Motif-based analysis, as defined in ISO 12085 (the CNOMO method), segments the surface profile into a series of motifs—discrete elements typically bounded by significant local peaks and valleys that characterize the repetitive surface pattern. Each motif is defined based on rules involving curvature, prominence, and spacing thresholds. The depth of a motif is the vertical distance between its peak and valley. Rdmd is calculated as the average of these depths over all detected motifs.
    * **Use Case:** Useful for characterizing surfaces with repeating patterns (e.g., machined textures), where average depth of local features is more informative than point-wise extremes.
    * *Formula concept: Rdmd = (1/N) * Σᵢ(depth of motifᵢ), for N detected motifs*

* **Rdmn**
    * **Name:** **Mean Height of Roughness Motifs**
    * **Description:** Rdmn is the average height of the individual roughness motifs, where the height of a motif typically refers to the peak value measured relative to a common datum such as the profile mean line. In motif-based analysis, these motifs reflect regular or quasi-regular geometric features like ridges or dimples. While Rdmd focuses on the vertical difference between the peak and valley of a single motif, Rdmn isolates the vertical position of the peak relative to a reference. This is particularly useful in applications involving contact surfaces or tribological performance.
    * **Use Case:** Helps assess the average prominence of surface structures, especially when predicting initial contact conditions or coating thickness variations.
    * *Formula concept: Rdmn = (1/N) * Σᵢ(height of motifᵢ peak), for N detected motifs*

* **Rdq (or RΔq)**
    * **Name:** **Root Mean Square Slope of the Roughness Profile**
    * **Description:** Rdq measures the root mean square value of the local slope of the roughness profile over the evaluation length. It provides information about the steepness or angularity of the surface texture.

* **Rdsk**
    * **Name:** **Reduced Skewness (Core Skewness of the Roughness Profile)**
    * **Description:** Rdsk is a non-standard but commonly used derived parameter that quantifies the **skewness of the core material portion** of a roughness profile. Unlike Rsk, which calculates skewness over the entire roughness profile, Rdsk focuses only on the central section of the **material ratio curve** (typically from 10% to 80% or 20% to 80% material ratios), filtering out extreme peaks and valleys.

      This helps provide a more stable and functionally relevant measurement of asymmetry in the surface topography, especially important for **contact, sealing, or wear applications**, where the central material is load-bearing.

    * **Use Case:** Rdsk is commonly used in advanced surface texture analysis to minimize the influence of surface anomalies, scratches, or rare outliers. It is particularly useful in automotive, tribological, and aerospace applications where surface behavior under load is critical.

    * *Formula concept:*
      Rdsk = (1 / Rq_c³) × (1/L_c) × ∫ (Z_c(x))³ dx
      where:
      - **Z_c(x)** is the profile height within the selected material ratio range (core)
      - **Rq_c** is the root mean square roughness within that same core
      - **L_c** is the length of the core region

    * **Interpretation:**
      - **Rdsk > 0**: More high peaks than deep valleys in the core
      - **Rdsk < 0**: More valleys than peaks in the core
      - **Rdsk ≈ 0**: Symmetric core distribution

* **Rhsc**
    * **Name:** **High Spot Count (Roughness Profile)**
    * **Description:** Rhsc measures the number of roughness profile peaks projecting above a predefined height level (or threshold) within the evaluation length.

* **Rk**
    * **Name:** **Core Roughness Depth** (ISO 13565-2)
    * **Description:** Rk represents the depth of the core (main working) part of the surface profile, determined from the Abbott-Firestone (material ratio) curve. It often relates to the wear behavior of a surface during its stable operational life.

* **Rk+vk** (*Note: Often written as Rk + Rvk*)
    * **Name:** *Sum of Core Roughness Depth and Reduced Valley Depth*
    * **Description:** This is not a single parameter but the sum of two parameters from the Rk family (ISO 13565-2): `Rk` (Core Roughness Depth) and `Rvk` (Reduced Valley Depth). `Rvk` relates to the depth of the deeper valleys below the core roughness. The sum `Rk + Rvk` indicates the depth from the mean line (or a reference) to the bottom of the primary valley structures relevant for fluid retention or debris entrapment.

* **Rku**
    * **Name:** **Kurtosis of the Roughness Profile**
    * **Description:** Rku measures the "peakedness" or sharpness of the amplitude distribution curve (histogram of profile heights) for the roughness profile.
        * `Rku > 3`: Indicates a spiky roughness profile.
        * `Rku < 3`: Indicates a bumpy roughness profile with flattened peaks/valleys.
        * `Rku = 3`: Indicates a Gaussian distribution of roughness heights.

* **Rm0, Rm2, Rm4**
    * **Name:** **Roughness Profile Spectral Moments**
    * **Description:** These are statistical moments calculated from the Power Spectral Density (PSD) of the *roughness* profile.
        * **Rm0:** Zeroth moment. Equal to the variance of the roughness profile heights (`Rq^2`).
        * **Rm2:** Second moment. Related to the variance of the roughness profile *slope*.
        * **Rm4:** Fourth moment. Related to the variance of the roughness profile *curvature*.

* **Rmax**
    * **Name:** **Maximum Roughness Height** *(Definition varies, often deprecated in favor of Rz or Rt in ISO)*
    * **Description:** Rmax typically refers to the largest single peak-to-valley height found within *any single sampling length* inside the overall evaluation length. *Caution: The definition of Rmax can differ significantly between standards (e.g., ISO 4287:1984, JIS B 0601, ASME B46.1). In modern ISO 4287 (1997 onwards), `Rz` often refers to the average maximum height over sampling lengths, and `Rt` is the total maximum height over the evaluation length. Rmax is less commonly used or has specific definitions in certain contexts.*

* **Rmq**
    * **Name:** **Root Mean Square Slope of the Roughness Profile**
    * **Description:** Rmq characterizes the average steepness of the roughness profile. It is the root mean square value of the slope (first derivative) of the roughness profile. Higher Rmq values indicate sharper texture features.
    * *Formula concept: Rmq = sqrt[(1/l) * ∫₀ˡ (dZ/dx)² dx]*

* **Rmr1**
    * **Name:** **Material Ratio at Level 1**
    * **Description:** Rmr1 is the percentage of the profile’s evaluation length at which the material ratio curve (Abbott–Firestone curve) intersects a height level corresponding to 5% of the bearing area. It reflects how much material is present near the peak of the surface.
    * *Derived from the material ratio curve based on a defined height level (typically at 5%).*

* **Rmr2**
    * **Name:** **Material Ratio at Level 2**
    * **Description:** Rmr2 is the material ratio corresponding to a deeper section of the surface (typically 50–80% bearing ratio), indicating how much of the profile is load-bearing at that level.
    * *Derived similarly to Rmr1 but at a different vertical slice of the Abbott–Firestone curve.*

* **Rp**
    * **Name:** **Maximum Profile Peak Height**
    * **Description:** Rp is the height of the highest peak of the roughness profile above the mean line within the evaluation length. It highlights the most prominent upward deviation.
    * *Rp = max(Z(x)) for Z(x) > 0*

* **Rpc**
    * **Name:** **Peak Count**
    * **Description:** Rpc counts the number of peaks per unit length on the roughness profile. It provides insight into surface density of asperities.
    * *Rpc = Number of peaks / Evaluation length*

* **Rpk**
    * **Name:** **Reduced Peak Height**
    * **Description:** Rpk represents the average height of the most prominent peaks above the core roughness area, usually derived from the upper 5% of the material ratio curve. It is used in bearing surface evaluations.
    * *Part of the Rk family parameters from the bearing area curve.*

* **Rpk***
    * **Name:** **Corrected Reduced Peak Height**
    * **Description:** Rpk* is a modified or corrected form of Rpk, adjusted to reduce sensitivity to noise or to fit specific standards. Not universally standardized; interpretation depends on measurement context.

* **Rpk/k**
    * **Name:** **Ratio of Reduced Peak Height to Core Roughness Depth**
    * **Description:** Rpk/k is a dimensionless ratio that compares the height of the peaks (Rpk) to the central core roughness depth (Rk), indicating how pronounced the peaks are relative to the main surface structure.

* **Rpk+k**
    * **Name:** **Sum of Reduced Peak Height and Core Roughness**
    * **Description:** Rpk+k quantifies the combined contribution of the peak region and the core roughness depth, often used as a combined indicator of initial contact behavior and surface load-bearing capacity.

* **Rpm**
    * **Name:** **Mean Height of Peaks**
    * **Description:** Rpm is the average height of all identified peaks within the profile relative to the mean line, typically used for characterizing general peak morphology.

* **Rpm/3z**
    * **Name:** **Normalized Mean Peak Height**
    * **Description:** Rpm/3z is the mean peak height normalized by three times the profile standard deviation (3σ). It helps compare peak prominence across surfaces with different roughness scales.

* **Rpm7**
    * **Name:** **Mean Height of 7 Highest Peaks**
    * **Description:** Rpm7 averages the height of the seven highest peaks above the mean line. It is used to evaluate surface peak outliers that may influence performance like wear or lubrication.

* **Rpm7l**
    * **Name:** **Mean Height of 7 Lowest Valleys**
    * **Description:** Rpm7l calculates the mean depth of the seven lowest valleys beneath the mean line. This metric helps quantify outlier valleys which might affect material retention or coating thickness.

* **Rpq**
    * **Name:** **Root Mean Square Peak Height**
    * **Description:** Rpq is the root mean square value of all peak heights above the mean line within the evaluation length. It gives a statistical measure of peak height variation.

* **Rq**
    * **Name:** **Root Mean Square Roughness**
    * **Description:** Rq is the root mean square of the roughness profile deviations from the mean line. It is a widely used parameter that reflects overall surface texture energy.
    * *Formula concept: Rq = sqrt[(1/l) * ∫₀ˡ Z(x)² dx]*

* **Rs**
    * **Name:** **Mean Spacing of Profile Irregularities**
    * **Description:** Rs is the mean spacing between adjacent local peaks or valleys (depending on definition) along the profile, providing information about surface texture frequency.

* **Rsk**
    * **Name:** **Skewness of the Roughness Profile**
    * **Description:** Rsk quantifies the asymmetry of the roughness profile. A positive Rsk means more high peaks, and a negative Rsk indicates deeper valleys.
    * *Formula concept: Rsk = (1/Rq³) * [(1/l) * ∫₀ˡ Z(x)³ dx]*

* **Rsm**
    * **Name:** **Mean Spacing of Profile Irregularities**
    * **Description:** Rsm is the average spacing between adjacent peaks and valleys along the profile, often used in evaluating wave-like features or textures.
    * *Measured peak-to-peak or valley-to-valley spacing averaged over the profile.*

* **Rt**
    * **Name:** **Total Height of the Roughness Profile**
    * **Description:** Rt is the vertical distance between the highest peak and the lowest valley in the roughness profile within the evaluation length. It represents the extreme height variation.
    * *Rt = Rp + Rv*

* **Rtwi**
    * **Name:** **Waviness Total Height**
    * **Description:** Rtwi is the total height of the waviness profile, capturing the distance from the highest peak to the lowest valley in the waviness component (filtered profile). It reflects macro-texture features.
    * *Measured over the waviness profile.*

* **Rv**
    * **Name:** **Maximum Profile Valley Depth**
    * **Description:** Rv is the maximum depth of the valleys below the mean line within the evaluation length. It characterizes the deepest groove or pit on the profile.
    * *Rv = min(Z(x)) for Z(x) < 0*

* **RVc**
    * **Name:** **Core Valley Depth**
    * **Description:** RVc is the average depth of valleys in the core roughness area. It helps assess the load-bearing capacity of the surface by measuring valley depths below the core region.

* **Rvk**
    * **Name:** **Reduced Valley Depth**
    * **Description:** Rvk is the average depth of the valleys below the core roughness region, usually taken from the lower 5% of the Abbott–Firestone curve. It indicates potential for lubricant retention.
    * *Part of the Rk parameter family.*

* **Rvk***
    * **Name:** **Corrected Reduced Valley Depth**
    * **Description:** Rvk* is a refined or corrected version of Rvk to improve measurement repeatability or match a particular standard or analysis context.

* **Rvk/k**
    * **Name:** **Ratio of Reduced Valley Depth to Core Roughness Depth**
    * **Description:** Rvk/k compares the depth of valleys (Rvk) to the core roughness (Rk), offering a normalized insight into valley prominence relative to overall surface roughness.

* **Rvm**
    * **Name:** **Mean Depth of Valleys**
    * **Description:** Rvm is the average depth of all identified valleys in the profile relative to the mean line. It helps characterize valley morphology beyond just the deepest point.

* **Rvo**
    * **Name:** **Oil Retention Volume Indicator**
    * **Description:** Rvo estimates the volume that valleys might hold oil or fluid. Often derived by integrating valley depths below the core roughness region.
    * *Related to Rvk and valley geometry.*

* **Rvq**
    * **Name:** **Root Mean Square Valley Depth**
    * **Description:** Rvq is the RMS value of the depths of all valleys below the mean line, indicating statistical variation of valley depths.

* **Ry**
    * **Name:** **Maximum Height of Profile**
    * **Description:** Ry is a legacy term equivalent to Rt in many standards, representing the maximum peak-to-valley height within the sampling length. Still used in some industrial contexts.

* **Rz1max**
    * **Name:** **Maximum Ten-Point Height**
    * **Description:** Rz1max measures the vertical distance between the five highest peaks and five deepest valleys within one sampling length, taking the maximum across all such segments.

* **RzDIN**
    * **Name:** **Ten-Point Height (DIN Standard)**
    * **Description:** RzDIN is the average height difference between the five highest peaks and five deepest valleys over the entire evaluation length, defined per DIN standards.

* **RzJIS**
    * **Name:** **Ten-Point Height (JIS Standard)**
    * **Description:** RzJIS is a variant of Rz defined by Japanese standards, often differing in how peaks and valleys are selected or averaged.

* **R3z**
    * **Name:** **Three-Zone Mean Roughness Height**
    * **Description:** R3z is the mean of three Rz values taken from consecutive sampling lengths, providing a smoothed and averaged ten-point height value.

* **Rx**
    * **Name:** **Extreme Roughness Height**
    * **Description:** Rx can denote the highest single peak-to-valley height across the entire surface profile, serving as a robustness measure for the worst-case surface deviation.
    * *May also be used generically as an extreme-value roughness metric.*

#### W: Waviness profile (after filtering to remove roughness)

* **Wa**
    * **Name:** **Arithmetic Mean Deviation of the Waviness Profile**
    * **Description:** Wa is the average of the absolute values of the waviness profile deviations from its mean line over the evaluation length. It is analogous to Ra but considers the long-wavelength components of the surface profile.

* **Wc**
    * **Name:** **Cutoff Wavelength for Waviness**
    * **Description:** Wc is the wavelength threshold used to separate waviness from roughness using a filter. It defines the spatial frequency below which features are classified as waviness.

* **Wcvx**
    * **Name:** **Waviness Profile Convexity**
    * **Description:** Wcvx quantifies the curvature or convexity of the waviness profile. It may be used in applications where profile curvature affects performance, such as sealing surfaces.

* **Wcvxl**
    * **Name:** **Long-Wavelength Convexity**
    * **Description:** Wcvxl isolates the convexity over longer wavelength components, further emphasizing the broad surface shape rather than fine texture.

* **Wcvxm**
    * **Name:** **Mean Convexity of Waviness Profile**
    * **Description:** Wcvxm represents the mean convex shape derived from the waviness profile, averaged across the evaluation length.

* **Wdq**
    * **Name:** **Root Mean Square Slope of the Waviness Profile**
    * **Description:** Wdq is the root mean square value of the slope of the waviness profile. It captures the general steepness of the broader form of the surface.

* **Weslp**
    * **Name:** **Waviness Edge Slope**
    * **Description:** Weslp measures the average slope of the waviness profile at the boundaries of the evaluation length. It's useful in assessing how the surface shape transitions near the edges.

* **Weslpl**
    * **Name:** **Lower Waviness Edge Slope**
    * **Description:** Weslpl is the average slope specifically at the lower end (start) of the evaluation length. It may help identify tapering or rising trends at the surface entry.

* **Wlslp**
    * **Name:** **Waviness Leading Slope**
    * **Description:** Wlslp characterizes the slope of the profile as it begins—used to assess initial surface rise or drop-off.

* **Wlslpl**
    * **Name:** **Lower Leading Slope of Waviness**
    * **Description:** Wlslpl isolates the slope feature at the very beginning of the evaluation profile, capturing the first major surface rise.

* **Wp**
    * **Name:** **Maximum Peak Height of Waviness Profile**
    * **Description:** Wp is the maximum height of the waviness profile above the mean waviness line, indicating the most prominent upward deviation in the long-wavelength component.

* **Wpc**
    * **Name:** **Waviness Peak Count**
    * **Description:** Wpc counts the number of peaks per unit length in the waviness profile, describing the frequency of broad surface undulations.

* **Wpl**
    * **Name:** **Mean Height of Waviness Peaks**
    * **Description:** Wpl is the average height of peaks in the waviness profile relative to the mean line, useful in evaluating periodic surface patterns.

* **Wpr**
    * **Name:** **Mean Peak Spacing in Waviness Profile**
    * **Description:** Wpr gives the average spacing between peaks in the waviness profile, helping quantify the wavelength characteristics of surface undulations.

* **Wq**
    * **Name:** **Root Mean Square Waviness**
    * **Description:** Wq is the RMS value of the waviness profile deviations from the mean waviness line, analogous to Rq but for long-wavelength components.

* **Ws**
    * **Name:** **Mean Spacing of Local Peaks in Waviness Profile**
    * **Description:** Ws is the average distance between local peaks in the waviness profile, providing a measure of the spacing of surface undulations after roughness filtering has been applied.

* **Wseg**
    * **Name:** **Waviness Segment Count**
    * **Description:** Wseg denotes the number of discrete segments or features identified in the waviness profile that meet defined criteria (e.g., peak-to-valley spacing, prominence). It can be used in pattern recognition or diagnostics of repetitive form deviations.

* **Wsegl**
    * **Name:** **Mean Length of Waviness Segments**
    * **Description:** Wsegl is the average length of the segments identified in the waviness profile. It is useful for assessing the characteristic scale of recurring waviness units or undulation cycles.

* **Wsm**
    * **Name:** **Mean Spacing of Profile Irregularities in Waviness Profile**
    * **Description:** Wsm is the mean spacing between peaks and valleys (zero crossings or curvature inflection points) in the waviness profile. It provides a general sense of surface undulation frequency.

* **Wt**
    * **Name:** **Total Height of the Waviness Profile**
    * **Description:** Wt is the vertical distance between the highest peak and the lowest valley in the waviness profile over the evaluation length. It is analogous to Rt but applies to the low-frequency portion of the surface.

* **Wtc**
    * **Name:** **Waviness Total Convexity**
    * **Description:** Wtc measures the total convex curvature or bowing of the waviness profile across the evaluation length. It may be derived by integrating local convex portions and is important in form error diagnostics.

* **Wv**
    * **Name:** **Maximum Valley Depth in Waviness Profile**
    * **Description:** Wv is the depth of the deepest valley in the waviness profile, measured from the mean waviness line. This metric highlights the most prominent concavity after roughness filtering.

* **Wvda**
    * **Name:** **Average Waviness Valley Depth**
    * **Description:** Wvda is the average of all valley depths found in the waviness profile, offering a statistical measure of concavity frequency and magnitude.

* **Wvdas**
    * **Name:** **Standard Deviation of Valley Depths in Waviness Profile**
    * **Description:** Wvdas quantifies the variability of valley depths within the waviness profile. High values may indicate irregular wear or distortion patterns.

* **Wvdc**
    * **Name:** **Cumulative Valley Depth**
    * **Description:** Wvdc is the sum of all valley depths across the waviness profile, potentially representing retained fluid volume or concavity over large areas.

* **Wvdd**
    * **Name:** **Difference Between Deepest and Shallowest Valleys**
    * **Description:** Wvdd measures the range in valley depths across the waviness profile, offering a spread metric for evaluating surface deformation consistency.

* **Wvddl**
    * **Name:** **Deepest Valley Length**
    * **Description:** Wvddl is the horizontal length of the deepest valley within the waviness profile. It can be important in tribological or sealing surface assessments.

* **Wvdm**
    * **Name:** **Mean Valley Depth Below Mean Line**
    * **Description:** Wvdm is the average depth of valleys that fall below the mean waviness line, excluding minor features or peaks.

* **Wvdmp**
    * **Name:** **Peak-Referenced Mean Valley Depth**
    * **Description:** Wvdmp measures the average valley depth referenced from the nearest adjacent peak rather than from the mean line, helpful in evaluating local pit geometry.

* **Wvoid**
    * **Name:** **Waviness Void Ratio**
    * **Description:** Wvoid estimates the proportion of the surface evaluation length composed of valleys or concave regions in the waviness profile. It may be relevant for predicting lubricant retention or void formation under load.

* **Wte**
    * **Name:** **Waviness Total Edge Deviation**
    * **Description:** Wte quantifies the deviation of the waviness profile near the edges of the evaluation length. It is useful in evaluating edge artifacts, form mismatch, or tapering trends.

* **Wx**
    * **Name:** **Extreme Waviness Height**
    * **Description:** Wx denotes the maximum peak-to-valley height in the waviness profile. Analogous to Rx, but specific to the filtered waviness component.

####  Bearing Ratios & Material Ratios

* **Pmr / Rmr**
    * **Name:** **Material Ratio (Primary / Roughness Profile)**
    * **Description:** Pmr (Primary Material Ratio) and Rmr (Roughness Material Ratio) represent the percentage of the profile that is solid material at a given height level, measured from the highest peak downward. These values are typically sampled at specific percentage levels (e.g., 0%, 5%, ..., 100%) and are used to generate the Abbott–Firestone Curve.
    * **Use Case:** Evaluates load-bearing capacity and wear behavior of surfaces.
    * *Formula concept: Rmr(c) = (material length at height c) / evaluation length × 100%*

* **tpi / tpa**
    * **Name:** **Intersection Points of the Material Ratio Curve**
    * **Description:** tpi (intersection on the Primary profile) and tpa (on the Roughness profile) are the heights at which a specified material ratio occurs. These are used inversely to Pmr/Rmr to find the vertical height at which a given percentage of the surface is composed of material.
    * **Use Case:** Commonly extracted for 10 levels (e.g., 0.1 to 1.0 in 0.1 increments) to characterize the shape of the material ratio curve.
    * *Formula concept: tpi(Pmr) = profile height at material ratio Pmr*

---

####  Htp Values (Peak Height at Material Ratio Level)

* **PHtp / RHtp**
    * **Name:** **Height at Material Ratio (Primary / Roughness Profile)**
    * **Description:** PHtp and RHtp represent the height of the primary or roughness profile at a defined material ratio level. Typically sampled at uniform intervals (e.g., every 10% from 0% to 100%), these values indicate the vertical location of material supporting regions and are derived from the inverse of the Abbott–Firestone curve.
    * **Use Case:** Valuable in tribological and contact analysis, where the support height for a given percentage of the profile is functionally relevant (e.g., for oil film retention or sealing).
    * *Formula concept: PHtp(x%) = height where material ratio curve = x%*

**Typical Sample Points:**
- Material Ratios: 0%, 10%, 20%, ..., 100%
- Values returned: Pmr10, Pmr20, ..., Pmr100 (or PHtp10, PHtp20, etc.)

**Profiles:**
- **Primary (P)**: Includes waviness and roughness
- **Roughness (R)**: Includes only short-wavelength features, after filtering waviness out

---

These values are often visualized together on **bearing ratio graphs** or **Abbott–Firestone curves**, where:
- X-axis = Material ratio (%)
- Y-axis = Profile height (μm)

###  Motif-Based Surface Texture Analysis (ISO 12085 Overview)

Motif analysis, sometimes referred to as the CNOMO method (from French standards), is a form of profile decomposition where the surface trace is segmented into individual, geometrically defined features—called **motifs**. This approach is different from classical roughness parameters that are based on filtering (e.g., Gaussian filters in ISO 4287). Instead, motifs are determined by:
- Identifying **inflection points** (where curvature changes)
- Finding **local maxima and minima**
- Applying rules to merge insignificant motifs based on height or spacing thresholds

**Advantages of Motif Analysis:**
- Captures **repeating texture units** like threads, machine marks, or etch pits
- Less sensitive to extreme values or measurement noise
- Intuitive for evaluating real-world shape features (e.g., in bearing surfaces or sealing zones)

**Typical Parameters from Motif Analysis Include:**
- **Sm**: Mean spacing of motifs
- **R**: Height of individual motifs
- **Rmr**: Material ratio at motif level
- **Rdmd** and **Rdmn**: Depth and height averages over all motifs

This technique is particularly useful for **functional characterization** where micro-geometry of repeating features affects wear, friction, sealing, or coating adhesion.

###  Filters

# Surface Texture Filters in Metrology

Standard UPR (undulations per revolution) ranges used for filters: 1-15 UPR, 1-50 UPR, 1-150 UPR, 1-500 UPR.
---

##  1. Gaussian Filter

The Gaussian filter is a linear low-pass filter used to separate roughness and waviness components.

**Equation:**

The filtered profile `z_f(x)` is the convolution of the profile `z(x)` with a Gaussian kernel `G(x)`:

```
z_f(x) = ∫ z(u) · G(x - u) du
G(x) = (1 / √(2π)λ_c) · exp( -x² / (2λ_c²) )
```

- `λ_c`: Cutoff wavelength (e.g., 0.8 mm for roughness, 2.5 mm for waviness)

---

##  2. Spline-Based Gaussian Filter (Adjustable Tension)

This filter uses spline smoothing with tension to approximate a Gaussian-like response, providing better control over boundary effects.

**Minimization functional:**

```
J[z_f] = Σ (z_i - z_f(x_i))² + α ∫ (d²z_f/dx²)² dx
```

- `α`: Tension parameter controlling stiffness of the spline

---

##  3. Valley Suppression Filter (ISO 13565-1)

Used for functional surfaces with deep valleys (e.g., plateau honing). Suppresses deep valleys to avoid biasing the mean line.

**Algorithm:**

1. Apply Gaussian filter
2. Replace points deviating below a threshold
3. Iterate until convergence

This is a **truncated Gaussian filter** that ignores excessive valleys during smoothing.

---

##  4. Robust Spline-Based Gaussian Filter

Improves upon spline-based filters by using robust regression to reduce sensitivity to outliers (e.g., scratches, dirt).

**Functional:**

```
Minimize: Σ ρ(z_i - z_f(x_i)) + α ∫ (d²z_f/dx²)² dx
```

- `ρ(e)`: Robust loss function (e.g., Huber or Tukey's biweight)
- Handles outliers more gracefully than squared-error loss

---

##  5. Morphological Closing Filter

A nonlinear filter that **removes valleys** (fills pits). Based on morphological operations using a circular structuring element.

**Operation:**

```
z_close(x) = dilation(z, B) followed by erosion(z, B)
```

Where `B` is a circular structuring element of specified radius.

---

## 6. Morphological Opening Filter

The dual of closing — **removes peaks** (cuts spikes).

**Operation:**

```
z_open(x) = erosion(z, B) followed by dilation(z, B)
```

---

## Summary Table

| Filter Name                        | Type     | Purpose                          | Notes                   |
|-----------------------------------|----------|----------------------------------|-------------------------|
| Gaussian Filter                   | Linear   | Separate roughness/waviness      | ISO standard            |
| Spline-Based Gaussian             | Linear   | Smooth with boundary control     | Tension-adjustable      |
| Valley Suppression (ISO 13565-1)  | Nonlinear| Remove deep valleys              | Iterative, asymmetric   |
| Robust Spline Gaussian            | Nonlinear| Reduce outlier influence         | Uses robust loss        |
| Morphological Closing             | Nonlinear| Suppress valleys                 | Dilation then erosion   |
| Morphological Opening             | Nonlinear| Suppress peaks                   | Erosion then dilation   |

###  Form Removal Methods

---

## 🔹 Instrument Reference (Mean Suppression)

**Definition**: Subtracts the arithmetic mean height from the profile.

**Mathematics**:
Given a discrete profile `z(x_i)` over `i = 1, ..., N`:

```
z̄ = (1 / N) * Σ z(x_i)
z_res(x_i) = z(x_i) - z̄
```

---

## 🔹 Least Squares Line

**Definition**: Fits a straight line `z(x) = a * x + b` to the profile via least squares and subtracts it.

**Minimization**:
Find `a`, `b` that minimize:

```
Σ [z(x_i) - (a * x_i + b)]²
```

**Form Removal**:

```
z_res(x_i) = z(x_i) - (a * x_i + b)
```

---

## 🔹 Least Squares Arc

**Definition**: Fits a circle `z(x) = sqrt(r² - (x - x₀)²) + z₀` via least squares.

**Minimization**:
Find center `(x₀, z₀)` and radius `r` that minimize:

```
Σ [z(x_i) - (z₀ ± sqrt(r² - (x_i - x₀)²))]²
```

**Residual**:

```
z_res(x_i) = z(x_i) - z_arc(x_i)
```

---

## 🔹 Fixed Radius

**Definition**: Subtracts a circle with a **user-specified radius** and origin (or fitted height).

**Form**:
If radius `r` and center `x₀` are fixed:

```
z_ref(x_i) = sqrt(r² - (x_i - x₀)²) + z₀
```

Then subtract `z_ref` as usual.

---

## 🔹 Least Squares Polynomial (User-Specified Order)

**Definition**: Fits a polynomial of order `n` via least squares:

```
z(x) = a₀ + a₁ x + a₂ x² + ... + aₙ xⁿ
```

**Minimization**:
Find coefficients `{a_k}` minimizing:

```
Σ [z(x_i) - Σ a_k x_i^k]²
```

Then subtract the fitted polynomial from the profile.

---

## 🔹 Spline Filter (Bandpass Waviness)

**Definition**: Fits a smoothing spline with user-defined cutoff wavelength(s), often used to extract waviness after form is removed.

**Objective Function**:

```
J[z_f] = Σ [z(x_i) - z_f(x_i)]² + α ∫ (d²z_f/dx²)² dx
```

- `α` controls the stiffness of the spline (linked to cutoff wavelength)

---

## 🔹 Asphere (User-Defined Coefficients)

**Definition**: Subtracts a user-defined even-order aspheric surface model:

```
z(r) = (r² / [R * (1 + sqrt(1 - (1 + k) * r² / R²))]) + Σ A_{2n} * r^{2n}
```

Where:
- `r` = radial distance
- `R` = vertex radius of curvature (optionally optimized)
- `k` = conic constant
- `A_{2n}` = aspheric coefficients

---

## 🔹 Free Form (User-Defined Surface)

**Definition**: Subtracts a user-supplied coordinate-based height map or mathematical model, optionally after applying a pre-filter (e.g., low-pass to smooth input).

**Generic Form**:

```
z_res(x_i) = z(x_i) - z_form(x_i)
```

`z_form(x)` may be an arbitrary function, lookup table, or interpolated mesh.

---

###  Surface Measurement Data Types

- SigmaSurf: `*.sig`
- Jenoptik/Hommelwerke: `*.pro`, `*.pip`, `*.asc`, `*.smd`, `*.hwp`, `*.hfm`
- Mitutoyo: `*.csv`, `*.dat`, `*.mes`
- Form Talysurf: `*.prf`, `*.fts`, `*.ruf`, `*.mod`
- Talysurf 10: `*.ten`
- Talysurf 6: `*.six`
- Talyrond: `*.str`
- Surtronic 3+: `*.stp`
- Renishaw: `*.xls`, `*.zpx`, `*.rsf`, `*.out`, `*.txt`
- Metrex: `*.prf`
- Mahr: `*.pcd`, `*.prf`, `*.txt`, `*.pra`, `*.s2p`
- TSK/Zeiss: `*.rst`, `*.txt`, `*.tx1`, `*.tx2`, `*.nc`, `*.ynz`
- Federal: `*.dir`
- Somicronic: `*.pro`, `*.smd`
- 2Pros: `*.hdr`
- BTI: `*.csv`
- BM ASCII: `*.pr`
- FRT X-axis Text: `*.txt`
- Precision Devices: `*.pdi`, `*.dat`, `*.nds`
- Standard Data: `*.sdf`
- Hexagon Coordinate: `*.yz`
- CellView Text: `*.txt`
- Jenoptik Evovis XML: `*.xml`
- Zygo Profile: `*.txt`
- CMM Data: `*.out`, `*.txt`
- OmniSurf Settings: `*.ini`
- 2-Column Radius: `*.sip`
- Insitutec: `*.ist`
- Multi-column SIG: `*.sigh`
