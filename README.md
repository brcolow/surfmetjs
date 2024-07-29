Four Circles for Determining Roundness (ISO 4291:1985 (current version as of 2020)):

least squares circle (LSC)
minimum zone circle (MZC)
minimum circumscribed circle (MCC)
maximum inscribed circle (MIC)

These can all be solved as non-linear minimization problems.

# Least Squares Circle (LSC)

https://github.com/AlliedToasters/circle-fit/blob/master/src/circle_fit/circle_fit.py
https://people.cas.uab.edu/~mosya/cl/LM.m
http://www.math.uab.edu/~chernov/cl/MATLABcircle.html
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

Here's a breakdown of why the term (r^2) is subtracted:

    The squared distance between a point and the circle center represents the ideal distance if the point lies perfectly on the circle.
    Subtracting r^2 from this term effectively removes the squared radius component from the ideal distance, focusing on minimizing the deviations from the ideal scenario.

Will use Levenberg–Marquardt algorithm to solve this.

The below three would need an adapted LM algorithm:

    Core functionality: The LM algorithm excels at iteratively refining estimates for parameters to minimize a non-linear objective function. This core functionality is applicable to all four circle fitting problems we discussed.
    Objective functions: MZC, MCC, and MIC all require minimizing or maximizing a function (radius) subject to geometric constraints expressed through inequalities. These constraints can be incorporated into the optimization process.

Implementation Details:

For MZC, MCC, and MIC, the objective function becomes the radius (r) to be minimized/maximized. The constraints can be handled in two ways:

    Penalty Function Approach: Introduce a penalty term in the objective function that heavily penalizes violations of the constraints.  The LM algorithm then minimizes the combined objective function, effectively pushing the solution towards configurations that satisfy the constraints.

    Projected Gradient Approach: Modify the update step within the LM algorithm to project the gradient direction onto a feasible space that satisfies the constraints.

Alternative Algorithms:

While the LM algorithm is a popular choice, other non-linear optimization techniques can also be used for these problems, such as:

    Sequential Quadratic Programming (SQP): This method can be more efficient for certain types of constraints.
    Simulated Annealing (SA): This is a stochastic approach that can be useful for complex problems with many local minima.

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

Minimize:

 r

Subject to:

min( (x_i - a)^2 + (y_i - b)^2 ) = r^2  for at least one data point (i)

where the minimization is taken over all data points.

Here, the objective is still to minimize the radius (r). However, the constraint ensures that at least one data point lies exactly on the circle (i.e., the minimum squared distance from a data point to the center is equal to the square of the radius).

# Maximum Inscribed Circle (MIC)

Maximize:

 r

Subject to:

(x_i - a)^2 + (y_i - b)^2 >= r^2  for all data points (i = 1 to N)

This formulation maximizes the radius (r) while guaranteeing that all data points are outside or on the circle through the squared distance constraint (reversed inequality compared to MZC).

These objective functions, along with the associated constraints, are used within non-linear optimization algorithms to find the circle parameters (a, b, r) that satisfy the specific geometric conditions for each case (MZC, MCC, MIC).

## Papers

The paper we will try and use for the implementation of MZC, MCC and MIC will be calvo2015:

Roque Calvo, Emilio Gómez,
Accurate evaluation of functional roundness from point coordinates,
Measurement,
Volume 73,
2015,
Pages 211-225,
ISSN 0263-2241,
https://doi.org/10.1016/j.measurement.2015.04.009.

