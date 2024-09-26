# Cubic interpolation spline
By default, it calculates the spline and its first two derivatives for the function x^3 - x.
The approximated function is changed manually (in the functions f, df and ddf of the Spline class).
When creating an object, the class constructor specifies the step, the left edge of the segment, the relaxation parameter, and the number of nodal points (the number of split segments + 1).
