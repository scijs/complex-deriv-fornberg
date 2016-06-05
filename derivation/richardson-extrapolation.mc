kill(all);

/*
 * Compute and confirm the coefficients and methodology of repeated Richardson Extrapolation
 * since the second extrapolation was left to be implicitly understood:
 */

/*
 * The series for which we want a32 and a64 --> 0
 */
A0: a0 + r0**32 * a32 + r0**64 * a64 + r0**96 * a96;
A1: a0 + r1**32 * a32 + r1**64 * a64 + r1**96 * a96;
A2: a0 + r2**32 * a32 + r2**64 * a64 + r2**96 * a96;

/*
 * First richardson extrapolation, subject to a32 = 0
 */
eq1: expand(A1 - (A1 - A0) * c1);
c1solve: rhs(solve(coeff(eq1, a32) = 0, c1)[1]);
extrap1: eq1, c1 = c1solve;
extrap1: expand(factor(extrap1));

/*
 * Second richardson extrapolation, subject to a32 = 0
 */
eq2: expand(A2 - (A2 - A1) * c2);
c2solve: rhs(solve(coeff(eq2, a32) = 0, c2)[1]);
extrap2: eq2, c2 = c2solve;
extrap2: expand(factor(extrap2));

/*
 * Third richardson extrapolation, subject to a64 = 0
 */
eq3: expand(extrap2 + (extrap2 - extrap1) * c3);
c3solve: rhs(solve(coeff(eq3, a64) = 0, c3)[1]);

/*
 * The final answer:
 */
ans: extrap2 + (extrap2 - extrap1) * c3solve;
