'use strict';

module.exports = complexDeriv;

var ndfft = require('ndfft');
var cmod = require('complex-modulus');
var EPSILON = require('number-epsilon');
var EPSILON_3_14 = Math.pow(EPSILON, 3 / 14);

// Premature optimization, probably, but in case this needs to be run many many times,
// avoid allocating arrays. JavaScript is single-threaded, so this is presumably safe:
var circleCache = new Array(5);
var _buf = [];
var _re = [];
var _im = [];
var _bsr = [[], [], []];
var _bsi = [[], [], []];
var _rs = [];
var convergenceCheckPoints = [[-0.4, 0.3], [0.7, 0.2], [0.02, -0.06]];

// Cache points on a unit circle since it's a bunch of sines and cosines that are nice
// not to need to repeat and since there are only five possible unit circles based on
// the number of points m.
function getCachedCirclePoints (rangeNum, m) {
  var i, theta, circle, re, im;
  circleCache[rangeNum] = circle = new Array(2);
  re = circle[0] = new Array(m);
  im = circle[1] = new Array(m);
  for (i = m - 1; i >= 0; i--) {
    theta = 2 * Math.PI * i / m;
    re[i] = Math.cos(theta);
    im[i] = Math.sin(theta);
  }

  return circle;
}

// Test for poor convergence based on three function evaluations. To avoid letting randomness
// enter the algorithm, I've used three fixed points, as defined above (guaranteed to be random
// every time...) This test evaluates the function at the three points and returns false if
// the relative error is greater than 1e-3.
function convergenceIsPoor (z0r, z0i, r, f, m, br, bi) {
  var i, fr, fi, znr, zni, zr, zi, tmp, errr, erri;
  var fEvalMax = 0;
  var fEval = _buf;
  var absErrMax = 0;

  for (i = 0; i < convergenceCheckPoints.length; i++) {
    zr = convergenceCheckPoints[i][0];
    zi = convergenceCheckPoints[i][1];

    if (f.length === 2) {
      // If use-supplied function outputs array:
      fEval = f(z0r + r * zr, z0i + r * zi);
    } else {
      // if user-supplied function writes to buffer:
      f(fEval, z0r + r * zr, z0i + r * zi);
    }

    // Compute the Taylor series approximation. Start with the first term:
    fr = br[0];
    fi = bi[0];

    // Initialize z^i accumulator with z:
    znr = zr;
    zni = zi;

    // On each loop iteration, multiply by an additional z, and add term: b[i] * z^i:
    for (i = 1; i < m; i++,
      tmp = znr,
      znr = tmp * zr - zni * zi,
      zni = zni * zr + tmp * zi
    ) {
      fr += br[i] * znr - bi[i] * zni;
      fi += bi[i] * znr + br[i] * zni;
    }

    // b (output of fft) is unnormalized, so divide by the number of points here:
    fr /= m;
    fi /= m;

    errr = fEval[0] - fr;
    erri = fEval[1] - fi;

    absErrMax = Math.max(absErrMax, cmod(errr, erri));
    fEvalMax = Math.max(fEvalMax, cmod(fEval[0], fEval[1]));
  }

  return absErrMax / fEvalMax > 1e-3;
}

function complexDeriv (out, f, n, a, b, options, status) {
  var m, i, rangeNum, dirChangeCnt, iter, maxIters, radiusGrowthFactor, isDegenerate, cp, crat, m1, m2, ef10, ef21, ef20;
  var cr, ci, circle, needsSmaller, prevNeedsSmaller, e1r, e1i, e2r, e2i, numSamples, bcr, bci, taylor, factorial, r;

  // Local references, continuing with potential premature optimization (maybe a couple ms
  // faster than reallocating every time, but the scatter in the benchmark is a little high
  // to be 100% certain). At any rate, this hopefully avoids some garbage collection if the
  // derivated is computed *repeatedly*.
  var buf = _buf;
  var fEvalr = _re;
  var fEvali = _im;
  var bsr = _bsr;
  var bsi = _bsi;
  var rs = _rs;

  // Handle variadic form in case output array is not provided:
  if (typeof out === 'function') {
    status = options;
    options = b;
    b = a;
    a = n;
    n = f;
    f = out;
    out = [[], []];
  } else {
    if (!Array.isArray(out[0])) {
      out[0] = [];
    }
    if (!Array.isArray(out[1])) {
      out[1] = [];
    }
  }

  // Get options and set defaults:
  options = options || {};
  r = options.r === undefined ? 0.6580924658 : options.r;
  maxIters = options.maxIterations === undefined ? 30 : options.maxIterations;
  taylor = options.taylor === undefined ? false : !!options.taylor;

  // Select the number of points to use based on Fornberg's heuristics:
  if (n <= 6) {
    m = 8;
    rangeNum = 0;
  } else if (n <= 12) {
    m = 16;
    rangeNum = 1;
  } else if (n <= 25) {
    m = 32;
    rangeNum = 2;
  } else if (n <= 51) {
    m = 64;
    rangeNum = 3;
  } else if (n <= 100) {
    m = 128;
    rangeNum = 4;
  } else {
    throw new Error('complex-deriv-fornberg: Number of derivatives requested (' + n + ') is greater than upper limit of 100');
  }

  // Create new re/im arrays if the size has shrunk. This is needed because ndfft reads the
  // length of the transform from the input, so if it's leftover from a larger run, re/im
  // will cause problems. Otherwise can avoid reallocating tons of little arrays by reusing.
  if (m < fEvalr.length) {
    _re = fEvalr = new Array(m);
    _im = fEvali = new Array(m);
  }

  // Get an interleaved list of re/im pairs on the unit circle:
  circle = circleCache[rangeNum];

  // If not already cached for this range number, then compute:
  if (!(circle = circleCache[rangeNum])) {
    circle = getCachedCirclePoints(rangeNum, m);
  }

  // Unpack the components, again premature opt...
  cr = circle[0];
  ci = circle[1];

  // Initial grow/shring factor for the circle:
  radiusGrowthFactor = 2;

  // A factor for testing against the targeted geometric progression of fourier coefficients:
  crat = 1 / Math.exp(Math.log(1e-4) / (m - 1));

  // Number of direction changes:
  dirChangeCnt = 0;

  // Degenerate if either the first half or the last half of the coefficients are essentially
  // zero compared to the other. This means we can't adapt the error:
  isDegenerate = false;

  // Number of accumulated samples for Richardson extrapolation:
  numSamples = 0;

  // Start iterating. The goal of this loops is to select a circle radius that yields a nice
  // geometric progression of the coefficients (which controls the error), and then to accumulate
  // *three* successive approximations as a function of the circle radius r so that we can
  // perform Richardson Extrapolation and zero out error terms, *grealy* improving the quality
  // of the approximation.
  for (iter = 0; iter < maxIters; iter++) {
    // Evaluate the function on this unit circle:
    for (i = m - 1; i >= 0; i--) {
      if (f.length === 2) {
        // If use-supplied function outputs array:
        buf = f(a + r * cr[i], b + r * ci[i]);
      } else {
        // if user-supplied function writes to buffer:
        f(buf, a + r * cr[i], b + r * ci[i]);
      }
      fEvalr[i] = buf[0];
      fEvali[i] = buf[1];
    }

    // Compute the fourier transform. This is the main step in converting evaluations around
    // the circle into coefficients b:
    ndfft(1, fEvalr, fEvali);

    // If we've changed twice or it's degenerate, then start storing:
    if (dirChangeCnt > 1 || isDegenerate) {
      // Store coefficients. These become our output once extrapolated.
      for (i = 0, cp = 1 / m; i <= n; i++, cp /= r) {
        bsr[numSamples][i] = fEvalr[i] * cp;
        bsi[numSamples][i] = fEvali[i] * cp;
      }

      // Store the current radius. This is the independent variable for richardson extrap.
      rs[numSamples] = r;

      numSamples++;

      if (numSamples === 3) break;
    }

    // If not degenerate, check for geometric progression in the fourier transform:
    if (!isDegenerate) {
      // Check the magnitude of the first half of the range:
      for (m1 = 0, cp = 1, i = 0; i < m / 2; i++, cp *= crat) {
        bcr = fEvalr[i] * cp;
        bci = fEvali[i] * cp;
        m1 = Math.max(cmod(bcr, bci), m1);
      }
      // Check the magnitude of the second half of the range:
      for (m2 = 0; i < m; i++, cp *= crat) {
        bcr = fEvalr[i] * cp;
        bci = fEvali[i] * cp;
        m2 = Math.max(cmod(bcr, bci), m2);
      }

      // If there's an extreme mismatch, then we can consider the geometric progression
      // degenerate, whether one way or the other, and just alternate directions instead
      // of trying to target a specific error bound (not ideal, but not a good reason to
      // fail catastrophically):
      isDegenerate = m1 / m2 < 1e-8 || m2 / m1 < 1e-8;
    }

    if (isDegenerate) {
      // If degenerate, then simply ignore error management and alternate bigger/smaller
      // so that we can perform richardson extrapolation:
      needsSmaller = iter % 2 === 0;
    } else {
      // Otherwise check the progression and if false then perform a convergence check:
      needsSmaller = m1 < m2 || convergenceIsPoor(a, b, r, f, m, fEvalr, fEvali);
    }

    if (iter > 0 && needsSmaller !== prevNeedsSmaller) {
      dirChangeCnt++;
    }

    // Once we've started changing directions, we've found our range so start taking the
    // square root of the growth factor so that richardson extrapolation is well-behaved:
    if (dirChangeCnt > 0) {
      radiusGrowthFactor = Math.sqrt(radiusGrowthFactor);
    }

    if (needsSmaller) {
      r /= radiusGrowthFactor;
    } else {
      r *= radiusGrowthFactor;
    }

    prevNeedsSmaller = needsSmaller;
  }

  // Begin Richardson Extrapolation. Presumably we have b[i]'s around three successive circles
  // and can now extrapolate those coefficients, zeroing out higher order error terms.
  //
  // Extrapolation factors for Richardson Extrapolation:
  ef10 = 1 / (1 - Math.pow(rs[0] / rs[1], m));
  ef21 = 1 / (1 - Math.pow(rs[1] / rs[2], m));
  ef20 = 1 / (1 - Math.pow(rs[0] / rs[2], m));

  if (status) {
    status.truncationError = status.truncationError || [];
    status.roundingError = status.roundingError || [];
  }

  // Extrapolate derivative terms, one at a time:
  //   Note: only output the first n as requested--which is usually somewhat fewer than the
  //   number of terms we've actually calculated at this point, but creates the smallest
  //   amount of confusion.
  for (i = 0, factorial = 1; i <= n; i++, factorial *= i) {
    // The first Richardson extrapolation:
    e1r = bsr[1][i] - (bsr[1][i] - bsr[0][i]) * ef10;
    e1i = bsi[1][i] - (bsi[1][i] - bsi[0][i]) * ef10;

    // Second Richardson extrapolation:
    e2r = bsr[2][i] - (bsr[2][i] - bsr[1][i]) * ef21;
    e2i = bsi[2][i] - (bsi[2][i] - bsi[1][i]) * ef21;

    // Combined Richardson extrapolation (see derivation/richardson.mc for maxima code that
    // verifies zeroing of error terms):
    //
    // Store the third update so that we can compute the truncation error, if needed
    var rr = (e2r - e1r) * ef20;
    var ri = (e2i - e1i) * ef20;

    // Third and final richardson extrapolation Update:
    out[0][i] = e2r - rr;
    out[1][i] = e2i - ri;

    // If this is *not* a taylor series approximation (this is the default behavior), then
    // multiply by the factorial term to convert taylor series --> actual derivatives.
    if (!taylor) {
      out[0][i] *= factorial;
      out[1][i] *= factorial;
    }

    // If status to be output, then compute the truncation error:
    if (status) {
      status.truncationError[i] = EPSILON_3_14 * cmod(rr, ri);

      status.roundingError[i] = EPSILON / Math.pow(rs[2], i) * Math.max(m1, m2);

      // Ditto on taylor --> derivative conversion via factorial:
      if (!taylor) {
        status.truncationError[i] *= factorial;
        status.roundingError[i] *= factorial;
      }
    }
  }

  if (status) {
    status.degenerate = isDegenerate;
    status.iterations = iter;
    status.finalRadius = r;
    status.failed = iter === maxIters;
  }

  // Output is also written in-place if storage was provided.
  return out;
}
