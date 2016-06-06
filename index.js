'use strict';

module.exports = complexDeriv;

var ndfft = require('ndfft');
var cmod = require('complex-modulus');
var EPSILON = require('number-epsilon');
var EPSILON_3_14 = Math.pow(EPSILON, 3 / 14);

// Premature optimization, probably, but in case this needs to be run many many times,
// avoid allocating arrays. JavaScript is single-threaded, so this is presumably safe:
var circles = new Array(5);
var _buf = [];
var _re = [];
var _im = [];
var _bsr = [[], [], []];
var _bsi = [[], [], []];
var _rs = [];
var refPts = [[-0.4, 0.3], [0.7, 0.2], [0.02, -0.06]];

// Cache points on a unit circle since it's a bunch of sines and cosines that are nice
// not to need to repeat and since there are only five possible unit circles based on
// the number of points m.
function getCircle (rangeNum, m) {
  var i, theta, circle, re, im;
  circles[rangeNum] = circle = new Array(2);
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
  var i, re, im, znr, zni, zr, zi, tmp, dzr, dzi;
  var fmax = 0;
  var buf = _buf;
  var dzmax = 0;

  // For each of the three evaluation points:
  for (i = 0; i < 3; i++) {
    // Compute the first exact value:
    zr = refPts[i][0];
    zi = refPts[i][1];

    // Evaluate the function:
    if (f.length === 2) {
      // If use-supplied function outputs array:
      buf = f(z0r + r * zr, z0i + r * zi);
    } else {
      // if user-supplied function writes to buffer:
      f(buf, z0r + r * zr, z0i + r * zi);
    }

    // Compute the Taylor series approximation. Start with the first term:
    re = br[0];
    im = bi[0];

    // Initialize z^i accumulatorwith z^1:
    znr = zr;
    zni = zi;

    // On each loop iteration, multiply by an additional z, and add term: b[i] * z^i:
    for (i = 1; i < m; i++,
      tmp = znr,
      znr = tmp * zr - zni * zi,
      zni = zni * zr + tmp * zi
    ) {
      re += br[i] * znr - bi[i] * zni;
      im += bi[i] * znr + br[i] * zni;
    }

    // Perform the normalization by m here (i.e. converting b -> a):
    re /= m;
    im /= m;

    // The difference between computed and expected value:
    dzr = buf[0] - re;
    dzi = buf[1] - im;

    // Maximum error:
    dzmax = Math.max(dzmax, cmod(dzr, dzi));

    // Maximum function value:
    fmax = Math.max(fmax, cmod(buf[0], buf[1]));
  }

  // if (dzmax / fmax > 1e-3) {
  //   console.log('failed test =', dzmax / fmax)
  // }

  // Convergence is poor if max error relative to max function value is large:
  return dzmax / fmax > 1e-3;
}

function complexDeriv (out, f, n, a, b, options, status) {
  var m, i, rangeNum, changeCnt, iter, maxIters, radiusFactor, isDegenerate, cp, crat, m1, m2, ef10, ef21, ef20;
  var cr, ci, circle, needsSmaller, prevDir, e1r, e1i, e2r, e2i, k, bcr, bci, taylor, factorial, r;

  // Local references, continuing with potential premature optimization (maybe a couple ms
  // faster than reallocating every time, but the scatter in the benchmark is a little high
  // to be 100% certain). At any rate, this hopefully avoids some garbage collection if the
  // derivated is computed *repeatedly*.
  var buf = _buf;
  var re = _re;
  var im = _im;
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
  if (m < re.length) {
    _re = re = new Array(m);
    _im = im = new Array(m);
  }

  // Get an interleaved list of re/im pairs on the unit circle:
  circle = circles[rangeNum];

  // If not already cached for this range number, then compute:
  if (!(circle = circles[rangeNum])) {
    circle = getCircle(rangeNum, m);
  }

  // Unpack the components, again premature opt...
  cr = circle[0];
  ci = circle[1];

  // Initial grow/shring factor for the circle:
  radiusFactor = 2;

  // A factor for testing against the targeted geometric progression of fourier coefficients:
  crat = 1 / Math.exp(Math.log(1e-4) / (m - 1));

  // Number of direction changes:
  changeCnt = 0;

  // Degenerate if either the first half or the last half of the coefficients are essentially
  // zero compared to the other. This means we can't adapt the error:
  isDegenerate = false;

  // Number of accumulated samples for Richardson extrapolation:
  k = 0;

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
      re[i] = buf[0];
      im[i] = buf[1];
    }

    // Compute the fourier transform. This is the main step in converting evaluations around
    // the circle into coefficients b:
    ndfft(1, re, im);

    // If we've changed twice or it's degenerate, then start storing:
    if (changeCnt > 1 || isDegenerate) {
      // Store coefficients. These become our output once extrapolated.
      for (i = 0, cp = 1 / m; i <= n; i++, cp /= r) {
        bsr[k][i] = re[i] * cp;
        bsi[k][i] = im[i] * cp;
      }

      // Store the current radius. This is the independent variable for richardson extrap.
      rs[k] = r;

      // Increment the storage counter:
      k++;

      // If we've accumulated three, then break;
      if (k === 3) break;
    }

    // If not degenerate, check for geometric progression in the fourier transform:
    if (!isDegenerate) {
      // Check the magnitude of the first half of the range:
      for (m1 = 0, cp = 1, i = 0; i < m / 2; i++, cp *= crat) {
        bcr = re[i] * cp;
        bci = im[i] * cp;
        m1 = Math.max(cmod(bcr, bci), m1);
      }
      // Check the magnitude of the second half of the range:
      for (m2 = 0; i < m; i++, cp *= crat) {
        bcr = re[i] * cp;
        bci = im[i] * cp;
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
      needsSmaller = m1 < m2 || convergenceIsPoor(a, b, r, f, m, re, im);
    }

    // If beyond the first iteration and direction changed, then increment counter:
    if (iter && needsSmaller !== prevDir) {
      changeCnt++;
    }

    // If changes have occured, start taking the square root of the factor. The effect is
    // that we first double/half the size, then once we're in a reasonable range start
    // changing it by a smaller factor that yields good numbers for Richardson extrapolation:
    if (changeCnt > 0) {
      radiusFactor = Math.sqrt(radiusFactor);
    }

    // If needs smaller, divide by radiusFactor, otherwise mult:
    if (needsSmaller) {
      r /= radiusFactor;
    } else {
      r *= radiusFactor;
    }

    // Store the previous direction so we can find out when it changes
    prevDir = needsSmaller;
  }

  // Begin Richardson Extrapolation. Presumably we have b[i]'s around three successive circles
  // and can now extrapolate those coefficients, zeroing out higher order error terms.
  //
  // Extrapolation factors for Richardson Extrapolation:
  ef10 = 1 / (1 - Math.pow(rs[0] / rs[1], m));
  ef21 = 1 / (1 - Math.pow(rs[1] / rs[2], m));
  ef20 = 1 / (1 - Math.pow(rs[0] / rs[2], m));

  // If we need to outpust status, allocate or get per-term error storage:
  if (status) {
    status.truncationError = status.truncationError || [];
    status.roundingError = status.roundingError || [];
  }

  // Extrapolate derivative terms, one at a time:
  //   Note: only output the first n as requested--which is usually somewhat fewer than the
  //   number of terms we've actually calculated at this point, but is what's probably expected
  //   by the user.
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
      // Truncation error from last richardson update and with ϵ^1/14 safety factor:
      status.truncationError[i] = EPSILON_3_14 * cmod(rr, ri);

      // Rounding error, from previously-computed factors max_i |b_i| / c_i
      status.roundingError[i] = EPSILON / Math.pow(rs[2], i) * Math.max(m1, m2);

      // If not taylor, then multiply the error terms by the factorial also:
      if (!taylor) {
        status.truncationError[i] *= factorial;
        status.roundingError[i] *= factorial;
      }
    }
  }

  // If outputting status, write some extra info
  if (status) {
    status.degenerate = isDegenerate;
    status.iterations = iter;
    status.finalRadius = r;
    status.failed = iter === maxIters;
  }

  // Return the final extrapolated result (which may also have been written to `out` as an
  // initial input, but nice to return anyway):
  return out;
}
