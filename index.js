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

function getCircle (range, m) {
  var i, theta, circle, re, im;
  circles[range] = circle = new Array(2);
  re = circle[0] = new Array(m);
  im = circle[1] = new Array(m);
  for (i = m - 1; i >= 0; i--) {
    theta = 2 * Math.PI * i / m;
    re[i] = Math.cos(theta);
    im[i] = Math.sin(theta);
  }

  return circle;
}

function convergenceIsPoor (z0r, z0i, r, f, m, br, bi) {
  var i, re, im, znr, zni, zr, zi, tmp, dzr, dzi;
  var fmax = 0;
  var buf = _buf;
  var dzmax = 0;

  for (i = 0; i < 3; i++) {
    // Compute the first exact value:
    zr = refPts[i][0];
    zi = refPts[i][1];

    // Evaluate the function:
    if (f.length === 2) {
      buf = f(z0r + r * zr, z0i + r * zi);
    } else {
      f(buf, z0r + r * zr, z0i + r * zi);
    }

    // Compute the Laurent series approximation:
    re = br[0];
    im = bi[0];

    // Initialize z^n:
    znr = zr;
    zni = zi;

    // On each loop iteration, multiply by an additional z:
    for (i = 1; i < m; i++,
      tmp = znr,
      znr = tmp * zr - zni * zi,
      zni = zni * zr + tmp * zi
    ) {
      re += br[i] * znr - bi[i] * zni;
      im += bi[i] * znr + br[i] * zni;
    }

    // Perform the normalization by m here:
    re /= m;
    im /= m;

    dzr = buf[0] - re;
    dzi = buf[1] - im;

    dzmax = Math.max(dzmax, cmod(dzr, dzi));
    fmax = Math.max(fmax, cmod(buf[0], buf[1]));
  }

  // if (dzmax / fmax > 1e-3) {
  //   console.log('failed test =', dzmax / fmax)
  // }

  return dzmax / fmax > 1e-3;
}

function complexDeriv (out, f, n, a, b, options, status) {
  var m, i, range, changes, iter, maxIters, fac, degen, cp, crat, m1, m2, ef10, ef21, ef20;
  var cr, ci, circle, sm, psm, e1r, e1i, e2r, e2i, k, bcr, bci, taylor, factorial, r;
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

  options = options || {};
  r = options.r === undefined ? 0.6580924658 : options.r;
  maxIters = options.maxIterations === undefined ? 30 : options.maxIterations;
  taylor = options.taylor === undefined ? false : !!options.taylor;

  // Select the number of points to use based on Fornberg's heuristics:
  if (n <= 6) {
    m = 8;
    range = 0;
  } else if (n <= 12) {
    m = 16;
    range = 1;
  } else if (n <= 25) {
    m = 32;
    range = 2;
  } else if (n <= 51) {
    m = 64;
    range = 3;
  } else if (n <= 100) {
    m = 128;
    range = 4;
  } else {
    throw new Error('complex-deriv-fornberg: Number of derivatives requested (' + n + ') is greater than upper limit of 100');
  }

  // ndfft reads the length, so if it's leftover from a larger run, this will cause problems
  // otherwise can avoid reallocating tons of little arrays:
  if (m < re.length) {
    _re = re = new Array(m);
    _im = im = new Array(m);
  }

  // Get an interleaved list of re/im pairs on the unit circle:
  circle = circles[range];
  if (!(circle = circles[range])) {
    circle = getCircle(range, m);
  }
  cr = circle[0];
  ci = circle[1];

  // Initial enlargement factor:
  fac = 2;

  // A factor for testing against the targeted geometric progression of fourier coefficients:
  crat = 1 / Math.exp(Math.log(1e-4) / (m - 1));

  // Number of direction chnages:
  changes = 0;

  // Degenerate if either the first half or the last half of the coefficients are essentially
  // zero compared to the other. This means we can't adapt the error:
  degen = false;

  // Number of accumulated samples for Richardson extrapolation:
  k = 0;

  for (iter = 0; iter < maxIters; iter++) {
    // Evaluate the function on this unit circle:
    for (i = m - 1; i >= 0; i--) {
      if (f.length === 2) {
        buf = f(a + r * cr[i], b + r * ci[i]);
      } else {
        f(buf, a + r * cr[i], b + r * ci[i]);
      }
      re[i] = buf[0];
      im[i] = buf[1];
    }

    // Compute the fourier transform:
    ndfft(1, re, im);

    // If we've changed twice or it's degenerate, then start storing:
    if (changes > 1 || degen) {
      for (i = 0, cp = 1 / m; i <= n; i++, cp /= r) {
        bsr[k][i] = re[i] * cp;
        bsi[k][i] = im[i] * cp;
      }

      rs[k] = r;

      k++;

      // If we've accumulated three, then break;
      if (k === 3) break;
    }

    // If not degenerate, check for geometric progression
    // in the fourier transform:
    if (!degen) {
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
      degen = m1 / m2 < 1e-8 || m2 / m1 < 1e-8;
    }

    if (degen) {
      // If degenerate, then simply ignore error management and alternate bigger/smaller
      // so that we can perform richardson extrapolation:
      sm = iter % 2 === 0;
    } else {
      // Otherwise check the progression and if false then perform a convergence check:
      sm = m1 < m2 || convergenceIsPoor(a, b, r, f, m, re, im);
    }

    // If beyond the first iteration and direction changed, then increment counter:
    if (iter && sm !== psm) {
      changes++;
    }

    // If changes have occured, start taking the square root of the factor:
    if (changes > 0) {
      fac = Math.sqrt(fac);
    }

    // If needs smaller, divide by fac, otherwise mult:
    if (sm) {
      r /= fac;
    } else {
      r *= fac;
    }

    psm = sm;
  }

  // Extrapolation factors for Richardson Extrapolation:
  ef10 = 1 / (1 - Math.pow(rs[0] / rs[1], m));
  ef21 = 1 / (1 - Math.pow(rs[1] / rs[2], m));
  ef20 = 1 / (1 - Math.pow(rs[0] / rs[2], m));

  if (status) {
    status.truncationError = status.truncationError || [];
    status.roundingError = status.roundingError || [];
  }

  // Note: only output the first n as requested:
  for (i = 0, factorial = 1; i <= n; i++, factorial *= i) {
    // The first Richardson extrapolation:
    e1r = bsr[1][i] - (bsr[1][i] - bsr[0][i]) * ef10;
    e1i = bsi[1][i] - (bsi[1][i] - bsi[0][i]) * ef10;

    // Second Richardson extrapolation:
    e2r = bsr[2][i] - (bsr[2][i] - bsr[1][i]) * ef21;
    e2i = bsi[2][i] - (bsi[2][i] - bsi[1][i]) * ef21;

    // Combined Richardson extrapolation:
    // Store the components for computation of truncation error:
    var rr = (e2r - e1r) * ef20;
    var ri = (e2i - e1i) * ef20;
    // Update:
    out[0][i] = e2r - rr;
    out[1][i] = e2i - ri;

    // By default, unless *not* the taylor series, then multiply by i! to compute the derivative:
    if (!taylor) {
      out[0][i] *= factorial;
      out[1][i] *= factorial;
    }

    if (status) {
      // Truncation error from last richardson update and with ϵ^1/14 safety facto:
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

  if (status) {
    status.degenerate = degen;
    status.iterations = iter;
    status.finalRadius = r;
    status.failed = iter === maxIters;
  }

  // Extrapolated result:
  return out;
}
