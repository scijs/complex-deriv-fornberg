'use strict';

var assert = require('chai').assert;
var deriv = require('../');
assert.almostEqual = require('./util/almost-equal');
var almostEqual = require('almost-equal');
var Complex = require('complex.js');
var cmod = require('complex-modulus');

var status = {};

describe('complex-deriv-fornberg', function () {
  describe('input format', function () {
    it('accepts a function that writes re/im in-place', function () {
      var ans = deriv(function (out, a, b) {
        var c = b * b + (1 - a) * (1 - a);
        out[0] = (1 - a) / c;
        out[1] = b / c;
      }, 5, 0, 0, {}, status);

      assert.almostEqual(ans, [[1, 1, 2, 6, 24, 120], [0, 0, 0, 0, 0, 0]], status.truncationError, status.roundingError);
    });

    it('accepts a function that returns re/im in an array', function () {
      var ans = deriv(function (a, b) {
        var c = b * b + (1 - a) * (1 - a);
        return [(1 - a) / c, b / c];
      }, 5, 0, 0, {}, status);

      assert.almostEqual(ans, [[1, 1, 2, 6, 24, 120], [0, 0, 0, 0, 0, 0]], status.truncationError, status.roundingError);
    });
  });

  describe('sin(x / 1e-6)', function () {
    var c = 1e-3;

    function f (out, a, b) {
      var s = Complex(a / c, b / c).sin();
      out[0] = s.re;
      out[1] = s.im;
      return out;
    }

    function fp (out, a, b) {
      var s = Complex(a / c, b / c).cos().div(c);
      out[0] = s.re;
      out[1] = s.im;
      return out;
    }

    function fpp (out, a, b) {
      var s = Complex(a / c, b / c).sin().div(c * c);
      out[0] = -s.re;
      out[1] = -s.im;
      return out;
    }
    function fppp (out, a, b) {
      var s = Complex(a / c, b / c).cos().div(c * c * c);
      out[0] = -s.re;
      out[1] = -s.im;
      return out;
    }

    function analytical (n, a, b) {
      var y0 = f([], a, b);
      var y1 = fp([], a, b);
      var y2 = fpp([], a, b);
      var y3 = fppp([], a, b);
      var ans = [[y0[0], y1[0], y2[0], y3[0]], [y0[1], y1[1], y2[1], y3[1]]];
      return ans;
    }

    // Spiral outwards and try a *lot* of different values:
    for (var rr = 1e-10, aarg = 0; rr < 0.1; rr *= 1.3, aarg++) {
      (function (r, arg) {
        var re = r * Math.cos(arg);
        var im = r * Math.sin(arg);
        it('at z = ' + re + ' + ' + im + ' * i', function () {
          var ans = deriv(f, 3, re, im, {r: 2e-3}, status);
          assert.isFalse(status.degenerate, 'is not degenerate');
          assert.almostEqual(ans, analytical(3, re, im), status.truncationError, status.roundingError, 11);
        });
      }(rr, aarg));
    }
  });

  describe('sin(x / 1e6)', function () {
    var c = 1e6;

    function f (out, a, b) {
      var s = Complex(a / c, b / c).sin();
      out[0] = s.re;
      out[1] = s.im;
      return out;
    }

    function fp (out, a, b) {
      var s = Complex(a / c, b / c).cos().div(c);
      out[0] = s.re;
      out[1] = s.im;
      return out;
    }

    function fpp (out, a, b) {
      var s = Complex(a / c, b / c).sin().div(c * c);
      out[0] = -s.re;
      out[1] = -s.im;
      return out;
    }
    function fppp (out, a, b) {
      var s = Complex(a / c, b / c).cos().div(c * c * c);
      out[0] = -s.re;
      out[1] = -s.im;
      return out;
    }

    function analytical (n, a, b) {
      var y0 = f([], a, b);
      var y1 = fp([], a, b);
      var y2 = fpp([], a, b);
      var y3 = fppp([], a, b);
      var ans = [[y0[0], y1[0], y2[0], y3[0]], [y0[1], y1[1], y2[1], y3[1]]];
      return ans;
    }

    // Spiral outwards and try a *lot* of different values:
    for (var rr = 1e-6, aarg = 0; rr < 1e7; rr *= 1.1, aarg++) {
      (function (r, arg) {
        var re = r * Math.cos(arg);
        var im = r * Math.sin(arg);
        it('at z = ' + re + ' + ' + im + ' * i', function () {
          var ans = deriv(f, 3, re, im, {}, status);
          assert.almostEqual(ans, analytical(3, re, im), status.truncationError, status.roundingError);
        });
      }(rr, aarg));
    }
  });

  describe('sin(x)', function () {
    function f (out, a, b) {
      var s = Complex(a, b).sin();
      out[0] = s.re;
      out[1] = s.im;
      return out;
    }

    function fp (out, a, b) {
      var s = Complex(a, b).cos();
      out[0] = s.re;
      out[1] = s.im;
      return out;
    }

    function fpp (out, a, b) {
      var s = Complex(a, b).sin();
      out[0] = -s.re;
      out[1] = -s.im;
      return out;
    }
    function fppp (out, a, b) {
      var s = Complex(a, b).cos();
      out[0] = -s.re;
      out[1] = -s.im;
      return out;
    }

    function analytical (n, a, b) {
      var i;
      var y0 = f([], a, b);
      var y1 = fp([], a, b);
      var y2 = fpp([], a, b);
      var y3 = fppp([], a, b);
      var ans = [[y0[0], y1[0], y2[0], y3[0]], [y0[1], y1[1], y2[1], y3[1]]];
      for (i = 4; i <= n; i++) {
        ans[0][i] = ans[0][i % 4];
        ans[1][i] = ans[1][i % 4];
      }
      return ans;
    }

    it('at z = 0', function () {
      var ans = deriv(f, 8, 0, 0, {}, status);
      assert.isFalse(status.degenerate, 'Is not degenerate');
      assert.almostEqual(ans, analytical(8, 0, 0), status.truncationError, status.roundingError);
    });

    // Spiral outwards and try a *lot* of different values:
    for (var rr = 1e-6, aarg = 0; rr < 1e2; rr *= 1.1, aarg++) {
      (function (r, arg) {
        var re = r * Math.cos(arg);
        var im = r * Math.sin(arg);
        it('at z = ' + re + ' + ' + im + ' * i', function () {
          var ans = deriv(f, 8, re, im, {}, status);
          assert.almostEqual(ans, analytical(8, re, im), status.truncationError, status.roundingError);
        });
      }(rr, aarg));
    }
  });

  describe('e^z', function () {
    function f (out, a, b) {
      var ex = Math.exp(a);
      out[0] = ex * Math.cos(b);
      out[1] = ex * Math.sin(b);
      return out;
    }

    function analytical (n, a, b) {
      var i;
      var ans = [[], []];
      for (i = 0; i <= n; i++) {
        var fp = f([], a, b);
        ans[0][i] = fp[0];
        ans[1][i] = fp[1];
      }
      return ans;
    }

    it('at z = 2 + i', function () {
      var ans = deriv(f, 1, 2, 1, {r: 1}, status);
      assert.isFalse(status.degenerate, 'Is not degenerate');
      assert.almostEqual(ans, analytical(1, 2, 1), status.truncationError, status.roundingError, 6);
    });

    // Spiral outwards and try a *lot* of different values:
    for (var rr = 1e-6, aarg = 0; rr < 1e2; rr *= 1.1, aarg++) {
      (function (r, arg) {
        var re = r * Math.cos(arg);
        var im = r * Math.sin(arg);
        it('at z = ' + re + ' + ' + im + ' * i', function () {
          var ans = deriv(f, 8, re, im, {}, status);
          assert.almostEqual(ans, analytical(8, re, im), status.truncationError, status.roundingError);
        });
      }(rr, aarg));
    }
  });

  describe('1 + z + z^2', function () {
    function f (out, a, b) {
      out[0] = 1 + a + a * a - b * b;
      out[1] = b * (2 * a + 1);
      return out;
    }

    function fp (out, a, b) {
      out[0] = 1 + 2 * a;
      out[1] = 2 * b;
      return out;
    }

    function fpp (out, a, b) {
      out[0] = 2;
      out[1] = 0;
      return out;
    }

    function analytical (n, a, b) {
      var i;
      var y0 = f([], a, b);
      var y1 = fp([], a, b);
      var y2 = fpp([], a, b);
      var ans = [[y0[0], y1[0], y2[0]], [y0[1], y1[1], y2[1]]];
      for (i = 3; i <= n; i++) {
        ans[0][i] = 0;
        ans[1][i] = 0;
      }
      return ans;
    }

    it('is degenerate', function () {
      deriv(f, 7, 2, 1, {}, status);
      assert.isTrue(status.degenerate, 'Is not degenerate');
    });

    it('at z = 2 + i', function () {
      var ans = deriv(f, 7, 2, 1, null, status);
      assert.almostEqual(ans, analytical(7, 2, 1), status.truncationError, status.roundingError);
    });
  });

  describe('1 / (1 - z):', function () {
    function f (out, a, b) {
      var a1 = 1 - a;
      var c = b * b + a1 * a1;
      out[0] = a1 / c;
      out[1] = b / c;
      return out;
    }

    function fp (out, a, b) {
      var c = Math.pow(b * b + a * a - 2 * a + 1, 2);
      out[0] = -(b - a + 1) * (b + a - 1) / c;
      out[1] = 2 * (1 - a) * b / c;
      return out;
    }

    function fpp (out, a, b) {
      var c = Math.pow(b * b + Math.pow(1 - a, 2), 3);
      out[0] = 2 * (Math.pow(1 - a, 3) - 3 * (1 - a) * b * b) / c;
      out[1] = -2 * (b * b * b - 3 * Math.pow(1 - a, 2) * b) / c;
      return out;
    }

    function analytical (a, b) {
      var y0 = f([], a, b);
      var y1 = fp([], a, b);
      var y2 = fpp([], a, b);
      return [[y0[0], y1[0], y2[0]], [y0[1], y1[1], y2[1]]];
    }

    it('outputs optional status', function () {
      assert.almostEqual(deriv(f, 2, 0, 0, {}, status), analytical(0, 0), status.truncationError, status.roundingError);
      assert.isFalse(status.degenerate, 'Is not degenerate');
      assert.equal(status.iterations, 7);
      assert(Math.abs(status.finalRadius - 0.25) < 0.125);
    });

    it('first 50 terms of the taylor expansion at z = 0', function () {
      var ans = deriv(f, 50, 0, 0, {taylor: true}, status);
      for (var i = 0; i < 50; i++) {
        assert(almostEqual(ans[0][i], 1, status.truncationError[i]));
        assert(almostEqual(ans[1][i], 0, status.truncationError[i]));
      }
    });

    it('computes rounding error at z = 0', function () {
      var computed = deriv(f, 8, 0, 0, {}, status);
      var actual = analytical(0, 0);

      for (var i = 0; i < 3; i++) {
        var actualError = cmod(computed[0][i] - actual[0][i], computed[1][i] - actual[1][i]);

        // Should be less than the truncation error (with safety factor of ϵ^1/14):
        assert(actualError < status.truncationError[i]);

        // Should be larger than the rounding error:
        assert(actualError > status.roundingError[i]);
      }
    });

    it('at z = 0', function () {
      assert.almostEqual(deriv(f, 2, 0, 0, null, status), analytical(0, 0), status.truncationError, status.roundingError);
    });

    it('at z = i', function () {
      assert.almostEqual(deriv(f, 2, 0, 1, null, status), analytical(0, 1), status.truncationError, status.roundingError);
    });

    it('at z = -i', function () {
      assert.almostEqual(deriv(f, 2, 0, -1, null, status), analytical(0, -1), status.truncationError, status.roundingError);
    });

    it('at z = 1000', function () {
      assert.almostEqual(deriv(f, 2, 1000, 0, null, status), analytical(1000, 0), status.truncationError, status.roundingError);
    });

    it('at z = 10000', function () {
      assert.almostEqual(deriv(f, 2, 10000, 0, {r: 10000}, status), analytical(10000, 0), status.truncationError, status.roundingError);
    });

    it('at z = 10000 + 10000j', function () {
      assert.almostEqual(deriv(f, 2, 10000, 10000, {r: 10000}, status), analytical(10000, 10000), status.truncationError, status.roundingError);
    });

    it('at z = 1 + i', function () {
      assert.almostEqual(deriv(f, 2, 1, 1, null, status), analytical(1, 1), status.truncationError, status.roundingError);
    });

    // Spiral outwards and try a *lot* of different values:
    for (var rr = 1e-6, aarg = 0; rr < 1e2; rr *= 1.1, aarg++) {
      (function (r, arg) {
        var re = r * Math.cos(arg);
        var im = r * Math.sin(arg);
        it('at z = ' + re + ' + ' + im + ' * i', function () {
          assert.almostEqual(deriv(f, 2, re, im, {}, status), analytical(re, im), status.truncationError, status.roundingError, 50);
          assert.isFalse(status.degenerate, 'Is not degenerate');
        });
      }(rr, aarg));
    }

    it('performs well', function () {
      var i;
      var n = 10000;
      for (i = 0; i < n; i++) {
        deriv(f, 2, 0, 0);
      }
      var t1 = Date.now();

      for (i = 0; i < n; i++) {
        deriv(f, 3, 0, 0);
      }

      var t2 = Date.now();
      var perEval = (t2 - t1) / n;

      console.log('Total time for ' + n + ' evaluations:', (t2 - t1) + ' ms');
      console.log('Avg time per evaluation', 1000 * perEval + ' µs');

      // No assertion here; just for reporting
    });
  });

  describe('log(z)', function () {
    var f = function (out, a, b) {
      var fz = Complex(a, b).log();
      out[0] = fz.re;
      out[1] = fz.im;
      return out;
    };

    function analytical (n, a, b) {
      var i;
      var ans = [[], []];
      var fz = f([], a, b);
      ans[0][0] = fz[0];
      ans[1][0] = fz[1];
      var fact = 1;
      for (i = 1; i <= n; i++, fact *= Math.max(1, i - 1)) {
        var fp = Complex(a, b).pow(-i).mul(fact * Math.pow(-1, i + 1));
        ans[0][i] = fp.re;
        ans[1][i] = fp.im;
      }
      return ans;
    }

    /* it('fails along the branch cut', function () {
      var ans = deriv(f, 20, -1, 1e-2, {r: 1e-2}, status);
      assert.isTrue(status.degenerate, 'Is degenerate');
      assert.almostEqual(ans, analytical(20, 1, 1), status.truncationError, status.roundingError, 6);
    }); */

    it('f^(0) ... f^(20) at at z = 1 + i with r0 = 1e5', function () {
      var ans = deriv(f, 20, 1, 1, {r: 10000}, status);
      assert.isFalse(status.degenerate, 'Is not degenerate');
      assert.almostEqual(ans, analytical(20, 1, 1), status.truncationError, status.roundingError, 6);
    });

    it('f^(0) ... f^(20) at at z = 1 + i', function () {
      var ans = deriv(f, 20, 1, 1, null, status);
      assert.isFalse(status.degenerate, 'Is not degenerate');
      assert.almostEqual(ans, analytical(20, 1, 1), status.truncationError, status.roundingError);
    });

    it('f^(0) ... f^(20) at at z = 0', function () {
      var ans = deriv(f, 20, 1, 1, null, status);
      assert.isFalse(status.degenerate, 'Is not degenerate');
      assert.almostEqual(ans, analytical(20, 1, 1), status.truncationError, status.roundingError);
    });
  });

  describe('e^z / (cos(z)^3 + sin(z)^3)', function () {
    var f = function (out, a, b) {
      var fz = Complex(Math.E).pow(a, b).div(Complex(a, b).cos().pow(3).add(Complex(a, b).sin().pow(3)));
      out[0] = fz.re;
      out[1] = fz.im;
    };

    it('f^(0) ... f^(11) at z = 0', function () {
      var ans = deriv(f, 11, 0, 0, null, status);
      assert.almostEqual(ans,
        [
          [1, 1, 4, 4, 28, -164, 64, -13376, 47248, -858224, 13829824, -112705856],
          [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
        ],
        status.truncationError, status.roundingError, 3
      );
    });

    it('f^(51) at z = 0', function () {
      var fp = deriv(f, 51, 0, 0, {}, status);
      assert.isFalse(status.degenerate, 'Is not degenerate');
      assert.almostEqual([[fp[0][50]], [fp[1][50]]], [[0.1464836745910329202251099956e70], [0]], status.truncationError[50], status.roundingError[50], 2);
    });
  });
});
