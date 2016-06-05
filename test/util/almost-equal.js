'use strict';

var assert = require('chai').assert;
var ae = require('almost-equal');
var cmod = require('complex-modulus');

function isNumeric (n) {
  return !isNaN(parseFloat(n)) && isFinite(n);
}

module.exports = function almostEqual (actual, expected, trunc, rounding, safety) {
  var i;
  var an = actual[0].length;
  var en = actual[0].length;
  safety = safety || 1;
  trunc = trunc === undefined ? ae.DBL_EPSILON : trunc;
  rounding = rounding === undefined ? ae.DBL_EPSILON : rounding;
  assert.equal(an, en, 'Actual length (' + an + ') not equal to expected length (' + en + ')');
  for (i = 0; i < en; i++) {
    assert(isNumeric(actual[0][i]), 'Actual is numeric');
    assert(isNumeric(actual[1][i]), 'Actual is numeric');
    assert(isNumeric(expected[0][i]), 'Expected is numeric');
    assert(isNumeric(expected[1][i]), 'Expected is numeric');
    var dre = actual[0][i] - expected[0][i];
    var dim = actual[1][i] - expected[1][i];
    var d = cmod(dre, dim);

    var tr = Array.isArray(trunc) ? trunc[i] : trunc;
    var rn = Array.isArray(rounding) ? rounding[i] : rounding;

    var bound = Math.max(tr, rn);

    assert(isNumeric(tr), 'Truncation error tolerance is numeric (' + tr + ')');
    assert(isNumeric(rn), 'Rounding error tolerance is numeric (' + rn + ')');

    // console.log('d / bound =', d / bound)

    if (d > bound * safety) {
      assert(false,
        'i = ' + i + ':\n\texpected (' + expected[0][i] + ' + ' + expected[1][i] + 'j) !~=' +
        '\n\tactual   (' + actual[0][i] + ' + ' + actual[1][i] + 'j)' +
        '\n\terr:         ' + (d) +
        '\n\ttrunc bound: ' + (tr) +
        '\n\tround bound: ' + (rn) +
        '\n\tratio:       ' + (d / bound)
      );
    }
  }
};
