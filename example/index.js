'use strict';

var deriv = require('../');

console.log(
  deriv(function (a, b) {
    var c = b * b + (1 - a) * (1 - a);
    return [(1 - a) / c, b / c];
  }, 5, 0, 0)
);
