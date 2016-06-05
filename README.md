# complex-deriv-fornberg

[![Build Status][travis-image]][travis-url] [![npm version][npm-image]][npm-url]  [![Dependency Status][daviddm-image]][daviddm-url] [![js-semistandard-style][semistandard-image]][semistandard-url]

> Compute the derivative of a complex analytic function using the method of Fornberg

## Introduction

This module uses the method of [Fornberg](#references) to compute the derivatives (with error bounds) of a complex anaytic function. In essence, it uses a Fourier Transform to invert function values into the Taylor Series coefficients, from which computation of the derivative is then trivial.

## Restrictions

It uses the coefficients themselves to control the truncation error, so the error will not be properly bounded for functions like low-order polynomials whose Taylor series coefficients are nearly zero.

## Example

To compute the first five derivatives of `1 / (1 - z)` at `z = 0`:

```javascript
var deriv = require('complex-deriv-fornberg');

function f(a, b) {
  var c = b * b + (1 - a) * (1 - a);
  return [(1 - a) / c, b / c];
}

deriv(f, 5, 0, 0)

// =>
// [ [ 1.0000000000000138,
//    1.0000000000000107,
//    2.0000000000000178,
//    6.000000000000274,
//    23.99999999999432,
//    120.00000000001907 ],
//  [ 3.260944019499375e-17,
//    -7.34322880255595e-17,
//    -1.0399525753674522e-15,
//    -1.9409005869118293e-14,
//    -1.3097221675086528e-13,
//    7.683410164000554e-13 ] ]

```

## Installation

```bash
$ npm install complex-deriv-fornberg
```

## API

#### `require('complex-deriv-fornberg')([output, ]f, n, a, b[, options[, status]])`
Compute the derivative of a complex analytic function `f` at `a + b * i`.

**Parameters**:
- `output` (optional). Optional array of arrays into which the output is written. If not provided, arrays will be allocated and returned.
- `f`: function of format `function([out, ]a, b)`, that evaluates the function at `a + b * i`, either `a` into `out[0]` and `b` into `out[1]` or simply returning `[a, b]`.
- `n`: Number of derivatives to compute where 0 represents the value of the function and `n` represents the nth derivative. Maximum number is 100.
- `a`: Real component of `z` at which to evaluate the derivatives
- `b`: Imaginary component of `z` at which to evaluate the derivatives
- `options`: Optional object of configuration parameters
  - `r` (default: `0.6580924658`): Initial radius at which to evaluate. For well-behaved functions, the computation should be insensitive to the initial radius to within about four orders of magnitude.
  - `maxIters` (default: `30`): Maximum number of iterations
  - `taylor`: (default: `false`): If false, output represents the derivatives of the function. If true, the output represents Taylor series coefficients, differing only in multiplication by a factorial.
- `status`: Optional object into which output information is written. Fields are:
  - `degenerate`: True if the algorithm was unable to bound the error
  - `iterations`: Number of iterations executed
  - `finalRadius`: Ending radius of the algorithm
  - `failed`: True if the maximum number of iterations was reached
  - `truncationError`: An array containing approximate bounds of the truncation error achieved for each component of the solution
  - `roundingError`: An array containing approximate bounds of the rounding error achieved for each component of the solution

**Returns**: Returns the real and imaginary components of the derivatives in arrays, i.e. `[[re1, re2, ...], [im1, im2, ...]]`, also writing the arrays to `output`, if provided.

## References

\[1\] Fornberg, B. (1981). [Numerical Differentiation of Analytic Functions](https://amath.colorado.edu/faculty/fornberg/Docs/ACM_81_1.pdf). ACM Transactions on Mathematical Software (TOMS), 7(4), 512–526. http://doi.org/10.1145/355972.355979

## License

&copy; 2016 Scijs Authors. MIT License.

## Authors

Ricky Reusser

[npm-image]: https://badge.fury.io/js/complex-deriv-fornberg.svg
[npm-url]: https://npmjs.org/package/complex-deriv-fornberg
[travis-image]: https://travis-ci.org/scijs/complex-deriv-fornberg.svg?branch=master
[travis-url]: https://travis-ci.org//complex-deriv-fornberg
[daviddm-image]: https://david-dm.org/scijs/complex-deriv-fornberg.svg?theme=shields.io
[daviddm-url]: https://david-dm.org//complex-deriv-fornberg
[semistandard-image]: https://img.shields.io/badge/code%20style-semistandard-brightgreen.svg?style=flat-square
[semistandard-url]: https://github.com/Flet/semistandard
