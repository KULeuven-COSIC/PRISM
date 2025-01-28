# PRISM

SageMath implementation of **PRISM: PRime degree ISogeny Mechanism**.

We give a proof of concept implementation of PRISM. In some cases, cleaner and more readable code is preferred over a fully optimized implementation.

The code, and in particular the ideal-to-isogeny translation algorithm, is based on the SQIsign2D-West SageMath implementation, which has been privately shared with us by the authors. This can be found in the folder `sqisign2d_west`.
The code to compute (2,2)-isogenies using theta coordinates is based on [ThetaIsogenies/two-isogenies](https://github.com/ThetaIsogenies/two-isogenies). It can be found in `theta_isogenies` and `theta_structures`.
The Kummer line code is based on
[FESTA-PKE/FESTA-SageMath](https://github.com/FESTA-PKE/FESTA-SageMath). It can be found in `montgomery_isogenies`.

## How to run

There are two different version of the protocol that can be run

- `sage --python test_id.py` runs the identification protocol `PRISM_id`;
- `sage --python test_sign.py` runs the signature `PRISM_sign.py`;

All these implementations use NIST Security Level I parameters.
