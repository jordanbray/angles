use crate::helpers::{bi, to_f64};
use nalgebra::RealField;
use num::BigRational;
use num_traits::identities::Zero;
use std::ops::{Add, Sub};

#[derive(Copy, Clone, PartialEq, PartialOrd, Debug, Default)]
pub struct Angle {
    clockwise: bool,
    units: Option<u64>,
}

impl Angle {
    pub fn pi_2() -> Angle {
        Angle {
            clockwise: false,
            units: Some(1u64 << 62),
        }
    }

    pub fn pi() -> Angle {
        Angle {
            clockwise: false,
            units: Some(1u64 << 63),
        }
    }

    pub fn two_pi() -> Angle {
        Angle {
            clockwise: false,
            units: None,
        }
    }

    pub fn units(&self) -> Option<u64> {
        self.units
    }

    pub fn clockwise(&self) -> bool {
        self.clockwise
    }

    pub fn pi_bigint() -> BigRational {
        bi(245850922) / bi(78256779)
    }

    /// Returns the cosine of (x * 1<<63) / pi with a very high precision
    ///
    /// This only works over the range [-pi/2 to pi/2]
    fn cos_partial(x: i64) -> f64 {
        assert!(x > -(1<<62) && x < (1<<62));

        let pi = 1i128 << 63;

        let s1 = 8;
        let s2 = s1 - 4;

        // 58 bits per thing
        let a = 1152921504606846976;
        let b = -5689439577989151082;
        let c = 4679376491554475694;
        let d = -1539453160513334483;
        let e = 271317744433296946;
        let f = -29753320046888181;
        let g = 2224647666372392;
        let h = -120639134191123;
        let i = 4959362087252;
        let j = -155841981042;

        let coeffs = vec![i, h, g, f, e, d, c, b, a];

        let mut sum: i128 = j;
        let x = x as i128;

        for i in 0..coeffs.len() {
            let t1 = sum * x;
            let t2 = t1 / pi;
            let t3 = t2 * x;
            let t4 = t3 / pi;
            sum = t4 + coeffs[i];
        }

        // Keep in mind 1<<60 is representable exactly with a float
        (sum as f64) / ((1i128<<60) as f64)
    }

    pub fn cos(&self) -> f64 {
        match (self.clockwise, self.units) {
            (_, None) => 1.0,
            (_, Some(0)) => 1.0,
            (false, Some(x)) => {
                let x = x as i64; // because of the symmetry, this always results in the same angle

                if x < -(1<<62) {
                    // x is in the range [-pi to pi/2)
                    -Angle::cos_partial(x.wrapping_add(i64::min_value()))
                } else if x > (1<<62) {
                    -Angle::cos_partial(x.wrapping_add(i64::min_value()))
                } else {
                    Angle::cos_partial(x)
                }
            },
            (true, Some(x)) => {
                let x = x as i64; // because of the symmetry, this always results in the same angle
                let x = -x;

                if x < -(1<<62) {
                    // x is in the range [-pi to pi/2)
                    -Angle::cos_partial(x.wrapping_add(i64::min_value()))
                } else if x > (1<<62) {
                    -Angle::cos_partial(x.wrapping_add(i64::min_value()))
                } else {
                    Angle::cos_partial(x)
                }
            }
        }
    }

    pub fn sin(&self) -> f64 {
        (*self - Angle::pi_2()).cos()
    }

    pub fn atan(slope: f64) -> Angle {
        slope.atan().into()
    }
}

impl From<f64> for Angle {
    fn from(real: f64) -> Angle {
        // should this be try_from, given NaN and Infinity for floats?
        let mut real = real;
        if real < f64::zero() {
            real = real.rem_euclid(f64::two_pi());
        } else if real > f64::zero() {
            real = -(-real.rem_euclid(f64::two_pi()));
        }

        if real == f64::zero() {
            Angle {
                clockwise: false,
                units: Some(0),
            }
        } else if real == f64::two_pi() {
            Angle {
                clockwise: false,
                units: None,
            }
        } else if real == -f64::two_pi() {
            Angle {
                clockwise: true,
                units: None,
            }
        } else if real > f64::zero() {
            let i = (real * ((Angle::pi().units.unwrap() as f64) / f64::pi())) as u64;
            Angle {
                clockwise: false,
                units: Some(i),
            }
        } else if real < f64::zero() {
            let i = ((-real) * ((Angle::pi().units.unwrap() as f64) / f64::pi())) as u64;
            Angle {
                clockwise: true,
                units: Some(i),
            }
        } else {
            panic!("NaN or Infinity passed into From<T: RealField> to Angle");
        }
    }
}

impl From<Angle> for f64 {
    fn from(angle: Angle) -> f64 {
        match (angle.units, angle.clockwise) {
            (None, true) => -f64::two_pi(),
            (None, false) => f64::two_pi(),
            (Some(0), _) => 0.0,
            (Some(x), clockwise) => {
                let pi = Angle::pi_bigint();
                let mut ratio = bi(x as i64) * -(pi.clone()) / bi(i64::min_value());
                if x > (1u64 << 63) {
                    // but, when converting we ended up with a negative value when we shouldn't have
                    ratio += bi(2) * pi;
                }

                if clockwise {
                    ratio *= bi(-1);
                }

                // ratio is now the value for real
                to_f64(ratio).unwrap()
            }
        }
    }
}

fn add_units(a: Option<u64>, b: Option<u64>) -> Option<u64> {
    match (a, b) {
        (None, None) => None,
        (Some(x), None) => Some(x),
        (None, Some(y)) => Some(y),
        (Some(x), Some(y)) => {
            if x == 0 && y == 0 {
                Some(0)
            } else if x.wrapping_add(y) == 0 {
                None
            } else {
                Some(x.wrapping_add(y))
            }
        }
    }
}

fn sub_units(a: Option<u64>, b: Option<u64>) -> (Option<u64>, bool) {
    match (a, b) {
        (None, None) => (None, false),
        (Some(x), None) => (Some(x), false),
        (None, Some(y)) => (Some(y), true),
        (Some(x), Some(y)) => {
            if x > y {
                (Some(x.wrapping_sub(y)), false)
            } else {
                (Some(y.wrapping_sub(x)), true)
            }
        }
    }
}

impl Add for Angle {
    type Output = Self;

    fn add(self, other: Self) -> Self {
        match (self.clockwise, other.clockwise) {
            (false, false) => Angle {
                units: add_units(self.units, other.units),
                clockwise: false,
            },
            (false, true) => {
                let (units, flip) = sub_units(self.units, other.units);
                Angle {
                    units: units,
                    clockwise: flip, // false ^ flip == flip
                }
            }
            (true, false) => {
                let (units, flip) = sub_units(self.units, other.units);
                Angle {
                    units: units,
                    clockwise: true ^ flip,
                }
            }
            (true, true) => Angle {
                units: add_units(self.units, other.units),
                clockwise: true,
            },
        }
    }
}

impl Sub for Angle {
    type Output = Self;

    fn sub(self, other: Self) -> Self {
        let to_add = Angle {
            units: other.units,
            clockwise: !other.clockwise,
        };

        self + to_add
    }
}
