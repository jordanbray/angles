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

    fn sin_partial(x: i64, clockwise: bool) -> BigRational {
        let n = 12;
        let pi = Angle::pi_bigint();
        let zero = bi(0);
        if x == 0 {
            zero
        } else {
            // Convert x into radians, but keep it as a rational, based on the computed ratio of pi
            let mut x = bi(x) * -pi / bi(i64::min_value());
            if clockwise {
                x *= bi(-1);
            }

            let mut sum = x.clone();
            let mut t = x.clone();
            let x_2 = x.clone() * x.clone();

            for i in 1..n {
                t = (t * bi(-1) * x_2.clone()) / (bi(2) * bi(i) * (bi(2) * bi(i) + bi(1)));
                sum = sum + t.clone();
            }
            sum
        }
    }

    fn cos_partial(x: i64) -> BigRational {
        let n = 13;
        let pi = Angle::pi_bigint();
        let one = bi(1);
        if x == 0 {
            one
        } else {
            let x = bi(x) * -pi / bi(i64::min_value());

            let mut sum = one.clone();
            let mut t = one.clone();
            let x_2 = x.clone() * x.clone();
            for i in 1..n {
                t = t * bi(-1) * x_2.clone() / bi(2 * i * (2 * i - 1));
                sum = sum + t.clone();
            }
            sum
        }
    }



    pub fn sin(&self) -> f64 {
        match (self.units, self.clockwise) {
            (None, _) => 0.0,                         // 360 degrees
            (Some(x), _) if x == 1 << 63 => 0.0,      // 180 degrees
            (Some(x), false) if x == 1 << 62 => 1.0,  // 90 degrees counter-clockwise
            (Some(x), true) if x == 1 << 62 => -1.0,  // 90 degrees clockwise
            (Some(x), false) if x == 3 << 62 => -1.0, // 270 degrees counter-clockwise
            (Some(x), true) if x == 3 << 62 => 1.0,   // 270 degrees clockwise
            (Some(x), clockwise) => {
                // any other value
                let pi_2 = 1i64 << 62;
                let minus_pi = i64::min_value();

                // we want to keep our angle within this range... somehow...
                // this will be our weapon
                // sin(x) = -sin(x + 180)
                //
                // We can just cast this to an i64 due to the way the units work out
                let x: i64 = x as i64;
                to_f64(if x < pi_2 {
                    -Angle::sin_partial(x - minus_pi, clockwise)
                } else if x > pi_2 {
                    -Angle::sin_partial(x + minus_pi, clockwise)
                } else {
                    Angle::sin_partial(x, clockwise)
                })
                .unwrap()
            }
        }
    }

    pub fn cos(&self) -> f64 {
        match self.units {
            None => 1.0,                     // 360 degrees
            Some(x) if x == 1 << 63 => -1.0, // 180 degrees
            Some(x) if x == 1 << 62 => 0.0,  // 90 degrees
            Some(x) if x == 3 << 62 => 0.0,  // 270 degrees
            Some(x) => {
                // any other value
                let pi_2 = 1i64 << 62;
                let minus_pi = i64::min_value();

                // we want to keep our angle within this range... somehow...
                // this will be our weapon
                // sin(x) = -sin(x + 180)
                //
                // We can just cast this to an i64 due to the way the units work out
                let x: i64 = x as i64;
                to_f64(if x < pi_2 {
                    -Angle::cos_partial(x - minus_pi)
                } else if x > pi_2 {
                    -Angle::cos_partial(x + minus_pi)
                } else {
                    Angle::cos_partial(x)
                })
                .unwrap()
            }
        }
    }

    pub fn tan(&self) -> f64 {
        self.sin() / self.cos()
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
