use nalgebra::{RealField, try_convert, convert};
use std::ops::{Add, Sub, Mul};
use serde_derive::{Serialize, Deserialize};

#[derive(Copy, Clone, PartialEq, PartialOrd, Debug, Default, Serialize, Deserialize)]
pub struct Angle {
    clockwise: bool,
    units: Option<u64>,
}

impl Angle {
    pub fn is_zero(&self) -> bool {
        self.units == Some(0)
    }

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

    /// Returns the cosine of (x * 1<<63) / pi with a very high precision
    ///
    /// This only works over the range [-pi/2 to pi/2]
    fn cos_partial(x: i64) -> f64 {
        assert!(x >= -(1<<62) && x <= (1<<62));

        let pi = 1i128 << 63;

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

    pub fn cosf64(&self) -> f64 {
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

    pub fn cos<T: RealField>(self) -> T {
        convert(self.cosf64())
    }

    pub fn sin(&self) -> f64 {
        (*self - Angle::pi_2()).cos()
    }

    pub fn tan(&self) -> f64 {
        self.sin() / self.cos::<f64>()
    }

    pub fn acos(x: f64) -> Angle {
        let mut min = 0;
        let mut max = 1<<63;

        let mut error = max - min;
        while error > 1 {
            let middle = (max + min) / 2;
            let cos: f64 = (Angle { units: Some(middle), clockwise: false }).cos();
            if cos == x {
                return Angle { units: Some(middle), clockwise: false };
            } else if cos < x {
                max = middle;
            } else {
                min = middle;
            }
            error = max - min;
        }

        Angle { units: Some((max + min) / 2), clockwise: false }
    }

    pub fn asin(x: f64) -> Angle {
        let mut min: i64 = -(1<<62);
        let mut max: i64 = 1<<62;

        while min + 1 < max {
            let middle = (max + min) / 2;
            let sin = (Angle { units: Some(middle as u64), clockwise: false }).sin();
            if sin == x {
                break;
            } else if sin < x {
                min = middle;
            } else {
                max = middle;
            }
        }

        let x = ((max + min) / 2) as u64;
        if x > (1<<63) {
            Angle { units: Some(u64::max_value() - x + 1), clockwise: true }
        } else {
            Angle { units: Some(x), clockwise: false }
        }
    }

    pub fn atan(x: f64) -> Angle {
        let mut min: i64 = -(1<<62);
        let mut max: i64 = 1<<62;

        while min + 1 < max {
            let middle = (max + min) / 2;
            let sin = (Angle { units: Some(middle as u64), clockwise: false }).tan();
            if sin == x {
                break;
            } else if sin < x {
                min = middle;
            } else {
                max = middle;
            }
        }

        let x = ((max + min) / 2) as u64;
        if x > (1<<63) {
            Angle { units: Some(u64::max_value() - x + 1), clockwise: true }
        } else {
            Angle { units: Some(x), clockwise: false }
        }
    }

    pub fn atan2(y: f64, x: f64) -> Option<Angle> {
        if x == 0.0 && y == 0.0 {
            None
        } else {
            let r = if x == 0.0 {
                y.abs()
            } else if y == 0.0 {
                x.abs() 
            } else {
                (x * x + y * y).sqrt()
            };

            if (r + x).abs() >= y.abs() { // x == 0, y != 0 or x and y both != 0
                Some(Angle::atan(y / (r + x)) * 2)
            } else if y.abs() > 0.0 { // x and y both != 0
                Some(Angle::atan((r - x) / y) * 2)
            } else if x < 0.0 { // x < 0 && y == 0
                Some(Angle::pi())
            } else { // x > 0 && y == 0
                Some(Angle::pi() * -1)
            }
        }
    }

    pub fn into<T: RealField>(self) -> T {
        match (self.units, self.clockwise) {
            (None, true) => -T::two_pi(),
            (None, false) => T::two_pi(),
            (Some(0), _) => T::zero(),
            (Some(x), clockwise) => {
                let x: T = convert::<f64, T>(x as f64) / convert::<f64, T>((1u64<<63) as f64);
                let x = x * T::pi();

                if !clockwise {
                    x
                } else {
                    -x
                }
            }
        }
    }
}

impl<T: RealField> From<T> for Angle {
    fn from(real: T) -> Angle {
        // should this be try_from, given NaN and Infinity for floats?
        let mut real = real;
        if real < T::zero() {
            real = -convert::<f64, T>(-f64::rem_euclid(try_convert(real).unwrap(), f64::two_pi()));
        } else if real > T::zero() {
            real = convert::<f64, T>(f64::rem_euclid(try_convert(real).unwrap(), f64::two_pi()));
        }

        if real == T::zero() {
            Angle {
                clockwise: false,
                units: Some(0),
            }
        } else if real == T::two_pi() {
            Angle {
                clockwise: false,
                units: None,
            }
        } else if real == -T::two_pi() {
            Angle {
                clockwise: true,
                units: None,
            }
        } else if real > T::zero() {
            let i = try_convert((real * T::from_u64(Angle::pi().units.unwrap()).unwrap()) / T::pi()).unwrap();
            Angle {
                clockwise: false,
                units: Some(i as u64),
            }
        } else if real < T::zero() {
            let i = try_convert((-real * T::from_u64(Angle::pi().units.unwrap()).unwrap()) / T::pi()).unwrap();
            Angle {
                clockwise: true,
                units: Some(i as u64),
            }
        } else {
            panic!("NaN or Infinity passed into From<T: RealField> to Angle");
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

impl Mul<i64> for Angle {
    type Output = Angle;

    fn mul(self, rhs: i64) -> Self {
        let flip_clockwise = rhs < 0;

        let rhs: u64 = if rhs > 0 {
            rhs as u64
        } else if rhs == i64::min_value() {
            1<<63
        } else {
            rhs.abs() as u64
        };

        match self.units {
            None => Angle { units: None, clockwise: self.clockwise ^ flip_clockwise },
            Some(x) => {
                let units = x * rhs;
                if units == 0 { // for example, pi * 2
                    Angle { units: None, clockwise: self.clockwise ^ flip_clockwise }
                } else {
                    Angle { units: Some(units), clockwise: self.clockwise ^ flip_clockwise }
                }
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
