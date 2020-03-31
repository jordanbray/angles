use nalgebra::{RealField, try_convert, convert};
use std::ops::{Add, AddAssign, Sub, SubAssign, Mul, MulAssign, Div, DivAssign, Neg};
use serde_derive::{Serialize, Deserialize};

/// Implements an angle structure where angles are stored as integers.
///
/// Also keeps track of if the angle is going clockwise or counter-clockwise.  It also stores
/// a full circle different from a 0 degree arc.  This gets rid of precision issues when
/// storing angles as `pi` can be represented exactly if you change the units.
#[derive(Copy, Clone, PartialEq, PartialOrd, Debug, Default, Serialize, Deserialize)]
pub struct Angle {
    clockwise: bool,
    units: Option<u64>,
}

impl Angle {
    /// Is this angle zero?
    ///
    /// ```
    /// use angles::Angle;
    ///
    /// assert!(Angle::zero().is_zero());
    /// assert!(!Angle::pi().is_zero());
    /// ```
    pub fn is_zero(&self) -> bool {
        self.units == Some(0)
    }

    /// Divide pi by some number, and return the result
    ///
    /// ```
    /// use angles::Angle;
    ///
    /// assert_eq!(Angle::pi_over(2) * 2, Angle::pi());
    /// ```
    pub fn pi_over(n: i64) -> Angle {
        Angle::pi() / n
    }

    /// Create a 0-degree angle
    ///
    /// ```
    /// use angles::Angle;
    ///
    /// assert_eq!(Angle::zero() * 2, Angle::zero());
    /// ```
    pub fn zero() -> Angle {
        Angle {
            clockwise: false,
            units: Some(0)
        }
    }

    /// Create an angle of pi over two counter-clockwise (for convenience).
    ///
    /// ```
    /// use angles::Angle;
    ///
    /// assert_eq!(Angle::pi_2(), Angle::pi_over(2));
    /// assert_eq!(Angle::pi_2() * 2, Angle::pi());
    /// ```
    pub fn pi_2() -> Angle {
        Angle {
            clockwise: false,
            units: Some(1u64 << 62),
        }
    }

    /// Create an angle of pi counter-clockwise.
    ///
    /// ```
    /// use angles::Angle;
    ///
    /// assert_eq!(Angle::pi() / 2, Angle::pi_2());
    /// ```
    pub fn pi() -> Angle {
        Angle {
            clockwise: false,
            units: Some(1u64 << 63),
        }
    }

    /// Create an angle that is `two_pi` radians counter-clockwise.
    ///
    /// ```
    /// use angles::Angle;
    ///
    /// assert_eq!(Angle::two_pi() * 2, Angle::pi() * 2);
    /// ```
    pub fn two_pi() -> Angle {
        Angle {
            clockwise: false,
            units: None,
        }
    }

    /// Get the number of units in this angle.  This will be some number between 0 and 2^64 - 1
    /// where a full circle is 2^64 units.  If this represents a full cirle, `None` will be
    /// returned.
    ///
    /// ```
    /// use angles::Angle;
    ///
    /// assert_eq!(Angle::zero().units(), Some(0));
    /// assert_eq!(Angle::two_pi().units(), None);
    /// ```
    pub fn units(&self) -> Option<u64> {
        self.units
    }

    /// Is this angle clockwise or counter-clockwise?
    ///
    /// ```
    /// use angles::Angle;
    ///
    /// assert_eq!(Angle::pi().clockwise(), false);
    /// assert_eq!((-Angle::pi()).clockwise(), true);
    /// ```
    pub fn clockwise(&self) -> bool {
        self.clockwise
    }

    /// Returns the cosine of (x * 1<<63) / pi with a very high precision
    ///
    /// This only works over the range [-pi/2 to pi/2]
    fn cos_partial<T: RealField>(x: i64) -> T {
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
        T::from_i64(sum as i64).unwrap() / T::from_i64(1i64<<60).unwrap()
    }

    /// Returns the cosine of the angle.
    ///
    /// ```
    /// use angles::Angle;
    ///
    /// assert_eq!(Angle::zero().cos::<f64>(), 1.0f64);
    /// assert_eq!(Angle::pi_2().cos::<f64>(), 0.0f64);
    /// assert_eq!(Angle::pi().cos::<f64>(), -1.0f64);
    /// assert_eq!((Angle::pi_2() * 3).cos::<f64>(), 0.0f64);
    /// assert_eq!(Angle::two_pi().cos::<f64>(), 1.0f64);
    /// ```
    pub fn cos<T: RealField>(&self) -> T {
        match (self.clockwise, self.units) {
            (_, None) => T::one(),
            (_, Some(0)) => T::one(),
            (false, Some(x)) => {
                let x = x as i64; // because of the symmetry, this always results in the same angle

                if x < -(1<<62) {
                    // x is in the range [-pi to pi/2)
                    -Angle::cos_partial::<T>(x.wrapping_add(i64::min_value()))
                } else if x > (1<<62) {
                    -Angle::cos_partial::<T>(x.wrapping_add(i64::min_value()))
                } else {
                    Angle::cos_partial(x)
                }
            },
            (true, Some(x)) => {
                let x = x as i64; // because of the symmetry, this always results in the same angle
                let x = x.wrapping_neg();

                if x < -(1<<62) {
                    // x is in the range [-pi to pi/2)
                    -Angle::cos_partial::<T>(x.wrapping_add(i64::min_value()))
                } else if x > (1<<62) {
                    -Angle::cos_partial::<T>(x.wrapping_add(i64::min_value()))
                } else {
                    Angle::cos_partial(x)
                }
            }
        }
    }

    /// Returns the cosine of the angle.
    ///
    /// ```
    /// use angles::Angle;
    ///
    /// assert_eq!(Angle::zero().sin::<f64>(), 0.0f64);
    /// assert_eq!(Angle::pi_2().sin::<f64>(), 1.0f64);
    /// assert_eq!(Angle::pi().sin::<f64>(), 0.0f64);
    /// assert_eq!((3 * Angle::pi_2()).sin::<f64>(), -1.0f64);
    /// assert_eq!(Angle::two_pi().sin::<f64>(), 0.0f64);
    /// ```
     pub fn sin<T: RealField>(&self) -> T {
        (*self - Angle::pi_2()).cos()
    }

    /// Returns the tangent of the angle.
    ///
    /// ```
    /// use angles::Angle;
    ///
    /// assert_eq!(Angle::pi_over(4).tan::<f64>(), 1.0f64);
    /// ```
    pub fn tan<T: RealField>(&self) -> T {
        self.sin::<T>() / self.cos::<T>()
    }


    /// Compute the arccos of the number.
    ///
    /// ```
    /// use angles::Angle;
    ///
    /// assert_eq!(Angle::acos(0.0f64), Angle::pi_over(2));
    /// ```
    pub fn acos<T: RealField>(x: T) -> Angle {
        let mut min = 0;
        let mut max = 1<<63;

        let mut error = max - min;
        while error > 1 {
            let middle = (max + min) / 2;
            let cos: T = (Angle { units: Some(middle), clockwise: false }).cos();
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

    /// Compute the arcsin of the number.
    ///
    /// ```
    /// use angles::Angle;
    ///
    /// assert_eq!(Angle::asin(0.5f64).sin::<f64>(), 0.5f64);
    /// ```
    pub fn asin<T: RealField>(x: T) -> Angle {
        let mut min: i64 = -(1<<62);
        let mut max: i64 = 1<<62;

        while min + 1 < max {
            let middle = (max + min) / 2;
            let sin: T = (Angle { units: Some(middle as u64), clockwise: false }).sin();
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

    /// Compute the arctan of the number.
    ///
    /// ```
    /// use angles::Angle;
    /// 
    /// assert_eq!(Angle::atan(1.0f64), Angle::pi_over(4));
    /// ```
    pub fn atan<T: RealField>(x: T) -> Angle {
        let mut min: i64 = -(1<<62);
        let mut max: i64 = 1<<62;

        while min + 1 < max {
            let middle = (max + min) / 2;
            let tan: T = (Angle { units: Some(middle as u64), clockwise: false }).tan();
            if tan == x {
                break;
            } else if tan < x {
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

    /// Compute the atan(y / x) but keeping the sign of y and x.
    ///
    /// ```
    /// use angles::Angle;
    ///
    /// assert_eq!(Angle::atan2(1.0f64, 0.0).unwrap(), Angle::pi_2());
    /// ```
    pub fn atan2<T: RealField>(y: T, x: T) -> Option<Angle> {
        if x == T::zero() && y == T::zero() {
            None
        } else {
            let r = if x == T::zero() {
                y.abs()
            } else if y == T::zero() {
                x.abs() 
            } else {
                (x * x + y * y).sqrt()
            };

            if (r + x).abs() >= y.abs() { // x == 0, y != 0 or x and y both != 0
                Some(Angle::atan(y / (r + x)) * 2)
            } else if y.abs() > T::zero() { // x and y both != 0
                Some(Angle::atan((r - x) / y) * 2)
            } else if x < T::zero() { // x < 0 && y == 0
                Some(Angle::pi())
            } else { // x > 0 && y == 0
                Some(Angle::pi() * -1)
            }
        }
    }

    /// Convert the angle into radians.  This can lose precision pretty easily
    ///
    /// ```
    /// use angles::Angle;
    /// use nalgebra::RealField;
    ///
    /// assert_eq!(Angle::pi().into::<f64>(), f64::pi());
    /// ```
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
    /// Convert from radians into an angle.
    ///
    /// ```
    /// use angles::Angle;
    /// use nalgebra::RealField;
    ///
    /// assert_eq!(Angle::from(f64::pi()), Angle::pi());
    /// ```
    fn from(real: T) -> Angle {
        // should this be try_from, given NaN and Infinity for floats?
        let mut real = real;
        if real < T::zero() {
            real = -convert::<f64, T>(f64::rem_euclid(-try_convert(real).unwrap(), f64::two_pi()));
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
            let i = try_convert(((-real) * T::from_u64(Angle::pi().units.unwrap()).unwrap()) / T::pi()).unwrap();
            Angle {
                clockwise: true,
                units: Some(i as u64),
            }
        } else {
            panic!("NaN or Infinity passed into From<T: RealField> to Angle");
        }
    }
}

#[test]
fn negative_is_positive() {
    assert_eq!(Angle::from(-2.0f64).units, Angle::from(2.0f64).units);
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
                let units = x.wrapping_mul(rhs);
                if units == 0 && x != 0 && rhs != 0 { // we wrapped to 0, so make this two pi, instead of 0 radians
                    Angle { units: None, clockwise: self.clockwise ^ flip_clockwise }
                } else {
                    Angle { units: Some(units), clockwise: self.clockwise ^ flip_clockwise }
                }
            }
        }
    }
}

impl Mul<Angle> for i64 {
    type Output = Angle;

    fn mul(self, rhs: Angle) -> Angle {
        rhs * self
    }
}

impl MulAssign<i64> for Angle {
    fn mul_assign(&mut self, rhs: i64) {
        *self = *self * rhs;
    }
}

impl Div<i64> for Angle {
    type Output = Angle;

    fn div(self, rhs: i64) -> Angle {
        let invert_clockwise = rhs < 0;
        let rhs = if rhs > 0 {
            rhs as u64
        } else if rhs == i64::min_value() {
            1u64<<63
        } else {
            (-rhs) as u64
        };

        match self.units {
            None => if rhs == 1 {
                Angle { units: self.units, clockwise: self.clockwise ^ invert_clockwise }
            } else {
                Angle { units: Some(((1i128<<64)/(rhs as i128)) as u64), clockwise: self.clockwise ^ invert_clockwise }
            },
            Some(x) => {
                Angle { units: Some(x / rhs), clockwise: self.clockwise ^ invert_clockwise }
            }
        }
    }
}

impl DivAssign<i64> for Angle {
    fn div_assign(&mut self, rhs: i64) {
        *self = *self / rhs;
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

impl AddAssign<Angle> for Angle {
    fn add_assign(&mut self, other: Self) {
        *self = *self + other;
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

impl SubAssign<Angle> for Angle {
    fn sub_assign(&mut self, other: Self) {
        *self = *self - other;
    }
}

impl Neg for Angle {
    type Output = Angle;

    fn neg(self) -> Self::Output {
        Angle {
            units: self.units,
            clockwise: !self.clockwise,
        }
    }
}
