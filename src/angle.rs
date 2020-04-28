use nalgebra::{convert, try_convert, RealField};
use serde_derive::{Deserialize, Serialize};
use std::ops::{Add, AddAssign, Div, DivAssign, Mul, MulAssign, Neg, Sub, SubAssign};

/// Implements an angle structure where angles are stored as integers.
///
/// Also keeps track of if the angle is going clockwise or counter-clockwise.  It also stores
/// a full circle different from a 0 degree arc.  This gets rid of precision issues when
/// storing angles as `pi` can be represented exactly if you change the units.
#[derive(Copy, Clone, Eq, PartialEq, Debug, Default, Serialize, Deserialize)]
pub struct Angle {
    clockwise: bool,
    units: Option<u64>,
}

use quickcheck::{Arbitrary, Gen};
use rand::Rng;

impl Arbitrary for Angle {
    fn arbitrary<G: Gen>(g: &mut G) -> Self {
        let clockwise = g.gen();
        let units = g.gen();
        Angle {
            clockwise: clockwise,
            units: if g.gen() { None } else { Some(units) },
        }
    }
}

impl Angle {
    /// Is this angle zero?
    ///
    /// ```
    /// use integer_angles::Angle;
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
    /// use integer_angles::Angle;
    ///
    /// assert_eq!(Angle::pi_over(2) * 2, Angle::pi());
    /// ```
    pub fn pi_over(n: i64) -> Angle {
        Angle::pi() / n
    }

    /// Don't change the direction of the angle, but possibly
    /// change the direction to get to that angle.
    ///
    /// ```
    /// use integer_angles::Angle;
    ///
    /// assert_eq!(Angle::pi_over(2).set_clockwise(true), -3 * Angle::pi_over(2));
    /// ```
    pub fn set_clockwise(&self, clockwise: bool) -> Angle {
        if clockwise && !self.clockwise {
            *self - Angle::two_pi()
        } else if !clockwise && self.clockwise {
            *self + Angle::two_pi()
        } else {
            *self
        }
    }

    /// Add or subtract an angle depending on the value of the `add` boolean.
    ///
    /// ```
    /// use integer_angles::Angle;
    ///
    /// assert_eq!(Angle::pi() + Angle::pi(), Angle::pi().add_or_subtract(true, Angle::pi()));
    /// assert_eq!(Angle::pi() - Angle::pi(), Angle::pi().add_or_subtract(false, Angle::pi()));
    /// ```
    pub fn add_or_subtract(&self, add: bool, other: Angle) -> Angle {
        if add {
            *self + other
        } else {
            *self - other
        }
    }

    /// Create a 0-degree angle
    ///
    /// ```
    /// use integer_angles::Angle;
    ///
    /// assert_eq!(Angle::zero() * 2, Angle::zero());
    /// ```
    pub fn zero() -> Angle {
        Angle {
            clockwise: false,
            units: Some(0),
        }
    }

    /// Determine if the angles are equal while allowing wrapping (meaning, 0 == 2pi).
    ///
    /// ```
    /// use integer_angles::Angle;
    ///
    /// assert!(Angle::zero().wrapping_eq(Angle::two_pi()));
    /// assert!(Angle::pi().wrapping_eq(-Angle::pi()));
    /// ```
    pub fn wrapping_eq(&self, other: Angle) -> bool {
        self.to_i64() == other.to_i64()
    }

    pub fn wrapping_approx_eq(&self, other: Angle, error: Angle) -> bool {
        let error = error.to_i64().abs();
        (self.to_i64() - other.to_i64()) < error
    }

    /// Get the length of this as an arc given a specific radius
    ///
    /// ```
    /// use integer_angles::Angle;
    /// use nalgebra::RealField;
    ///
    /// assert_eq!(Angle::zero().length(1.0f64), 0.0f64);
    /// assert_eq!(Angle::two_pi().length(1.0f64), f64::two_pi());
    /// assert_eq!(Angle::pi_2().length(1.0f64), f64::frac_pi_2());
    /// ```
    pub fn length<T: RealField>(&self, radius: T) -> T {
        if let Some(x) = self.units {
            radius * convert::<f64, T>(x as f64) / convert::<f64, T>((1u64 << 63) as f64) * T::pi()
        } else {
            radius * T::two_pi()
        }
    }

    /// Create an angle of pi over two counter-clockwise (for convenience).
    ///
    /// ```
    /// use integer_angles::Angle;
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
    /// use integer_angles::Angle;
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
    /// use integer_angles::Angle;
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
    /// use integer_angles::Angle;
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
    /// use integer_angles::Angle;
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
        assert!(x >= -(1 << 62) && x <= (1 << 62));

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
            sum *= x; // this order of operations is 'blessed' in that it prevents
            sum /= pi; // overflow on the i128 value, while keeping enough bits of
            sum *= x; // precision to convert to a float
            sum /= pi;
            sum += coeffs[i];
        }

        // Keep in mind 1<<60 is representable exactly with a float
        T::from_i64(sum as i64).unwrap() / T::from_i64(1i64 << 60).unwrap()
    }

    /// Returns the angle as an i64 (for doing math), and wrap 2pi to 0.  If clockwise, return
    /// the negative of the result.
    ///
    /// ```
    /// use integer_angles::Angle;
    ///
    /// assert_eq!(Angle::pi().to_i64(), i64::min_value());
    /// assert_eq!(Angle::two_pi().to_i64(), 0);
    /// assert_eq!(-Angle::pi_2().to_i64(), i64::min_value() / 2);
    /// assert_eq!(Angle::pi_2().to_i64(), -(i64::min_value() / 2));
    /// ```
    pub fn to_i64(&self) -> i64 {
        (self.units.unwrap_or(0) as i64).wrapping_mul(if self.clockwise { -1 } else { 1 })
    }

    /// Returns the cosine of the angle.
    ///
    /// ```
    /// use integer_angles::Angle;
    ///
    /// assert_eq!(Angle::zero().cos::<f64>(), 1.0f64);
    /// assert_eq!(Angle::pi_2().cos::<f64>(), 0.0f64);
    /// assert_eq!(Angle::pi().cos::<f64>(), -1.0f64);
    /// assert_eq!((Angle::pi_2() * 3).cos::<f64>(), 0.0f64);
    /// assert_eq!(Angle::two_pi().cos::<f64>(), 1.0f64);
    /// ```
    pub fn cos<T: RealField>(&self) -> T {
        let x = self.to_i64();
        if x < -(1 << 62) {
            // x is in the range [-pi to pi/2)
            -Angle::cos_partial::<T>(x.wrapping_add(i64::min_value()))
        } else if x > (1 << 62) {
            -Angle::cos_partial::<T>(x.wrapping_add(i64::min_value()))
        } else {
            Angle::cos_partial(x)
        }
    }

    /// Returns the cosine of the angle.
    ///
    /// ```
    /// use integer_angles::Angle;
    ///
    /// assert_eq!(Angle::zero().sin::<f64>(), 0.0f64);
    /// assert_eq!(Angle::pi_2().sin::<f64>(), 1.0f64);
    /// assert_eq!(Angle::pi().sin::<f64>(), 0.0f64);
    /// assert_eq!((3 * Angle::pi_2()).sin::<f64>(), -1.0f64);
    /// assert_eq!(Angle::two_pi().sin::<f64>(), 0.0f64);
    /// ```
    pub fn sin<T: RealField>(&self) -> T {
        (*self - Angle::pi_2()).cos() // no need to redo the polynomial, as the subtraction of pi_2 is lossless
    }

    /// Returns the tangent of the angle.
    ///
    /// ```
    /// use integer_angles::Angle;
    ///
    /// assert_eq!(Angle::pi_over(4).tan::<f64>(), 1.0f64);
    /// ```
    pub fn tan<T: RealField>(&self) -> T {
        self.sin::<T>() / self.cos::<T>() // here... yeah, I should make a new polynomial
    }

    /// Compute the arccos of the number.
    ///
    /// ```
    /// use integer_angles::Angle;
    ///
    /// assert_eq!(Angle::acos(0.0f64), Angle::pi_over(2));
    /// ```
    pub fn acos<T: RealField>(x: T) -> Angle {
        let mut min = 0;
        let mut max = 1 << 63;

        let mut error = max - min;
        while error > 1 {
            let middle = (max + min) / 2;
            let cos: T = (Angle {
                units: Some(middle),
                clockwise: false,
            })
            .cos();
            if cos == x {
                return Angle {
                    units: Some(middle),
                    clockwise: false,
                };
            } else if cos < x {
                max = middle;
            } else {
                min = middle;
            }
            error = max - min;
        }

        Angle {
            units: Some((max + min) / 2),
            clockwise: false,
        }
    }

    /// Compute the arcsin of the number.
    ///
    /// ```
    /// use integer_angles::Angle;
    ///
    /// assert_eq!(Angle::asin(0.5f64).sin::<f64>(), 0.5f64);
    /// ```
    pub fn asin<T: RealField>(x: T) -> Angle {
        let mut min: i64 = -(1 << 62);
        let mut max: i64 = 1 << 62;

        while min + 1 < max {
            let middle = (max + min) / 2;
            let sin: T = (Angle {
                units: Some(middle as u64),
                clockwise: false,
            })
            .sin();
            if sin == x {
                break;
            } else if sin < x {
                min = middle;
            } else {
                max = middle;
            }
        }

        let x = ((max + min) / 2) as u64;
        if x > (1 << 63) {
            Angle {
                units: Some(u64::max_value() - x + 1),
                clockwise: true,
            }
        } else {
            Angle {
                units: Some(x),
                clockwise: false,
            }
        }
    }

    /// Compute the arctan of the number.
    ///
    /// ```
    /// use integer_angles::Angle;
    ///
    /// assert_eq!(Angle::atan(1.0f64), Angle::pi_over(4));
    /// ```
    pub fn atan<T: RealField>(x: T) -> Angle {
        let mut min: i64 = -(1 << 62);
        let mut max: i64 = 1 << 62;

        while min + 1 < max {
            let middle = (max + min) / 2;
            let tan: T = (Angle {
                units: Some(middle as u64),
                clockwise: false,
            })
            .tan();
            if tan == x {
                break;
            } else if tan < x {
                min = middle;
            } else {
                max = middle;
            }
        }

        let x = ((max + min) / 2) as u64;
        if x > (1 << 63) {
            Angle {
                units: Some(u64::max_value() - x + 1),
                clockwise: true,
            }
        } else {
            Angle {
                units: Some(x),
                clockwise: false,
            }
        }
    }

    /// Compute the atan(y / x) but keeping the sign of y and x.
    ///
    /// Returns None if both y and x are 0.
    ///
    /// ```
    /// use integer_angles::Angle;
    ///
    /// assert_eq!(Angle::atan2(1.0f64, 0.0).unwrap(), Angle::pi_2());
    /// assert_eq!(Angle::atan2(0.0, 0.0f64), None);
    /// ```
    pub fn atan2<T: RealField>(y: T, x: T) -> Option<Angle> {
        if x == T::zero() && y == T::zero() {
            None // I strongly disagree with IEEE on the correct result here.
                 // This is clearly the correct result.
        } else {
            let r = if x == T::zero() {
                y.abs()
            } else if y == T::zero() {
                x.abs()
            } else {
                x.hypot(y)
            };

            if (r + x).abs() >= y.abs() {
                // x == 0, y != 0 or x and y both != 0
                Some(Angle::atan(y / (r + x)) * 2)
            } else if y.abs() > T::zero() {
                // x and y both != 0
                Some(Angle::atan((r - x) / y) * 2)
            } else if x < T::zero() {
                // x < 0 && y == 0
                Some(Angle::pi())
            } else {
                // x > 0 && y == 0
                Some(Angle::pi() * -1)
            }
        }
    }

    /// Convert the angle into radians.  This can lose precision pretty easily
    ///
    /// ```
    /// use integer_angles::Angle;
    /// use nalgebra::RealField;
    ///
    /// assert_eq!(Angle::pi().radians::<f64>(), f64::pi());
    /// ```
    pub fn radians<T: RealField>(self) -> T {
        match (self.units, self.clockwise) {
            (None, true) => -T::two_pi(),
            (None, false) => T::two_pi(),
            (Some(0), _) => T::zero(),
            (Some(x), clockwise) => {
                let x: T = convert::<f64, T>(x as f64) / convert::<f64, T>((1u64 << 63) as f64);
                let x = x * T::pi();

                if !clockwise {
                    x
                } else {
                    -x
                }
            }
        }
    }

    pub fn degrees<T: RealField>(self) -> T {
        match (self.units, self.clockwise) {
            (None, true) => -T::two_pi(),
            (None, false) => T::two_pi(),
            (Some(0), _) => T::zero(),
            (Some(x), clockwise) => {
                let x: T = convert::<f64, T>(x as f64) / convert::<f64, T>((1u64 << 63) as f64);
                let x = x * convert::<f64, T>(180.0f64);

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
    /// Convert from radians into an angle.  Panics if NaN or infinity is received
    ///
    /// ```
    /// use integer_angles::Angle;
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
            // This is technically lossy, but I don't care
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
            let i =
                try_convert((real * T::from_u64(Angle::pi().units.unwrap()).unwrap()) / T::pi())
                    .unwrap();
            Angle {
                clockwise: false,
                units: Some(i as u64),
            }
        } else if real < T::zero() {
            let i =
                try_convert(((-real) * T::from_u64(Angle::pi().units.unwrap()).unwrap()) / T::pi())
                    .unwrap();
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

impl Mul<i64> for Angle {
    type Output = Angle;

    fn mul(self, rhs: i64) -> Self {
        let flip_clockwise = rhs < 0;

        let rhs: u64 = if rhs > 0 {
            rhs as u64
        } else if rhs == i64::min_value() {
            1 << 63
        } else {
            rhs.abs() as u64
        };

        match self.units {
            None => Angle {
                units: None,
                clockwise: self.clockwise ^ flip_clockwise,
            },
            Some(x) => {
                let units = x.wrapping_mul(rhs);
                if units == 0 && x != 0 && rhs != 0 {
                    // we wrapped to 0, so make this two pi, instead of 0 radians
                    Angle {
                        units: None,
                        clockwise: self.clockwise ^ flip_clockwise,
                    }
                } else {
                    Angle {
                        units: Some(units),
                        clockwise: self.clockwise ^ flip_clockwise,
                    }
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
            1u64 << 63
        } else {
            (-rhs) as u64
        };

        match self.units {
            None => {
                if rhs == 1 {
                    Angle {
                        units: self.units,
                        clockwise: self.clockwise ^ invert_clockwise,
                    }
                } else {
                    Angle {
                        units: Some(((1i128 << 64) / (rhs as i128)) as u64),
                        clockwise: self.clockwise ^ invert_clockwise,
                    }
                }
            }
            Some(x) => Angle {
                units: Some(x / rhs),
                clockwise: self.clockwise ^ invert_clockwise,
            },
        }
    }
}

impl DivAssign<i64> for Angle {
    fn div_assign(&mut self, rhs: i64) {
        *self = *self / rhs;
    }
}

fn add(units1: Option<u64>, units2: Option<u64>) -> Option<u64> {
    match (units1, units2) {
        (None, None) => None,
        (None, Some(x)) => Some(x),
        (Some(x), None) => Some(x),
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

fn sub(units1: Option<u64>, units2: Option<u64>) -> Option<u64> {
    match (units1, units2) {
        (None, None) => Some(0),
        (None, Some(x)) => {
            if x == 0 {
                None
            } else {
                Some(u64::max_value().wrapping_sub(x - 1))
            }
        }
        (Some(x), None) => {
            if x == 0 {
                None
            } else {
                Some(u64::max_value().wrapping_sub(x - 1))
            }
        }
        (Some(x), Some(y)) => {
            if x > y {
                Some(x - y)
            } else {
                Some(y - x)
            }
        }
    }
}

impl Add for Angle {
    type Output = Self;

    fn add(self, other: Self) -> Self {
        if self.clockwise == other.clockwise {
            Angle {
                units: add(self.units, other.units),
                clockwise: self.clockwise,
            }
        } else {
            // either x - y or -x + y == y - x == -(x - y).  So, the result (x - y) is useful
            let result = sub(self.units, other.units);
            let flip_sign = if self.units.is_none() {
                false
            } else if other.units.is_none() {
                true
            } else {
                self.units.unwrap() < other.units.unwrap()
            };

            Angle {
                units: result,
                clockwise: self.clockwise ^ flip_sign,
            }
        }
    }
}

use std::cmp::Ordering;

impl Ord for Angle {
    fn cmp(&self, other: &Self) -> Ordering {
        if self == other {
            Ordering::Equal
        } else if !self.clockwise && other.clockwise {
            Ordering::Greater
        } else if self.clockwise && !other.clockwise {
            Ordering::Less
        } else if self.clockwise && other.clockwise {
            if self.units.is_none() {
                Ordering::Less
            } else if other.units.is_none() {
                Ordering::Greater
            } else if self.units.unwrap() < other.units.unwrap() {
                Ordering::Greater
            } else {
                Ordering::Less
            }
        } else {
            if self.units.is_none() {
                Ordering::Greater
            } else if other.units.is_none() {
                Ordering::Less
            } else if self.units.unwrap() < other.units.unwrap() {
                Ordering::Less
            } else {
                Ordering::Greater
            }
        }
    }
}

impl PartialOrd for Angle {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
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
