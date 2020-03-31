//! # Angles Done With Integers
//!
//! ```
//! use integer_angles::Angle;
//!
//! assert_eq!(Angle::pi_2().cos::<f64>(), 0.0f64);
//! ```
//!
//! Here we go, down the rabbit hole of floating-point instability and all sorts of crazy problems
//! that come with representing angles within computers.  The goal of this library is to solve the
//! following problems:
//!
//! * If you have multiple angles, and you add them together, the result you get should be exactly
//! correct.
//! * If you add multiple angles together and end up with a full circle, that should be exactly a
//! full circle.
//! * If you do trigonometry of some multiple of `pi` radians, you should end up with the exact
//! answer.
//! * Keep track of the difference between a `0` radian angle, and a `2 pi` radians angle.
//! * Keep track of if the angle is going clockwise or counter-clockwise starting at the positive x
//! axis.
//! * Do not allow the user to represent an angle outside the range [`-2 pi` to `2 pi`]
//!
//! The way this library does it's magic is the following:
//!
//! * Stores the angle in units of `[0..2**64)` where each unit is `1/(2**64)`th of a circle.
//! * This means that adding and subtracting angles (with wrapping) will always be correct, and
//!   always within the specified range. (No more range reduction!)
//! * This also means that you can (inside the library) *cast* an angle from `u64` to `i64` and
//!   end up with the same angle.
//! * Set a flag for a full circle, and allow units to be `0` for a `0` degree angle.
//! * This also means, for example, `pi` radians is exactly equal to `1<<63` units in this library.
//! * Keep track of the clockwise/counterclockwise-ness of the angle using a separate flag.
//! * Solves the Chebyshev to compute the sin/cos/tan using the new units (with more precision
//!   than the standard library).
//! * Uses a binary search (at the moment) to compute asin/acos/atan/atan2.
//!
//! Caveats:
//! * This library is slower than using an f64 (about 10 times slower to compute `cos`.  You've
//!   gotta wait a whole 80 ns to get the result!).
//! * ... Probably other things.

mod angle;

#[cfg(test)]
mod tests;

#[cfg(test)]
#[macro_use(quickcheck)]
extern crate quickcheck_macros;

pub use angle::Angle;
