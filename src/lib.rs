mod angle;
mod helpers;
#[cfg(test)]
mod tests;
#[cfg(test)]
#[macro_use(quickcheck)]
extern crate quickcheck_macros;

pub use angle::Angle;
