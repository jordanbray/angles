use crate::angle::Angle;
use crate::helpers::to_f64;
use nalgebra::RealField;
use num_traits::identities::Zero;
use std::convert::*;

#[test]
fn test_factorial_thing() {
    println!("{:?}", to_f64(Angle::pi_bigint()).unwrap());
    println!("3.141592653589793238462643383279502884197169399375105820974");
    println!("{:?}", Angle::pi_bigint());
}

#[quickcheck]
fn test_sin(angle: f64) {
    let sin = Angle::from(angle).sin();
    if angle.sin() != sin {
        println!("angle : {:?}", angle);
        println!("mine  : {:?}", sin);
        println!("theirs: {:?}", angle.sin());
        assert!((angle.sin() - sin).abs() < 0.000000000001);
    }
}

#[quickcheck]
fn test_cos(angle: f64) {
    let cos = Angle::from(angle).cos();
    if angle.cos() != cos {
        println!("angle : {:?}", angle);
        println!("mine  : {:?}", cos);
        println!("theirs: {:?}", angle.cos());
        assert!((angle.cos() - cos).abs() < 0.000000000001);
    }
}

#[quickcheck]
fn test_tan(angle: f64) {
    let tan = Angle::from(angle).tan();
    if angle.tan() != tan {
        println!("angle : {:?}", angle);
        println!("mine  : {:?}", tan);
        println!("theirs: {:?}", angle.tan());
        assert!((angle.tan() - tan).abs() < 0.00000001);
    }
}

#[quickcheck]
fn test_angle(angle: f64) {
    let angle = if angle < f64::zero() {
        angle.rem_euclid(f64::two_pi())
    } else if angle > f64::zero() {
        -(-angle.rem_euclid(f64::two_pi()))
    } else {
        angle
    };
    let s = Angle::from(angle);
    let pass = (f64::from(s) - angle).abs() < 0.00000000000001;
    if !pass {
        println!("start:  {:?}", angle);
        println!("end   : {:?}", f64::from(s));
        println!("middle: {:?}", s);
        assert!(pass);
    }
}

#[quickcheck]
fn test_add(angle1: f64, angle2: f64) {
    let angle = angle1 + angle2;
    let angle = if angle < f64::zero() {
        angle.rem_euclid(f64::two_pi())
    } else if angle > f64::zero() {
        -(-angle.rem_euclid(f64::two_pi()))
    } else {
        angle
    };

    let a1 = Angle::from(angle1);
    let a2 = Angle::from(angle2);

    let result = a1 + a2;

    let result_f64 = f64::from(result);
    let pass = (result_f64 - angle).abs() < 0.0000000000001;
    if !pass {
        println!("angle1 : {:?}", angle1);
        println!("angle2 : {:?}", angle2);
        println!("correct: {:?}", angle);
        println!("mine   : {:?}", result_f64);
        assert!(pass);
    }
}
