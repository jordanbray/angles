use crate::angle::Angle;
use nalgebra::RealField;
use quickcheck::{Arbitrary, Gen};
use rand::Rng;

#[quickcheck]
fn test_sin(angle: f64) {
    let sin = Angle::from(angle).sin();
    assert!((angle.sin() - sin).abs() < 0.00000000000001);
}

#[quickcheck]
fn test_cos(angle: f64) {
    let cos = Angle::from(angle).cos::<f64>();
    assert!((angle.cos() - cos).abs() < 0.00000000000001);
}

#[quickcheck]
fn test_tan(angle: f64) {
    let tan = Angle::from(angle).tan();
    let pass = (angle.tan() - tan).abs() < 0.000001;

    if !pass {
        if (angle.tan() > 1000.0  && tan > 1000.0) ||
           (angle.tan() < -1000.0 && tan < -1000.0) { // honestly, don't do tangent of these kinds of angles...
            return;
        }

        println!("angle : {:?}", angle);
        println!("mine  : {:?}", tan);
        println!("theirs: {:?}", angle.tan());
        assert!(pass);
    }
}

#[derive(Copy, Clone, PartialEq, PartialOrd, Debug, Default)]
struct SmallFloat {
    x: f64
}

impl Arbitrary for SmallFloat {
    fn arbitrary<G: Gen>(g: &mut G) -> Self {
        SmallFloat { x: g.gen_range(-1.0, 1.0) }
    }
}

#[quickcheck]
fn test_acos(x: SmallFloat) {
    let x = x.x;
    let angle: f64 = Angle::acos(x).into();
    let other = x.acos();

    if (angle - other).abs() >= 0.0000000000001 {
        println!("Angle: {}\nMine  : {}\nTheirs: {}\n", x, angle, other);
        assert!(false);

    }
}

#[quickcheck]
fn test_asin(x: SmallFloat) {
    let x = x.x;
    let angle: f64 = Angle::asin(x).into();
    let other = x.asin();

    if (angle - other).abs() >= 0.0000000000001 {
        println!("Angle: {}\nMine  : {}\nTheirs: {}\n", x, angle, other);
        assert!(false);
    }
}

#[quickcheck]
fn test_atan(x: f64) {
    let angle: f64 = Angle::atan(x).into();
    let other = x.atan();

    if (angle - other).abs() >= 0.0000000000001 {
        println!("Angle: {}\nMine  : {}\nTheirs: {}\n", x, angle, other);
        assert!(false);
    }
}

#[quickcheck]
fn test_atan2(x: SmallFloat, y: SmallFloat) {
    let x = x.x;
    let y = y.x;
    if x == 0.0 && y == 0.0 {
        return;
    }

    let angle: f64 = Angle::atan2(y, x).unwrap().into();
    let other = f64::atan2(y, x);
    if (angle - other).abs() >= 0.000000000000001 {
        println!("x: {}, y: {}\nMine  : {}\nTheirs: {}\n", x, y, angle, other);
        assert!(false);
    }
}


#[quickcheck]
fn test_angle(angle: f64) {
    let angle = if angle < 0.0 {
        angle.rem_euclid(f64::two_pi())
    } else if angle > 0.0 {
        -(-angle.rem_euclid(f64::two_pi()))
    } else {
        angle
    };
    let s = Angle::from(angle);
    let pass = (s.into::<f64>() - angle).abs() < 0.00000000000001;
    if !pass {
        println!("start:  {:?}", angle);
        println!("end   : {:?}", s.into::<f64>());
        println!("middle: {:?}", s);
        assert!(pass);
    }
}

#[quickcheck]
fn test_add(angle1: f64, angle2: f64) {
    let angle = angle1 + angle2;
    let angle = if angle < 0.0 {
        angle.rem_euclid(f64::two_pi())
    } else if angle > 0.0 {
        -(-angle.rem_euclid(f64::two_pi()))
    } else {
        angle
    };

    let a1 = Angle::from(angle1);
    let a2 = Angle::from(angle2);

    let result = a1 + a2;

    let result_f64 = result.into::<f64>();
    let pass = (result_f64 - angle).abs() < 0.0000000000001;
    if !pass {
        println!("angle1 : {:?}", angle1);
        println!("angle2 : {:?}", angle2);
        println!("correct: {:?}", angle);
        println!("mine   : {:?}", result_f64);
        assert!(pass);
    }
}

