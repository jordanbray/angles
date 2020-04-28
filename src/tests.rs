use crate::angle::Angle;
use nalgebra::RealField;
use quickcheck::{Arbitrary, Gen};
use rand::Rng;

#[quickcheck]
fn test_sin(angle: Angle) {
    let sin: f64 = angle.sin();
    assert!((angle.radians::<f64>().sin() - sin).abs() < 0.00000000000001);
}

#[quickcheck]
fn test_cos(angle: Angle) {
    let cos: f64 = angle.cos::<f64>();
    assert!((angle.radians::<f64>().cos() - cos).abs() < 0.00000000000001);
}

#[quickcheck]
fn test_tan(angle: Angle) {
    let tan: f64 = angle.tan();
    let pass = (angle.radians::<f64>().tan() - tan).abs() < 0.000001;

    if !pass {
        if (angle.radians::<f64>().tan() > 1000.0 && tan > 1000.0)
            || (angle.radians::<f64>().tan() < -1000.0 && tan < -1000.0)
        {
            // honestly, don't do tangent of these kinds of angles...
            return;
        }

        println!("angle : {:?}", angle);
        println!("mine  : {:?}", tan);
        println!("theirs: {:?}", angle.radians::<f64>().tan());
        assert!(pass);
    }
}

#[derive(Copy, Clone, PartialEq, PartialOrd, Debug, Default)]
struct SmallFloat {
    x: f64,
}

impl Arbitrary for SmallFloat {
    fn arbitrary<G: Gen>(g: &mut G) -> Self {
        SmallFloat {
            x: g.gen_range(-1.0, 1.0),
        }
    }
}

#[quickcheck]
fn test_acos(x: SmallFloat) {
    let x = x.x;
    let angle: f64 = Angle::acos(x).radians();
    let other = x.acos();

    if (angle - other).abs() >= 0.0000000000001 {
        println!("Angle: {}\nMine  : {}\nTheirs: {}\n", x, angle, other);
        assert!(false);
    }
}

#[quickcheck]
fn test_asin(x: SmallFloat) {
    let x = x.x;
    let angle: f64 = Angle::asin(x).radians();
    let other = x.asin();

    if (angle - other).abs() >= 0.0000000000001 {
        println!("Angle: {}\nMine  : {}\nTheirs: {}\n", x, angle, other);
        assert!(false);
    }
}

#[quickcheck]
fn test_atan(x: f64) {
    let angle: f64 = Angle::atan(x).radians();
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

    let angle: f64 = Angle::atan2(y, x).unwrap().radians();
    let other = f64::atan2(y, x);
    if (angle - other).abs() >= 0.000000000000001 {
        println!("x: {}, y: {}\nMine  : {}\nTheirs: {}\n", x, y, angle, other);
        assert!(false);
    }
}

#[quickcheck]
fn test_angle(angle: f64) {
    let angle = if angle < 0.0 {
        -(-angle).rem_euclid(f64::two_pi())
    } else if angle > 0.0 {
        angle.rem_euclid(f64::two_pi())
    } else {
        angle
    };
    let s = Angle::from(angle);
    let pass = (s.radians::<f64>() - angle).abs() < 0.00000000000001;
    if !pass {
        println!("start:  {:?}", angle);
        println!("end   : {:?}", s.radians::<f64>());
        println!("middle: {:?}", s);
        assert!(pass);
    }
}

fn range_reduce(f: f64) -> f64 {
    if f < 0.0 {
        -(-f).rem_euclid(f64::two_pi())
    } else if f > 0.0 {
        f.rem_euclid(f64::two_pi())
    } else {
        f
    }
}

#[quickcheck]
fn test_add(a1: Angle, a2: Angle) {
    let angle1: f64 = a1.radians();
    let angle2: f64 = a2.radians();
    let angle = if a1 == Angle::two_pi() && a2 == Angle::two_pi() {
        f64::two_pi()
    } else if a1 == -Angle::two_pi() && a2 == -Angle::two_pi() {
        -f64::two_pi()
    } else {
        range_reduce(angle1 + angle2)
    };

    let result = a1 + a2;

    let result_f64 = result.radians::<f64>();
    let pass = (result_f64 - angle).abs() < 0.0000000000001;
    if !pass {
        println!("angle1 : {:?}", angle1);
        println!("angle2 : {:?}", angle2);
        println!("correct: {:?}", angle);
        println!("mine   : {:?}", result_f64);
        assert!(pass);
    }
}

#[quickcheck]
fn test_simple_ordering(a1: Angle, a2: Angle) {
    let angle1: f64 = a1.radians();
    let angle2: f64 = a2.radians();

    assert_eq!(angle1 == angle2, a1 == a2);
    assert_eq!(angle1 < angle2, a1 < a2);
    assert_eq!(angle1 > angle2, a1 > a2);
}

#[quickcheck]
fn test_ord(a1: Angle, a2: Angle, a3: Angle) {
    if a1 < a2 && a2 < a3 {
        assert!(a1 < a3);
    }
    if a1 == a2 && a2 == a3 {
        assert!(a1 == a3);
    }

    if a1 > a2 && a2 > a3 {
        assert!(a1 > a3);
    }
}
