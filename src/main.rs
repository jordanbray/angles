mod helpers;

use crate::helpers::to_f64;
use num::{bigint::Sign, rational::Ratio, BigInt, BigRational, FromPrimitive};
use std::str::FromStr;

fn nu(x: &str) -> BigRational {
    BigRational::new(BigInt::from_str(x).unwrap(), BigInt::from_str("1").unwrap())
}

fn r(num: i64, denom: u64) -> BigRational {
    BigRational::new(BigInt::from_i64(num).unwrap(), BigInt::from_u64(denom).unwrap())
}

fn cos(x: u64) -> f64 {
    //let a_j_denom = 661656681561401158;
    //let d_g_denom = 11248163586543819688;
    //let f_i_denom = 3749387862181273229;
    let pi = 1u64 << 63;

    // let a = r(-201303869774232317,  a_j_denom);
    let a = r(-2806138795432098259, pi);
    //let b = r(-7280320382166766806, b_c_denom);
    let b = r(-8954675523107343155, pi);
    //let c = r(2270997891633034674,  b_c_denom);
    let c = r(2793290427581747328, pi);
    let d = r(-268325731830568274, pi);
    //let d = r(-327230834238204960,  d_g_denom);
    //let e = r(193335643301312,      138866217117824934);
    let e = r(12841183429369569, pi);

    //let f = r(-50229282213546 * 3, f_i_denom);
    let f = r(-370687195373421, pi);
    //let g = r(8754017027998,               d_g_denom);
    let g = r(7178198933982, pi);
    // let h = r(-1400277267,          129289236626940456);
    let h = r(-99894457770, pi);
    //let i = r(422594474,            f_i_denom);
    let i = r(1039568643, pi);
    //let j = r(-682350,              a_j_denom);
    let j = r(-9511834, pi);

    let coeffs = vec![a, b, c, d, e, f, g, h, i, j];

    let x = Ratio::new(
        BigInt::from_u64(x).unwrap(),
        BigInt::from_u64(9223372036854775808).unwrap(),
    );

    let two = Ratio::new(BigInt::from_i8(2).unwrap(), BigInt::from_i8(1).unwrap());

    let mut p2 = x.clone();
    let mut p1 = two.clone() * p2.clone() * x.clone()
        - Ratio::new(BigInt::from_i8(1).unwrap(), BigInt::from_i8(1).unwrap());
    let mut sum = Ratio::new(BigInt::from_i8(0).unwrap(), BigInt::from_i8(1).unwrap());
    for i in 0..coeffs.len() {
        if i == 0 {
            sum += coeffs[i].clone() // is same a sum += ... * 1
        } else {
            sum += coeffs[i].clone() * p1.clone();

            // advance two terms
            p2 = p1.clone() * two.clone() * x.clone() - p2;
            p1 = p2.clone() * two.clone() * x.clone() - p1;
        }
    }
    to_f64(sum).unwrap()
}

fn main() {
    for arg in std::env::args().skip(1) {
        println!("cos({}) = {}", nu(&arg), cos(u64::from_str(&arg).unwrap()));
    }
}

#[test]
fn cos_does_not_lose_precision() {
    let questions = vec![
        0,
        144115188075855872,
        288230376151711744,
        432345564227567616,
        576460752303423488,
        720575940379279360,
        864691128455135232,
        1008806316530991104,
        1152921504606846976,
        1297036692682702848,
        1441151880758558720,
        1585267068834414592,
        1729382256910270464,
        1873497444986126336,
        2017612633061982208,
        2161727821137838080,
        2305843009213693952,
        2449958197289549824,
        2594073385365405696,
        2738188573441261568,
        2882303761517117440,
        3026418949592973312,
        3170534137668829184,
        3314649325744685056,
        3458764513820540928,
        3602879701896396800,
        3746994889972252672,
        3891110078048108544,
        4035225266123964416,
        4179340454199820288,
        4323455642275676160,
        4467570830351532032,
        4611686018427387904,
    ];
    let answers = vec![
        1.0,
        0.9987954562051724,
        //  0.99879545620517239271
        0.9951847266721969,
        //  0.99518472667219688624
        0.989176509964781,
        //  0.98917650996478097345
        0.9807852804032304,
        //  0.98078528040323044912
        0.970031253194544,
        //  0.97003125319454399260
        0.9569403357322088,
        //  0.95694033573220886493
        0.9415440651830208,
        //  0.94154406518302077841
        0.9238795325112867,
        //  0.92387953251128675612
        0.9039892931234433,
        //  0.90398929312344333158
        0.881921264348355,
        //  0.88192126434835502971
        0.8577286100002721,
        //  0.85772861000027206990
        0.8314696123025452,
        //  0.83146961230254523707
        0.8032075314806449,
        //  0.80320753148064490980
        0.773010453362737,
        //  0.77301045336273696081
        0.7409511253549591,
        //  0.74095112535495909117
        0.7071067811865476,
        //  0.70710678118654752440
        0.6715589548470184,
        //  0.67155895484701840062
        0.6343932841636455,
        //  0.63439328416364549821
        0.5956993044924334,
        //  0.59569930449243334346
        0.5555702330196022,
        //  0.55557023301960222474
        0.5141027441932218,
        //  0.51410274419322172659
        0.47139673682599764,
        //  0.47139673682599764855
        0.4275550934302821,
        //  0.42755509343028209432
        0.3826834323650898,
        //  0.38268343236508977172
        0.33688985339222005,
        //  0.33688985339222005068
        0.2902846772544624,
        //  0.29028467725446236763
        0.2429801799032639,
        //  0.24298017990326388994
        0.19509032201612828,
        //  0.19509032201612826784 (exactly equal) 7 pi / 16
        0.14673047445536175,
        //  0.14673047445536175165
        0.0980171403295606,
        //  0.09801714032956060199 15 / 32
        0.049067674327418015,
        //  0.04906767432741801425
        0.0,
    ];
    for i in 0..33 {
        assert_eq!(cos(questions[i]), answers[i]);
    }
}
