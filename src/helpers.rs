use num::{bigint::BigInt, BigRational, FromPrimitive, ToPrimitive};
use num_integer::Integer;
use num_traits::identities::Zero;
use num_traits::sign::Signed;

/// Convert a BigRational to an f64
pub fn to_f64(n: BigRational) -> Option<f64> {
    assert_eq!(
        std::f64::RADIX,
        2,
        "only floating point implementations with radix 2 are supported"
    );

    // Inclusive upper bound to the range of exactly-representable ints in an f64.
    const MAX_EXACT_INT: u64 = 1u64 << std::f64::MANTISSA_DIGITS;

    let numer: BigInt = n.numer().abs();
    let denom: BigInt = n.denom().abs();
    let sign: BigInt = n.numer().signum() * n.denom().signum();
    let flo_sign = if sign.is_negative() { -1.0 } else { 1.0 };

    if numer.is_zero() {
        return Some(0.0 * flo_sign);
    }

    // Fast track: both sides can losslessly be converted to f64s. In this case, letting the
    // FPU do the job is faster and easier. In any other case, converting to f64s may lead
    // to an inexact result: https://stackoverflow.com/questions/56641441/.
    if let (Some(n), Some(d)) = (numer.to_u64(), denom.to_u64()) {
        if n <= MAX_EXACT_INT && d <= MAX_EXACT_INT {
            return Some(flo_sign * n.to_f64().unwrap() / d.to_f64().unwrap());
        }
    }

    // Otherwise, the goal is to obtain a quotient with at least 55 bits. 53 of these bits will
    // be used as the mantissa of the resulting float, and the remaining two are for rounding.
    // There's an error of up to 1 on the number of resulting bits, so we may get either 55 or
    // 56 bits.
    let mut numer = numer;
    let mut denom = denom;
    let (is_diff_positive, absolute_diff) = match numer.bits().checked_sub(denom.bits()) {
        Some(diff) => (true, diff),
        None => (false, denom.bits() - numer.bits()),
    };

    // Filter out overflows and underflows. After this step, the signed difference fits in an
    // isize.
    if is_diff_positive && absolute_diff > std::f64::MAX_EXP as usize {
        return Some(std::f64::INFINITY * flo_sign);
    }
    if !is_diff_positive
        && absolute_diff > -std::f64::MIN_EXP as usize + std::f64::MANTISSA_DIGITS as usize + 1
    {
        return Some(0.0 * flo_sign);
    }
    let diff = if is_diff_positive {
        absolute_diff.to_isize().unwrap()
    } else {
        -absolute_diff.to_isize().unwrap()
    };

    // Shift is chosen so that the quotient will have 55 or 56 bits. The exception is if the
    // quotient is going to be subnormal, in which case it may have fewer bits.
    let shift: isize =
        std::cmp::max(diff, std::f64::MIN_EXP as isize) - std::f64::MANTISSA_DIGITS as isize - 2;
    if shift >= 0 {
        denom <<= shift as usize
    } else {
        numer <<= -shift as usize
    };

    let (quotient, remainder) = numer.div_rem(&denom);

    // This is guaranteed to fit since we've set up quotient to be at most 56 bits.
    let mut quotient = quotient.to_u64().unwrap();
    let n_rounding_bits = {
        let quotient_bits = 64 - quotient.leading_zeros() as isize;
        let subnormal_bits = std::f64::MIN_EXP as isize - shift;
        std::cmp::max(quotient_bits, subnormal_bits) - std::f64::MANTISSA_DIGITS as isize
    } as usize;
    debug_assert!(n_rounding_bits == 2 || n_rounding_bits == 3);
    let rounding_bit_mask = (1u64 << n_rounding_bits) - 1;

    // Round to 53 bits with round-to-even. For rounding, we need to take into account both
    // our rounding bits and the division's remainder.
    let ls_bit = quotient & (1u64 << n_rounding_bits) != 0;
    let ms_rounding_bit = quotient & (1u64 << (n_rounding_bits - 1)) != 0;
    let ls_rounding_bits = quotient & (rounding_bit_mask >> 1) != 0;
    if ms_rounding_bit && (ls_bit || ls_rounding_bits || !remainder.is_zero()) {
        quotient += 1u64 << n_rounding_bits;
    }
    quotient &= !rounding_bit_mask;

    // The quotient is guaranteed to be exactly representable as it's now 53 bits + 2 or 3
    // trailing zeros, so there is no risk of a rounding error here.
    let q_float = quotient as f64;
    Some(q_float * f64::exp2(shift as f64) * flo_sign)
}

/// Convert an i64 into a BigRational
pub fn bi(i: i64) -> BigRational {
    BigRational::from_integer(FromPrimitive::from_i64(i).unwrap())
}
