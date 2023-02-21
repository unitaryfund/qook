extern crate qook;
use qook::qrack_simulator;
use std::f64::consts::FRAC_1_SQRT_2;
fn main() {
    // *ALL* QrackSimulator methods can report exception in the underlying C++, hence "unwrap()."
    let q_reg = qrack_simulator::QrackSimulator::new(12).unwrap();
    for i in 0..11 {
        q_reg.h(i).unwrap();
    }
    // Test mtrx array parameter
    q_reg.mtrx(&[FRAC_1_SQRT_2, 0.0, FRAC_1_SQRT_2, 0.0, FRAC_1_SQRT_2, 0.0, -FRAC_1_SQRT_2, 0.0], 11).unwrap();
    println!("Probability of highest bit for |1> state: {}", q_reg.prob(11).unwrap());
    println!("Random number less than 2^12: {}", q_reg.m_all().unwrap());
}
