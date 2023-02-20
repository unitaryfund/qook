extern crate qook;
use qook::qrack_simulator;
fn main() {
    let q_reg = qrack_simulator::QrackSimulator::new(12).unwrap();
    for i in 0..12 {
        q_reg.h(i).unwrap();
    }
    println!("Random number less than 2^12: {}", q_reg.m_all().unwrap());
}
