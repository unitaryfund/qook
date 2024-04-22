// (C) Daniel Strano and the Qrack contributors 2017-2023. All rights reserved.
//
// Use of this source code is governed by an MIT-style license that can be
// found in the LICENSE file or at https://opensource.org/licenses/MIT.

use std::ffi::CString;

use qrack_error::QrackError;
use qrack_simulator::QrackSimulator;
use qrack_system;

pub struct QrackCircuit {
    // Class that exposes the QNeuron class of Qrack
    //
    // This model of a "quantum neuron" is based on the concept of a "uniformly controlled"
    // rotation of a single output qubit around the Pauli Y axis, and has been developed by
    // others. In our case, the primary relevant gate could also be called a
    // single-qubit-target multiplexer.
    //
    // (See https://arxiv.org/abs/quant-ph/0407010 for an introduction to "uniformly controlled
    // gates.)
    //
    // QrackNeuron is meant to be interchangeable with a single classical neuron, as in
    // conventional neural net software. It differs from classical neurons in conventional
    // neural nets, in that the "synaptic cleft" is modelled as a single qubit. Hence, this
    // neuron can train and predict in superposition.
    //
    // Attributes:
    //     cid(u64): Corresponding circuit id.
    cid: u64
}

impl Clone for QrackCircuit {
    fn clone(&self) -> Self {
        let cid;
        unsafe {
            cid = qrack_system::init_qcircuit_clone(self.cid);
        }
        Self{
            cid
        }
    }
}

impl Drop for QrackCircuit {
    fn drop(&mut self) {
        unsafe {
            qrack_system::destroy_qcircuit(self.cid);
        }
    }
}

impl QrackCircuit {
    // constructors
    pub fn new() -> Self {
        let cid;
        unsafe {
            cid = qrack_system::init_qcircuit(false, false);
        }
        Self{cid}
    }

    pub fn get_qubit_count(&self) -> u64 {
        // Get count of qubits in circuit
        unsafe {
            qrack_system::get_qcircuit_qubit_count(self.cid)
        }
    }

    pub fn inverse(&self) -> QrackCircuit {
        let cid;
        unsafe {
            cid = qrack_system::qcircuit_inverse(self.cid);
        }
        Self{
            cid
        }
    }

    pub fn past_light_cone(&self, q: Vec<u64>) -> QrackCircuit {
        let cid;
        let mut _q = q.to_vec();
        unsafe {
            cid = qrack_system::qcircuit_past_light_cone(self.cid, _q.len() as u64, _q.as_mut_ptr());
        }
        Self{
            cid
        }
    }

    pub fn swap(&self, q1: u64, q2: u64) -> () {
        // Add a 'Swap' gate to the circuit
        //
        // Args:
        //     q1: qubit index #1
        //     q2: qubit index #2
        unsafe {
            qrack_system::qcircuit_swap(self.cid, q1, q2)
        }
    }

    pub fn mtrx(&self, m: &[f64;8], q: u64) -> () {
        // Operation from matrix.
        //
        // Applies arbitrary operation defined by the given matrix.
        //
        // Args:
        //     m(&[f64;8]): row-major complex list representing the operator.
        //     q(u64): the qubit number on which the gate is applied to.
        let mut _m = [m[0], m[1], m[2], m[3], m[4], m[5], m[6], m[7]];
        unsafe {
            qrack_system::qcircuit_append_1qb(self.cid, _m.as_mut_ptr(), q)
        }
    }

    pub fn ucmtrx(&self, c: Vec<u64>, m: &[f64;8], q: u64, p: u64) -> () {
        // Multi-controlled arbitrary operator with arbitrary controls
        //
        // If all control qubits match 'p' permutation by bit order, then the arbitrary
        // operation by parameters is applied to the target qubit.
        //
        // Args:
        //     c(Vec<u64>): list of controlled qubits
        //     m(&[f64;8]): row-major complex list representing the operator.
        //     q(u64): target qubit
        //     p(u64): permutation of list of control qubits
        let mut _m = [m[0], m[1], m[2], m[3], m[4], m[5], m[6], m[7]];
        let mut _c = c.to_vec();
        unsafe {
            qrack_system::qcircuit_append_mc(self.cid, _m.as_mut_ptr(), _c.len() as u64, _c.as_mut_ptr(), q, p);
        }
    }

    pub fn run(&self, qsim: &QrackSimulator) -> Result<(), QrackError> {
        // Run circuit on simulator
        //
        // Run the encoded circuit on a specific simulator. The
        // result will remain in this simulator.
        //
        // Args:
        //     qsim(&QrackSimulator): QrackSimulator on which to run circuit
        // Raises:
        //     RuntimeError: QrackCircuit raised an exception.
        unsafe {
            qrack_system::qcircuit_run(self.cid, qsim.get_sid())
        }
        qsim.check_error()
    }

    pub fn out_to_file(&self, filename: &str) -> () {
        // Output optimized circuit to file
        //
        // Outputs the (optimized) circuit to a file named
        // according to the "filename" parameter.
        //
        // Args:
        //     filename: Name of file
        unsafe {
            qrack_system::qcircuit_out_to_file(self.cid, CString::new(filename).unwrap().into_bytes_with_nul().as_mut_ptr() as *mut i8)
        }
    }

    pub fn in_from_file(&self, filename: &str) -> () {
        // Read in optimized circuit from file
        //
        // Reads in an (optimized) circuit from a file named
        // according to the "filename" parameter.
        //
        // Args:
        //     filename: Name of file
        unsafe {
            qrack_system::qcircuit_in_from_file(self.cid, CString::new(filename).unwrap().into_bytes_with_nul().as_mut_ptr() as *mut i8)
        }
    }
}
