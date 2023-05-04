// (C) Daniel Strano and the Qrack contributors 2017-2023. All rights reserved.
//
// Use of this source code is governed by an MIT-style license that can be
// found in the LICENSE file or at https://opensource.org/licenses/MIT.

use neuron_activation_fn::NeuronActivationFn;
use qrack_error::QrackError;
use qrack_simulator::QrackSimulator;
use qrack_system;

pub struct QrackNeuron<'a> {
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
    //     nid(i64): Corresponding neuron id.
    //     simulator(QrackSimulator): Simulator for neuron
    //     controls(Vec<u64>): Neuron input qubits
    //     target(u64): Neuron output qubit
    //     activation_fn(NeuronActivationFn): Neuron activation function
    //     alpha(f64): Activation function parameter (if used)
    //     tolerance(f64): Rounding tolerance
    //     amp_count(u64): Count of amplitudes in training space
    nid: u64,
    simulator: &'a QrackSimulator,
    controls: Vec<u64>,
    target: u64,
    activation_fn: NeuronActivationFn,
    alpha: f64,
    tolerance: f64,
    amp_count: u64
}

impl Clone for QrackNeuron<'_> {
    fn clone(&self) -> Self {
        let nid;
        unsafe {
            nid = qrack_system::clone_qneuron(self.nid);
        }
        Self{
            nid,
            simulator: self.simulator,
            controls: self.controls.clone(),
            target: self.target,
            activation_fn: self.activation_fn.clone(),
            alpha: self.alpha,
            tolerance: self.tolerance,
            amp_count: self.amp_count
        }
    }
}

impl Drop for QrackNeuron<'_> {
    fn drop(&mut self) {
        unsafe {
            qrack_system::destroy_qneuron(self.nid);
        }
    }
}

impl QrackNeuron<'_> {
    // private functions
    fn get_error(&self) -> i32 {
        unsafe {
            qrack_system::get_error(self.nid)
        }
    }
    fn check_error(&self) -> Result<(), QrackError> {
        if self.get_error() != 0 {
            return Err(QrackError{});
        }
        return Ok(());
    }
    
    // constructors
    pub fn new<'a: 'b, 'b>(
        sim: &'a QrackSimulator,
        ctrls: Vec<u64>,
        trgt: u64,
        act_fn: NeuronActivationFn,
        a: f64,
        t: f64
   ) -> Result<QrackNeuron<'b>, QrackError> {
        let nid;
        let mut _controls = ctrls.to_vec();
        unsafe {
            nid = qrack_system::init_qneuron(
                sim.get_sid(),
                _controls.len() as u64,
                _controls.as_mut_ptr() as *mut u64,
                trgt,
                act_fn.clone() as u64,
                a,
                t
            );
            if qrack_system::get_error(nid) != 0 {
                return Err(QrackError{});
            }
        }
        let amp_cnt = 1 << (ctrls.len() + 1);
        Ok(QrackNeuron{
            nid,
            simulator: sim,
            controls: ctrls,
            target: trgt,
            activation_fn: act_fn,
            alpha: a,
            tolerance: t,
            amp_count: amp_cnt
        })
    }

    pub fn set_angles(&self, a: Vec<f32>) -> Result<(), QrackError> {
        // Directly sets the neuron parameters.
        //
        // Set all synaptic parameters of the neuron directly, by a list
        // enumerated over the integer permutations of input qubits.
        //
        // Args:
        //     a(Vec<f64>): List of input permutation angles
        //
        // Raises:
        //     RuntimeError: QrackNeuron C++ library raised an exception.

        let mut _a = a.to_vec();
        unsafe {
            qrack_system::set_qneuron_angles(self.nid, _a.as_mut_ptr());
        }
        self.check_error()
    }

    pub fn get_angles(&self) -> Result<Vec<f32>, QrackError> {
        // Directly gets the neuron parameters.
        //
        // Get all synaptic parameters of the neuron directly, as a list
        // enumerated over the integer permutations of input qubits.
        //
        // Raises:
        //     RuntimeError: QrackNeuron C++ library raised an exception.

        let mut result = vec![0.0;self.amp_count as usize];
        unsafe {
            qrack_system::get_qneuron_angles(self.nid, result.as_mut_ptr())
        }
        if self.get_error() != 0 {
            return Err(QrackError{});
        }
        Ok(result)
    }

    pub fn set_alpha(mut self, a: f64) -> Result<(), QrackError> {
        // Set the neuron 'alpha' parameter.
        //
        // To enable nonlinear activation, `QrackNeuron` has an 'alpha'
        // parameter that is applied as a power to its angles, before
        // learning and prediction. This makes the activation function
        // sharper (or less sharp).
        //
        // Raises:
        //     RuntimeError: QrackNeuron C++ library raised an exception.

        self.alpha = a;
        unsafe {
            qrack_system::set_qneuron_alpha(self.nid, a);
        }
        self.check_error()
    }

    pub fn set_activation_fn(mut self, f: NeuronActivationFn) -> Result<(), QrackError> {
        // Sets the activation function of this QrackNeuron
        //
        // Nonlinear activation functions can be important to neural net
        // applications, like DNN. The available activation functions are
        // enumerated in `NeuronActivationFn`. 
        //
        // Raises:
        //     RuntimeError: QrackNeuron C++ library raised an exception.

        self.activation_fn = f.clone();
        unsafe {
            qrack_system::set_qneuron_activation_fn(self.nid, f as u64);
        }
        self.check_error()
    }

    pub fn predict(&self, e: bool, r: bool) -> Result<f64, QrackError> {
        // Predict based on training
        //
        // "Predict" the anticipated output, based on input and training.
        // By default, "predict()" will initialize the output qubit as by
        // resetting to |0> and then acting a Hadamard gate. From that
        // state, the method amends the output qubit upon the basis of
        // the state of its input qubits, applying a rotation around
        // Pauli Y axis according to the angle learned for the input.
        //
        // Args:
        //     e(bool): If False, predict the opposite
        //     r(bool): If True, start by resetting the output to 50/50
        //
        // Raises:
        //     RuntimeError: QrackNeuron C++ library raised an exception.

        let result:f64;
        unsafe {
            result = qrack_system::qneuron_predict(self.nid, e, r);
        }
        if self.get_error() != 0 {
            return Err(QrackError{});
        }
        Ok(result)
    }

    pub fn unpredict(&self, e: bool) -> Result<f64, QrackError> {
        // Uncompute a prediction
        //
        // Uncompute a 'prediction' of the anticipated output, based on
        // input and training.
        //
        // Args:
        //     e(bool): If False, unpredict the opposite
        //
        // Raises:
        //     RuntimeError: QrackNeuron C++ library raised an exception.

        let result:f64;
        unsafe {
            result = qrack_system::qneuron_unpredict(self.nid, e);
        }
        if self.get_error() != 0 {
            return Err(QrackError{});
        }
        Ok(result)
    }

    pub fn learn_cycle(&self, e: bool) -> Result<(), QrackError> {
        // Run a learning cycle
        //
        // A learning cycle consists of predicting a result, saving the
        // classical outcome, and uncomputing the prediction.
        //
        // Args:
        //     e(bool): If False, predict the opposite
        //
        // Raises:
        //     RuntimeError: QrackNeuron C++ library raised an exception.

        unsafe {
            qrack_system::qneuron_learn_cycle(self.nid, e);
        }
        self.check_error()
    }

    pub fn learn(&self, eta: f64, e: bool, r: bool) -> Result<(), QrackError> {
        // Learn from current qubit state
        //
        // "Learn" to associate current inputs with output. Based on
        // input qubit states and volatility 'eta,' the input state
        // synaptic parameter is updated to prefer the "e" ("expected")
        // output.
        //
        // Args:
        //     eta(double): Training volatility, 0 to 1
        //     e(bool): If False, predict the opposite
        //     r(bool): If True, start by resetting the output to 50/50
        //
        // Raises:
        //     RuntimeError: QrackNeuron C++ library raised an exception.

        unsafe {
            qrack_system::qneuron_learn(self.nid, eta, e, r);
        }
        self.check_error()
    }

    pub fn learn_permutation(&self, eta: f64, e: bool, r: bool) -> Result<(), QrackError> {
        // Learn from current classical state
        //
        // Learn to associate current inputs with output, under the
        // assumption that the inputs and outputs are "classical."
        // Based on input qubit states and volatility 'eta,' the input
        // state angle is updated to prefer the "e" ("expected") output.
        //
        // Args:
        //     eta(double): Training volatility, 0 to 1
        //     e(bool): If False, predict the opposite
        //     r(bool): If True, start by resetting the output to 50/50
        //
        // Raises:
        //     RuntimeError: QrackNeuron C++ library raised an exception.

        unsafe {
            qrack_system::qneuron_learn_permutation(self.nid, eta, e, r);
        }
        self.check_error()
    }
}
