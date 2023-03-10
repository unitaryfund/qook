// (C) Daniel Strano and the Qrack contributors 2017-2023. All rights reserved.
//
// Pauli operators are specified for "b" (or "basis") parameters.
//
// Use of this source code is governed by an MIT-style license that can be
// found in the LICENSE file or at https://opensource.org/licenses/MIT.

#[derive(Debug)]
pub struct QrackError;

impl std::fmt::Display for QrackError {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        write!(f, "QrackSimulator C++ library raised exception, (or Vector arguments have unexpected lengths)")
    }
}
