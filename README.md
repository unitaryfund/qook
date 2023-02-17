# qook
Pure language standard, "safe" Rust bindings for the pure C++11/OpenCL Qrack quantum computer simulator library

(**Qook** is just pure Qrack.)

To use this package, it's helpful to be familiar with [unitaryfund/qrack](https://github.com/unitaryfund/qrack). Users gain **much** more control over options by building **unitaryfund/qrack** and [unitaryfund/qook](https://github.com/unitaryfund/qook) from source. For advanced users, building from source is the intended primary method of **Qook** package distribution and use.

Import and instantiate [`QrackSimulator`](https://github.com/unitaryfund/qook/blob/main/src/qrack_simulator.rs) instances. This simulator can perform arbitrary single qubit and controlled-single-qubit gates, as well as other specific gates like `SWAP`.

Any 2x2 bit operator matrix is represented by a list of 4 `complex` floating point numbers, in [**row-major order**](https://en.wikipedia.org/wiki/Row-_and_column-major_order).

Single and array "`b`" parameters represent [**Pauli operator bases**](https://en.wikipedia.org/wiki/Pauli_matrices). They are specified according to the enumeration of the [`Pauli`](https://github.com/unitaryfund/qook/blob/main/src/pauli.rs) class.

`MC[x]` and `MAC[x]` methods are controlled single bit gates, with as many control qubits as you specify via the Rust vector `c` argument. `MCX` is multiply-controlled Pauli X, and `MACX` is "anti-"controlled Pauli X, i.e. "anti-control" activates the gate if all control bits are specifically **off**, as opposed to **on**.

The package installation directory contains a `qrack_cl_precompile` folder with executables for supported platforms, to compile OpenCL kernels once, beforehand, avoiding the need to recompile "just-in-time" every time that you load this package in a Python interpreter. If you no longer want to use precompiled kernels, or if precompilation fails, just delete the `~/.qrack` directory, or the equivalent `.qrack` sub-directory in the user home folder of your operating system.

To load the required **unitaryfund/qrack** libraries from a different location, set the `QOOK_SHARED_LIB_PATH` environment variable. If an absolute path is set for this variable, it will be treated as an absolute path. If a relative path is used, the **qook** package `qrack_system` directory is the root.

The API is meant to exactly mirror (Python-based) PyQrack. See [https://pyqrack.readthedocs.io/en/latest/](https://pyqrack.readthedocs.io/en/latest/) for an API reference.

Please feel welcome to open an issue, if you'd like help. 😃

**For their work on PyQrack, special thanks go to Zeeshan Ahmed, for bug fixes and design suggestions, Ashish Panigrahi, for documentation and design suggestions, WingCode, for documentation, and to the broader community of Qrack contributors, for years of happy Qracking! You rock!**
