# NLOPT (SLSQP-Only) with Integrated Parallel QR

This repository contains an SLSQP-only version of the NLOPT library, enhanced with a parallel QR factorization module. All necessary code for SLSQP and the integrated parallel QR factorization is located under `NLOPT/`; there is no separate build required for the `ParallelQR` portion.

> **Note**
> Because this repository only includes the SLSQP part of NLOPT, you will only be able to compile and use the SLSQP algorithm from the NLOPT suite. If you need other algorithms, you should use the complete NLOPT repository or library.

---

## Contents

-   **NLOPT/**
    Contains the modified SLSQP algorithm, adapted to use parallel QR factorization. It also includes the CMake configuration that automatically integrates the Parallel QR routines.

-   **ParallelQR/**
    The code for parallel QR factorization is stored here. **No separate build step is needed**; it is automatically compiled through the NLOPT CMake configuration.

---
## Building and Installing

1.  **Clone this repository:**

    ```bash
    git clone https://github.com/PDCRL/ParSQP.git
    cd ParSQP
    ```

2.  **Configure and Build NLOPT :**

    ```bash
    cd NLOPT
    mkdir build && cd build
    cmake -DCMAKE_INSTALL_PREFIX=/path/to/nlopt_install ..
    make
    make install
    ```
    Replace `/path/to/nlopt_install` with your desired installation directory. Make sure the folder in the path is already exist. 

    After installation, you should find `include`, `lib`, and possibly `bin` folders under `nlopt_install`.

    For more detailed information, visit https://nlopt.readthedocs.io/en/latest/NLopt_Installation/

4.  **ParallelQR Directory**

    The `ParallelQR` folder **does not require a separate build step**. All necessary code is already included in the main CMake configuration of `NLOPT`.

---
