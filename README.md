# Aeroservoelasticity - Politecnico di Milano

This repository contains the MATLAB implementations for the Aeroservoelasticity course workshops at Politecnico di Milano. The projects focus on flutter prediction, stability analysis of rotary wings, and non-linear aeroelastic phenomena.

### Workshop 1: Goland Wing Flutter & Control Effectiveness
* **Objective:** Developed a solver for the flutter problem in the Goland Wing model.
* **Methods:** Implemented standard *k* and *p-k* methods to determine the flutter velocity and damping loops.
* **Extended Analysis:** Studied the aeroelastic influence of adding a rigid control surface, analyzing the shifts in flutter velocity, divergence, and control reversal boundaries.

### Workshop 2: Puma Helicopter Blade Stability
* **Objective:** Aeroelastic study of the pitch-flap coupling in a Puma helicopter blade.
* **Analysis:** Generated comprehensive stability maps to analyze the effects of varying the pitch link stiffness and the blade's Center of Gravity (CG) offset. Evaluated the impact of modal truncation (varying the number of modes) and different aerodynamic assumptions (e.g., vectorial Theodorsen function).
* **Implementation:** The MATLAB code is highly optimized, leveraging orthogonal Legendre polynomials, pre-computed structural matrices, and parallel computing (`parfor`) to drastically reduce the execution time of heavy iterative parameter sweeps.

### Workshop 3: Limit Cycle Oscillations (LCO)
* **Team:** Co-developed with Alberto Rivero García and Alejandro Rivera Míguez.
* **Objective:** Non-linear aeroelastic analysis of the Goland Wing model. 
* **Focus:** Investigated the Limit Cycle Oscillations (LCO) triggered by a structural non-linearity, specifically simulating the presence of free-play in the control surface hinge.


