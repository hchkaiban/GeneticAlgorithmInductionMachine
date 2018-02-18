# GeneticAlgorithmInductionMachine
Genetic Algorithm for off-line identification of induction machine's parameters

Accurate parameters are critical for the control of the machine. Classical methods like open load or blocked rotor tests are intrusive.
The proposed method based on a genetic algorithm only requires standard measurements as input in stationary (permanent regime) operation.

- GA_IMparamIdent_ISSN 1330-3651_Paper.pdf: Theory and equations
- IM_GA_model.slx: Matlab Simulink simulation model
- IM_Param.m: real machine parameters for validation
- IM_AG_Fit_Sfunc.c: S-function in C code of the genetic algorithm
- IMSimulation_data.ods: simulation measurements with required inputs

Results:
Proof of concept validated.
Accuracy shall be be improved with more accurate measurements (in particular cos Phi) and further operating points. 

![Identified Parameters](https://github.com/hchkaiban/GeneticAlgorithmInductionMachine/blob/master/AG_Scope.png)
 
