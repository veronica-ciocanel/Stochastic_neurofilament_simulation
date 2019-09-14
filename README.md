# Stochastic_neurofilament_simulation

Sample Matlab code for simulating a stochastic model of neurofilament transport through internodes and nodes of Ranvier in axonal constrictions.

This code accompanies a preprint on "A Mechanism for Neurofilament Transport Acceleration through Nodes of Ranvier". 

The main code nf_stochastic.m sets up kinetic and spatial parameters and carries out a sample simulation of neurofilament kinetics through an axonal segment and can be easily modified to simulate a nodal constriction. The function update_nfs.m carries out the update in each neurofilament's location and state at each time interval.

Other functions include initialize.m, which generates the initial spatial locations of the neurofilament centers and their kinetic states; content.m, which calculates the total neurofilament content in space; and contenton.m, which calculates the on-track neurofilament content.
