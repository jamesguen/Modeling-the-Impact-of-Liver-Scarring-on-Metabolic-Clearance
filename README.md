# Modeling-the-Impact-of-Liver-Scarring-on-Metabolic-Clearance
A Physiologically Based Pharmacokinetic (PBPK) model developed in MATLAB to simulate and quantify the non-linear clearance kinetics of ethanol and its highly toxic intermediate, acetaldehyde, across varying degrees of liver health.

Developed as part of the BME 121: Problem-Based Learning (PBL) curriculum within the Department of Biomedical Engineering at the University of California, Irvine.

## Overview
This script simulates a multi-compartment, non-steady-state mass balance model of human alcohol metabolism. It evaluates how progressive liver damage (fibrosis/cirrhosis) slows metabolic pathways, causing a non-proportional, prolonged accumulation of toxic metabolic byproducts in both the liver tissue and the circulatory system.

When alcohol is consumed, it is processed primarily in the liver. Ethanol is converted into highly toxic acetaldehyde via alcohol dehydrogenase (ADH), which is further metabolized into the less harmful acetate via aldehyde dehydrogenase (ALDH).

To model this system, the script solves a coupled non-linear system of four ordinary differential equations (ODEs) using MATLAB's built-in ode45 solver. To simulate varying stages of liver disease, a dimensionless scaling parameter, **$\alpha$**, is introduced:

**$\alpha$ = 1.0:** Healthy liver tissue with 100% metabolic capacity.
**$\alpha$ < 1.0:** Progressively scarred or damaged liver tissue, representing a fractional reduction in normal enzyme efficiency (Vmax).

## Model Architecture

The code models a well-mixed, lumped system across two main physiological pools (Blood and Liver) across four distinct state variables:

* **$E_B(t)$**: Blood Ethanol Concentration (mM)
* **$E_L(t)$**: Liver Ethanol Concentration (mM)
* **$A_B(t)$**: Blood Acetaldehyde Concentration (mM)
* **$A_L(t)$**: Liver Acetaldehyde Concentration (mM)

## Kinetic Equations

The reaction rates for ethanol metabolism ($v_{ADH}$) and acetaldehyde metabolism ($v_{ALDH}$) utilize non-steady-state Michaelis-Menten kinetics scaled by the liver health fraction ($\alpha$):

$$v_{ADH} = \alpha \cdot V_{max,E} \cdot \left(\frac{E_L}{K_{m,E} + E_L}\right)$$

$$v_{ALDH} = \alpha \cdot V_{max,A} \cdot \left(\frac{A_L}{K_{m,A} + A_L}\right)$$

The differential dynamics governing the entire system are solved simultaneously using the following system of equations:

$$\frac{dE_B}{dt} = \frac{Q}{V_B}(E_L - E_B) + \frac{u(t)}{V_B}$$

$$\frac{dE_L}{dt} = \frac{Q}{V_L}(E_B - E_L) - v_{ADH}$$

$$\frac{dA_B}{dt} = \frac{Q}{V_B}(A_L - A_B)$$

$$\frac{dA_L}{dt} = \frac{Q}{V_L}(A_B - A_L) + v_{ADH} - v_{ALDH}$$

Where $u(t)$ represents the alcohol input function over a specified drinking window.

## Getting Started

Prerequisites
MATLAB (R2016b or newer recommended)
Optimization Toolbox / Global Optimization Toolbox (Not strictly required, standard ode45 is built into base MATLAB)

## Running the Simulation
Clone this repository or copy the code file directly into your MATLAB working directory.
Ensure you have a helper function for hex2rgb in your path, or replace instances of hex2rgb('#hex') with standard MATLAB RGB triplets (e.g., [0.1, 0.1, 0.5]) if the custom plotting colors throw an error.
Execute the script from the MATLAB Command Window:
***run('liver_clearance_sim.m')***

## Expected Outputs
The script initializes a custom dark-themed GUI dashboard featuring a 2x2 subplot visualizer grid. It loops through 5 distinct $\alpha$ variations ($\alpha$ = 1.0, 0.9, 0.75, 0.5, 0.3) and color-codes each curve accordingly to plot concentration profiles over time:

Top Left: Blood Ethanol Concentration E_B(t) over a 6-hour window.
Top Right: Liver Ethanol Concentration E_L(t).
Bottom Left: Blood Acetaldehyde Concentration A_B(t).
Bottom Right: Liver Acetaldehyde Concentration A_L(t).

Key Analytical Takeaway: As $\alpha$ approaches 0.3 (advanced cirrhosis/severe scarring), the metabolism shifts into zero-order saturation limits. The liver curves begin to directly mirror the body blood pool, showing an inability to clear metabolic volumes effectively, resulting in a dangerous time-extension of toxic acetaldehyde exposure.

## Credits
Code Implementation: James Guentert
Project Team Members: Abigail Marciniak, James Guentert, Kiyan Nguyen Abel, Satya Govindarajan
Course Instructor: Professor Anna Grosberg, University of California, Irvine (UCI)
