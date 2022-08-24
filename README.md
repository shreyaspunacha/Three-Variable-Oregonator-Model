# Three Variable Oregonator Model

## Introduction

This project contains solutions to three-variable Oregonator model written using Finite Difference Explicit scheme and C programming language. The code is used to generate pinned spiral waves in Belousov-Zhabotinsky chemical reaction and eliminate them using low-voltage electric field.

The three-variable Oregonator model is as follows:

$$
\frac{\partial u}{\partial t} = \frac{1}{\epsilon} (qw - uw + u - u^{2}) + D_{u} \nabla^{2} u
$$

$$
\frac{\partial u}{\partial t} = u - v + D_{v} \nabla^{2} v) - K_{v} E \frac{\partial v}{\partial x}
$$

$$
\frac{\partial w}{\partial t} = \frac{1}{\epsilon^{`}} (-qw - uw + fv) + D_{w} \nabla^{2} w - K_{w} E \frac{\partial w}{\partial x}
$$


The above model equations and the parameters are chosen from the paper titled "Forced parallel drift of spiral waves in the Belousov-Zhabotinsky reaction" by "Bernd Schmidt and Stefan C. Müller" VOLUME 55, NUMBER 4, APRIL 1997, Physical Review E.

# Instructions To Run The Code.

From the Terminal, navigate to the folder containing the codes and type

* make
* ./run
* python Plot.py

The make file is works on Mac and Linux systems. 
