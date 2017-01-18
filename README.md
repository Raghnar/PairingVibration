# PairingVibration
Code for the calculation of Pairing Vibration with phenomenological Random Phase Approximation (RPA). Eventually reads the Atomic Mass Evaluation Database for an automatic usage.

Code is not (and never will be) complete, releasing as process, not as product. For suggestions, contact me.

Andrea Idini, Univ. of Surrey


##How to Run:
Use Makefile to compile.
Change nucleus and proton/neutron vibration in Input\_WS.in. 

If you use the option 2, the program will scan the Atomic Mass Evaluation 2012 for masses and fix nuclear levels with a standard Bohr and Mottelson Wood Saxon, enjoy the rest.

If you use the option 1, you can specify levels and other inputs following the format in the example.


##Output:
Factor.dat and RPA.dat, for X and Y RPA amplitudes and energy of the phonons (in two formats).
Dispersion.dat, dispersion relation of RPA.


##Physics:
If you don't have much idea what a Pairing Vibration and what a Random Phase Approximation (RPA) is, have a torough look at:

- Physics of Atomic Nuclei 77, 941: http://link.springer.com/article/10.1134/S106377881407014X
- Physica Scripta, Volume 91, Number 6:  http://iopscience.iop.org/article/10.1088/0031-8949/91/6/063012/meta
- Phys. Rev. C 87, 054321, http://journals.aps.org/prc/abstract/10.1103/PhysRevC.87.054321

and refs. therein.


###License
The Software is provided "as is" according to the GNUv3 License. 
I do not offer support.
However, feel free to drop me an email with the 
scientific/didactic/business/life problem you want to discuss, I usually answer.


Read file LICENSE for more information on GNUv3.
