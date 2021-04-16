# **Dimsym**

## Geometric and algebraic techniques for differential equations

### With modelling applications

## Symmetry Determination and Linear Differential Equation Package

**Dimsym** is a program primarily for the determination of symmetries of differential equations. It also can be used to compute symmetries of distributions of vector fields or differential forms on finite dimensional manifolds, symmetries of geometric objects (e.g. isometries), and also to solve linear partial differential equations.

To use its primary function the user specifies a system of ordinary and/or partial differential equations and the type of symmetry to be found (Lie point, Lie-Backlund or some user-provided ansatz).

**Dimsym** then produces the corresponding determining equations - a system of linear partial differential equations for the generator of the generic symmetry. It proceeds to solve these equations, reporting any special conditions required to produce a solution. Finally, **Dimsym** gives the generators of the symmetry group, which may of course be infinite dimensional.

The program allows the user to compute Lie brackets, vector derivatives and so on, and it has an interface with the **REDUCE** package _EXCALC_ so that all the machinery of calculus on manifolds can be utilised from within the program. Its use can be interactive or batch, and there are extensive tracing options.

If you don't have **REDUCE**, or if you work on a stand alone IBM compatible, you may wish to use [Alan Head's LIE program](http://archives.math.utk.edu/software/msdos/adv.diff.equations/lie/.html).

## **Dimsym** and related files

- [Source Code](https://github.com/reduce-algebra/dimsym/tree/master/src/)
- [Examples](https://github.com/reduce-algebra/dimsym/tree/master/examples/)
- [Manual](https://github.com/reduce-algebra/dimsym/tree/master/doc/)

  - It is **highly recommended** that you read the manual and work through the example files.

## Getting Started

- **Read the** **_Dimsym_** **user's manual.**

Given that you want to find symmetries of your differential equation, it is assumed that you know only the basics of **REDUCE**, such as how to run **REDUCE** on your machine, and how to compile and load **Dimsym**, which may vary between different implementations.

More advanced usage of **Dimsym** such as performing manual manipulations of the determining equations if needed will require greater familiarity with **REDUCE**, while the user who wants to push the program to its limits may need to be familiar with the algorithms involved.

To learn the basics, chapter 2 of the **Dimsym** manual features an example to work through. Chapter 17 consists of more examples files to guide you. Read chapter 3 on terminology and chapter 4 on how to best enter your equations and you will be ready to use the program seriously.

## Acknowledgement

The development of this program was funded in part by the Australian Research Committee, the Department of Mathematics at La Trobe University, the Department of Mathematics at the University of Wollongong, and the CSIRO's Division of Material Science. The author also acknowledges an Australian Postgraduate Research Award.

The authors particularly wish to thank the following beta-testers of the program: Phil Broadbridge, Ted Fackerell, Greg Reid, and Willy Sarlet. Special thanks go to Alan Head, the author of [LIE](http://archives.math.utk.edu/software/msdos/adv.diff.equations/lie/.html), who has been a continuing source of ideas and assistance. We also thank Charles Wright for his assistance in integrating his **REDUCE** _ODESOLVE_ package into **Dimsym**.

## Authors

- James Sherring
- Geoff Prince
- Michael Jerie

## Homepage

- <https://www.latrobe.edu.au/mathematics-and-statistics/research/geometric-and-algebraic-techniques-for-differential-equations/dimsym>
