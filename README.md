<a name="top"></a>

# WenOOF [![GitHub tag](https://img.shields.io/github/tag/Fotran-FOSS-Programmers/WenOOF.svg)]() [![Join the chat at https://gitter.im/Fotran-FOSS-Programmers/WenOOF](https://badges.gitter.im/Join%20Chat.svg)](https://gitter.im/Fotran-FOSS-Programmers/WenOOF?utm_source=badge&utm_medium=badge&utm_campaign=pr-badge&utm_content=badge)

[![License](https://img.shields.io/badge/license-GNU%20GeneraL%20Public%20License%20v3,%20GPLv3-blue.svg)]()
[![License](https://img.shields.io/badge/license-BSD2-red.svg)]()
[![License](https://img.shields.io/badge/license-BSD3-red.svg)]()
[![License](https://img.shields.io/badge/license-MIT-red.svg)]()

[![Status](https://img.shields.io/badge/status-alpha-orange.svg)]()

### WenOOF, WENO interpolation Object Oriented Fortran library

- WenOOF is a pure Fortran (KISS) library for computing WENO interpolations;
- WenOOF is Fortran 2003+ standard compliant;
- WenOOF is OOP designed;
- WenOOF is a Free, Open Source Project.

#### Table of Contents

+ [What is WenOOF?](#what-is-wenoof?)
	+ [What is WENO?](#what-is-weno?)
+ [Main features](#main-features)
+ [Status](#status)
+ [Copyrights](#copyrights)
+ [Documentation](#documentation)

#### Issues

[![GitHub issues](https://img.shields.io/github/issues/Fotran-FOSS-Programmers/WenOOF.svg)]()
[![Ready in backlog](https://badge.waffle.io/Fotran-FOSS-Programmers/WenOOF.png?label=ready&title=Ready)](https://waffle.io/Fotran-FOSS-Programmers/WenOOF)
[![In Progress](https://badge.waffle.io/Fotran-FOSS-Programmers/WenOOF.png?label=in%20progress&title=In%20Progress)](https://waffle.io/Fotran-FOSS-Programmers/WenOOF)
[![Open bugs](https://badge.waffle.io/Fotran-FOSS-Programmers/WenOOF.png?label=bug&title=Open%20Bugs)](https://waffle.io/Fotran-FOSS-Programmers/WenOOF)

## What is WenOOF?

Modern Fortran standards (2003+) have introduced support for Object Oriented Programming. Exploiting new features like Abstract Data Type (ADT) is now possible to develop a KISS library for computing Weighted Essentially Non-Oscillatory (WENO) interpolation on ADT making the development of new numerical schemes faster, easier and clearer.

### What is WENO?

Go to [Top](#top)

## Main features

WenOOF is aimed to be a KISS-pure-Fortran library for computing WENO interpolation, it being:

+ [ ] Pure Fortran implementation;
+ [ ] KISS and user-friendly:
  + [ ] simple API;
  + [ ] easy building and porting on heterogeneous architectures;
+ [ ] comprehensive:
  + [ ] central schemes;
  + [ ] upwind biased schemes;
  + [ ] hybrid schemes;
+ [ ] efficient:
  + [ ] high scalability on parallel architectures:
    + [ ] support for shared memory multi/many cores architecture;
    + [ ] support for distributed memory cluster;
    + [ ] support for GPGPU/accelerators device;
+ [ ] well documented:
  + [ ] clear documentation of schemes implementations;
  + [ ] complete API reference;
  + [ ] comprehensive wiki:
+ [ ] collaborative developed;
+ [ ] FOSS licensed;

Any feature request is welcome.

Go to [Top](#top)

## Status

WenOOF project is just started. Nothing has been done. We are searching for Fortraners enthusiast joining our team.

## Copyrights

WenOOF is an open source project, it is distributed under a multi-licensing system:

+ for FOSS projects:
  - [GPL v3](http://www.gnu.org/licenses/gpl-3.0.html);
+ for closed source/commercial projects:
  - [BSD 2-Clause](http://opensource.org/licenses/BSD-2-Clause);
  - [BSD 3-Clause](http://opensource.org/licenses/BSD-3-Clause);
  - [MIT](http://opensource.org/licenses/MIT).

Anyone is interest to use, to develop or to contribute to WenOOF is welcome, feel free to select the license that best matches your soul!

More details can be found on [wiki](https://github.com/Fotran-FOSS-Programmers/WenOOF/wiki/Copyrights).

Go to [Top](#top)

## Documentation

Besides this README file the WenOOF documentation is contained into its own [wiki](https://github.com/Fotran-FOSS-Programmers/WenOOF/wiki). Detailed documentation of the API is contained into the [GitHub Pages](http://Fotran-FOSS-Programmers.github.io/WenOOF/index.html) that can also be created locally by means of [ford tool](https://github.com/cmacmackin/ford).

Go to [Top](#top)
