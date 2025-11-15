# Getting Started

## Installation

### Quick Installation Instructions

First, clone the Rhea repository onto your machine, then change to the new
directory `rhea` with the cloned code:

```sh
git clone https://bitbucket.org/ccgo/rhea.git rhea
cd rhea
```

Set up Rhea's submodules:

```sh
git submodule init
git submodule update
```

Prepare Rhea for compilation with the `bootstrap` script (this needs to be done
only once after a fresh `git clone`):

```sh
./bootstrap
```

Now run the usual configure and make cycle:

```sh
./configure
make
```

The configure script will try to find necessary tools in your path.  When you
run configure you may optionally use the `--program-prefix=PREFIX` argument to
change the default installation directory.
