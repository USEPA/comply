# Fortran version of COMPLY code

### Compiling code on Linux system

Instructions for compiling "Level 1" on a linux system:

```
make -f makefile1
```

Instructions for compiling "Level 2" on a linux system:

```
make -f makefile2
```

### Code application on Linux system 

Instructions for running "Level 1" on a linux system:

First, edit the input file "input.level1.txt" to reflect the radionuclides, phase state (S, L, G), and amount for comparison to the possession limit table. 


```
./level1
```

Instructions for running "Level 2" on a linux system:

First, edit the input file "input.level2.txt" to reflect the radionuclides and emissions amounts. 

```
./level2
```
