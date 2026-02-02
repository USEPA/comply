# Fortran version of COMPLY-r code

***

### Compiling code on Linux system

Instructions for compiling "Level 4" on a linux system:

```
make -f makefile4r
```



***

### Code application on Linux system 



Instructions for running "Level 4" on a linux system:

First, edit the input file "input.level4r.txt" to reflect the stack characteristics and emissions amounts. 

Note, the wind rose and downwind receptor distance information is combined for each stack on the input file!

In the example below, NNE is the sector wind is blowing FROM, the Freq is the wind frequency from that sector, the speed (m/s) is the wind speed from that sector, the "15000" value represents a downwind receptor distance of 15,000 m and that receptor is located in the SSW sector (wind blowing TO this sector).

Dir,Freq,Speed(m/s)|NNE,0.081,4.62,15000.,SSW

Enter the direction at the command line below to run the model. 
```
./level4r
```

Output will automatically be generated into a file called output.level4r.txt.

