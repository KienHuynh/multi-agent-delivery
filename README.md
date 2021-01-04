# Package delivery with a group of agents

This is the repository for the following problem:

Given n agents of different speeds, a package(s) from a source(s) and a destination(s), find a way for the agents to deliver the package to the destination as fast as possible.

Two agents can exchange a package. Agents can have different speed.

### ===========HOW TO COMPILE AND RUN===========

Dependencies:

* FLTK 1.3.5
  * All include files are in ndrones/includes/
  * All necessary lib files for linking are in ndrones/lib/
  * Dependencies needed to compile:
    * fltkd.lib
    * fltkformsd.lib
    * fltkgld.lib
    * fltkimagesd.lib
    * fltkjpegd.lib
    * fltkpngd.lib
    * fltkzlibd.lib
* cxxopts
  * All include files are in ndrones/includes
  * No linking needed

To run the program in the command line:
```
ndrones.exe -g [0|1] -i [input_path] -o [output_path]

# -g: turning on GUI (1), or off (0)
# input_path: full or relative path to input file
# output_path: full or relative path to output file you want to write
```
To test if the program is running correctly, run it with the ndrones/ndrones/testcase/input0.txt

Input file: the format of an input file is described in ndrones/ndrones/inputDescription.txt

Output file: same as input file.

### ===========PROBLEM DETAILS===========
There are two main problems:

1 There are multiple agents, but there is only 1 destination to consider. These problems are depicted as 
```
TWODIM EUCLID DISCRETE SINGLE_ID
```
on the top of each input file.

  1.1 It is possible to have more than one sources here, each serves the same purpose, i.e. the agents can pick up the package from any of those sources.

2 There are more than one destination to consider. This will then be on top of the input file:
```
TWODIM EUCLID DISCRETE
```

Sources and targets can be points:
```
SINGLE_POINT
```
or polygons (regional):
```
POLY
```

Sampling method:

The algorithm work on a discrete set of points. All of the handoffs occur on these points. The recommended sampling method is
```
APGRID
x1 y1 sx
x2 y2 sy
```
Which create a grid from (x1, y1) to (x2, y2) in which any two adjacent points are spaced by sx/sy along the x/y direction.

### ===========VISUALIZATION (GUI MODE)===========
![](viz_aid/figs/gui/animation1.gif?raw=true)

The main GUI:
![](viz_aid/figs/gui/1.jpg?raw=true  =100x100)
