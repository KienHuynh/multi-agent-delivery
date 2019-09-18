# Package delivery with a group of agents

This is the repository for the following problem:

Given n drones of different speeds, 1 package and 1 destination, find a way for the drones to deliver the package to the destination as fast as possible.

Two drones can exchange a package. Drones can have different speed or different range (future implementation).

Dependencies:

* FLTK 1.3.5
* cxxopts

To run the program in the command line:
```
ndrones.exe -g [0|1] -i [input_path] -o [output_path]

# -g: turning on GUI
# input_path: full or relative path to input file
# output_path: full or relative path to output file you want to write
```
