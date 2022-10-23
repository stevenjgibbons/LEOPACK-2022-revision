In this example, there are two stages:

First we have to run the itfvf program on the input file
itfvf.input to generate the inhomogeneous temperature vectors:

itfvf < itfvf.input 

generates the files 

example_aHEATFLUX.ints
example_aHEATFLUX.vecs
example_aHEATFLUX.xarr

Second, we run the time-stepping code using these vectors.

o2ibtctsc2 < example_a.input 


