# SimuCell
Brain tumor growth Monte-Carlo simulator

*NEW: We now provide here, together with the source code, a flowchart of the program, an executable file (SimuCellWinx64.exe) compiled for Windows x64 and a batch file (RunExample_Command_lines.bat) containing two command lines running the simulation with two sets of parameters. One set implements clonal-dynamic based only on inherited variation of cell-cycle length, the other set adds to the former simulation a clonal-based competition (with an alienWeight=3 and lineageWeight=0.1).


A series of questions were raised by the reviewers following the first submission of the manuscript reporting the results of SimuCell simulations. We think that the questions and the relative answers can be relevant and helpful for all potential readers, therefore we decided to include here the point by point list of answers 

Q-how large is the simulated region?

A- The space in which the simulation takes place was designed to be virtually unconstrained, since it is a cube with a 6 cm-long side. However, given that all the simulation parameters were tuned to reproduce real-world, simulated gliomas never exceeded 7.3 mm of diameter, in agreement with experimental gliomas. In particular, without implementing competition, simulated gliomas have a diameter of 6.9 ± 0.5 mm. Implementing competition, simulated gliomas have a diameter of 5.2 ± 0.06 mm.

Q-are cells modelled as point particles in a continuous space or as e.g. spheres? Or was a discrete space used?

A- Cells were modelled as point particles in a continuous space. Maximum local density is constrained by the link between cell density and apoptosis rate. 

Q-do cells interact with each other, at what distance is this is possible and what are the consequences?

A- Cells interact each other by reciprocally modifying their fitness (apoptosis probability). In the simulation, we explored a range of maximum interaction distances from 50 μm to 300 μm. 

Q-is it possible for cells to undergo multiple events during the considered time units of 1 hour?

A- At each time unit, cells can perform one or more of the following actions: migrate, undergo mitosis, quiescence, or apoptosis. The migration is always performed but the displacement might occasionally be zero. No action can be repeated in the same time unit. 
 
Q-how is migration modelled?

A- Each hour X,Y and Z coordinates of the cells are updated by adding a normally distributed random number to each of them. The mean value of the random added number was implemented as a simulation parameter and tuned to match experimental in-vivo speed (Supplementary figure 5c). 

Q-how is competition implemented?

A- As we mentioned before, cells interact each other by reciprocally modifying their fitness (cell death probability). Higher is the number of neighbouring cells, higher is the probability of cell death. When clonal-based competition is implemented, the neighbours of a given cell C are accounted in two different class: cells belonging to the same clone of C and cells belonging to different clones. The cell death probability is then calculated giving different weight to the two numbers: and neighbours belonging to the same clone as C have a lighter effect than neighbours belonging to the same clone. The actual weights were implemented as simulation parameters and their ratio have been systematically explored in a range from 1 (no clonal-based competition) to 30. 

Q- a list of model parameters should be provided and their default values and ranges studied.

A- The default values of the parameters are listed in the first lines of the code. The list of parameters for the two most relevant simulations are collected in the folder "Simulation_parameter_list". The studied ranges for the parameters are summarized in the file " Explored_ranges-Parameters.txtt".
