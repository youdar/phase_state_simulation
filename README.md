# mix_scale_simulation

The code in this repository demonstrates the use of a hybrid simulation method, stochastic (Gillespie’s method) & deterministic (base of the equation describing the process), for biologically complex phenomena, where different processes happening on a very different time scales. Typical analysis of such complex systems or processes are broken down to region of behavior, simplifying the mathematical equations describing the process by using the different scales involved. 
For example, in Alzheimer’s disease the disease develops over many years while the aggregation of Aβ peptide, that might be an important player in the disease, might be affected by processes in the cell, on the scale of fraction of a second.

This approach promotes the idea of keeping the process intact, exploring it via clever simulation strategies.

This code include tools to create a “phase-state diagrams” allowing to understand how parameters that are known at one state might look like at a different state. For instance in Alzheimer’s disease related experiment in mice, both the concentration and the rate of production of Aβ are very different than those at humans. The phase-state diagrams allows to explore specific parameters and a particular state, and see what are, for example, equivalent concentration and rate of production of Aβ in humans. 

