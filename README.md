# KMC_GUI
Kinetic Monte Carlo simulation with Qt5 GUI

Manual for Kinetic Monte Carlo (KMC) simulator 

This program is not finished and can be full of bugs. If you find any bugs not listed here, please send e-mail to LauriSikk@gmail.com
Current buglist:

1) File menu -> exit does nothing
2) Dragging-dropping links between nodes can be complicated when multiple links already exist - move nodes to isolated space and try again (this is caused by link objects being above node objects on canvas)

Fast tutorial.


Toolbar is composed of reactions (green) and species (blue) node buttons. Left mouse click places the currently selected tool.
Left click and drag moves nodes, right click opens menu to delete/edit/select the node. Each node has 2 lines of text: 1st line
is the name of the reaction or species; second line is concentration [mol/l] for speces and rate constant for reactions. [1/sec] for monomolecular reaction, [1/M*sec] for bimolecular reaction. Right click - edit allows to change the name and rate constant/concentration values (names have to be unique!). Species and reaction nodes have plugs - gray triangles to connect them to other reactions or species.
Connections between reactions and species are made by clicking and dragging from one plug to another.

Files can be saved, loaded etc., extension is *.kmc

To edit KMC parameters open Edit menu -> KMC parameters.
timestep for data storage: concentrations and species numbers will be saved using this interval (for 0.1s, time series will be 0 s, 0.1 s, 0.2 s, ...)
max simulation time: time when the simulation will be stopped
total number of starting molecules: number of species present at t=0s. Larger number reduces the errors from random number generation but also 
increases simulation time. For simple systems less than 10000 is enough
number of repeats: this number of simulations will be done and results will be averaged. For simple systems, more than 10 is not needed. Larger number
will greatly increase overall computational time.
Simulation volume: informational number to indicate what is the current volume of simulation "box" - calculated from total number of starting molecules and concentration of all species.

To run the simulation, first save the file and then select Edit-> Run
After simulation is finished, 2 csv files will be generated: filename_concentration.csv and filename_population.csv These hold the data of finished simulation (concentrations and species count) in comma separated values format (so you can easily load it to excel or somewhere else).

To plot the results, select Edit-> Plot results. This will open new window with matplotlib graph and a list of all species. To include or exclude species to graph select or deselt them from the list and click show button. To zoom in or change plot appearance use matplotlib toolbar.

If you load a file and the population and concentration csv files are present, you should be able to plot results without running simulations again.



