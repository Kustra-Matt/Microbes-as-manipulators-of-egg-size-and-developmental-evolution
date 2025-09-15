Description of code and data files affiliated with Kustra and Carrier 2025 "*Microbes as manipulators of egg size and developmental evolution*." https://journals.asm.org/doi/10.1128/mbio.03655-24

For questions, contact: Matthew Kustra: [matthewckustra@gmail.com](mailto:matthewckustra@gmail.com) 

To properly run the **analysis** code, first unzip both **"Code.zip"** and **"Data.zip."** Then, make sure the Code and Data files are in the directories described below. This is most easily accomplished by moving all code files in the ***Analysis*** directory of the unzipped **"Code.zip"** to the directory above the unzipped **"Data.zip."**

Below is a description of all the code and data files broken down by the directories in which code and data files should be organized in.

**Running Sims** and **Analysis** files uploaded to **Zenodo**.

***Running_Sims***: Folder that contains all the files needed to run simulations. Only occurs in **"Code.zip"** and does not require any data.

* **"QGPop_Simulation.jl"**: Julia script that was used for our quantitative genetic population simulations. Code is heavily commented. Session information for packages and versions are given in **"QGPop_Simulation_Session_Info.txt"** within the ***SessionInfo*** folder.
* **"Grid_Simulation.R"**: R script that was used to run the invasion grid model. Code is heavily commented. Session information for packages and version are given in **"Grid_Simulation_session_info.txt"** within the ***SessionInfo*** folder.

***Analysis***: Folder that contains all R code to generate figures/summarize data.

* **"Graph_Grid.R"**: R script to generate figures and analyses related to invasion grid analysis. Specifically makes Figs. 1B,C; Fig S1B; Fig. S3, S4. Uses files in the ***Data*** folder. Session information for packages and R version are given in **"Graph_Grid_session_info.txt"** within the ***SessionInfo*** folder.
* **"FigureS2.R"**: R script to generate figures and analyses related to figures S2 of the invasion grid analysis. Specifically makes Fig S2. Uses files in the ***Data*** folder. Session information for packages and R version are given in **"Graph_Grid_session_info.txt"** within the ***SessionInfo*** folder.
* **"Graph_QGPop.R**": R script to generate figures from quantitative genetic population model. Specifically makes: Fig.2,3, and Fig S6,7,8,9,10,11. Uses files in the ***Data*** folder. Session information for packages and R version are given in **"QGPop_graphing_session_info.txt"** within the ***SessionInfo*** folder.
* **"Paramater_graphing.R**": R script to generate graphs from some of the equations. Specifically Figures S1A and S5.Session information for packages and R version are given in **"Parm_graphing_session_info.txt"** within the ***SessionInfo*** folder.

***Data***: Folder with all the data generated from the model that were used in analysis/graphing.

* "**fem_grid.csv**": Data file of simulation output from feminization grid based simulations. This is used by "**Graph_Grid.R**". Here are the descriptions of the column variables:
  * ***DIRF***: Direction of egg size evolution (larger,stable,smaller, either, NA) from the egg size given in the "Egg_Size" column. "larger" means larger egg sizes should evolve. "stable" means the current egg size (i.e., resident egg size) has the highest fitness. "smaller" means smaller egg sizes should evolve. "Either" means both larger and smaller egg sizes have higher fitness than the current egg size. NA means there was an error in the calculation
  * ***RESWF***: Fitness of resident/current/focal egg size (the egg size given in the Egg_size column). Fitness is measured as the number  of offspring produced (zygotes per μL)
  * ***LWF***: Fitness of mutant with an egg size that is 10 uM smaller than the resident egg size (i.e. egg size in "Egg_Size" column. Fitness is measured as the number  of offspring produced (zygotes per μL)
  * ***HWF***: Fitness of mutant with an egg size that is 10 uM larger than the resident egg size (i.e. egg size in "Egg_Size" column. Fitness is measured as the number  of offspring produced (zygotes per μL)
  * ***Egg_Size***: Egg cross sectional area (mm2)
  * ***Egg_Size2***: Egg size diameter (um)
  * ***Larval_m***: Larval mortality
  * ***Egg_number*** Density of eggs produced. (eggs per uL)
  * ***Feminization*** Feminization rate. The proportion of individuals that become female
  * ***Time*** Number of generations to reach equilibrium
  * ***Density*** Stable density (individuals/m^2)
  * ***Sex_ratio*** stable sex ratio
  * ***Error*** Says whether there was an error
  * ***nif0*** Density of females m^2
  * ***nim0*** Density of males m^2
  * ***cost*** Enhanced growth rate/cost
  * ***MK*** Male killing Rate. The proportion of males killed
  * ***b*** *B* parameter. Specifically, B is the parameter determining the importance of egg size on survival, where larger values of B equate to egg size being more important to offspring survival 
* "**MK_grid.csv**": Data file of simulation output from Male killing grid based simulations. This is used by "**Graph_Grid.R**". Here are the descriptions of the column variables:
  * ***DIRF***: Direction of egg size evolution (larger,stable,smaller, either, NA) from the egg size given in the "Egg_Size" column. "larger" means larger egg sizes should evolve. "stable" means the current egg size (i.e., resident egg size) has the highest fitness. "smaller" means smaller egg sizes should evolve. "Either" means both larger and smaller egg sizes have higher fitness than the current egg size. NA means there was an error in the calculation
  * ***RESWF***: Fitness of resident/current/focal egg size (the egg size given in the Egg_size column). Fitness is measured as the number  of offspring produced (zygotes per μL)
  * ***LWF***: Fitness of mutant with an egg size that is 10 uM smaller than the resident egg size (i.e. egg size in "Egg_Size" column. Fitness is measured as the number  of offspring produced (zygotes per μL)
  * ***HWF***: Fitness of mutant with an egg size that is 10 uM larger than the resident egg size (i.e. egg size in "Egg_Size" column. Fitness is measured as the number  of offspring produced (zygotes per μL)
  * ***Egg_Size***: Egg cross sectional area (mm2)
  * ***Egg_Size2***: Egg size diameter (um)
  * ***Larval_m***: Larval mortality
  * ***Egg_number*** Density of eggs produced
  * ***Feminization*** Feminization rate. The proportion of individuals that become female
  * ***Time*** Number of generations to reach equilibrium
  * ***Density*** Stable density (individuals/m^2)
  * ***Sex_ratio*** stable sex ratio
  * ***Error*** Says whether there was an error
  * ***nif0*** Density of females m^2
  * ***nim0*** Density of males m^2
  * ***cost*** Enhanced growth rate/cost. Proportion of increase/decrease in fecundity
  * ***MK*** Male killing Rate. The proportion of males killed
  * ***b*** *B* parameter. Specifically,B is the parameter determining the importance of egg size on survival, where larger values of B equate to egg size being more important to offspring survival 
* "**FemQGPop.csv**": Data file of simulation output from quantitative genetic model. This is used by "**Graph_QGPop.R**". The column names and descriptions are below:
  * ***V1*** column rows (please ignore).
  * ***Generation*** Generation of the simulation
  * ***EggSize*** Egg diameter uM
  * ***Density*** Population density (individuals/m^2)
  * ***SexRatio*** Proportion of population male
  * ***MK*** Feminization Rate. The proportion of individuals that become female
  * ***GR*** Enhanced growth rate. Proportion of increase in fecundity
  * ***B***: *B* parameter. Specifically, B is the parameter determining the importance of egg size on survival, where larger values of B equate to egg size being more important to offspring survival 
  * ***V_b***: V parameter. The V parameter alters the proportional influence that microbial abundance—a proxy for the functional capacity of that microbial population—has on the offspring survival-egg size relationship, with smaller V values representing a greater influence on that relationship
* "**Fig2QG.csv**": Data file of simulation output from quantitative genetic model used for figure 2. This is used by "**Graph_QGPop.R**". The column names and descriptions are below:
  * ***Generation*** Generation of the simulation
  * ***EggSize*** Egg diameter uM
  * ***Density*** Population density (individuals/m^2)
  * ***SexRatio*** Proportion of population male
  * ***MK*** Feminization Rate (The proportion of individuals that become female), or male killing rate (proportion of males that are killed). Depends on Cat2 column
  * ***GR*** Enhanced growth rate. Proportion of increase in fecundity
  * ***B***: *B* parameter. Specifically, B is the parameter determining the importance of egg size on survival, where larger values of B equate to egg size being more important to offspring survival 
  * ***V_b***: V parameter. The V parameter alters the proportional influence that microbial abundance—a proxy for the functional capacity of that microbial population—has on the offspring survival-egg size relationship, with smaller V values representing a greater influence on that relationship.
  * ***Cat***: Full description of reproductive manipulation and enhanced growth
  * ***Cat2***: What reproductive manipulation
  * ***GRCat***: Whether enhanced growth is present
* "**FigS6lossQG.csv**": Data file of simulation output from quantitative genetic model used for figure S5. This is used by "**Graph_QGPop.R**". The column names and descriptions are below:
  * ***Generation*** Generation of the simulation
  * ***EggSize*** Egg diameter uM
  * ***Density*** Population density (individuals/m^2)
  * ***SexRatio*** Proportion of population male
  * ***MK*** Feminization Rate (The proportion of individuals that become female), or male killing rate (proportion of males that are killed). Depends on Cat2 column 
  * ***GR*** Enhanced growth rate. Proportion of increase in fecundity
  * ***B***: *B* parameter. Specifically, B is the parameter determining the importance of egg size on survival, where larger values of B equate to egg size being more important to offspring survival 
  * ***V_b***: V parameter. The V parameter alters the proportional influence that microbial abundance—a proxy for the functional capacity of that microbial population—has on the offspring survival-egg size relationship, with smaller V values representing a greater influence on that relationship.
  * ***Cat***: Full description of reproductive manipulation and enhanced growth
* "**FigS7GQ.csv**": Data file of simulation output from quantitative genetic model used for figure S6. This is used by "**Graph_QGPop.R**". The column names and descriptions are below:
  * ***Generation*** Generation of the simulation
  * ***EggSize*** Egg diameter uM
  * ***Density*** Population density (individuals/m^2)
  * ***SexRatio*** Proportion of population male
  * ***MK*** Feminization Rate (The proportion of individuals that become female), or male killing rate (proportion of males that are killed). Depends on Cat2 column
  * ***GR*** Enhanced growth rate. Proportion of increase in fecundity
  * ***B***: *B* parameter. Specifically, B is the parameter determining the importance of egg size on survival, where larger values of B equate to egg size being more important to offspring survival 
  * ***V_b***: V parameter. The V parameter alters the proportional influence that microbial abundance—a proxy for the functional capacity of that microbial population—has on the offspring survival-egg size relationship, with smaller V values representing a greater influence on that relationship
* "**MK_QGPop.csv**": Data file of simulation output from male killling quantitative genetic model. This is used by "**Graph_QGPop.R**". The column names and descriptions are below:
  * ***V1*** Row number (ignore)
  * ***Generation*** Generation of the simulation
  * ***EggSize*** Egg diameter uM
  * ***Density*** Population density (individuals/m^2)
  * ***SexRatio*** Proportion of population male
  * ***MK*** male killing rate (proportion of males that are killed)
  * ***GR*** Enhanced growth rate. Proportion of increase in fecundity
  * ***B***: *B* parameter. Specifically, B is the parameter determining the importance of egg size on survival, where larger values of B equate to egg size being more important to offspring survival 
  * ***V_b***: V parameter. The V parameter alters the proportional influence that microbial abundance—a proxy for the functional capacity of that microbial population—has on the offspring survival-egg size relationship, with smaller V values representing a greater influence on that relationship
