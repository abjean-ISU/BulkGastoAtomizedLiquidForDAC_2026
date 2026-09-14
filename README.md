This repository contains the data and code relevent to the manuscript Techno-economic and life cycle assessment of a bulk gas-to-atomized liquid approach for physisorptive direct air capture: performance characterization and research pathways submitted to Energy and Environmental Science

Structure Folders:
Droplet Tower Modeling Files: Modeling files used to predict mass transfer and liquid holdup in BGAL contactors. User-input parameters such as density, Henry's constant, and diffusivity are sorbent/gas-dependent. PD-gqmom_0_1.m was shared by Alberto Passalacqua and should be cited: Fox, R. O., Laurent, F., & Passalacqua, A. (2023). The generalized quadrature method of moments. Journal of Aerosol Science, 167, 106096. https://doi.org/10.1016/j.jaerosci.2022.106096. Droplet towers were updated from 2026 paper submitted to AIChE Journal to process a flowing gas stream with changing composition instead of queicent gas phase.

IL_Aspen: contains Aspen plus .bkp file for IL system model used for MCDA and the partially optimized model. All stream results, model results, and .txt results for each process simulation are included.

LCA: contains OpenLCA results from project simulations, Inventory assessments generated from Aspen Plus results, OpenLCA JSON-LD files for custom processes and flows, and plotting .py files. 

MCDA Files: BGAL Settings.xlsx contains the setpoints used to define the Aspen Plus adsorption and desorption excel blocks for each process simulation. MCDA_SA.xlsx contains MCDA computations and process simulation results; python file used to generate MCDA figures. 

ORCA files: contains inputs, outputs, and molecular structure files used to generate COSMOSAC property parameters for Aspen Plus. More information on ORCA can be found at: https://orcaforum.kofo.mpg.de/ and Neese, F.; "Software update: the ORCA program system, version 4.0" WIREs Comput Mol Sci 2017, e1327. doi: 10.1002/wcms.1327; also contains .py file used to process ORCA simulation outputs into Aspen Plus COSMOSAC property parameters format.

TEA Files: contains TEA computation files used to compute CAPEX and OPEX for DAC systems; also contains results files for all SAs and .py files generated for figure generation

Water_Aspen: contains Aspen plus .bkp file for water system model used for MCDA and the partially optimized model. All stream results, model results, and .txt results for each process simulation are included.

