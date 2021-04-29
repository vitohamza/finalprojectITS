This is electrical engineering bachelor thesis at Institut Teknologi Sepuluh Nopember. The title of the thesis is " OPTIMAL ALLOCATION AND SIZING OF ENERGY
STORAGE CONSIDERING LOAD SHEDDING "

The optimal power flow is taken using metaheuristic AI method called differential evolution. The code was written using MATLAB 2020b. Using MATPOWER 6.0 and Differential Evolution Algorithm taken from Yarpiz

Run the main program via InputbattDE. **Core files** are: 
1. ESInv (function for solving energy storage sizing)
2. DCOPF_SOLVER (function for solving the DC OPF)
3. rundcopf 
4. mpoption

**Source code flow**
rundcopf --> dcopf_solver --> opf_execute --> opf --> runopf --> rundcopf 

rundcopf --> runopf --> opf --> opf_execute --> dcopf_solver
