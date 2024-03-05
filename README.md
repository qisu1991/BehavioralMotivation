# BehavioralMotivation

DESCRIPTION
-----------
These are Julia programs (version 1.5.3) and Matlab programs (version 2023b) designed for the production of a paper entitled "Evolutionary dynamics of behavioral motivations."

FILE
-----
<pre>
Figure3.jl      obtain the frequency of cooperation for a range of need threshold \alpha and four types of games and six types of networks (see Figure 3 in the main text)
                provide detailed instructions for each function
Figure4.mlx     obtain the frequency of cooperation for a range of need threshold \alpha, number of levels L, and benefit factor \omega (see Figure 4 in the main text)
Figure5.jl      obtain the fixation probability of motivation A for a range of benefits b, on random regular networks and scale-free networks (see Figure 5 in the main text)
Figure6.jl      obtain the evolutionary trajectories of behavioral motivations and cooperation frequencies (see Figure 6 in the main text)
</pre>
                        
INSTALLATION
------------
<pre>
*) Download the open-source software Julia 1.5.3 or a newer version from https://julialang.org/downloads/.
*) Install Julia on your local PC. The expected installation time is approximately five minutes.
*) Download the open-source software VS Code 1.77.3 or a newer version from https://code.visualstudio.com/download.
*) Install VS Code on your local PC. The expected installation time is approximately ten minutes.
*) Install the Julia VS Code extension from https://code.visualstudio.com/docs/languages/julia. The expected installation time is approximately five minutes.
*) Run VS Code, add Figure3.jl, and install a solver package: 
   1) Type the following command in TEMINAL, and then press ENTER:
      julia
   2) Type the following command, and then press ENTER:  
      using Pkg; Pkg.add("IterativeSolvers");
*) Click "Run and Debug" to execute the program.
*) The expected outputs are
   48 data files, corresponding to the cooperation frequency for 48 interaction scenarios, i.e. two behavioral motivations, four types of games, and six types of networks (see the dots in Figure 3 in the main text). The expected run time is approximately forty minutes;
*) Add Figure5.jl, and click "Run and Debug" to execute the program. The expected outputs are
   4 data files, corresponding to the fixation probability of motivation A for 4 interaction scenarios, namely two behavioral motivations and two types of networks (see the dots in Figure 5 in the main text). The total run time is approximately 120 hours. The time can be reduced by distributing the computations across multiple computers.
*) Add Figure6.jl, and click "Run and Debug" to execute the program. The expected outputs are
   6 data files, corresponding to the evolutionary trajectories of motivation intensity \lambda, need threshold \alpha, and the cooperation frequency xC, for interaction scenario b/c<(b/c)* and for b/c>(b/c)* (see the dots in Figure 6 in the main text). The total run time is approximately 120 hours.
*) Run Matlab, add Figure4.mlx and click "Run" to execute the program. The expected outputs are
   a figure showing the frequency of cooperation for a range of need threshold \alpha, number of levels L, and benefit factor \omega (see Figure 4 in the main text). The expected run time is approximately five minutes.
</pre>

QUESTIONS
---------
For any questions about this program, please contact
Dr. Qi Su, Email: qisu@sjtu.edu.cn.
