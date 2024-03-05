# BehavioralMotivation

DESCRIPTION
-----------
These are Julia programs (version 1.5.3) and Matlab programs (version 2023b) designed for the production of a paper entitled "Evolutionary dynamics of behavioral motivations."

FILE
-----
<pre>
Figure3.jl      obtain the frequency of cooperation for a range of need threshold \alpha for four types of games and six types of networks (see Figure 3 in the main text)
                provide detailed instructions for each function
Figure4.mlx     obtain the frequency of cooperation for a range of need threshold \alpha, number of levels L, and benefit factor w (see Figure 4 in the main text)
Figure5.jl      obtain the fixation probability of motivation A for a range of benefits b, on random regular networks and scale-free networks (see Figure 5 in the main text)
Figure6.jl      obtain the 
</pre>
                        
INSTALLATION
------------
<pre>
*) Download the open-source software Julia 1.5.3 or a newer version from https://julialang.org/downloads/.
*) Install Julia on your local PC. The expected installation time is approximately five minutes.
*) Download the open-source software VS Code 1.77.3 or a newer version from https://code.visualstudio.com/download.
*) Install VS Code on your local PC. The expected installation time is approximately ten minutes.
*) Install the Julia VS Code extension from https://code.visualstudio.com/docs/languages/julia. The expected installation time is approximately five minutes.
*) Run VS Code, add DynamicNetworks.jl, and install a solver package: 
   1) Type the following command in TEMINAL, and then press ENTER:
      julia
   2) Type the following command, and then press ENTER:  
      using Pkg; Pkg.add("IterativeSolvers");
*) Click "Run and Debug" to execute the program.
*) The expected outputs are
   1) Two critical benefit-to-cost ratios of examples from Figure 2a in the main text (by analytical computations), one with \alpha=0.5, the other with \alpha=0.8, corresponding to the two vertical lines in Figure 2b. The expected run time is approximately three minutes;
   2) A set of critical benefit-to-cost ratios presented in the inset of Figure 3a in the main text (by analytical computations). The expected run time is approximately five minutes.
</pre>

QUESTIONS
---------
For any question about this program, please contact
Dr. Qi Su, Email: qisu@sjtu.edu.cn.
