# Welcome to Reachability-based Trajectory Safeguard
## What does this do?
This method uses the dynamics model and Reachability computation to ensures the safety of a decision-making agent(human or RL agent). We take advantage of parameterized trajectories and "adjust" the parameter selected by the decision-making agent to a guaranteed safe parameter close by.

### Paper:
Please cite our [paper](https://arxiv.org/abs/2011.08421) as 

Shao, Y. S., Chen, C., Kousik, S., & Vasudevan, R. (2020). Reachability-based Trajectory Safeguard (RTS): A Safe and Fast Reinforcement Learning Safety Layer for Continuous Control. arXiv preprint arXiv:2011.08421.

#### Abstract:
Reinforcement Learning (RL) algorithms have achieved remarkable performance in decision making and control tasks due to their ability to reason about long-term, cumulative reward using trial and error. However, during RL training, applying this trial-and-error approach to real-world robots operating in safety critical environment may lead to collisions. 
To address this challenge, this paper proposes a Reachability-based Trajectory Safeguard (RTS), which leverages trajectory parameterization and reachability analysis to ensure safety during training and testing.
This method ensures an agent with continuous action space can be trained from scratch safely in real-time.
By ensuring safety with RTS, this paper demonstrates that the proposed algorithm is not only safe, but can achieve a higher reward in a considerably shorter training time when compared to RTD, RTS with a discrete action space, and a baseline RL algorithm.
## Questions and Bugs
Please contact Yifei Shao(syifei) for questions regarding Car or Drone example, and Chao Chen(joecc) for questions regarding the cartpole example. All emails end with @umich.edu

## Dependencies
Step1: Install MATLAB 2020a. Since its RL toolbox is a bit inflexible and so modify MATLABIntallPath/toolbox/rl/rl/+rl/+env/MATLABEnvironment.m to have a the IsDone flag do a little more than what it does now: Change Line 243 from 'if isdone' to 'if abs(isdone - 1) < 0.1 || abs(isdone - 3) < 0.1 || abs(isdone - 4) < 0.1 || abs(isdone - 5) < 0.1'. Then restart MATLAB.

Step2: Clone all repositories and checking out to the correct branch

Step3: add all to MATLAB path. You should be good to go!

Sanity Check: run run_highway_testing and use the arrow keys on the figure to drive the car around, it should edit your inputs so that it never crashes.

### All required:
[RTD](https://github.com/ramvasudevan/RTD) 

[RTD_tutorial](https://github.com/skousik/RTD_tutorial) 

[simulator](https://github.com/skousik/simulator)

[CORA](https://tumcps.github.io/CORA/) checkout to commit 484c54e0d7990312741fddde5a9c9309d3e8808c

[zono_RTD_turtlebot_example](https://github.com/pdholmes/zono_RTD_turtlebot_example)

MATLAB_2020a + Just install all Toolboxes

### Drone:
[quadrotor_RTD](https://github.com/skousik/quadrotor_RTD)

edits on the repo: Change bounds 

## How to use this repo?
### Level 1 Observe Result: 
Use common_evaluation.m to see the training plots of the three examples for different methods, also use that file to tally up experiment random simulation result. To visualize how each agent performs, use run_xxx_eval.m and load different agents in agent&exp to look at the behavior of different agents.

### Level 2 Do evaluation :
Run run_xxx_eval.m till completion and save the experience to observe how good it is

### Level 3 Do training:
Run run_xxx_training with plot_sim_flag turned off, so it automatically uses parallel pool. WIth 16 parpool workers, Car training takes about 10 hours, Drone 2 hours, and Cartpole in no time.

### Level 4 Do offline FRS computation:
Car: Run gen_frs_idea5.m to get the FRS file. You may wish to clean it up using clean_up_FRS.m

Drone: The FRS was computed in the depended quadrotor_RTD repository

Cartpole: run gen_cartpole_frs.m, documentation under construction.

### Different Modes: 
Change S.safety_layer = 'Z' for proposed method, 'Z' with S.discrete_flag = 1 for discrete version of proposed method, 'N' for No safety, 'R' with HLP = [] for reward optimizing RTD, 'R' with HLP defined for RTD.

## Common Bugs
### Most episodes are very short
Make sure you have modified the rl toolbox isdone flag, and the start location is not already in collision



