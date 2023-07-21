Multistate Model (MSM) Transition Estimation

We used a parallel computing framework for estimating the transitions in MSM shown in Figure 1 of the main paper. Due to the computational demands, initial disease transitions are run on separate nodes, while the tasks for fitting later transitions, where appropriate, are grouped together on a single node.


R Code Description

hpc_code_01 to hpc_code_05: These codes fit the transitions from 'none' to the first diagnosis.

hpc_code_1 to hpc_code_5: These codes fit the transitions following the diagnosis of the first condition.

hpc_code_cp3: This code fits all the transitions with two pre-existing conditions.

hpc_code_cp4: This code fits all the transitions with three pre-existing conditions.

hpc_code_345: This code, fits transitions from CKD+MH+HF; it is run separately due to alterations in the code.

hpc_code_cp5: This code fits all the transitions with four pre-existing conditions.


Slurm Script Description

cprd_MM_xx slurm scripts correspond to the above R codes and are used to submit the jobs to the High-Performance Computing (HPC) environment.
