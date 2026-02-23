# Meeting 2/13/26
## 1. Plans:
### 1.1 Diversity of Genes expression

A 0-parameters fit of the observed diversity as a function of the number of cells using the empirical expression probability of individual genes confirm that there are no correlation between the expression of different genes.

### 1.2 Fitting of dose response curves

It is hard to really infer 2 parameters / curves, as a result we try to fit $EC_{50}$ the half max activation concentration. It is not working perfectly, but I don't think it just don't work.

-> What about the amplitude

### 1.3 Mutual Information at fixed concentration

This is different from previous approaches, where the agent tries to measure the both the ligand identity AND concentration.

-> Computing the empirical mutual information shows that on the pool of ligand proposed, the array of the octopus is overall pretty bad....
-> The question of the heteromerization, in this case, is really not relevant. Comparing the homo-pentamers arrays, the difference is already massive.
-> We see that we need more than $10^4$ ligands with 5 homo-pentamers to feel the difference.
-> Adding the concentration dimension, would justify more heteromerization.

## 2. What's next ?

-> It seems that some of these components are either not part of the octopus environment or he just doesn't care.
    -> Need to model the environment
    -> Classify the chemical in familly ?
-> Looking for optimal array to specific optimization strategy might reveal why so many homomers do not respond.