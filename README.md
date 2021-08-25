New Version！

[![INFORMS Journal on Computing Logo](https://INFORMSJoC.github.io/logos/INFORMS_Journal_on_Computing_Header.jpg)](https://pubsonline.informs.org/journal/ijoc)

This archive is distributed in association with the [INFORMS Journal on
Computing](https://pubsonline.informs.org/journal/ijoc) under the [MIT License](LICENSE).

The data in this repository is an archive of the data that were used in the research reported on in the paper [An Interior-Point Differentiable Path-Following Method to Compute Stationary Equilibria in Stochastic Games](link) by Chuangyin Dang, P.Jean-Jacques Herings and Peixuan Li.

To cite this data, please cite the [research article](link) and the data itself.

Below is the BibTex for citing this version of the data.

```
@article{Data.IJOC.link,
  author =        {Chuangyin Dang, P.Jean-Jacques Herings and Peixuan Li},
  publisher =     {INFORMS Journal on Computing},
  title =         {Data for An Interior-Point Differentiable Path-Following
Method to Compute Stationary Equilibria in Stochastic Games},
  year =          {2021},
  doi =           {doi},
  url =           {https://github.com/INFORMSJoC/link},
}  
```

## Content
The IPM works in MatLab software. Almost all numerical experiments are randomly generated, and the experimental results can be found in the numerical section of the paper. This repository gives the code files so that one can reproduce these results by running the code.

This repository includes:

1. The code in the folder **CoASLTP** is for comparing the proposed interior-point difierentiable path-following method (**IPM**) and the **ASLTP**, where the file [ycsgse.m](CoASLTP/ycsgse.m) is the main program of the **IPM** and [dltpsgse.m](CoASLTP/dltpsgse.m) is the main program of the **ASLTP**. The code in this folder has been used in Section 4.1.
2. The code in the folder **CoPathsolver** is for comparing the proposed **IPM** and the **path solver**, where the file [trysg.m](CoPathsolver/trysg.m) is the main program. The code in this folder has been used in Sections 4.2 and 4.3.
3. The folder **Bargaining** includes the code for computing a solution to the bargaining model, which has been presented in Section 4.4.

We next illustrate how to associate the code files with the numerical results(e.g., tables and figures) presented in the paper.
1. By running the files [exm1.m](CoASLTP/exm1.m), [exm2.m](CoASLTP/exm2.m), [exm3.m](CoASLTP/exm3.m),  [exm4.m](CoASLTP/exm4.m), [exm5.m](CoASLTP/exm5.m) in the folder **CoASLTP**, one can obtain the results for the five fundamental examples in Section 4.1.
2. By running the file [inputs225.m](CoASLTP/inputs225.m) in the folder **CoASLTP**, one can get a stationary equilibrium for a stochastic game with two players, two states and five actions for each player in each state. The computational costs of the **IPM** and **ASLTP** for solving this instance are obtained as well. By repeatedly running [inputs225.m](CoASLTP/inputs225.m) for ten times, one can obtain the average computational costs for both methods, which are shown in the first row of Table 1 and Table 2. Through changing the parameters *n*; *d*; *m*; *pd0* in [inputs225.m](CoASLTP/inputs225.m), we can attain various instances. The average computational costs for solving these stochastic games are shown in Table 1, Table 2, and Table 4.
3. By implementing the file  [se225.m](CoPathsolver/se225.m) in the folder **CoPathsolver**, one can get the comparison results between the proposed **IPM** and the **path solver** for computing a stationary equilibrium in a randomly generated stochastic game with two players, two states and five actions. Similarly, by changing the parameters *n*; *d*;*m*, we attain various stochastic games with difierent scales. The comparison results are included in Table 3. By running the file [r1.m](CoPathsolver/r1.m), one may obtain the success rates of the two methods for 100 randomly generated stochastic games, which are recorded in Figure 5.
4. By implementing the file [bargaining.m](Bargaining/bargaining.m) in the folder **Bargaining**, one can get Figure 6, which shows a solution to the presented bargaining model.
