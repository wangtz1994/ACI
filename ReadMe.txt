This is the R code for simulations of ACI method described in our ICML'20 paper "Cost-effectively identifying causal effects when only response variable is observable". It is developed by Tian-Zuo Wang (URL: www.lamda.nju.edu.cn/wangtz, Email: wangtz@lamda.nju.edu.cn). If you have some questions, please feel free to contact me.

The code contain four files.

distribution_query learning: the main program of the simulations for 100 rounds.

distribution_max_main: the main program of the ACI methods

distribution_function: some relevant functions

generate_graph: generate random causal graph and set the response variable

Notice: In the step to estimate the causal effects, our code just uses a linear model to estimate E[Y|X,Z] in our code. If you want to use this code in the non-linear data, please try to use a "stronger" regression model. The linear model is not suitable for the non-linear case.