# SGII
In this work, we used the Systemic Gene Importance Index (SGII) to systematically identify essential lncRNAs in mouse and human genome based on the lncRNA-protein heterogeneous interaction network. The ```'data'``` folder stores the data needed in the experiment, including mouse and human data. The ```'src'``` folder contains the code used in the method SGII. The ```'result'``` folder contains the ```'human'``` folder, ```'mouse'``` folder and ```'performance'``` folder. The ```'human'``` folder and  ```'mouse'``` folder include the score results of lncRNAs using BC, CC, DC, EC and GIC methods, and the ```'performance'``` folder contains the performance results obtained by combining different methods.

# Run
To run our supplied program, you need to configure the python 3 environment.

Run our code with the following process:
1. Calculate the BC, CC, DC, and EC scores of lncRNAs:
```
python src/cal_centrality_score.py {organism} 
```
2. Calculate the GIC scores of lncRNAs:
```
python src/cal_GIC_score.py {organism}
```
3. Obtain various centrality scores of lncRNAs satisfying sequence length requirements:
```
get_centrial_score_valid.py {organism}
```
4. Generate essential lncRNAs result set:
```
python src/get_essentialLncRNAs.py {organism} {K} {T} {z}
```
5. Calculate the performance of different combinations of methods:
+ only GIC
```
python cal_performance_GIC.py {human/mouse/mouseHomologousOfHuman}
```
+ one centrality combined with GIC
```
python cal_performance_oneCentrality+GIC.py {human/mouse/mouseHomologousOfHuman} {centrality}
```
+ two centralities combined with GIC
```
python cal_performance_twoCentrality+GIC.py {human/mouse/mouseHomologousOfHuman} {centrality_1} {centrality_2}
```
+ four centralities combined with GIC
```
python cal_performance_fourCentrality+GIC.py {human/mouse/mouseHomologousOfHuman}
```

Parameters involved in commands:
+ organism: mouse or human
+ K: the threshold of the centrality method
+ T: GIC method threshold
+ z: Degree threshold
+ centrality: centrality method, include BC, CC, DC, EC
