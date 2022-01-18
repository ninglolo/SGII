# SGII
In this work, we used the Systemic Gene importance Index (SGII) to systematically identify essential lncRNAs in mouse and human genome based on the lncRNA-protein heterogeneous interaction network. The ```'data'``` folder stores the data needed in the experiment, including mouse and human data. The ```'src'``` folder contains the code used in the method SGII. The ```'result'``` folder holds the results for mouse and human, including the score results of lncRNAs using BC, CC, DC, EC and GIC methods.

# Run
To run our supplied program, you need to configure the python 3 environment.

Run our code with the following command:
+ Calculate the BC, CC, DC, and EC scores of lncRNAs:
```
python src/cal_centrality_score.py {organism} 
```
+ Calculate the GIC scores of lncRNAs:
```
python src/cal_GIC_score.py {organism}
```
+ Generate essential lncRNAs result set:
```
python src/get_essentialLncRNAs.py {organism} {K} {T} {z}
```

Parameters involved in commands:
+ organism: mouse or human
+ K: the threshold of the centrality method
+ T: GIC method threshold
+ z: Degree threshold
