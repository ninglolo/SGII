# SGII
In this work, we used the Systemic Gene Importance Index (SGII) to systematically identify essential lncRNAs in mouse and human genome based on the lncRNA-protein heterogeneous interaction network. The ```'data'``` folder stores the data needed in the experiment, including mouse and human data. The ```'src'``` folder contains the code used in the method SGII. The ```'result'``` folder contains the ```'human'``` folder, ```'mouse'``` folder and ```'performance'``` folder. The ```'human'``` folder and  ```'mouse'``` folder include scores of lncRNAs using BC, CC, DC, EC and GIC methods, respectively, and also include essential lncRNAs predicted by SGII under the optimal parameter combination. The ```'performance'``` folder contains the performance results obtained by combining different methods.

# Run
To run our supplied program, you need to configure the python 3 environment.

Run our code with the following steps:
1. Calculate the BC, CC, DC, and EC scores of lncRNAs:
```
python src/cal_centrality_score.py {organism} 
```
In this step, we we run the code using files LPI.csv and LPPI.csv, which are stored in the ```'data'``` folder. The results of running the code are stored in the ```'result'``` folder, named BC_score_allLncRNAs.csv, CC_score_allLncRNAs.csv, DC_score_allLncRNAs.csv and EC_score_allLncRNAs.csv, respectively.

2. Calculate the GIC scores of lncRNAs:
```
python src/cal_GIC_score.py {organism}
```
In this step, we run the code using files ncName_ncID_transID.csv, eng.csv and transcripts_seq.fasta, which are stored in the ```'data'``` folder. The result from running the code are stored in the ```'result'``` folder, named GIC_score.csv.

3. Obtain various centrality scores of lncRNAs satisfying sequence length requirements:
```
get_centrial_score_valid.py {organism}
```
In this step, we aim to obtain the centrality scores of lncRNAs that meet the length requirements. Results BC_score_allLncRNAs.csv, CC_score_allLncRNAs.csv, DC_score_allLncRNAs.csv and EC_score_allLncRNAs.csv of Step 1 and GIC_score.csv of Step 2 are needed, and these results are saved in the ```'result'``` folder. Finally, We also save the results BC_score.csv, CC_score.csv, DC_score.csv and EC_score.csv of the code in this step in the ```'result'``` folder.

4. Generate essential lncRNAs result set:
```
python src/get_essentialLncRNAs.py {organism} {K} {T} {Z}
```
In this step, we predict the essential lncRNAs according to the values of parameters K, T and Z. When running the code, we need to use the file LPI.csv, results BC_score.csv, CC_score.csv, DC_score.csv and EC_score.csv obtained in step 3 and result GIC_score.csv obtained in step 2, where LPI.csv is saved in the ```'data'``` folder, and BC_score.csv, CC_score.csv, DC_score.csv, EC_score.csv and GIC_score.csv are saved in the ```'result'``` folder. Finally, we save the predicted results in the ```'result'``` folder, named essentialLncRNAs.csv.

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
In this step, we run the above codes using files BC_score.csv, CC_score.csv, DC_score.csv, EC_score.csv and GIC_score.csv, which are saved in the ```'result'``` folder. In addition, when calculating the prediction performance of human, we need the file esslnc_ homo.csv, which is saved in the ```'data'``` folder. Finally, we save the performance results in the folder ```'result/ performance'```.


Parameters involved in commands:
+ organism: mouse or human
+ K: the threshold of the centrality method
+ T: GIC method threshold
+ Z: Degree threshold
+ centrality: one centrality method of BC, CC, DC, EC
+ centrality_1/centrality_2: one centrality method of BC, CC, DC, EC, centrality_1 and centrality_2 have different values

When we calculated the performance of different combinations of methods, we calculated not only the performance of humans and mouse, but also the mouse homologous to human.

The detailed explanations of the above files are saved in the ```'data'``` and ```'result'``` folders respectively.

