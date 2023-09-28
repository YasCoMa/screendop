# ScreenDOP - Screening of strategies for disease outcome prediction

Pipeline that performs a screening of disease outcome prediction evaluations following two strategies, one of them combines protein interactions data with normalized table of genes by samples of a gene expression experiment for some genetic disease in order to predict the outcome of this disease for the patients (each sample). The other strategy uses gene set enrichment with kegg pathways to build the features that will enter the classifier. In both strategies, it is possible to set up configurations to test the better parameter combination. The pipeline receives a json configuration file allowing the test of many experiments by execution. Here, we present examples for two scenarios: Leukaemia and Ovarian cancer datasets.

## Summary

The data preparation pipeline contains tasks for two distinct scenarios: [leukaemia](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE425) that contains microarray data for 119 patients and [ovarian](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE140082) cancer that contains next generation sequencing data for 380 patients.

The disease outcome prediction pipeline offers two strategies for this task:

**Graph kernel method**: It starts generating personalized networks for each patient using the interactome file provided and generate the patient network checking if each PPI of the interactome has both proteins up regulated or down regulated according to the gene expression table provided. The first step generate a set of graphs for the patients that are evaluated with 4 distinct kernels for graph classification, which are: Linear kernel between edge histograms, Linear kernel between vertex histograms and the Weisfeiler lehman. These kernels functions calculate a similarity matrix for the graphs and then this matrix is used by the support vector machine classifier. Then the predictions are delivered to the last task that exports a report with the accuracy reached by each kernel. It allows some customizations about the network parameters to be used, such as the DEG cutoff to determine up and down regulated based on the log2 fold change, which will determine the topology and the labels distribution in the specific sample graphs. It is also possible customize the type of node/edge attributes passed to the kernel function, which may be only label, only weight or both.

**GSEA based pathway scores method**: This method is faster and do not rely on tensor inputs such as the previous method. It uses geneset enrichment analysis on the pathways from KEGG 2021 of Human, and uses the scores of the pathways found enriched for the samples to build the numerical features matrix, that is then delivered to the AdaBoost classifier. The user may choose balance the dataset using oversampling strategy provided by SMOTE.

## Usage Instructions
### Preparation:
1. ````git clone https://github.com/YasCoMa/screendop.git````
2. ````cd screendop````
3. Decompress screening_ovarian/raw_expression_table.tsv.tar.xz
4. Create conda environment to handle dependencies: ````conda env create -f drugresponse_env.yml````
5. ````conda activate drugresponse_env````
6. Setup an environment variable named "path_workflow_screendop" with the full path to this workflow folder

### Data preparation - File ````data_preparation_for_pipeline.py```` :

#### Files decompression

- Decompress data_preparation/lekaemia.tar.xz
- Decompress data_preparation/ovarian/GSE140082_data.tar.xz
    - Put the decompressed file GSE140082_series_matrix.txt in data_preparation/ovarian/
    
#### Pipeline parameters

- __-rt__ or __--running_type__ <br>
	Use to prepare data for the desired scenario: <br>
	1 - Run with Leukaemia data <br>
	2 - Run with Ovarian cancer data

#### Running modes examples

1. Run for Leukaemia data: <br>
````python3 data_preparation_for_pipeline.py -rt 1 ```` 

In this case, you must have [R](https://www.r-project.org/) installed and also the library [limma](https://bioconductor.org/packages/release/bioc/html/limma.html), it is used to determine DEGs from microarray data. For this dataset, the files are already prepared in the folder.

2. Run for Ovarian cancer data: <br>
````python3 data_preparation_for_pipeline.py -rt 2 ```` 

In this case, you must have [R](https://www.r-project.org/) installed and also the library [DESeq](https://bioconductor.org/packages/release/bioc/html/DESeq.html), because this scenario treats next generation sequencing data

### Disease outcome prediction execution - File ````main.py````:

#### Pipeline parameters

- __-rt__ or __--running_step__ <br>
	Use to prepare data for the desired scenario: <br>
	1 - Run graph kernel method <br>
	2 - Run gsea based pathway scores method

- __-cf__ or __--configuration_file__ <br>
	File with the expression values for the genes by sample/patient in tsv format<br>
	
	Example of this file: config.json
		
#### Input configuration file

- Configuration file keys (see also the example in config.json):
    - **folder** (mandatory for both methods): working directory
    - **identifier**: project identifier to be used in the result files
    - **mask_expression_table** (mandatory for both methods): Gene expression values file with the result of the fold change normalized value of a certain gene for each sample, already pruned by the significance (p-value). 
    - **raw_expression_table** (mandatory for both methods): Raw gene expression values already normalized following the method pf preference of the user.
    - **labels_file** (mandatory for both methods): File with the prognosis label for each sample
    - **deg_cutoff_up**: Cutoff value to determine up regulated gene. Default value is 1.
    - **deg_cutoff_down**: Cutoff value to determine down regulated gene. Default value is -1.
    - **nodes_enrichment**: Node attributes to be used in the screening evaluation. It may be a list combining the options "label", "weight" or "all". Examples: ["all", "weight"], ["label"], ["label", "weight"]. Default value is ["all"].
    - **edges_enrichment**: Edge attributes to be used in the screening evaluation. It may be a list combining the options "label", "weight" or "all". Examples: ["all", "weight"], ["label"], ["label", "weight"]. Default value is ["all"].
    - **flag_balance**: Flag to indicate whether the user wants to balance the samples in each outcome class, by SMOTE oversampling. Values may be false or true. Default value is false.

#### Running modes examples
1. Running disease outcome prediction by graph kernel method: <br>
	````python3 main.py -rt 1 -cf config.json````

2. Running disease outcome prediction by gsea enriched network method: <br>
	````python3 main.py -rt 2 -cf config.json````

## Reference
Martins, Y. C. (2023). Multi-task analysis of gene expression data on cancer public datasets. medRxiv, 2023-09.

## Bug Report
Please, use the [Issue](https://github.com/YasCoMa/screendop/issues) tab to report any bug.
