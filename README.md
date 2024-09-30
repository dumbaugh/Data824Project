Shiny Project -- Data Visualization and Acquisition
Author: David Umbaugh
Date Created: 2021-07-29


.•´¯`•._.•´¯`•._.•´¯`•._.•´¯`•._.•´¯`•._.•´¯`•._APAP Single Cell Project _.•´¯`•._.•´¯`•.__.•´¯`•._.•´¯`•.__.•´¯`•._.•´¯`•.



About the data:

Acetaminophen overdose is the most common cause of liver failure in the United States, resulting in over 10,000 hospitalizations annually.
To understand how different cell types within the liver respond to acetaminophen hepatotoxicity, single cell RNA sequencing experiments have been
performed. The data in the application comes from three different single cell RNA sequencing experiments:

1. Umbaugh DS, Ramachandran A, Jaeschke H. Spatial Reconstruction of the Early Hepatic Transcriptomic Landscape after an Acetaminophen Overdose using Single-cell RNA Sequencing. Toxicol Sci. 2021 May 12:kfab052. doi: 10.1093/toxsci/kfab052.
2. Walesky, C.M., Kolb, K.E., Winston, C.L. et al. Functional compensation precedes recovery of tissue mass following acute liver injury. Nat Commun 11, 5785 (2020). https://doi.org/10.1038/s41467-020-19558-3
3. Kolodziejczyk, A.A., Federici, S., Zmora, N. et al. Acute liver failure is regulated by MYC- and microbiome-dependent programs. Nat Med 26, 1899–1911 (2020). https://doi.org/10.1038/s41591-020-1102-2

All three experiments were carried out in fasted, 8 week old, C57Bl6J mice. The Umbaugh and Walesky data treated mice with 300 mg/kg acetaminophen while Kolodziejczyk treated mice with 500 mg/kg APAP.
These three experiments are complementary as they examined different phases of the acetaminophen hepatotoxicity timecourse, an early injury phase <6h after acetaminophen exposure, injury phase,
and recovery phase (>48h). Additionally, dose dependent responses can be assessed as the Kolodziejczyk paper examined the injury phase at the higher dose (500 mg/kg). The injury peaks around 24h after APAP exposure,
while the liver is functionally recovered by 96h post APAP.

Preprocessing of the data:

Standard preprocessing steps were employed to remove low confidence cells and low quality data. Then the 3 datasets were integrated and batch corrected using Seurat v 4.02. All the datasets were 
manually annotated by cell type and underwent cluster analysis. The two primary annotation classifications are the "biological replicate" class and "celltype" class. The biological replicate represents the different
experimental groups, control, APAP2h (mice sacrificed 2h after APAP exposure), APAP6h, APAP20h_HighDose, APAP24h, APAP48h, APAP96h. The celltype class represents the parenchymal and non-parenchymal cell populations 
found in the liver (hepatocytes, immune cells, neutrophils, endothelial cells and hepatic stellate cells, the hepatocytes are divided into their lobular location, e.g. pericentral hepatocytes or periportal hepatocytes,
as they have found to be distinct with specialized functions dependent on anatomical location).   

What the application can do:

The application takes advantage of Shinys navBarPage() function and grid layout system. The primary interactions occur across three tabs, each offering different visualization options for the data.

Tab 1, "Overview: Data Visualization and Acquisition", provides a simple overview of the dataset, links to the respective papers, and an example image of what the application can produce.

Tab 2, "Visualization by UMAP and violin plot" allows the user to get a global view of the single cell gene expression data and the distribution of particular genes across cell types or biological replicates.
Essentially, the user can explore the data first by selecting how they want the groups defined, either by 'BiologicalReplicate' or 'CellType2'. After selecting the class of interest, the user can enter any gene to visualize
the expression level of that gene at the global level (UMAP representation) and at a sub-global level (Violin plot). However, the user may be interested in understanding how a gene changes in a specific cell type (such as neutrophils),
along the acetaminophen toxicity timecourse (represented by 'BiologicalReplicate'). The user can "Choose an identity class to split by" which will allow for this sort of examination. Finally, the user can choose a cell subpopulation to
examine using a checkbox input.

Tab 3, "Dotplot comparing gene expression across groups" allows the user to iput a list of genes separated by a comma, so they can see the proportion of cells in an identity class that expresses the user-inputed genes, as well as the expression
levels. Again, there is an option for the user to select a different identity class, so they can assess how the gene(s) change with respect to cell type or throughout the APAP toxicity timecourse.

Tab 4, "Feature scatter comparing two different genes" allows the user to enter two different genes, to see the relationship between the expression levels of those genes (e.g. are they directly proportional, inversely, no relationship etc.).
Again, the user has control over the group classifications, either by cell type or biological replicate. The user can also define which cell population(s) to include in this analysis via a checkbox input. For example, some genes are largely
restricted to certain cell subpopulations. For example, PP hepatocytes are high in Cyp2f2 expression, low in Cyp2e1 expression, while PC hepatocytes are high in Cyp2e1/low Cyp2f2. Other cell subpopulations in the liver express very low amounts of
either gene, therefore, its most interesting to only include PC and PP hepatocytes in the analysis. 

Representative images of the output:

NOTE: The single cell RNA seq RDS object had to be downsampled for upload to Github (1400 cells per group). The images were collected on the RDS object downsampled to 3000 cells per group. The overall structure of the data was maintained.

"image1_tab1": Introduces the data, no user inputs involved

"image2_tab2": Cyp2e1 was inputted as the gene, the identity class biological replicate was selected and split by 'celltype2'. This changes the violin plot such that the expression level for Cyp2e1 can be assessed for each biological replicate
(e.g. control, APAP2h, APAP6h etc.) split by the cell types. For example, Control is then split into the 6 defined cell types. So for example it is clear that Cyp2e1 was predominantly expressed in hepatocytes, particularly pericentral hepatocytes,
and that there is a loss of Cyp2e1 expression in PC and PP hepatocytes at 24h following APAP exposure. This is due the necrotic cell death of PC hepatocytes following APAP toxicity. Additionally, I selected pericentral and periportal hepatocytes
in the check grid, which then only graphs the PC and PP hepatocytes in the UMAP.

"image3_tab3": Here I input a list of genes that are marker genes for each of the candidate cell types. For example, Cyp2e1, Cyp1a2 and Rgn are PC hepatocyte markers, while Cxcr2 is a neutrophil marker. The dotplot highlights the average expression
and the number of cells in each group (in this case 'CellType2' e.g. the cells).

"image4_tab4": This tab allows for the comparison of two genes to understand the correlative relationship between them. For example, Cyp2e1 is a pericentral hepatocyte marker, while Cyp2f2 is a periportal hepatocyte marker, therefore these two genes
should yield a negative correlation, however if I include all the cell types in the liver I won't be able to assess the relationship of these two genes in the relevant cell populations of interest. Therefore, I first chose the 'celltype2'
identity, then chose only pericentral and periportal hepatocytes to graph. Correspondingly, I got a correlation of as -0.72 in the center title region of the graph. If I were to include all the cell types in this analysis I get a correlation of 0.23 
(not shown -- try it yourself!).
