# Data preparation for PopAmaranth
This is a guide for data prepration for integrating into JBrowse:
The code for generating the calculations and transforming the files ready for tracks is available [here](https://github.com/cropevolution/amaranthGDB/tree/master/code). 


## 1. Calculate population genetics 

In this step, [ANGSD](https://www.popgen.dk/angsd/index.php/Thetas,Tajima,Neutrality_tests) is used for calculating summary statistics. The input files are a `text file with the location of the bam files` and the `reference genome .fasta file`. 
Note that while in the _site allele frequency calculation_ step that the option _anc_ was used, even though an ancestral genome was used. This is solved in the next subtep (_site allele frequency calculation_)
 
 
Use the script `code/Neutrality_tests.sh`

``` bash
#!/bin/bash -l
#SBATCH --cpus-per-task=4
#SBATCH --mem=42gb
#SBATCH --time=16:00:00
#SBATCH --account=UniKoeln
#SBATCH --mail-user=jgoncal1@uni-koeln.de
#SBATCH --error /scratch/jgoncal1/logs/errors/%j
#SBATCH -D /scratch/jgoncal1/outputs/
#SBATCH -o /scratch/jgoncal1/logs/%j

source /home/jgoncal1/.bashrc
module load gnu/4.8.2

BAM_FILE="example_list.txt" # file contianig list with BAM files
POPULATION="test_popultatoion" # Population being analysed
NIND=3 # minimum number of individuals with data to be included in the analysis (recomended at least 1/3 of the total number of samples)

## Substep 1. site allele frequency calculation

``` bash
~/angsd/angsd -out $POPULATION \
-bam /projects/ag-stetter/jdias/data/acc_list/$BAM_FILE \
-doSaf 1 \
-GL 2 \
-P 24 \
-anc /projects/ag-stetter/jdias/data/genomes/Ahypochondriacus_459_v2.0.fa \
-remove_bads 1 \
-minMapQ 30 \
-minQ 20 \
-minInd $NIND


## Substep 2. calculte SFS, use fold in case of not knwoing the ancestral state

 ~/angsd/misc/realSFS  $POPULATION.saf.idx \
 -P 4 \
 -fold 1 > test_script$POPULATION.sfs

 ## calcualtion of theta per site

 ~/angsd/misc/realSFS \
 saf2theta $POPULATION.saf.idx \
 -sfs $POPULATION.sfs \
 -outname $POPULATION  \
 -fold 1 \
 -P 4


 ## Substep 3. calcualtion of statistics from theta estimations. (If folded site frequency spectrum was used, only tD, tP and Tajima are correctly estimated)

 ~/angsd/misc/thetaStat \
 do_stat $POPULATION.thetas.idx \
 -win 5000 \
 -step 5000 \
 -outnames $POPULATION.thetasWindow.gz
```

**FST calculation**

```bash
#first calculate per pop saf for each populatoin
~/angsd/angsd -b list1  -anc hg19ancNoChr.fa -out pop1 -dosaf 1 -gl 1
~/angsd/angsd -b list2  -anc hg19ancNoChr.fa -out pop2 -dosaf 1 -gl 1
#calculate the 2dsfs prior
 ~/angsd/misc/realSFS pop1.saf.idx pop2.saf.idx >pop1.pop2.ml
#prepare the fst  
 ~/angsd/misc/realSFSfst index pop1.saf.idx pop2.saf.idx -sfs pop1.pop2.ml -fstout here
#get the global estimate
 ~/angsd/misc/realSFSS fst stats here.fst.idx 
-> FST.Unweight:0.069395 Fst.Weight:0.042349
# get the window calucaltion
 ~/angsd/misc/realSFS fst stats2 here.fst.idx -win 50000 -step 10000 >slidingwindow

```


## 2. Tracks preparation
### 2.1 Convert .thetasWindow.gz to .bedgraph
From the output obtained from step 1, the files need to be converted into a _.bedgraph_ format. This is done on the first substep (by population):


```{bash}
#!/bin/bash

STATS_FILE=confirmation_quitensis.thetasWindow.gz.pestPG
POPULATION=quitensis

awk '{print $1, $2, $4, $5, $9, $14}' $STATS_FILE |  sed '1d;$d' | sed 's/(/\t/g' | sed 's/)/\t/g' | sed 's/,/\t/g' | awk '{print $7, $5, $6, $8, $9, $10, $11}' | sed 's/ /\t/g' | sort -k1,1 -k2,2n > $POPULATION_tests.bedgraph
```

### 2.2 load neutrality tests bedgrraph file for per site calculation of Watterson and nucleotide diversity
`code/Stats_track_preparation.Rmd`

Now, the outputs from the previous analysis can be restructured. 
Use `R`:

``` R
##### packages and necessary libraries
library(data.table)
library(RcppCNPy)
library("factoextra")
library('ggfortify')
library('cowplot')
library('dplyr')

nice_layout<-   theme_cowplot()+
    panel_border(color = "grey85", size = 1, linetype = 1,
  remove = FALSE, "black") 



  popualation_tests<-fread('/Users/josedias/mount/scratch/jgoncal1/outputs/caudatus_tests.bedgraph') %>%
          dplyr::rename(CHR=V1, StartPos=V2, EndPos=V3, tW=V4, tP=V5, Tajima=V6, nSites=V7) %>%
  mutate(Watterson=tW/nSites, pi=tP/nSites, species = "caudatus")
  
popualation_tests %>%
  select(CHR, StartPos, EndPos, pi) %>%
  fwrite("population_pi.bedgraph")

popualation_tests %>%
  select(CHR, StartPos, EndPos, tW) %>%
  fwrite("population_tW.bedgraph")


popualation_tests %>%
  select(CHR, StartPos, EndPos, Tajima) %>%
  fwrite("population_Tajima.bedgraph")

```
### 2.3  Convert begraph R output to BigWig for adding quantitative trakcs on the Genome Browser

``` bash
fa2size reference.fasta reference.sizes -detailed
bedGraphToBigWig in.bedGraph chrom.sizes out.bw
```

**Bonus**: Plot summary staitstics:  plot Watterson, pi and Tajima D

``` R
all_species.filtered <- full_join(caudatus_tests, cruentus_tests) %>%
  full_join(quitensis_tests) %>%
  full_join(hybridus_tests) %>%
  full_join(hypochondriacus_tests) %>%
  full_join(quitensis_tests) %>%
    filter(nSites>250)


all_species.filtered %>%
    group_by(species) %>%
summarise_if(is.numeric, mean)
    

  plot_Tajima_all_species_filtered<-            ggplot(all_species.filtered, aes(x=Tajima)) +
geom_density(aes(color = species, fill=species), alpha=0.4) +
nice_layout +

           labs(title = "Tajima")
      
      
             plot_pi_all_species_filtered<- ggplot(all_species.filtered, aes(x=pi)) +
geom_density(aes(color = species, fill=species), alpha=0.4) +
nice_layout +
        labs(title = "pi")

                   
             plot_Watterson_all_species_filtered<-     ggplot(all_species.filtered, aes(x=Watterson)) +
geom_density(aes(color = species, fill=species), alpha=0.4) +
nice_layout      +     
           labs(title = "Watterson")
           
         
           
title.filtered <- ggdraw() + 
draw_label(
"All species nsites >250",
fontface = 'bold',
x = 0,
hjust = 0
) 

plot_grid(plot_Tajima_all_species_filtered, plot_pi_all_species_filtered, plot_Watterson_all_species_filtered, labels="AUTO", title.filtered)
```

## 3. Add data to the Browser
Lastly, the files need to be added to location where the browser is hosted.

inside the `popamabrowser/` folder, there are the configuration files `jbrowse.conf`  `jbrowse_conf.json`. This is the **Dangerous Zone** as these files control the setup of `jbrowse`, including plugins.  In principle, they do not need to be altered. They also have the indication of where the data for the tracks is located. In the case of popamaranth this are:
- `dataRoot = data`
- `include  = {dataRoot}/trackList.json`

This means that data files for the tracks should be added to the `data` folder. 
It also means that for every file added, the information about this file **needs** to be added to the `/trackList.json` file.

Below there is an example for one data point.
The `bicolor_pivot` was mannualy added for each track to force the global mean yellow line. This allows the global mean being present, but the tracks will autoscale to the current view.
`urlTemplate` is where the file name should be added. In this example, `hybridus_pi.bw`. Each data point needs to between {} and separated by comma.

``` json
      {
         "autoscale" : "local",
         "bicolor_pivot" : "0.019854152683326895",
         "key" : "nucleotide diversity (π) A. hybridus",
         "label" : "nucleotide diversity (π) A. hybridus",
         "metadata" : {
            "Diversity" : "π",
            "Name" : "nucleotide diversity (π) A. hybridus",
            "Species" : "A. hybridus",
            "category" : "Diversity",
            "Track description" : "Nucleotide diversity (Nei, M. and W.-H. Li,  1979). Average nucleotide diversity (π) in 5kb windows. The yellow linerepresent the genome-wide meanπ. Values above the mean are represented in blue, values below the genome-wide mean are red. Shading in light and darker grey represents 1 standard deviationand 2 standard deviations from the mean, respectively. Maximum and minimum scale is adjusted for the local region in display."
         },
         "noFill" : true,
         "storeClass" : "JBrowse/Store/SeqFeature/BigWig",
         "style" : {
            "height" : "200",
            "neg_color" : "red",
            "pos_color" : "navy",
"origin_color": "yellow",
            "variance_band_color" : "rgba(0,0,0,0.04)"
         },
         "type" : "JBrowse/View/Track/Wiggle/XYPlot",
         "urlTemplate" : "hybridus_pi.bw",
         "variance_band" : true
      },
      
```

