# BamCleave
Working to extend Nigel's software to work on groups of cells more than single cells.


Extended existing software to split .bam files (gene mutation information files) by groups of cells, based on the unique identifier of a cell. 
Afterwards, using existing software such as IGv, gene mutations can be observed as whether they only occured in a type of cells. It is very useful when separating the gene mutations by immune cells and cancerous cells. 

The input consists of one .bam file and a mapping file (cellID to groupID). The output consists of a .bam file corresponding to each group of cells. 

The work here is based on previously developed software by dr. Nigel Dyer, University of Warick and it was extended during my undergraduate research project 2018-2019.
