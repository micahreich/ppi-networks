# Computational methods for protein function prediction in protein-protein interaction networks
Micah Reich, Sarah Fisher, Yoona Choi, Leon Xie (mreich@cmu.edu, skfisher@andrew.cmu.edu, yoonac@andrew.cmu.edu, leonx@andrew.cmu.edu)

Final project for 02-251 (Spring 2022)

## Abstract
Proteins are an integral component of biological processes, serving as the tools with which cells replicate, transcribe genomes, regulate cell function, and catalyze the reactions necessary for a cell to survive. Given that the range of possible functions a protein can serve is highly diverse, researchers seek to narrow the range of functions that an non-annotated protein may serve. The existence of experimentally-found protein-protein interaction networks opens the door to graph-based approaches for function annotation of novel proteins. In this paper, we describe three such approaches: a majority neighbor approach, a network flow approach, and a neighbor sequence alignment approach. We find that Functional Flow scores best on measures of accuracy and specificity across different radii; we also note that Sequence Alignment performs similarly at small radii, while Majority Approach performs well at smaller radii though suffers in measures of specificity as `r` increases.

## Introduction
Protein function prediction is an important problem in biology. With the advent of next-generation sequencing techniques, the amount of biological data available is at an all-time high. Moreover, computational techniques on genomic data have allowed researchers to identify new proteins. Traditional experimental techniques to research protein function are insufficient to keep up with this growth. In light of this, computational methods offer a scalable solution. 

Being able to predict protein functions allows for better understanding of key biological mechanisms; for example, understanding the functions of proteins involved in disease mechanisms can propel drug development, and identifying homologous proteins in different species can aid in evolutionary tracing. 

In the current literature, a plethora of methods for function annotation exist. These include analysis of sequence similarity, gene expression, structure, phylogeny, and more, employing graph analysis, deep learning, and traditional biological experimentation. We focus on graph theoretic methods within protein-protein interaction networks for the purpose of function annotation.

## Final Paper
To read the full paper, please visit [this link](https://github.com/micahreich/ppi-networks/blob/main/02251fnl_paper-3.pdf).
