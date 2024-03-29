# Meta-analysis-simulation-supplementary-code

NOTES ON R CODE

I would like to apologise in advance to anyone reading 
this code - I am not an experienced programmer and wrote it 
all over an extended period, with lots of alterations to fix bugs, 
improve efficiency, and sometimes to alter functionality. As such
there are likely sections which are redundant, inefficient, and 
with confusing variable names or comments.

There are multiple sets of R files. The first set (NS or LS as first characters) were submitted to a 
supercomputing cluster to generate the data. The second set of files (with Analysis in title)
generated variables for analysis from these data and collected them into a 
single document for all MD results and one for all LOR results. Copies of these documents
are also included as csv files for those who don't wish to rerun the simulation but want to use the results.
The data included in the csvs was examined in multiple ways, but the figures in the article (along with several others)
can be generated using the csvs with the files which have SimpleDiagrams in the title.
Finally the looped diagrams from the supplementary material can be generated using the files which have 
LoopDiagram in the title and the csvs.

Small alterations to the code have been made to aid comprehension, but the results
should be exactly the same as we found.

There are several sections of code that I have copied from various sources, including 
stackexchange and at least one function from the metafor package, which
was cited in the main paper.

Specific issues of note before running these scripts are:

1. There are passing references in the code to another type of publication bias
and a more agressive form of outcome reporting bias. These were experimental 
and not eventually analysed.

2. The analysis code is designed for a Windows computer, and will
search for file directories and create files. It will need adapting before
it works correctly in full.

3. The first set of scripts create very large files. This was chosen so that 
forms of analysis which require all the simulated values could be performed,
but if you are only interested in summary results for each condition, then 
by adding some of the code from the second set of files you could reduce the
exported file sizes dramatically. To reduce file size and increasing loading speed
RDS file types were used, which should be noted before use.

4. The first set of scripts are written to be performed in a parallel manner. 
It may well work in a situation where it cannot use parallel resources 
but it will be substantially less efficient.

5. Code was initially run using R 3.0.2, then checked to ensure function in R 4.0.3
