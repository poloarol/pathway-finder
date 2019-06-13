# pathway-finder


A Genomic Mining Tool which connects to the GenBank database via the NCBI Utils API. 
Using a queried gene, it builds a genomic pathway +/- bp around it and performs a blast to
find identical genes with the provided similarity ratio (Uses the Edit/Levenshtein ratio).
Finally it outputs the list of GenBank files (abc.gb) containing the various pathways.

# Installation

To install this program, PIP can be used as follows:

pip install pathway-finder


# Accessing the data

For now all the program only supports request over the E-UTILS web interface of blast. The 
programs bottleneck is the online blast. With an API key, 10 request can be done before pausing
the script for 3 seconds, but when lacking it only 3 request can be done.
