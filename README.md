# LSD2: LEAST-SQUARES METHODS TO ESTIMATE RATES AND DATES FROM SERIAL PHYLOGENIES

## Compile LSD2:
	Type 'make' from the folder 'src', you will have the executable file 'lsd2' in the 'src' folder.
	Note that C++ compiler and library support for the ISO C++ 2011 is required to compile the program from the code source. 
     
## Run LSD:

If you want to use the interface, type ./lsd2 in the terminal without parameters.
Otherwise, type './lsd2 <parameters>'  where the list of parameters usage can be obtained by './lsd2 -h'
The input tree file is required and should be specified by option -i. 
	
If the input date is provided, the program estimates will the absolute dates. 
	
	
The input date file is necessary to estimate absolute dates and can be specified by option -d. 
The input date file should contain the date of all tips and possiblly some internal nodes if known. 
If some tip dates are missing, the program just uses the subtree containing all defined date tips & nodes for the estimation. 
The missing tip dates would be inferred at the end using the estimated rate & dates.
In order to have unique solution, at least two different precise values of dates should be 
given. A tree with all tips having the same date and no further date 
information on internal nodes will not be able to infer absolute dates. 

If all tips have the same date, you can estimate relative dates using options -a and -z to 
specify the root date and tip date. 
	
Option -c is recommended to take into account the temporal constraints (date of a node >= date of its ancestors). 
It should be noticed that LSD always assumes an increasing-time order from root to tips, i.e the date of a node is smaller than that of its children. If your data has the reverse order, the simplest way is to take the negation of the
input date, and take the negation again of the output date to obtain your expected results.

Further options can be specified, see './lsd2 -h' for more details.
    
## Input files:


### Input_tree_file

Input tree(s) in newick format are compulsory. A tree can be either binary or polytomy. The input
file must contain one tree per line:

    ((A:0.12,D:0.12):0.3,(B:0.3,C:0.5):0.4);

    ((A:0.12,B:0.3):0.7,(C:0.5,D:0.8):0.1);

### Input_date_file

An input date file is optional. If it's not provided then the program estimates the relative dates by
assuming all tips have the same date (1 by default), and the root has date 0 by default.

Suppose that we have an input ((A:0.12,D:0.12):0.3,(B:0.3,C:0.5):0.4); then an example of
input date file can be as follows:

    5			# number of temporal constraints
    A 1999		# the date of A is 1999
    B 2000		# the date of B is 2000
    C l(1990)		# the date of C is at least 1990
    D b(1998,2000)	# the date of D is between 1998 and 2000
    mrca(A,B,C) u(2000)	# the date of the most recent ancestor of A,B, and C is at most 2000
    
You can also define the labels for internal nodes and use them to define their dates. 
For example you have an input tree: ((A:0.12,D:0.12)n1:0.3,(B:0.3,C:0.5)n2:0.4)root; 
then an input date file can be as follows:

    5
    A 2000
    n1 l(2001)
    C b(2001,2004)
    n2 u(2003)
    root b(1998,1999)

### given_rate_file

If the rates are known and you want to use it to infer the dates, then you can 
give them in a file. The file should have each rate per line which corresponds 
to each tree in the Input_tree_file, for example:

	0.0068	
	0.0052


### Example of Outgroup_file format:

	2
	outgroup1
	outgroup2

If there are more than 1 outgroups, than they must be monophyletic in the input trees.

### Example of Partition_file: 

You can partition the branch trees into several subsets that you know each subset
has a different rate. 
Suppose that we have a tree `((A:0.12,D:0.12)n1:0.3,((B:0.3,C:0.5)n2:0.4,(E:0.5,(F:0.2,G:0.3)n3:0.33)n4:0.22)n5:0.2)root;` then an example for Partition_file can be as follows:

    group1 {n1} {n5 n4}
    group2 {n3}

Each line defines a list of subtrees whose branches are supposed to have the same substitution rate. Each subtree is defined between {}: the first node is the root of the subtree and the following nodes (if there any) define its tips. 
If there's not any tip defined, then the subtree is extended down to the tips of the full tree. Hence, {n1} defines the subtree rooted at the node n1; and {n5 n4} defines the subtree rooted at n5 that has one tip as n4 and other tips as the ones of the full trees (here are B,C). 
As a consequence, in this example, the branches will be partitioned into 3 groups such that each group has a different rate: 

    (1) (n1,A), (n1,D), (n5,n4), (n5,n2), (n2,B), (n2,C); 
    (2) (n3,F), (n3,G); 
    (3) the remaining branches of the tree. 
    
Note that if the internal nodes don't have labels, then they can be defined by mrca of at least two tips, for example n1 is mrca(A,D)

## Some examples of command lines:

* for rooted tree, constrained mode, and using variances

    `./lsd2 -i rootedtree_file -d date_file -c -v 1`
    
* for rooted tree, constrained mode, using variances, remove outlier tips with k=3 in Tukey's fences

    `./lsd2 -i rootedtree_file -d date_file -c -v 1 -e 3`

* for rooted tree, constrained mode, using variances, using partition file 

    `./lsd2 -i rootedtree_file -d date_file -c -v 1 -p parition_file`

* for rooted tree, constrained mode, re-estimate the root position around the given root

    `./lsd2 -i rootedtree_file -d date_file -c -r l`

* similar to the previous example, but calculate confidence intervals from 100 simulated trees 
(sequence length is required via option -s to calculate confidence intervals. The program will use the min of sequence length and 1000 to generate branch lengths of simulated trees.)

    `./lsd2 -i rootedtree_file -d date_file -c -r l -f 100 -s 1700`

* for unrooted tree without outgroups, without constraints, estimate the root position

    `./lsd2 -i unrootedtree_file -d date_file -c -r a`

* for unrooted tree with outgroups, constrained mode, using variances from the estimated branch lengths (run LSD twice), 
remove outgroups to obtain the root

    `./lsd2 -i unrootedtree_file -d date_file -g outgroup_file -c -v 2 -s 1000`

* similar to the previous example, but keep outgroups in the tree, just estimate the root position defined by the outgroups

    `./lsd2 -i unrootedtree_file -d date_file -g outgroup_file -k -c -v 2 -s 1000`

* for rooted tree, constrained mode, and using given rates to estimate dates

    `./lsd2 -i rootedtree_file -d date_file -w given_rate_file -c` 

* for rooted tree, estimating relative dates with date root=0 and date of all leaves=1, under constrained mode

    `./lsd2 -i tree_file -c -a 0 -z 1`


## Output files: 

    .result : contain the estimated rates, root date, possibly confidente intervals, outlier tips and the value of the objective function.

    .nexus : trees in nexus format which contain information about the dates of internal nodes (named date), branch lengths, and the confidence intervals (named CI) if option -f was used.
    
    .date.nexus : trees in nexus format where branch lengths are rescaled to time unit by multiplying with the estimated rate. 

## Citation
If you use this software, please cite: “ Fast dating using least-squares criteria and algorithms”, T-H. To, M. Jung, S. Lycett, O. Gascuel, Syst Biol. 2016 Jan;65(1):82-97.
