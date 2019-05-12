# lsd-0.3beta
LSD: LEAST-SQUARES METHODS TO ESTIMATE RATES AND DATES FROM SERIAL PHYLOGENIES - v0.3beta by Thu-Hien TO

If you use this software, please cite: “ Fast dating using least-squares criteria and algorithms”, T-H. To, M. Jung, S. Lycett, O. Gascuel, Syst Biol. 2016 Jan;65(1):82-97.


How To Compile LSD:

     Use C++ compiler and library support for the ISO C++ 2011 to compile the program from the code source. From the folder src, type 'make':
     
How to Run LSD:

	After compiling the program, you have the lsd executable file in the src folder, then run it from a terminal. You can also use the binary files in the bin folder but sometimes it might not work on some computers that have different configuration as the machine compiled them. From the directory that contains the executable file:
	
		if you want to use the interface, type ./lsd
		
		if you want to use the command line, type ./lsd -i <"input_tree_file"> -d <"input_date_file"> (to estimate absolute dates)
		
		                                       or ./lsd -i <"input_tree_file"> -a root_date -z leaves_date (to estimate relative dates)
		                                       
			further options can be specified (-f -o -c -v -n -r -g -p ...). Type "./lsd -h" for help. Option -c is recommended to take into account the temporal constraints (date of a node >= date of its ancestors).

## Some examples of input files:


### Example of Input_tree_file format 

Newick format, can be either binary or polytomy, each tree per line:

    ((A:0.12,D:0.12):0.3,(B:0.3,C:0.5):0.4);

    ((A:0.12,B:0.3):0.7,(C:0.5,D:0.8):0.1);

### Example of Input_date_file format 

It's not necessary to give the temporal constraints for all tips. Suppose that we have an input ((A:0.12,D:0.12):0.3,(B:0.3,C:0.5):0.4); then we can have an input date for example as follows:

    5			# number of temporal constraints
    A 1999		# the date of A is 1999
    B 2000		# the date of B is 2000
    C l(1990)		# the date of C is at least 1990
    D b(1998,2000)	# the date of D is between 1998 and 2000
    mrca(A,B,C) u(2000)	# the date of the most recent ancestor of A,B, and C is at most 2000
    
You can also define the labels for internal nodes and use them to define their temporal constraints. For example you have an input tree: ((A:0.12,D:0.12)n1:0.3,(B:0.3,C:0.5)n2:0.4)root; then you can have an input date file as follows:

    5
    A 2000
    n1 l(2001)
    C b(2001,2004)
    n2 u(2003)
    root b(1998,1999)

### Example of given_rate_file format 

Each rate per line which corresponds to each tree in the Input_tree_file:

	0.0068	
	0.0052


### Example of Outgroup_file format:

	2
	outgroup1
	outgroup2

### Example of Partition_file: 

Suppose that we have a tree ((A:0.12,D:0.12)n1:0.3,((B:0.3,C:0.5)n2:0.4,(E:0.5,(F:0.2,G:0.3)n3:0.33)n4:0.22)n5:0.2)root; then an example for Partition_file can be as follows:

    group1 {n1} {n5 n4}
    group2 {n3}

Each line defines a list of subtrees that are supposed to have the same substitution rate. Each subtree is defined between {}: the first term is the root of the subtree and the following terms (if there any) define its tips. 
If there's not any tip defined, then the subtree is extended down to the tips of the full tree. Hence, {n1} defines the subtree rooted at the node n1; and {n5 n4} defines the subtree rooted at n5 that has one tip as n4 and other tips as the ones of the full trees (here are B,C). 
As a consequence, in this example, the branches will be partitioned into 3 groups such that each group has a different rate: (1) (n1,A), (n1,D), (n5,n4), (n5,n2), (n2,B), (n2,C); (2) (n3,F), (n3,G); (3) the remaining branches of the tree. Note that if the internal nodes don't have labels, then they can be defined by mrca of at least two tips, for example n1 is mrca(A,D)

## Some examples of command lines:

### for rooted tree, constrained mode, and using variances

    ./lsd -i rootedtree_file -d date_file -c -v 1

### for rooted tree, constrained mode, using variances, using partition file 

(note that sequence length is required via option -s to calculate variances)

    ./lsd -i rootedtree_file -d date_file -c -v 1 -s 500 -p parition_file

### for rooted tree, constrained mode, re-estimate the root position around the given root

    ./lsd -i rootedtree_file -d date_file -c -r l

### similar to the previous example, but calculate confidence intervals from 100 simulated trees 

(note that sequence length must be specified by option -s for calculating confidence intervals)

    ./lsd -i rootedtree_file -d date_file -c -r l -f 100 -s 1700

### for unrooted tree without outgroups, without constraints, estimate the root position

    ./lsd -i unrootedtree_file -d date_file -c -r a

### for unrooted tree with outgroups, constrained mode, using variances from the estimated branch lengths (run LSD twice), remove outgroups to obtain the root

    ./lsd -i unrootedtree_file -d date_file -g outgroup_file -c -v 2 -s 1000

### similar to the previous example, but keep outgroups in the tree, just estimate the root position defined by the outgroups

    ./lsd -i unrootedtree_file -d date_file -g outgroup_file -k -c -v 2 -s 1000

### for rooted tree, constrained mode, and using given rates to estimate dates

    ./lsd -i rootedtree_file -d date_file -w given_rate_file -c 

### for rooted tree, estimating relative dates with date root=0 and date of all leaves=1, with constraint

    ./lsd -i tree_file -c -a 0 -z 1


## Output files: 

    .result : contain the estimated rates, root date and the value of the objective function.

    .newick : trees in newick format with the new branch length (re-estimated by the program).

    .date.newick : trees in newick format where branch lengths are measured rescaled to time unit by multiplying with the estimated rate. 

    .nexus : trees in nexus format which contain information about the dates of internal nodes (named date), branch lengths, and the confidence intervals (named CI) if option -f was used.
