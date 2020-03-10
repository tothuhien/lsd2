# LSD2: LEAST-SQUARES METHODS TO ESTIMATE RATES AND DATES FROM PHYLOGENIES

## Compile/install LSD2:

### Compile from source:

Type *make* from the folder *src*, you will have the executable file *lsd2* in the same place.
Note that C++ compiler and library support for the ISO C++ 2011 is required to compile the program from sources. 
     
### Install via Homebrew:

Mac/Linux users can install lsd2 via Homebrew as follows:

`brew install brewsci/bio/lsd2`
     
## Run LSD2:

If you want to use the interface, type *./lsd2* without parameters in the terminal from the folder containing the executable file.
Otherwise, type *./lsd2 parameters*  where the list of parameters can be obtained by *./lsd2 -h*.

The input tree file is required and should be specified by option -i. 
	
The input date file is necessary to estimate absolute dates and can be specified by option -d. 
The input date file should contain the date of all tips and possiblly some internal nodes if known. 
If some tip dates are missing, the program just uses the subtree containing all defined date tips & nodes for estimating the rate. 
The missing tip dates would be inferred at the end using the estimated rate & dates.
In order to have unique solution, at least two nodes should have different given dates. 
A tree where all tips having the same date and no further date information on internal nodes will not be able to infer absolute dates. 
In this case, you can estimate relative dates using options -a and -z to specify the root date and tip date. 
	
Option -c is recommended to take into account the temporal constraints (date of a node >= date of its ancestors). 
It should be noticed that LSD2 always assumes an increasing-time order from root to tips, i.e the date of a node is smaller than that of its children. If your data has the reverse order, the simplest way is to take the negation of the
input date, and take the negation again of the output date to obtain your expected results.

Further options can be specified, see *./lsd2 -h* for more details.
    
## Input files format


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

### Given rate file

If the rates are known and you want to use it to infer the dates, then you can 
give them in a file. The file should have each rate per line which corresponds 
to each tree in the Input_tree_file, for example:

	0.0068	
	0.0052


### Outgroup file

	2
	outgroup1
	outgroup2

If there are more than 1 outgroups, than they must be monophyletic in the input trees.

### Partition file

You can partition the branch trees into several subsets that you know each subset
has a different rate. 

Suppose that we have a tree `((A:0.12,D:0.12)n1:0.3,((B:0.3,C:0.5)n2:0.4,(E:0.5,(F:0.2,G:0.3)n3:0.33)n4:0.22)n5:0.2)root;` 

<img src="examples/fig/ex.png" width="400"/>

then an example for partition file can be as follows:

    group1 1 {n1} {n5 n4}
    group2 1 {n3}

Each line defines a list of subtrees whose branches are supposed to have the same substitution rate. It starts by the name of the group (`group1`), then the prior proportion (`1`) of the group rate compared to the main rate. This is just a starting value, and the proportion will be estimated; giving an appropriate value helps to converge faster.
Each subtree is then defined between {}: the first node is the root of the subtree and the following nodes (if there any) define its tips. If the first node is a tip label then it takes the mrca of all tips as the root of the subtree.
If there's only root and not any tip defined, then the subtree is extended down to the tips of the full tree. Hence, {n1} defines the subtree rooted at the node n1; and {n5 n4} defines the subtree rooted at n5 that has one tip as n4 and other tips as the ones of the full trees (here are B,C). 
As a consequence, in this example, the branches will be partitioned into 3 groups such that each group has a different rate: 

- group1: (n1,A), (n1,D), (n5,n4), (n5,n2), (n2,B), (n2,C); 
- group2: (n3,F), (n3,G); 
- group0: the remaining branches of the tree. 
    
Note that if the internal nodes don't have labels, then they can be defined by mrca of at least two tips, for example n1 is mrca(A,D)

## Using variances

Variance is used to penalize long branch lengths. The variance formula of each branch v_i is proprtion to (b_i + b), where b (specified by option -b) is the pseudo constant added to adjust the dependency of variances to branch lengths. This parameter is a positive number < 1, and by defaul is the median branch length. It could be adjusted based on how much your input tree is relaxed. The smaller it is, the more variances are linear to branch lengths, which is more appropriate for strict clock tree. The bigger it is the less dependent of branch lengths on variances, which may be better for relaxed tree. Also, option `-v 2` may only be useful with strict clock tree.

## Some examples of command lines:

* If the input tree is rooted:

    - You want to estimate rate & dates under temporal constrained mode, using variances:

    `./lsd2 -i rootedtree_file -d date_file -c -v 1`

	- You want to remove outlier nodes with Zscore threshold 3:

    `./lsd2 -i rootedtree_file -d date_file -c -v 1 -e 3`

	- You know the tree partition where each part should have a different rate:

    `./lsd2 -i rootedtree_file -d date_file -c -v 1 -p parition_file`

	- You want to re-estimate the root position locally around the given root

    `./lsd2 -i rootedtree_file -d date_file -c -r l`

	- You want to calculate confidence intervals from 100 simulated trees. The sequence length used to build your tree was 1000, and you'd like to apply a lognormal relaxed clock of standard deviation 0.4 to the simulated branch lengths.

    `./lsd2 -i rootedtree_file -d date_file -c -r l -f 100 -s 1000 -q 0.4`
    
    (To calculate confidence intervals, the sequence length is required via option -s. The program generates simulated branch lengths using Poisson distributions whose mean equal to the estimated ones multiplied with sequence length. In addition, a lognormal relaxed clock is also applied to the branch lengths. This ditribution has mean 1 and standard deviation settable by users with option -q, by default is 0.2; 0 means strict clock. The bigger q is, the more your tree is relaxed and the bigger confidence intervals you should get).

	- If all tips are supposed to have the same date, you can still estimate the rate but only relative dates.
	
	`./lsd2 -i tree_file -c`

* If the input tree is unrooted, you should either specify outgroups or use option -r to estimate the root position.
	
	- If you don't have any outgroup and you want to estimate the root position:

    `./lsd2 -i unrootedtree_file -d date_file -c -r a -v 1`

	- If you have a list of outgroups and want to use them for rooting (also remove them):

    `./lsd2 -i unrootedtree_file -d date_file -g outgroup_file -c -v 1`
    
    - If you want to keep the outgroups in the tree, just estimate the root position on the branch that defined by the outgroups:
    
    `./lsd2 -i unrootedtree_file -d date_file -g outgroup_file -k -c -v 1`
    
## Output files: 

*.result* : contain the estimated rates, root date, possibly confidence intervals, outlier tips and the value of the objective function.

*.nexus* : trees in nexus format which contain information about the dates of internal nodes, branch lengths, and the confidence intervals (CI) if option -f was used.
    
*.date.nexus* : similar to *.nexus* trees, but branch lengths are rescaled to time unit by multiplying with the estimated rate. 

*.nwk* : trees in newick format 

## Citation
If you use this software, please cite: “ Fast dating using least-squares criteria and algorithms”, T-H. To, M. Jung, S. Lycett, O. Gascuel, Syst Biol. 2016 Jan;65(1):82-97.
