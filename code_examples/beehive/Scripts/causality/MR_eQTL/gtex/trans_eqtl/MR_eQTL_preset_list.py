#!/usr/bin/env python

import pandas as pd
import json
import numpy as np
from sys import argv, stderr
import json
import argparse
import subprocess
from collections import defaultdict

parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter, \
description=r"""

compute_t_MR.py uses Mendelian randomization to compute a 
statistic, here called t_MR, representative of the association
of a dependent variable (y) with independent variable (x) given
their association with a instrumental variable (z).

It takes as input z->y association results (e.g. trans-eQTL results),
z->x association results (e.g. best cis-eQTL results), 
dataframes of X, Y, and Z, and the pre-computed variance 
of the Mendelian randomization effect size (x->y) not yet scaled
by the difference in estimated y and true y.

The basic computation performed is as follows (in LaTeX notation):

\begin{eqnarray}
t_{MR} &=& \frac {\beta_{MR}} {var(\beta_{MR})} \\
&\text{where:}& \\
\beta_{MR} &=& \frac { \hat\beta_{y,z}} {\hat\beta_{x,z}} \\
var(\beta_{MR}) &=& \sigma^2 (x^T P_z x)^{-1} \\
\sigma^2 &=& \frac {(y - x \beta_{MR})^T (y - x \beta_{MR})} {n - \nu} \\
P_z &=& z(z^Tz)^{-1}z^T \\
\end{eqnarray}

""")

parser.add_argument("-i", "--trans_association_results", 
                    required = True,
                    dest="trans_association_results",
                    help='''
trans_association_results are header-less and space-delimited in the format:
y_ID z_ID beta_yz

NOTE: Exclude all cis results in trans_association_results.
Also may be best to exclude "local" results in trans_association_results,
that is, consider maintaining a buffer between cis 
results (e.g. < 100 Kb) and trans results (e.g. > 1-5 Mb) to make 
sure that trans results are not simply long distance cis results.
                   ''')
parser.add_argument("-X", "--expression_cis_dataframe", 
                    required = True,
                    dest="X",
                    help='''
X is a tab-delimited dataframe of the independent variable X
with header of sample names (e.g. individuals) and with row names
corresponding to each independent variable X (e.g. genes).
In the case of expression, each cell of the dataframe should 
be a floating point value.
                    ''')
parser.add_argument("-Y", "--expression_trans_dataframe",
                    required = True,
                    dest="Y",\
                    help='''
Y is a tab-delimited dataframe of the dependent variable Y
with header of sample names (e.g. individuals) and with row names
corresponding to each dependent variable Y (e.g. genes).
In the case of expression, each cell of the dataframe should 
be a floating point value.
                   ''')
parser.add_argument("-Z", "--genotypes_dataframe", 
                    required = True,
                    dest="Z",\
                    help='''
Z is a tab-delimited dataframe of the instrumental variable Z
with header of sample names (e.g. individuals) and with row names
corresponding to each instrumental variable Z (e.g. SNP ID).
In the case of genotypes, each cell of the dataframe should be 0,1,2
for the number of minor alleles for a SNP in a given individual.
                   ''')
parser.add_argument("-b", "--z_to_x_top_beta", 
                    required = True,
                    dest="z_to_x_top_beta",
                    help='''
z_to_x_top_beta is a header-less space-delimited dataframe 
of z->x associations in the format z_ID x_ID beta_xz, e.g.:
rs999 gene500 1.5.
                   ''')
parser.add_argument("-v", "--unscaled_var_beta_MR", 
                    required = True,
                    dest="unscaled_var_beta_MR",\
                    help='''
unscaled_var_beta_MR is a file containing a numpy array formed
using the script compute_unscaled_var_beta_MR.py and is 
of dimension |Z| by |X|.
                    ''')
parser.add_argument("-t", "--t_MRs_results_out", 
                    dest="t_MRs_results_out",\
                    help='''
All computed t_MRs are saved to the file t_MRs_results_out
as a numpy array.
                    ''')
parser.add_argument("-o", "--MR_association_results_out", 
                    dest="MR_association_results_out",\
                    help='''
Mendelian randomization results are saved to the file
MR_association_results_out in the tab-delimited header-less
format of x_ID y_ID t_MR.
                   ''')

args = parser.parse_args()

X = pd.read_csv(args.X, index_col=0, sep='\t')
Y = pd.read_csv(args.Y, index_col=0, sep='\t')
Z = pd.read_csv(args.Z, index_col=0, sep='\t')

assert(X.shape[1] == Z.shape[1])
assert(list(X.columns) == list(Z.columns))

X = np.array(X)
Z = np.array(Z)

inner = np.diag(np.diag(np.inner(Z,Z))) 
inner_inv = np.linalg.inv(inner)
unscaled_var_beta_MR = np.dot(np.dot(X,Z.T)**2, inner_inv).T

X = pd.read_csv(args.X, index_col=0, sep='\t')
Z = pd.read_csv(args.Z, index_col=0, sep='\t')

X_expr_dict = {gene:np.array(X.ix[gene]) for gene in list(X.index)}
Y_expr_dict = {gene:np.array(Y.ix[gene]) for gene in list(Y.index)}

z_to_x_top_beta = pd.read_csv(args.z_to_x_top_beta, header=None, delim_whitespace=True)

# v provided:
v = 3
n = X.shape[1]

MR_association_results_out = open(args.MR_association_results_out, 'w')

def compute_t_MRs(filename):
    
    total_num_lines = float(subprocess.check_output('wc -l < %s'%(filename),shell=True))
    t_MRs = []
    
    with open(filename, 'r') as f:
        for i, line in enumerate(f):
            # track progress
            if i % 100000 == 0:
                stderr.write('%s %0.1f%% complete\n'%(filename, 100 * i / total_num_lines)) 
            
            y_ID, z_ID, beta_yz = line.strip().split()
            beta_yz = float(beta_yz)            
            beta_xz, x_ID = z_to_x_top_beta[z_ID]
            
            i = row_dict[z_ID]                               
            j = col_dict[x_ID]
            
            y = Y_expr_dict[y_ID]
            x = X_expr_dict[x_ID]
            
            beta_MR = beta_yz / float(beta_xz)
            resid = (y - x * beta_MR)
            sigma_squared = np.inner(resid, resid) / ( n - v )
            var_beta_MR = sigma_squared / unscaled_var_beta_MR[i,j]
            t_MR = beta_MR**2 / var_beta_MR
            t_MRs.append(t_MR)
            MR_association_results_out.write('\t'.join([x_ID, y_ID, str(t_MR)]) + '\n')
            
    return(np.array(t_MRs))

t_MRs = compute_t_MRs(args.trans_association_results)

np.save(args.t_MRs_results_out, t_MRs)
MR_association_results_out.close()

