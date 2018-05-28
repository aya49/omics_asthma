BEGIN {FS = "\t"} # separator
{printf "%s%c", $1, "\t"} #get first column
{for (i = col; i <= NF; i += totalcol) printf ("%s%c", $i, i + totalcol <= NF ? "\t" : "\n");}
