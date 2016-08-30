# Bins_BedGraph_files

Needs Python 3.5

Input: bedGraph files or bedGraph files that are compressed (tar.gz).

Output: .csv files with same size bins. Column Names: CHR, Start, End, Value

USE:
-c <int> [number of cpus]
-i <char> [The name of the folder that the bedGraph files are located]
-o <char> [The name of the folder that the output will be saved]
-s <int> [The bin size]
-t <int> [The number of character that you want to remove at the and of the initial file name]

NEED:
-i

Default Values:
-c 1
-o out
-s 200
-t 0

Example:
python bedGraph_rebin.py -c 5 -i test -o out -s 100 -t 16
