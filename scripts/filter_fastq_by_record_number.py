import sys
import numpy as np

def filter_fastq(index_file):
    """
    Filter according to a file of integers. Each integer is the RECORD NUMBER, 
        not the line number.
    """
    with open(index_file, "r") as f_idx:
        
        # Read the entire index as a set
        index = set([int(x.strip()) for x in f_idx])
                
        # Iterate over the infile, if record number is in index, print
        line_number = 0
        for line in sys.stdin:
            if line_number / 4 in index:
                # Put together the line and the following 3
                record = "".join([
                    line,
                    sys.stdin.readline(),
                    sys.stdin.readline(),
                    sys.stdin.readline()])
                yield(record)
            else:
                # Move three lines
                sys.stdin.readline()
                sys.stdin.readline()
                sys.stdin.readline()
            line_number += 4



def write_fastq(iterator):
    """
    Just print whatever comes from the iterator.
    """
    for record in iterator:
        sys.stdout.write(record)



if __name__ == "__main__":
    
    usage = "filter_fastq_by_record_number.py indexes.txt < in.fastq > out.fastq"
    
    if len(sys.argv) != 2:
        sys.stderr("ERROR! Incorrect number of inputs.\n" + usage )
    
    index_file = sys.argv[1]
    
    iterator = filter_fastq(index_file)
    
    write_fastq(iterator)
