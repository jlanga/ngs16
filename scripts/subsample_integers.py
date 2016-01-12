#!/usr/bin/env python3

import numpy as np
import sys



def subsample_integers(set_size= 1000, n_subsamples= 10, subsample_sizes= 1000, seed= 1):
    """
    Make a sample without replacement of a set of size `set_size`.
    `subsample_sizes` should be a list of integers that denote the size of the 
        size of each subsample.
    Note that all subsamples are disjoint.
    """
    # Check inputs, check sizes of the subsamples and set sizes....
    # TODO
    
    np.random.seed(seed)
    
    subsamples = np.repeat(subsample_sizes, n_subsamples)
    
    # Create a subsample of the original set
    unsorted_array = np.random.choice(
        np.arange(set_size, dtype= np.uint32),
        size = sum(subsamples),
        replace= False
    )
    
    # Split the subsample, and throw away the last element of the np.split 
    # object (it contains the remaining parts)
    unsorted_samples = np.split(
        unsorted_array,
        np.cumsum(subsamples)[0:len(subsamples)]
    )
    
    # Sort each of the arrays
    sorted_samples = list(map(np.sort, unsorted_samples))
    
    return sorted_samples



def write_indexes(subsamples, prefix="./unnamed"):
    for i in range(len(subsamples) - 1 ):
        filename = prefix + "_%d.idx" % (i + 1)
        with open(filename, "w") as f:
            f.write("\n".join(map(str, subsamples[i])))



if __name__ == "__main__":
    
    usage = "subsample_integers.py tag out_dir set_size n_parts part_size seed"
    example_usage = "integer_sampling.py test 1000 10 100 42"
    
    if len(sys.argv) != 6:
        sys.exit(
            "ERROR! Wrong number of parameters\n" + 
                usage + "\n" +
                example_usage)
        
    prefix    = sys.argv[1]
    set_size  = int(sys.argv[2])
    n_parts   = int(sys.argv[3])
    part_size = int(sys.argv[4])
    seed      = int(sys.argv[5])
    
    subsamples = subsample_integers(
        set_size,
        n_parts, 
        part_size,
        seed
    )
    
    write_indexes(
        subsamples,
        prefix
    )

