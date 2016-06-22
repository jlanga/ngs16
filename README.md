---
title: Which tissues do we have to sequence to get more info per $?
author: Jorge Langa
email: jorgeeliseo.langa@ehu.eus, jorge.langa.arranz@gmail.com
---

To start, you should start creating a virtualenv with python3.4 (not tested with python2.7):

```
cd drer_ngs
virtualenv --python=python3.4 .
```

Numpy is required for the subsampling steps:

```
pip install numpy
```

The sofware that I'm using is:

- [trimmomatic-0.35]()

- [hisat2-2.0.1]()

- [samtools-0.1.19]()

- [bwa-0.7.19]()

To install them you can run `scripst/install_software.sh`. The executables will be in `./bin`.


