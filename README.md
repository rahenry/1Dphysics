This is a body of code used to manage simulations using the tenpy library. It also includes many other functions such as plotting utilities, HPC script generation/submission, exact diagonlisation code, functions to compute various exact results, and exceptional point detection. It is largely undocumented and has been uploaded for the use of other students.

I recommend using (mini)conda to create an environment with
`conda create --name <name> python pip`
`conda activate <name>`
`conda install --channel=conda-forge physics-tenpy`

A simulation or body of simulations is specified by a file in the data directory, which can be run using
`sh run <name>`
