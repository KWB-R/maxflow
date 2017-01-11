# Maxflow

##Requirements
- Installation of Flopy release 3.2.5  and
- Manual replacement of **"mfmnw2.py"** file (stored under: 
**"YOURPYTHONPATH\\Lib\\site-packages\\flopy\\modflow"**) with the following one: [https://github.com/modflowpy/flopy/blob/383438122fdb97623e9666fb11c8d9b84c9563d2/flopy/modflow/mfmnw2.py](https://github.com/modflowpy/flopy/blob/383438122fdb97623e9666fb11c8d9b84c9563d2/flopy/modflow/mfmnw2.py)). 


***If the second step is not performed Flopy will not be able to write the 
"wellfield.mnw2" file and the modelling will not be possible !!!***