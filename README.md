# statfish
This python library aims to make spatial transcriptomic analysis more efficient and more approachable. Several basic implementations for different processes of analysis are provided through a modular structure. The library focuses on customization of individual pipelines for a wide range of analysis tasks. 

This repository contains the raw code from my statfish Capstone project. To use the code, clone the repository and install dependencies listed in dependencies.txt. I have been using command prompts through mac terminal to execute the code, but usage should be relatively the same for others. cd into the directory containing the statfishpy root folder, this is where you can execute programs. 

Example usage: 
I have provided an example script in statfishpy/utils called demo_moran_i, which computes moran_i values for all genes in a dataset. Run the command "python -m statfishpy.utils.demo_moran_i" to execute this program. Note, you will have to change the filepath in the read_h5ad function of this program and you may have to adjust coordinate/volume keys. Currently the code only supports h5ad files. 

From this demo program you can infer usage of other tools in the library. Most process calls follow the same format:

```python
sf.spatial.build_spatial_graph(     #choose the process
    adata,                          #direct the function to your AnnData object
    x_key="center_x",               #provide variable keys as necessary
    y_key="center_y",
    k=30                            #adjust any other process parameters
)


