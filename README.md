# banded_visualization
Software for banded matrix visualization of binary data sets

## Example ##

View all the options:

```
python banded_visualization.py --help
```

Running with default parameters and showing resulting images to a popup window:

```bash
python -m banded_visualization <mode> <path-to-folder-with-data>
```
Running the included `tweets` example showing the clusters:

```bash
python -m banded_visualization clusters path-to-folder-tweets 
```
Running the included `tweets` example showing the rules and outputting the results into a folder:

```bash
python -m banded_visualization clusters path-to-folder-tweets -o path-to-output-folder/name
```
This results in the program to create files `name_cluster_1.png, name_cluster_2.png, ..., name_cluster_k.png` in the folder`path-to-output-folder`
