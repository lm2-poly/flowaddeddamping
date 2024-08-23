# FlowAddedDamping
flowaddeddamping is a custom script to allow the evaluation of flow-added damping on a hydrofoil using NASTRAN. The related article detailing the methodology will be submitted soon.

## Installation
Use the requirements.txt file to install the required packages using pip. A virtual environment using anaconda or miniconda is recommended.

```bash
conda create --name flowaddeddamping python=3.9
conda activate flowaddeddamping
pip install -r requirements.txt
```

Once the requirements installed, the main.py script should be functional using the virtual environment created.

## Usage
The main.py file is the main script to use this program. The necessary instructions are available in the main.py file to explain the use of the script. 
### Important: If NASTRAN's location is not recognized, add "nastran_location = ..." with your nastranw.exe location in the run method of the model class. (line 169 in the original main.py script)

Currently, only standalone straight hydrofoils are supported, but add-ons will be made to consider cascades, camber and non-straight hydrofoils. The hydrofoils can be simulated for static analysis under a gravity load (mainly for debugging), vacuum modal analysis, uncoupled vibro-acoustic analysis (mainly for debugging), coupled vibro-acoustic analysis, aeroelastic analysis and hydroelastic analysis.

The test cases considered in the article are given in the TestFiles directory. If other test cases are desired, simply add a .dat file detailing the hydrofoil's profile following the same format as presented.

The plots, similar to those found in the article, are output after the simulations using the main.py script.

To obtain detailed information about the analysis performed, use Simcenter 3D to import the .op2 file. In Simcenter, use the "Import Simulation Input File" feature, select "Simcenter Nastran", select the "OP2 Binary" radio button and open the corresponding .op2 file in your files.

## Support
For any questions related to the presented code, please refer to danick.lamoureux@polymtl.ca and clement.audefroy@polymtl.ca.

## Authors and acknowlegementds
The current authors of the present codes are Danick Lamoureux and Clément Audefroy, under the supervision of Prof. Frédérick Gosselin from Polytechnique Montréal and Prof. Sébastien Houde from Université Laval.
