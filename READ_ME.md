# Code for the thesis "Towards Less Greedy Quantum Coalition Structure Generation in Induced Subgraph Games"

This project contains the code for the master's thesis entitled 
"Towards Less Greedy Quantum Coalition Structure Generation in Induced Subgraph Games"
which was written by Daniëlle Schuman.

Most parts of the code in this project have been written by Daniëlle Schuman as well, with a few exceptions:
- Some methods in this project have been written by, or are adapted from, Code written by Jonas Nüßlein. These Methods are marked accordingly.
- The code in the folder `classical/best-partition-master` has been downloaded from `https://github.com/Alexander-Belyi/best-partition`. It has been written by Alexander Belyi as part of his paper "[Subnetwork Constraints for Tighter Upper Bounds and Exact Solution of the Clique Partitioning Problem](https://arxiv.org/abs/2110.05627)" and we introduced only a very small number of modifications to it which are marked accordingly.
- For some methods, especially a large part of those contained in the file `plotting.py`, we used the Large-Language Model ChatGPT (https://chat.openai.com/auth/login) to generate boilerplate code for basic functionalities like the reading and writing of files, the parsing of lists, dicts and strings or the generation of matplotlib plots. We subsequently checked this code carefully to ensure its correct functionality, and in the vast majority of the cases then adapted it to enhance its functionality (fixing bugs where necessary), improve readability regarding semantics and at times improve efficiency.


## Installing and running the code

### Quantum(-inspired) code

To run the quantum and quantum-inspired code in this project, first install the requirements using

```
pip install -r requirements.txt
```
Subsequently, install the minorminer package using
```
pip install minorminer==0.2.13
```
and ignore the resulting dependency conflict.

If you want to run the version of the code using the D-Wave quantum annealer, you furthermore need to
- configure the UQO middleware according to the instructions in `secrets_folder/configure_uqo.py`.
- create a file called `dwave_token.py` in the `secrets_folder`. This file needs to contain one line with `TOKEN = '...'` where ... is the token of your D-Wave Leap account. You can create a free D-Wave Leap account at https://cloud.dwavesys.com/leap/signup/. (Our algorithm will not use any of your QPU runtime, it just needs cloud access to up-to-date information about the quantum annealer's hardware graph.)

Our algorithms can then be run using:
```
python main.py
```
This can be passed the options:
- `--seed`: pass an integer to use as a custom seed for the simulator. If no seed is passed, the code will run all seeds from our experiments of the respective solver.
- `--solver`: The solver to run the algorithms with. Options are `"qbsolv"`, `"qaoa"`or `"dwave"`. The first one is the default.

Notice that the algorithms have long runtimes on either solver.

### Classical baseline
To install run the Belyi's code, which we used as a baseline, on our data, perform the following steps (mostly taken from Belyi's READ_ME):

1. Install [CMake](https://cmake.org/download/);
2. Install [CBC](https://github.com/coin-or/Cbc) as a solver. To do that, follow the instructions from the [official page](https://github.com/coin-or/Cbc). On Mac, Belyi et al. recommend using [Homebrew](https://brew.sh/), on ubuntu/debian use ```sudo apt install coinor-libcbc-dev```. Subsequently, adapt `classical/best-partition-master/cmake/Modules/FindCBC.cmake`, by replacing the paths to the different necessary components of CBC with those they have on your computer.
3. You can also (or additionally) use CPLEX if you have it (we did for our experiments). CMake should find it if installed in the default location, otherwise provide the path in `CPLEX_ROOT`.
4. Now, delete all contents from `classical/best-partition-master/build`. Then navigate to this folder on your command line and run CMake: 
```bash
cmake ..
```
5. Run ```make``` on your commandline;
6. Run built program: ```./BestPartition```.


## Our experiments
Notice that while the current configuration of `main.py` should be usable to reproduce our experiments, we used several older versions of the code with hardcoded configurations (especially regarding what graphs to run when) to produce the results in the thesis. A lot of them can be found in the other files in this project whose name starts with `main_`.

The results from our experiments can be found in the folder `results` in this project, with the respective command line outputs in `run_txt`. Visualizations of our experiments can be found in `plots`.