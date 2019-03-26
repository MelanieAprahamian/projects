**TUTORIAL**

The following is a tutorial for utilizing the prediction script (TargetID.py) described in the publication:
 Kim, SS, Aprahamian, ML, Lindert, S. Improving inverse docking target identification with Z‐score selection. Chem Biol Drug Des. 2019; 00: 1– 12. https://doi.org/10.1111/cbdd.13453 .

Text presented in fixed width are commands to be entered into the terminal.

*Example input and output files are included in the folder Tutorial.*

Tutorial: Prediction of Target/Ligand Pairs from Docking Scores

This tutorial provides all the information and steps necessary to run and analyze docking results from an inverse docking experiment to predict the correct target/ligand pairs. To utilize this script, you must have first performed an inverse docking experiment and compiled all of the docking scores for each potential target/ligand pair. 

1.	Prepare your working directory that will contain all of the input and output files along with the python script. Download the Tutorial directory included with the Supporting Information and cd into it:

`> cd Tutorial`

2.	For this tutorial you can use either the included example input file (example_input.txt) or generate your own. This file should contain 3 white space separated columns consisting of the target name, ligand name, and docking score. The file should contain as many lines as there are docking scores. If the docking software used does not always provide scores for poses with very high scores, please just omit the missing results from list.

3.	Execute the script using the following command specifying your input file using the -f flag:

`> python TargetID.py -f example_input.txt`

4.	The script will evaluate the combined z-score for each entry in the input file and print the predicted target/ligand pairs to the screen. Printed to the screen will be something like the following:

`> Predicted target/ligand pairs:
target_A ligand_A
target_B ligand_B`

Please also see the file example_output.txt for another example of what the output should look like.

5.	 If you wish to have the scripts predictions written to a file rather than just to the terminal window, please use the following command instead:

`> python TargetID.py -f example_input.txt > output.txt`

6.	OPTIONAL: If you are interested in more than just the target/ligand pairings, the -verbose flag can be included in the command line to additionally print the calculated average and standard deviations of the docking scores across both the targets and the ligands (used in the Z-score calculations) along with the evaluated combined Z-score for each target/ligand pair. An example of this output is included in the file example_output_verbose.py. The exact command to obtain the full output is:

`> python TargetID.py -f example_input.txt -verbose`
