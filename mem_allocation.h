
/*
This header file written by Sabyasachi Sutradhar, Yale university (sabyasachi.Sutradhar@yale.edu) 2022, as a part of
Neuonal morphogenesis main code.
This code simulates branching morphogensis of Drosophila class-IV dendritic arbor. branches are generates from the origin (0,0)
and then the tips go through three state dynamics G-P-S with different rates and also new branches spawn from an existing branchesrandomly
with a rate.
This header executes the dynamic memory aloocation of relevant variables.
Copyright @ Sabyasachi Sutradhar
*/

/////////////////////////////////
Dendrite_Length=dvector(1,Max_Tips);
Xbase=dvector(1,Max_Tips);
Ybase=dvector(1,Max_Tips);
Xtip=dvector(1,Max_Tips);
Ytip=dvector(1,Max_Tips);
Dendrite_Direction=dvector(1,Max_Tips);
Deletion_index=ivector(1,Max_Tips);
state=ivector(1,Max_Tips);
v=dvector(1,Max_Tips);
