%%%GREEN FINANCE AND AMBIGUITY%%%%%
%%%MARCO CARLI%%%%%

The codes are solved in Dynare 5.5.
In order to retrieve all results under ambiguity of the papers, add the following code at line 311 to the file dyn_first_order_solver in the matlab folder in your Dynare installation: 

amb_idx = find(strcmp(DynareModel.endo_names, 'amb'));
amb_state_vec = ismember(dr.state_var,amb_idx);
amb_state_idx = find(amb_state_vec);
sh_ag_idx = find(strcmp(DynareModel.exo_names, 'sh_ag'));
dr.ghx(:,amb_state_idx) = dr.ghx(:,amb_state_idx)+dr.ghu(:,sh_ag_idx);


After this addition, baseline results can be retrieved running the .mod file ambig_envir_fin.
All policy exercises are run in the .mod file ambig_envir_fin_opt by activating one of the rules in lines 114-127 at a time.

Figures are stored in the Graphs folder and can be reproduced running the figures .m files.
Moment matching and data files are stored in the empirics folder.