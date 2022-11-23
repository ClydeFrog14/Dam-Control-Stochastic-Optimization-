sets

         s " Scenarios " /low, medium, high/  ;




alias(s,s1) ;



parameters

         firstmonth(s) " First months expected scenarios "/ low -50, medium 100, high 250 /

         secondmonth(s) " Second months expected scenarios "/ low 50, medium 150, high 350 /

         probability(s) " Probabilities of scenarios "/ low 0.25, medium 0.50, high 0.25/

;



scalars

         initial_water_level  "Water Level at the beginning "/ 0 /

         dam_limit "When the water level exceeds 250, there will be a flood " /250/

         upper_release "Maximum amount of releasable water per month " /200/

         penalty_cost "Penalty cost of flood per mm of water level  " /10000/

         importation_cost "Importation cost of 1 mm water level  " /5000/

         small_m  "This will be used in binary constraints (If the flood amount or importation of water amount are under 5 mm, we can neglect the penalty cost of flood or importation cost of water per mm of water level)" /5/

         big_M "Big_M will be used in binary constraints  " /99999/

         ;


positive variables


x(s)          " Water level (mm) at the end of first month "

x1(s,s1)      " Water level (mm) at the end of second month "

y(s)          " Water released (mm) at the end of first month "

y1(s,s1)      " Water released (mm) at the end of second month "

z(s)          " Water imported (mm) at the end of first month "

z1(s,s1)      " Water imported (mm) at the end of second month "

fw(s)         " Water overflowed (mm) at the end of first month "

fw1(s,s1)     " Water overflowed (mm) at the end of second month "

;


binary variables



alfa(s)     " if flood happens at the end of first month, alfa is 1. otherwise 0 "

alfa1(s,s1) " if flood happens at the end of second month, alfa1 is 1. otherwise 0 "

beta(s)     " if importation of water happens at the end of first month, beta is 1. otherwise 0 "

beta1(s,s1)  " if importation of water happens at the end of second month, beta1 is 1. otherwise 0 "

;


variable

obj_1  " First objective functions value "

obj_2  " Second objective functions value "

;

equations

         obj_function_1  " Minimize expected costs "

         obj_function_2  " Minimize the probability of violating the maximum and minimum water levels "

         importation_constraint_lower(s) " Binary constraint for beta (if importation of water happens at the end of first month, beta is 1. otherwise 0)"

         importation_constraint_1_lower(s,s1)  " Binary constraint for beta1 (if importation of water happens at the end of second month, beta1 is 1. otherwise 0)"

         importation_constraint_upper(s)   " Binary constraint for beta (if importation of water happens at the end of first month, beta is 1. otherwise 0)"

         importation_constraint_1_upper(s,s1)   " Binary constraint for beta1 (if importation of water happens at the end of second month, beta1 is 1. otherwise 0)"

         flood_constraint_lower(s)  " Binary constraint for alfa (if flood happens at the end of first month, alfa is 1. otherwise 0)"

         flood_constraint_1_lower(s,s1)  " Binary constraint for alfa1 (if flood happens at the end of first month, alfa1 is 1. otherwise 0)"

         flood_constraint_upper(s) " Binary constraint for alfa (if flood happens at the end of first month, alfa is 1. otherwise 0)"

         flood_constraint_1_upper(s,s1) " Binary constraint for alfa1 (if flood happens at the end of first month, alfa1 is 1. otherwise 0)"

         dam_constraint(s)  " This equation restricts the water level for the first month "

         dam_constraint_1(s,s1)  " This equation restricts the water level for the second month "

         release_constraint(s)  " Constraint of maximum releasable water level  for the first month "

         release_constraint_1(s,s1)  " Constraint of maximum releasable water level  for the second month"

         equilibrium_constraint(s)   " Equations that connect water level at the beginning and at the end of first month (with the effects of flood, importation of water and scenario)"

         equilibrium_constraint_1(s,s1) " Equations that connect water level at the end of first month and at the end of second month (with the effects of flood, importation of water and scenario)"

;



obj_function_1.. obj_1 =e= penalty_cost *( sum(s, probability(s) * fw(s)) + sum((s, s1), probability(s) * probability(s1) * fw1(s,s1)   )) + importation_cost * ( sum(s, probability(s) * z(s)) + sum((s, s1), probability(s) * probability(s1) * z1(s,s1)   ))   ;

obj_function_2.. obj_2 =e= sum( s, probability(s) * alfa(s) ) + sum( (s,s1), probability(s) * probability(s1) * alfa1(s,s1) ) + sum( s, probability(s) * beta(s) ) + sum( (s,s1), probability(s) * probability(s1) * beta1(s,s1) ) ;

option limrow= 1000;


flood_constraint_lower(s).. alfa(s) * small_m =l= fw(s)  ;

flood_constraint_1_lower(s,s1).. alfa1(s,s1) * small_m =l= fw1(s,s1)    ;

option limrow= 1000;


flood_constraint_upper(s).. big_M * alfa(s) =g= fw(s)         ;

flood_constraint_1_upper(s,s1)..  big_M * alfa1(s,s1) =g= fw1(s,s1)    ;

option limrow= 1000;


importation_constraint_lower(s)..  small_m * beta(s) =l= z(s)    ;

importation_constraint_1_lower(s,s1).. small_m * beta1(s,s1) =l= z1(s,s1)   ;

option limrow= 1000;



importation_constraint_upper(s).. big_M * beta(s) =g= z(s);

importation_constraint_1_upper(s,s1).. big_M * beta1(s,s1) =g= z1(s,s1);

option limrow= 1000;


dam_constraint(s).. x(s) =l= dam_limit ;

dam_constraint_1(s,s1).. x1(s,s1) =l= dam_limit ;

option limrow= 1000;


release_constraint(s)..     y(s) =l= upper_release     ;

release_constraint_1(s,s1)..   y1(s,s1) =l= upper_release     ;

option limrow= 1000;


equilibrium_constraint(s).. x(s) =e= initial_water_level + firstmonth(s) - y(s) + z(s) - fw(s) ;

equilibrium_constraint_1(s,s1).. x1(s,s1) =e= x(s1) + secondmonth(s) - y1(s,s1) + z1(s,s1) - fw1(s,s1) ;

option limrow= 1000;





model dam_control_obj_1 /obj_function_1, dam_constraint, dam_constraint_1, release_constraint, release_constraint_1, equilibrium_constraint, equilibrium_constraint_1 /;

solve dam_control_obj_1 using lp minimizing obj_1;


model dam_control_obj_2 /obj_function_2, dam_constraint, dam_constraint_1, release_constraint, release_constraint_1, equilibrium_constraint, equilibrium_constraint_1, importation_constraint_lower, importation_constraint_1_lower, importation_constraint_upper, importation_constraint_1_upper, flood_constraint_lower, flood_constraint_1_lower, flood_constraint_upper, flood_constraint_1_upper/;

solve dam_control_obj_2 using mip minimizing obj_2;
