###############################################################################################
########                    MD simulation of a 2D cell model                        ###########
########                      Effects of mutating one cell                          ###########
########                   Script by J. R. Bordin - Oct 2022                        ########### 
########                         Powered by ESPRESSO                                ###########
###############################################################################################
set systemTime [clock seconds]
puts "The initial time is: [clock format $systemTime -format %H:%M:%S]"
puts "The initial date is: [clock format $systemTime -format %D]"


set pi [expr 4.*atan (1.0) ];#set value of pi
t_random seed [ pid ] ;#starts the random number generator
set sigma 1.0 ;#distance unit, in micrometers


#### four main parameters
set kstr 10.0 ;#coil strength that controls the cell deformation
set force 5.00 ;#constant force pulling in the x-direction
set phi 0.60 ;#cell packing fraction = Acell/Atot
set eps 2.00 ;#adhesion interaction parameter

#Set the cell radius, number of and area occupied by the  cells
set radius 3.0 ;#cell radius in units of sigma 
set n_cell 101 ;#total number of cells
set Acell [expr $pi*($radius+0.5)*($radius+0.5)] ;#are of one cell
set Acell [expr $n_cell*$Acell] ;#area occupied by the cells
set Atot [expr $Acell/$phi]

#Set the simulation box
set box_ly [expr pow($Atot/2.,0.5)]  ;# in units of sigma
set box_lx [expr 2.*$box_ly] ;# in units of sigma
set box_lz 50.0
setmd box_l $box_lx $box_ly $box_lz;#defines the initial simulation box
setmd periodic 1 1 0 ;#PBC in x and y directions

cellsystem domain_decomposition

#set the temperature and pressure 
set temp 1.00
set gamma_0 1.0 ;#same dumping parameter from Langevin/DPD thermostat
setmd skin 0.4;#for thermostat
#sets the thermosthat
thermostat set langevin $temp $gamma_0 ;# set the Langevin thermosthat

#set the integration data
set time_step 0.005
setmd time_step $time_step; #MD time step
#warm up steps
set integ_steps_warm 1000
#equilibration steps
set integ_steps_equil 100000
#production cycles
set integ_steps 200
set total_cycles 2000
set run_time [expr $integ_steps*$total_cycles*$time_step]

#each species needs a number to "define" his type. here I set a number to each species.
set cell 0 ;#the polymer ring
set ghost 1 ;#the central bead
set cell2 2 ;#non adhesive cell
set ghost2 3 ;#the central bead
set insert 4 ;#for insertion try



######################defining the bonded interactions


#bond interaction for two neighbour cell beads
set sig 1.0 ;#size of each bead
set ring 21
set lambda [expr 1.0*$sig]
inter $ring harmonic 100.0  $lambda  [expr 5.*$lambda]

#bond with the central bead
set rigid 22
inter $rigid harmonic 100.0 $radius [expr 2*$radius]

###for the mutated cell
#bond interaction for two neighbour cell beads
set sig 1.0 ;#size of each bead
set ring2 212
set lambda [expr 1.0*$sig]
inter $ring2 harmonic $kstr $lambda  [expr 5.*$lambda]

#bond with the central bead
set rigid2 222
inter $rigid2 harmonic $kstr $radius [expr 2*$radius]


#################################### defining the non-bonded interactions

set WCA_cut [expr pow(2,1./6.)] ;#repulsive WCA interaction
set r_cut 2.5

inter $cell $cell lennard-jones $eps [expr 1.*$lambda] [expr 1.*$lambda*$r_cut] auto

inter $cell2 $cell2 lennard-jones $eps [expr 1.*$lambda] [expr 1.*$lambda*$WCA_cut] auto

inter $cell $cell2 lennard-jones $eps [expr 1.5*$lambda] [expr 1.5*$lambda*$WCA_cut] auto

inter $ghost $ghost lennard-jones [expr 3.*$eps]  $radius [expr $radius*$WCA_cut] auto


#placing the ring-like cell#placing the ring-like cell
set nions_per_ring  [expr int(((2.*$pi*$radius)/$lambda))]
#set nadhesive [int($fad*$nions_per_ring)]
set teta  [expr 2.0*$pi/($nions_per_ring-1)]
set ns 0
for {set cell_counter 0} { $cell_counter < $n_cell } {incr cell_counter 1} {
    puts "$cell_counter"
    #insert the i-th ring
    if {$cell_counter == 0} {
	set posy0 [expr $box_ly*0.5]
	set posx0 [expr $box_lx*0.5]
	set posz0 [expr $box_lz*0.5]
	part $ns pos $posx0 $posy0 $posz0 q 0 type $ghost mass 1.0 fix 0 0 1;# first include the central bead for each ring
	part $ns mol $cell_counter
    }
    if {$cell_counter > 0} {
	set overlap 1
	set try 0
	while { $overlap == 1} {
	    incr try 1
	    if {$try > 100000} {
		puts "more than 100000 attempts were made to include the cell $cell_counter"
		puts "try to put fewer cells... bye bye"
		exit
	    }
	    set posy0 [expr $box_ly*[t_random]]
	    set posx0 [expr $box_lx*[t_random]]
	    set posz0 [expr $box_lz*0.5]
	    part $ns pos $posx0 $posy0 $posz0 q 0 type $ghost mass 1.0 fix 0 0 1;# first include the central bead for each ring
	    part $ns mol $cell_counter
	    set overlap_test2 [analyze mindist $ghost2 $ghost]	  
	    set overlap_test [analyze mindist $ghost $ghost]
	    if {$overlap_test > [expr 2.*$radius] && $overlap_test2 > [expr 2.*$radius]} {
		set overlap 0
	    }
	}
    }
    set gid $ns
    set ns [expr $ns+1]
    for {set j 1 } { $j < [expr $nions_per_ring] } {incr j 1} {
	#insert the ring
	set posz [expr $box_lz*0.5]
	set posy [expr $posy0 +($radius)*cos($j*$teta)]
	set posx [expr $posx0 +($radius)*sin($j*$teta)]
	if {$cell_counter == 0} {
	    part $ns pos $posx $posy $posz q 0 type $cell mass 1.0 fix 0 0 1
	    part $ns mol $cell_counter	
	    part $ns bond $rigid $gid ;#bond the bead with the inner monomer
	    if { $j > 1 } {
		part $ns bond $ring [expr $ns-1]
	    }
	}
	if {$cell_counter > 0} {
	    part $ns pos $posx $posy $posz q 0 type $cell mass 1.0 fix 0 0 1
	    part $ns mol $cell_counter	
	    part $ns bond $rigid $gid ;#bond the bead with the inner monomer
	    if { $j > 1 } {
		part $ns bond $ring [expr $ns-1]
	    }
	}
	set ns [expr $ns+1]
    }
}

for { set j 0 } { $j < [setmd n_part] } {incr j [expr $nions_per_ring]} {
    part [expr $j+1] bond $ring [expr $nions_per_ring+$j-1]
}
set n_cell_part [setmd n_part]
part auto_exclusions $nions_per_ring 

set systemTime [clock seconds]
puts "End of insertion step: [clock format $systemTime -format %H:%M:%S]"

analyze set chains 0 $n_cell [expr $nions_per_ring+1]

set Nmolecules [setmd n_part]

######################### Initial warm steps with caped interaction  #############################################
inter forcecap 200.0
integrate $integ_steps_warm
inter forcecap 0
inter $ghost $ghost lennard-jones $eps 0.0 $WCA_cut auto

###################### End of Warm steps #############################################

#define a random number for each pressure
t_random seed [ pid ]
#open the files
set id "phi-$phi-k-$kstr-f-$force-eps-$eps"
set xyzfile [open "snap-$id.xyz" "w"] ;#open the xyz file
set trajfile [open "traj-$id.xyz" "w"] ;#open the CM trajectory file
set velfile [open "vel-$id.dat" "w"] ;#open the CM velocity file
puts $velfile [expr ($total_cycles)]
set obsfile [open "obs-$id.dat" "w"] ;#open the observable file
######################### Equilibration steps #############################################
integrate $integ_steps_equil
set equil_time [setmd time]
part 0 ext_force $force 0. 0.
###################### End of Equilibration steps #############################################

puts $obsfile "\#Langevin Dynamics simulation of 2D cells in a solvent"
puts $obsfile "\#Powered by Espresso! - tcl script by J. R. Bordin - April 22"
puts $obsfile "\#---------------------------------------------------------------------"
puts $obsfile "\# $n_cell cells"
puts $obsfile "\# Total number of particles is [setmd n_part] "
puts $obsfile "\# Temperature = $temp and damping parameter = $gamma_0 and packing fraction = $phi"

set agreg_file [open "agreg-$id.dat" "w"]
puts $agreg_file "\#--Time -- N_aggreg -- <n_mol> ---- err(+-) - N_max(aggreg) - N_min(aggreg)"

puts $xyzfile "   $n_cell_part "
puts $xyzfile "    Atoms   "
for {set nion 0} { $nion < [setmd n_part] } { incr nion } {
    set x1  [expr $sigma*[lindex [part $nion print folded_position] 0]]
    set y1  [expr $sigma*[lindex [part $nion print folded_position] 1]]
    set z1  [expr $sigma*[lindex [part $nion print folded_position] 2]]
    if { [lindex [part $nion ] 6] == $cell } {
	puts $xyzfile " O   $x1   $y1   $z1"
    } 
    if { [lindex [part $nion ] 6] == $cell2 } {
	puts $xyzfile " N   $x1   $y1   $z1"
    }
    if { [lindex [part $nion ] 6] == $ghost } {
	puts $xyzfile " H  $x1   $y1   $z1"
    }
    if { [lindex [part $nion ] 6] == $ghost2 } {
	puts $xyzfile " P  $x1   $y1   $z1"
	set vx1  [expr [lindex [part $nion print v] 0]]
	set vy1  [expr [lindex [part $nion print v] 1]]
	set v2 [expr $vx1*$vx1+$vy1*$vy1]
	puts $velfile "[expr (pow($v2,1./2.))]"
	flush $velfile
    }
    
    flush $xyzfile
}

set x00 [expr $sigma*[lindex [part 0 print pos] 0]]
set y00 [expr $sigma*[lindex [part 0 print pos] 1]]
set z00 [expr $sigma*[lindex [part 0 print pos] 2]]
puts $trajfile " [expr ([setmd time] - $equil_time)]   0. 0."
flush $trajfile


###################main simulation steps - production of results################################################
#first, mutate the first cell

part 0 type $ghost2
part 0 bond delete

for { set i 1 } { $i < $nions_per_ring } {incr i } {
    part $i bond delete
    part $i type $cell2
    part $i bond $rigid2 0
    if { $i > 1 } {
	part $i bond $ring2 [expr $i-1]
    }
}
part [expr $nions_per_ring-1] bond $ring2 1

set vel [observable new particle_velocities type [list $ghost]]
set vacf [correlation new obs1 $vel corr_operation scalar_product tau_max [expr 2000.*$time_step] dt [expr $time_step]]
correlation $vacf autoupdate start
set vel2 [observable new particle_velocities type [list $ghost2]]
set vacf2 [correlation new obs1 $vel2 corr_operation scalar_product tau_max [expr 2000.*$time_step] dt [expr $time_step]]
correlation $vacf2 autoupdate start

set stress [observable new stress_tensor]

for {set i 1} { $i <= $total_cycles} { incr i} {
    
    integrate $integ_steps
    analyze append
    set energy [analyze energy];#evaluates the energy
    set u [lindex [lindex $energy 0] 1] ;#potential energy
    set k [lindex [lindex $energy 1] 1] ;#kinectic energy
    set pxx [expr [lindex [observable $stress print] 0]] ;#xx component of the stress tensor 
    set pyy [expr [lindex [observable $stress print] 4]] ;#yy component of the stress tensor
    set pxy [expr 0.5*($pxx+$pyy)] ;#xy component of the stress tensor
    puts $obsfile [format "%.3e %.5e %.5e %.5e " [expr ([setmd time] - $equil_time)]  $pxy  [expr $u/[setmd n_part]]  [expr $k/[setmd n_part]]   ] ;#write the quantities to the obs file
    flush $obsfile

    set agreg [analyze aggregation 2.5 0 $n_cell]
    puts $agreg_file [format "%.3e %.5e %.5e %.5e %.5e %.5e" [setmd time] [lindex $agreg 9] [lindex $agreg 5]  [lindex $agreg 7] [lindex $agreg 1] [lindex $agreg 3] ]
    flush $agreg_file
    
    puts $xyzfile "   $n_cell_part "
    puts $xyzfile "    Atoms   "
    for {set nion 0} { $nion < [setmd n_part] } { incr nion } {
	set x1  [expr $sigma*[lindex [part $nion print folded_position] 0]]
	set y1  [expr $sigma*[lindex [part $nion print folded_position] 1]]
	set z1  [expr $sigma*[lindex [part $nion print folded_position] 2]]
	if { [lindex [part $nion ] 6] == $cell } {
	    puts $xyzfile " O   $x1   $y1   $z1"
	} 
	if { [lindex [part $nion ] 6] == $cell2 } {
	    puts $xyzfile " N   $x1   $y1   $z1"
	    }
	if { [lindex [part $nion ] 6] == $ghost } {
	    puts $xyzfile " H  $x1   $y1   $z1"
	}
	if { [lindex [part $nion ] 6] == $ghost2 } {
	    puts $xyzfile " P  $x1   $y1   $z1"
	    set vx1  [expr [lindex [part $nion print v] 0]]
	    set vy1  [expr [lindex [part $nion print v] 1]]
	    set v2 [expr $vx1*$vx1+$vy1*$vy1]
	    puts $velfile "[expr (pow($v2,1./2.))]"
	    flush $velfile
	    set x11 [expr $sigma*[lindex [part $nion print pos] 0]]
	    set y11 [expr $sigma*[lindex [part $nion print pos] 1]]
	    set z11 [expr $sigma*[lindex [part $nion print pos] 2]]
	    puts $trajfile " [expr ([setmd time] - $equil_time)]   [expr $x11-$x00]   [expr $y11-$y00]"
	flush $trajfile
	}
	
	flush $xyzfile
    }
}

################################## end of main steps

correlation $vacf write_to_file "vacf-cor-$id.dat"
correlation $vacf2 write_to_file "vacf2-cor-$id.dat"

set systemTime [clock seconds]
puts "Production end time is: [clock format $systemTime -format %H:%M:%S]"
puts "Production end date is: [clock format $systemTime -format %D]"


#close $vacf_file
close $xyzfile
close $velfile
close $obsfile
close $agreg_file
#################end of program######################
exit


