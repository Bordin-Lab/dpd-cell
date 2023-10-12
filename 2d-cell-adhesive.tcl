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


#### two main parameters: spring constant and force pulling the mutated cell
set kstr 2.0 ;#coil strength for the bond between the central ghost bead and the ring. It controls the cell deformation
set force 5.00
#Set the simulation box
set box_lx 180.0  ;# in units of sigma
set box_ly 180.0 ;# in units of sigma 
set box_lz 180.0 
setmd box_l $box_lx $box_ly $box_lz;#defines the initial simulation box
setmd periodic 1 1 0 ;#PBC in x and y directions
set Atot [expr $box_lx*$box_ly]

#Set the cell radius and number of cells
set radius 10.0 ;#cell radius in units of sigma 
set n_cell 26 ;#total number of cells
set Acell [expr $pi*($radius+0.5)*($radius+0.5)] ;#are of one cell
set Acell [expr $n_cell*$Acell] ;#area occupied by the cells



cellsystem domain_decomposition

#set the temperature and pressure 
set temp 0.50
set press 0.01
set gamma_0 1.0 ;#same dumping parameter from Langevin/DPD thermostat
set gamma_v 0.001 ;#coupling parameter for the volume control
set piston_mass 0.0001 ;#mass of the piston that will control the pressure
setmd skin 0.4;#for thermostat
#sets the thermosthat
thermostat set npt_isotropic $temp $gamma_0 $gamma_v ;# set the NPT simulation
integrate set npt_isotropic $press $piston_mass 1 1 0

#set the integration data
set time_step 0.01
setmd time_step $time_step; #MD time step
#warm up steps
set integ_steps_warm 5000
#equilibration steps
set integ_steps_equil 5000
#production cycles
set integ_steps 500
set total_cycles 1000
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
inter $ring harmonic 20.0  $lambda  [expr 5*$lambda]

#bond with the central bead
set rigid 22
inter $rigid harmonic 20.0 $radius [expr 5.*$radius]

###for the mutated cell
#bond interaction for two neighbour cell beads
set sig 1.0 ;#size of each bead
set ring2 212
set lambda [expr 1.0*$sig]
inter $ring2 harmonic $kstr $lambda  [expr 5.*$lambda]

#bond with the central bead
set rigid2 222
inter $rigid2 harmonic $kstr $radius [expr 2.*$radius]


#################################### defining the non-bonded interactions

set WCA_cut [expr pow(2,1./6.)] ;#repulsive WCA interaction
set r_cut 2.5
set eps 0.50 ;#interaction parameter

inter $cell $cell lennard-jones $eps [expr 1.*$lambda] [expr 1.*$lambda*$r_cut] auto

inter $cell2 $cell2 lennard-jones $eps [expr 1.*$lambda] [expr 1.*$lambda*$WCA_cut] auto

inter $cell $cell2 lennard-jones $eps [expr 1.5*$lambda] [expr 1.5*$lambda*$WCA_cut] auto

#placing the ring-like cell
set nions_per_ring  [expr int(((2.*$pi*$radius)/$lambda))]
#set nadhesive [int($fad*$nions_per_ring)]
set teta  [expr 2.0*$pi/($nions_per_ring-1)]
set ns 0
for {set cell_counter 0} { $cell_counter < $n_cell } {incr cell_counter 1} {
    #insert the i-th ring
    if {$cell_counter == 0} {
	set posy0 [expr $box_ly*[t_random]]
	set posx0 [expr $box_lx*[t_random]]
	set posz0 [expr $box_lz*0.5]
	part $ns pos $posx0 $posy0 $posz0 q 0 type $ghost2 mass 1.0 fix 0 0 1;# first include the central bead for each ring
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
	    if {$overlap_test > [expr 2.*$radius+0.5] && $overlap_test2 > [expr 2.*$radius+0.5]} {
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
	    part $ns pos $posx $posy $posz q 0 type $cell2 mass 1.0 fix 0 0 1
	    part $ns mol $cell_counter	
	    part $ns bond $rigid2 $gid ;#bond the bead with the inner monomer
	    if { $j > 1 } {
		part $ns bond $ring2 [expr $ns-1]
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
part auto_exclusions 10

set n_cell_part [setmd n_part]
set n_part [setmd n_part]

set systemTime [clock seconds]
puts "End of insertion step: [clock format $systemTime -format %H:%M:%S]"

analyze set chains 0 $n_cell $nions_per_ring

set Nmolecules [setmd n_part]

######################### Initial warm steps with caped interaction  #############################################
inter forcecap 200.0
integrate $integ_steps_warm
inter forcecap 0
###################### End of Warm steps #############################################

#define a random number for each pressure
t_random seed [ pid ]
#open the files
set id "k-$kstr-f-$force"
set xyzfile [open "snap-$id.xyz" "w"] ;#open the xyz file
set trajfile [open "traj-$id.xyz" "w"] ;#open the CM trajectory file
set velfile [open "vel-$id.dat" "w"] ;#open the CM velocity file
puts $velfile [expr ($total_cycles)]
set obsfile [open "obs-$id.dat" "w"] ;#open the observable file
set rhofile [open "rho-$id.dat" "w"] ;#open the rho x p file
set rg_file [open "rg-$id.dat" "w"];#open the rg x p file
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
puts $obsfile "\# Temperature = $temp and damping parameter = $gamma_0 and pressure = $press"

set vel [observable new particle_velocities type [list $ghost]]
set vacf [correlation new obs1 $vel corr_operation scalar_product tau_max [expr 2000.*$time_step] dt [expr $time_step]]
correlation $vacf autoupdate start
set vel2 [observable new particle_velocities type [list $ghost2]]
set vacf2 [correlation new obs1 $vel2 corr_operation scalar_product tau_max [expr 2000.*$time_step] dt [expr $time_step]]
correlation $vacf2 autoupdate start

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
	puts $trajfile " [expr ([setmd time] - $equil_time)]   $x1   $y1   $z1"
	flush $trajfile
    }
    
    flush $xyzfile
}

###################main simulation steps - production of results################################################
set rho 0.0
set rgsum 0.0
for {set i 1} { $i <= $total_cycles} { incr i} {
    
    integrate $integ_steps
    
    analyze append
    set rg [analyze rg 0 $n_cell $nions_per_ring]
    set rg2 [expr [lindex $rg 0]]
    set rgsum [expr $rgsum+$rg2]
    set lxnow [expr [lindex [setmd box_l] 0]]
    set lynow [expr [lindex [setmd box_l] 1]]
    set rho [expr $rho + ($n_cell/($lxnow*$lynow))]
    puts $obsfile " [expr ([setmd time] - $equil_time)] [expr $n_cell/($lxnow*$lynow)] $rg2"
    flush $obsfile
    
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
	    puts $trajfile " [expr ([setmd time] - $equil_time)]   $x1   $y1   $z1"
	    flush $trajfile
	}
	
	flush $xyzfile
    }
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
puts $rhofile " $press [expr $rho/$total_cycles] "
flush $rhofile
puts $rg_file " $press [expr $rgsum/$total_cycles] "
flush $rg_file
#################end of program######################
exit


