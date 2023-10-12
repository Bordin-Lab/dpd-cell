###############################################################################################
########                    MD simulation of a 2D cell model                        ###########
########                    Scrit by J. R. Bordin - Apr 2022                        ########### 
########                         Powered by ESPRESSO                                ###########
###############################################################################################
set systemTime [clock seconds]
#puts "The initial time is: [clock format $systemTime -format %H:%M:%S]"
#puts "The initial date is: [clock format $systemTime -format %D]"


set pi [expr 4.*atan (1.0) ];#set value of pi
t_random seed [ pid ]
set sigma 1.0 ;#distance unit, in micrometers
set tau0 1.0

#Set the simulation box
set box_lx 80.0  ;# in units of sigma
set box_ly 80.0 ;# in units of sigma 
set box_lz 80.0
setmd box_l $box_lx $box_ly $box_lz;#defines the initial simulation box
setmd periodic 1 1 0 ;#PBC in x and y directions
set Atot [expr $box_lx*$box_ly]


#Set the cell radius and number of cells
set radius 3.0 ;#cell radius in units of sigma = 10 microns
set n_cell 100
set Acell [expr $pi*($radius+0.5)*($radius+0.5)] ;#are of one cell
set Acell [expr $n_cell*$Acell] ;#area occupied by the cells
set kstr 10.0 ;#coil strength for the bond between the central ghost bead and the ring. It controls the cell deformation

cellsystem domain_decomposition


#set the temperature range
set temp 1.00
set p_ini 0.31; #initial pressure
set p_fin 0.50; #final pressure
set dp 0.01 ;#temperature step
set Np [expr (($p_fin-$p_ini)/$dp)+1]
set gamma_0 1.0 ;#same dumping parameter from Langevin/DPD thermostat
set gamma_v 0.001 ;#coupling parameter for the volume control
set piston_mass 0.0001 ;#mass of the piston that will control the pressure
set cut 1.00
setmd skin 0.4;#for thermostat
  

#set the integration data
set time_step 0.01
setmd time_step $time_step;#MD time step
#warm up
set integ_steps_warm 50000
#equil steps
set integ_steps_equil 0000
#production cycles
set integ_steps 500
set total_cycles 100
set run_time [expr $integ_steps*$total_cycles*$time_step]



#each species needs a number to "define" his type. here I set a number to each species.
set cell 0
set ghost 1
set insert 2



######################defining the bonded interactions


#bond interaction for two neighbour cell beads
set sig 1.0 ;#size of each bead
set ring 21
set lambda 1.0
inter $ring harmonic 100.0 $lambda  [expr 3.*$lambda]

#bond with the central bead
set rigid 22
inter $rigid harmonic $kstr $radius [expr 3.*$radius]

#################################### defining the non-bonded interactions

set WCA_cut [expr pow(2,1./6.)] ;#repulsive WCA interaction
set atract_cut [expr 2.25] ;#atractive 12-6 LJ interaction (cylinder-tail)
#repulsive interactions
set eps 1.0 ;#interaction parameter

inter $cell $cell lennard-jones 1.0 1.0 $WCA_cut auto

#placing the ring-like cell
set nions_per_ring  [expr int(((2.*$pi*$radius)/$lambda))]
set teta  [expr 2.0*$pi/($nions_per_ring-1)]
set innerangle 41   ;#inner angular interaction to keep the ring structure
set k_angle 100
inter $innerangle angle $k_angle $teta
#puts "Total number of cell particles is [expr $n_cell*$nions_per_ring] in $n_cell cells"
set ns 0
for {set cell_counter 0} { $cell_counter < $n_cell } {incr cell_counter 1} {
    puts "$cell_counter"
    #insert the i-th ring
    if {$cell_counter == 0} {
	set posy0 [expr $box_ly*[t_random]]
	set posx0 [expr $box_lx*[t_random]]
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
	    set overlap_test [analyze mindist $ghost $ghost]
#	    puts "$overlap_test"
	    if {$overlap_test > [expr 2.*$radius+0.5]} {
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
	part $ns pos $posx $posy $posz q 0 type $cell mass 1.0 fix 0 0 1
	part $ns mol $cell_counter	
	part $ns bond $rigid $gid ;#bond the bead with the inner monomer
	if { $j > 1 } {
	    part $ns bond $ring [expr $ns-1]
	    part $ns exclude [expr $ns-1]
	    part $gid bond $innerangle [expr $ns-1] $ns;#bond two adjacent beads with the inner monomer by a angular constraint
	}
	set ns [expr $ns+1]
    }
}

for {set j 0 } { $j < [setmd n_part] } {incr j [expr $nions_per_ring]} {
#    puts " $nions_per_ring [expr $j+1] [expr $nions_per_ring+$j-1]"
    part [expr $j+1] bond $ring [expr $nions_per_ring+$j-1]
    part $j bond $innerangle [expr $j+1] [expr $j+1]
    part [expr $j+1] exclude [expr $nions_per_ring+$j-1]
}

set n_cell_part [setmd n_part]
set n_part [setmd n_part]

#set systemTime [clock seconds]
#puts "End of insertion step: [clock format $systemTime -format %H:%M:%S]"

analyze set chains 0 $n_cell $nions_per_ring

set Nmolecules [setmd n_part]
set id "k-$kstr"
set rhofile [open "rho-$id.dat" "w"] ;#open the rho x p file
set rg_file [open "rg-$id.dat" "w"];#open the rg x p file

######################### Initial warm steps #############################################
thermostat set npt_isotropic $temp $gamma_0 $gamma_v ;# set the NPT simulation
integrate set npt_isotropic $p_ini $piston_mass 1 1 0
inter forcecap 200.0
integrate $integ_steps_warm
inter forcecap 0
###################### End of Warm steps #############################################


for { set Pcount 0 } { $Pcount < $Np} {incr Pcount 1} {
    #define a random number for each temperature
    t_random seed [ pid ]
    #define the pressure for this simulation
    set press  [format %1.4g [expr $p_ini+$Pcount*$dp]]
    #sets the thermosthat
   # thermostat langevin $temp $gamma_0
    thermostat set npt_isotropic $temp $gamma_0 $gamma_v ;# set the NPT simulation
    integrate set npt_isotropic $press $piston_mass 1 1 0
    #id for the files
    set id "k-$kstr-p-$press"
    set xyzfile [open "snap-$id.xyz" "w"] ;#open the xyz file
    set trajfile [open "traj-$id.xyz" "w"] ;#open the CM trajectory file
    set velfile [open "vel-$id.dat" "w"] ;#open the CM velocity file
    puts $velfile [expr ($n_cell*$total_cycles)]
    set obsfile [open "obs-$id.dat" "w"] ;#open the observable file

    ######################### Equilibration steps #############################################
    integrate $integ_steps_equil
    set equil_time [setmd time]
    #print the initial positions to the trajectory file
    for {set nion 0} { $nion < [setmd n_part] } { incr nion } {
	set x1  [expr $sigma*[lindex [part $nion print folded_position] 0]]
	set y1  [expr $sigma*[lindex [part $nion print folded_position] 1]]
	set z1  [expr $sigma*[lindex [part $nion print folded_position] 2]]
	if { [lindex [part $nion ] 6] == $ghost } {
	    puts $trajfile " [expr ([setmd time] - $equil_time)]   $x1   $y1   $z1"
	    flush $trajfile
	}
    } 
    ###################### End of Equilibration steps #############################################
        
    puts $obsfile "\#Langevin Dynamics simulation of 2D cells in a solvent"
    puts $obsfile "\#Powered by Espresso! - tcl script by J. R. Bordin - April 22"
    puts $obsfile "\#---------------------------------------------------------------------"
    puts $obsfile "\# $n_cell cells"
    puts $obsfile "\# Total number of particles is [setmd n_part] "
    puts $obsfile "\# Temperature = $temp and damping parameter = $gamma_0 and pressure = $press"
    
   # set pos [observable new particle_positions type [list $ghost]]
   # set msd [correlation new obs1 $pos corr_operation square_distance_componentwise tau_lin 16 tau_max $run_time  dt [expr $time_step]  compress1 discard1]
    
   # correlation $msd autoupdate start
    
    set vel [observable new particle_velocities type [list $ghost]]
    set vacf [correlation new obs1 $vel corr_operation scalar_product tau_max [expr 2000.*$time_step] dt [expr $time_step]]
    correlation $vacf autoupdate start
    
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
	    if { [lindex [part $nion ] 6] == $ghost } {
		puts $xyzfile " H  $x1   $y1   $z1"
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
    ################################## end of main steps
    
 #   correlation $msd write_to_file "msd-per-cell-$id.dat"
    correlation $vacf write_to_file "vacf-cor-$id.dat"
    
    set systemTime [clock seconds]
   # puts "Production end time is: [clock format $systemTime -format %H:%M:%S]"
   # puts "Production end date is: [clock format $systemTime -format %D]"
    
    
    #close $vacf_file
    close $xyzfile
    close $velfile
    close $obsfile
    puts $rhofile " $press [expr $rho/$total_cycles] "
    flush $rhofile
    puts $rg_file " $press [expr $rgsum/$total_cycles] "
    flush $rg_file
}
#################end of program######################
exit


