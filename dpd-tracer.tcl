###############################################################################################
########        MD simulation of a fluid based in a tabulated potential file        ###########
########                    Scrit by J. R. Bordin - Apr 2020                        ########### 
########                         Powered by ESPRESSO                                ###########
###############################################################################################

#set basic parameters
set sig0 1.0
set eps0  1.0

#set number of particles and molecules
set Ncol 800 ;#real numer of colloids is half of this
set Ntracer 200
set Npart [expr $Ncol+$Ntracer]
#set the temperature range
set temp_ini 0.05; #initial temperature
set temp_fin 0.40; #final temperature
set dtemp 0.05 ;#temperature step
set Ntemp [expr (($temp_fin-$temp_ini)/$dtemp)+1]

#barostat and thermostat data
set press 6.0 ;#pressure for the barostat
set gamma_0 1.0 ;#same dumping parameter from Langevin/DPD thermostat
set gamma_v 0.0002 ;#coupling parameter for the volume control
set piston_mass 0.001 ;#mass of the piston that will control the pressure
setmd skin 0.4;#for thermostat
 
set box_lx 30.0 ;# size in the x-direction
set box_ly $box_lx ;# size in the y-direction
set box_lz $box_lx ;# size in the z-direction 
setmd box_l $box_lx $box_lx $box_lx;#define the simulation box
cellsystem domain_decomposition ;#neighbor list
setmd periodic 1 1 1 ;#pbc in all directions


#set the integration data
set time_step 0.01
setmd time_step $time_step;#MD time step
#capped interaction
set integ_steps_cap 0
set cap 5.0
set cap_max 500.0
#warm up
set integ_steps_warm 25000
#production cycles
set integ_steps 250
set total_cycles 200
set avg_cycles  1

#defines pi
set pi [expr 4.*atan (1.0)]

#the species needs a number to "define" his type. here I set a number to each species.
set hc 0
set ss 1
set tracer 2
							
#################################### defining the non-bonded interactions
set rc1 1.12
set rc2 2.0
set fmax1 50.0
set fmax2 10.0
inter $hc $hc hat $fmax1 $rc1
inter $ss $ss hat $fmax2 $rc2
inter $tracer $tracer hat $fmax1 $rc1
inter $tracer $hc hat $fmax1 $rc1
inter $tracer $ss hat $fmax2 $rc2

#################################### defining the bonded interactions

set shake 10
#inter $shake rigid_bond 0.0001 0.01 0.01
inter $shake harmonic 100. 0.0
########################## insert the particles
for { set i 0 } { $i < $Ncol } {incr i 2} {
    set posx [expr $box_lx*[t_random]]
    set posy [expr $box_ly*[t_random]]
    set posz [expr $box_lz*[t_random]]
    part $i pos $posx $posy $posz type $hc
    part [expr $i+1] pos $posx $posy $posz type $ss
    part $i bond $shake [expr $i+1] 
}

for { set i $Ncol } {$i < $Npart} {incr i 1 } {
    set posx [expr $box_lx*[t_random]]
    set posy [expr $box_ly*[t_random]]
    set posz [expr $box_lz*[t_random]]
    part $i pos $posx $posy $posz type $tracer
}

##########start the isothermal loop

for { set Tcount 0 } { $Tcount < $Ntemp} {incr Tcount 1} {
    #define a random number for each temperature
    t_random seed [ pid ]
    #define the temperature for this simulation
    set temp  [format %1.4g [expr $temp_fin-$Tcount*$dtemp]]
    #sets the thermosthat
    #thermostat dpd $temp $gamma_0 3.0  
    #thermostat langevin $temp $gamma_0  
    thermostat inter_dpd $temp
    inter $hc $hc inter_dpd $gamma_0 $rc1 1 1 1 1
    inter $ss $ss inter_dpd $gamma_0 $rc2 1 1 1 1
    inter $tracer $tracer inter_dpd $gamma_0 $rc1 1 1 1 1
    inter $tracer $hc inter_dpd $gamma_0 $rc1 1 1 1 1
    inter $tracer $ss inter_dpd $gamma_0 $rc2 1 1 1 1
    
    set id "P-$press-T-$temp"
    set xyzfile [open "snapshot-$id.xyz" "w"]
    
    if { $Tcount == 0 } {
	##################start of MD warm up steps##############################################
	set deg_free 3
	#warm up steps - with capsulated force to relax the system in the NVT ensemble
	while {$cap < $cap_max } {
	    inter forcecap $cap
	    integrate $integ_steps_cap
	    set cap [expr $cap*5.0]
	    set act_min_dist [analyze mindist]
	}
	inter forcecap 0.0
    }
    #########set the barostat 
    thermostat set npt_isotropic $temp $gamma_0 $gamma_v ;# set the NPT simulation
    integrate set npt_isotropic $press $piston_mass 
    
    #warm up steps in NPT ensemble
    integrate $integ_steps_warm
    set equil_time [setmd time]
    
    ######################### End of warm up steps#############################################
    set obs_file [open "obs-$id.dat" "w"] ;#open the observables file
    puts $obs_file "\#MD simulation of tabulated potential"
    puts $obs_file "\#Powered by Espresso! - tcl script by J. R. Bordin - Dec 2022"
    puts $obs_file "\#---------------------------------------------------------------------"
    puts $obs_file "\#All physical quantities are in reduced units"
    puts $obs_file "\#Pressure controled with Andersen barostat"
    puts $obs_file "\#pressure = $press, volume damping parameter = $gamma_v and piston mass = $piston_mass "
    puts $obs_file "\#Temperature controled with Langevin thermostat"
    puts $obs_file "\# Temperature = $temp and damping parameter = $gamma_0"
    puts $obs_file "\#Time ## Density  ## Volume ## Potential Energy ## Entalpy ## Pressure ## Temperature"
    set lxnow [expr [lindex [setmd box_l] 0]]
    set lynow [expr [lindex [setmd box_l] 1]]
    set lznow [expr [lindex [setmd box_l] 2]]
    set energy [analyze energy];#evaluates the energy
    set press_now  [ setmd npt_p_inst] ;#evaluate the energy using the virial theorem
    set u [expr [lindex [lindex $energy 0] 1]] ;#potential energy
    set k [expr  [lindex [lindex $energy 1] 1]] ;#kinectic energy
    set volume [expr ($lxnow*$lynow*$lznow)]
    set rho [expr ($Npart/$volume)]
    set vmol [expr ($volume/$Npart)]
    puts $obs_file [format "%.4e %3.4f %3.4f %3.4f %3.4f %3.4f %3.4f" [expr [setmd time] - $equil_time] $rho  $vmol [expr $u/$Npart] [expr (($u+ $press_now*$volume)/$Npart)] $press $temp ] ;#write the quantities to the obs file
    
    set control_file [open "control-$id.dat" "w"] ;#open the control file
    puts $control_file "\#Time ## Pressure  ## Kinetic Energy "
    
    #set stress [observable new stress_tensor];# to evaluate the stress tensor. here I check the pxy parte of the tensor, that should be equal to the pressure 
    set u 0.0
    set k 0.0
    set vmol 0.0
    set rho 0.0
    set press_now 0.0
    ###################main simulation steps - production of results################################################
    for {set i 1} { $i <= $total_cycles} { incr i} {
	
	set lxnow [expr [lindex [setmd box_l] 0]]
	set lynow [expr [lindex [setmd box_l] 1]]
	set lznow [expr [lindex [setmd box_l] 2]]
	set energy [analyze energy];#evaluates the energy
	set press_now [expr $press_now + [ setmd npt_p_inst] ];#evaluate the energy using the virial theorem
	set u [expr $u + [lindex [lindex $energy 0] 1]] ;#potential energy
	set k [expr $k + [lindex [lindex $energy 1] 1]] ;#kinectic energy
	set volume [expr ($lxnow*$lynow*$lznow)]
	set rho [expr $rho + ($Npart/$volume)]
	set vmol [expr $vmol + ($volume/$Npart)]
	
	if { $i % $avg_cycles == 0 } {
	    set u [expr $u/$avg_cycles]
	    set k [expr $k/$avg_cycles]
	    set rho [expr $rho/$avg_cycles]
	    set vmol [expr $vmol/$avg_cycles]
	    set press_now [expr $press_now/$avg_cycles]
	    puts $obs_file [format "%.4e %3.4f %3.4f %3.4f %3.4f %3.4f %3.4f" [expr [setmd time] - $equil_time] $rho  $vmol [expr $u/$Npart] [expr (($u+ $press_now*$volume)/$Npart)] $press $temp ] ;#write the quantities to the obs file
	    
	    puts $control_file [format "%.3e %.5e %.5e "  [expr [setmd time] - $equil_time] $press_now [expr $k/$Npart] ]
	    set u 0.0
	    set k 0.0
	    set vmol 0.0
	    set rho 0.0
	    set press_now 0.0
	    
	    puts $xyzfile "   [setmd n_part] "
	    puts $xyzfile "    Atoms   "
	    
	    for {set nion 0} { $nion < [setmd n_part] } { incr nion } {
		set x1  [expr [lindex [part $nion print folded_position] 0]]
		set y1  [expr [lindex [part $nion print folded_position] 1]]
		set z1  [expr [lindex [part $nion print folded_position] 2]]
		if { [lindex [part $nion ] 6] == $hc } {
		    puts $xyzfile " C   $x1   $y1   $z1"
		} 
		if { [lindex [part $nion ] 6] == $ss } {
		    puts $xyzfile " S   $x1   $y1   $z1"
		} 
		if { [lindex [part $nion ] 6] == $tracer } {
		    puts $xyzfile " T   $x1   $y1   $z1"
		} 
		
	    }
	    flush $xyzfile
	}
    	integrate $integ_steps
    }
    close $xyzfile   
    close $obs_file
    close $control_file 
    thermostat off
}  

#################end of program######################
exit


