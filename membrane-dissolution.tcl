###############################################################################################
########                 MD simulation of a bilayer membrane with surfactant        ###########
########    Lipid model inspired in Cooke, Kremer and Deserno, PRE 011506 (2005)    ###########
########                    Scrit by J. R. Bordin - Dez 2019                        ########### 
########                         Powered by ESPRESSO                                ###########
###############################################################################################

set do_config "no" ;#write the configurations for further analysis
set dpd "no" ;#set yes to use DPD as a thermostat

set eps_ll 1.0     ;#lipid-lipid tail interaction in units of kT, with T = 300K
set eps_ss 1.0     ;#surfactant-surfactat tail interaction in units of kT, with T = 300K
set eps_ls 1.5     ;#surfactant-lipid tail interaction in units of kT, with T = 300K
set eps_ratio [expr $eps_ls] ;#ratio for the competition
set eps0 1.0      ;#for WCA repulsive interactions in units of kT, with T = 300K

set n_beadslipid 3;#number of monomers in each lipid molecule
set n_beadssurf 3;#number of monomers in each surfactant molecule
set csurf 10.0 ;#concentration of surfactant, in mmol/L
set cfluid 50.0 ;#fluid concentration in mmol/L
set id "ratio-$eps_ratio-csurf-$csurf-cfluid-$cfluid" ;#id for the files
set pi [expr 4.*atan (1.0) ];#set value of pi

#set simulation parameters
set temp 1.0; #temperature in reduces units, equal to 300 K
set gamma 1.0;#damping parameter
set press 0.024 ;#pressure = 1,01 x 10^5 Pa, using eps* = kT, with T = 300K, 
set gamma_v 0.0001 ;#coupling parameter for the volume control
set piston_mass 0.5 ;#mass of the piston that will control the pressure

#set the initial box size - it will vary, since we are fixing the pressure in the equilibration
set box_lx 20.0
set box_ly 20.0
set box_lz 20.0
setmd box_l $box_lx $box_ly $box_lz;#defines the initial simulation box
set volume [expr $box_lx* $box_ly* $box_lz];#initial volume of the simulation box

#set the integration data
set time_step 0.0025
setmd time_step $time_step;#MD time step
#capped interaction
set integ_steps_cap 50
set cap_cycles 10
set cap 50.0
#warm up in NPT ensemble
set integ_steps_warmNPT 0
#warm up in NVT ensemble
set integ_steps_warmNVT 0
#production cycles
set integ_steps 500
set total_cycles 500
set run_time [expr $integ_steps*$total_cycles*$time_step]


cellsystem domain_decomposition
setmd skin 0.4;#for thermostat
setmd periodic 1 1 0 ;#pbc in all directions
t_random seed [ pid ] ;#start the random number

#each species needs a number to "define" his type. here I set a number to each species.
set headlipid 0
set taillipid 1
set headsurf 2
set tailsurf 3
set fluid1 4
set fluid2 5
set wall 6

#defined the size (diameter) of each species in reduces units, with sigma* = 1 nm
set sig_headlipid 0.95
set sig_taillipid 1.0
set sig_headsurf 1.00
set sig_tailsurf 1.0
set sig_fluid 1.0

###########################INTERACTIONS

#define the properties of the bonded interaction between lipids/surfactants
set fene 10
inter $fene fene 30.0 1.5

#bend interaction for the lipid/surfactant rigidity
set bend 20
inter $bend harmonic 10.0 4.0


#set the non-bonded interactions interactions
set WCA_cut [expr pow(2,1./6.)] ;#repulsive WCA interaction
set atract_cut [expr 2.25] ;#atractive 12-6 LJ interaction (cylinder-taillipid)
set wc 1.5 ;#for the tail-tail interaction


inter $headlipid $headlipid lennard-jones $eps0 $sig_headlipid [expr $WCA_cut*$sig_headlipid] auto
inter $headsurf $headsurf lennard-jones $eps0 $sig_headlipid [expr $WCA_cut*$sig_headlipid] auto

inter $taillipid $taillipid lj-cos2 $eps_ll $sig_taillipid 0. [expr $wc*$sig_taillipid] 
inter $tailsurf $tailsurf lj-cos2  $eps_ss $sig_tailsurf 0. [expr $wc*$sig_tailsurf] 
inter $tailsurf $taillipid lj-cos2 $eps_ls [expr 0.5*($sig_taillipid+$sig_headsurf)] 0. [expr 2.*$wc]
inter $headsurf $headlipid lennard-jones $eps0 $sig_headlipid [expr $WCA_cut*$sig_headlipid] auto
inter $taillipid $headlipid lennard-jones $eps0 [expr 0.5*($sig_taillipid+$sig_headlipid)] [expr 0.5*($sig_taillipid+$sig_headsurf)*$WCA_cut] auto
inter $tailsurf $headlipid lennard-jones $eps0 [expr 0.5*($sig_tailsurf+$sig_headlipid)] [expr 0.5*($sig_taillipid+$sig_headsurf)*$WCA_cut] auto
inter $tailsurf $headsurf lennard-jones $eps0 [expr 0.5*($sig_tailsurf+$sig_headsurf)] [expr 0.5*($sig_taillipid+$sig_headsurf)*$WCA_cut] auto
inter $taillipid $headsurf lennard-jones $eps0 [expr 0.5*($sig_taillipid+$sig_headsurf)] [expr 0.5*($sig_taillipid+$sig_headsurf)*$WCA_cut] auto
inter $fluid1 $fluid1 lennard-jones $eps0 $sig_fluid [expr $WCA_cut*$sig_fluid] auto
inter $headlipid $fluid1 lennard-jones $eps0 [expr 0.5*($sig_fluid+$sig_headlipid)] [expr 0.5*($sig_fluid+$sig_headlipid)*$WCA_cut] auto
inter $headsurf $fluid1 lennard-jones $eps0 [expr 0.5*($sig_fluid+$sig_headsurf)] [expr 0.5*($sig_fluid+$sig_headsurf)*$WCA_cut] auto
inter $taillipid $fluid1 lennard-jones $eps0 [expr 0.5*($sig_fluid+$sig_taillipid)] [expr 0.5*($sig_fluid+$sig_taillipid)*$WCA_cut] auto
inter $tailsurf $fluid1 lennard-jones $eps0 [expr 0.5*($sig_fluid+$sig_tailsurf)] [expr 0.5*($sig_fluid+$sig_taillipid)*$WCA_cut] auto
inter $fluid2 $fluid2 lennard-jones $eps0 $sig_fluid [expr $WCA_cut*$sig_fluid] auto
inter $headlipid $fluid2 lennard-jones $eps0 [expr 0.5*($sig_fluid+$sig_headlipid)] [expr 0.5*($sig_fluid+$sig_headlipid)*$WCA_cut] auto
inter $headsurf $fluid2 lennard-jones $eps0 [expr 0.5*($sig_fluid+$sig_headsurf)] [expr 0.5*($sig_fluid+$sig_headsurf)*$WCA_cut] auto
inter $taillipid $fluid2 lennard-jones $eps0 [expr 0.5*($sig_fluid+$sig_taillipid)] [expr 0.5*($sig_fluid+$sig_taillipid)*$WCA_cut] auto
inter $tailsurf $fluid2 lennard-jones $eps0 [expr 0.5*($sig_fluid+$sig_tailsurf)] [expr 0.5*($sig_fluid+$sig_taillipid)*$WCA_cut] auto

for { set i $headlipid } { $i <= $wall } { incr i 1 } {
    inter $i $wall lennard-jones $eps0 1.0 $WCA_cut auto
}


######################################################## INSERT THE PARTICLES

#placing the lipids in a square mesh with size grid
set n_lip 0 ;#number of lipid molecules
set grid 1.05
set xgrid [expr int(($box_lx)/$grid)]
set ygrid [expr int(($box_ly)/$grid)]
#first include the upper layer

for { set xx 0 } { $xx < $xgrid } {incr xx 1} {
    for { set yy 0 } { $yy < $ygrid } {incr yy 1} {
	set posx [expr ($xx*$grid)]
	set posy [expr ($yy*$grid)]
	set posz [expr $box_lz/2.+$n_beadslipid]
	set i [expr [setmd n_part]]
	set j [expr $i+1]
	set k [expr $i+2]
	
	part $i pos $posx $posy $posz q 0 type $headlipid mass 1.0
	part $j pos $posx $posy [expr $posz-1.] q 0 type $taillipid mass 1.0	    
	part $k pos $posx $posy [expr $posz-2.] q 0 type $taillipid mass 1.0
	
	part $i bond $fene $j
	part $j bond $fene $k
	part $i bond $bend $k
	
	set n_lip [expr $n_lip+1]
    }
}
#now include the down layer
for { set xx 0 } { $xx < $xgrid } {incr xx 1} {
    for { set yy 0 } { $yy < $ygrid } {incr yy 1} {
	set posx [expr ($xx*$grid)]
	set posy [expr ($yy*$grid)]
	set posz [expr $box_lz/2.-($n_beadslipid)]
	set i [expr [setmd n_part]]
	set j [expr $i+1]
	set k [expr $i+2]
	part $i pos $posx $posy $posz q 0 type $headlipid mass 1.0
	part $j pos $posx $posy [expr $posz+1.] q 0 type $taillipid mass 1.0
	part $k pos $posx $posy [expr $posz+2.] q 0 type $taillipid mass 1.0
	
	part $i bond $fene $j
	part $j bond $fene $k
	part $i bond $bend $k
	
	set n_lip [expr $n_lip+1]   
    }
}

#DPD has long range hydrodynamics effects, Langevin local hydrodynamics effects
if {$dpd == "yes"} {
    thermostat dpd $temp $gamma 5.0 ;#uses DPD as thermostat
}  else {
    thermostat langevin $temp $gamma ;#used Langevin thermostat
}
#set the barostat for the equilibration steps
thermostat set npt_isotropic $temp $gamma $gamma_v ;# set the NPT simulation
integrate set npt_isotropic $press $piston_mass 1 1 0

################## membrane equilibration ##############################################
set deg_free 3
#warm up steps - with capsulated force to relax the system
#performs integ_steps_cap each cycle - 100 x 100 steps should be enough
for { set i 0 } { $i < $cap_cycles } { incr i } {
#    puts "$i $cap"
    inter forcecap $cap
    integrate $integ_steps_cap
    set cap [expr $cap*1.1]
}
inter forcecap 0
#warm up steps to reach equilibrium using NPT ensemble - at leat 2x10^5 steps, only one cycle
#can be set zero for a test run
integrate $integ_steps_warmNPT

################## End of NPT equilibraton steps#############################################

#size of the box after equilibration
set lxnow [expr [lindex [setmd box_l] 0]]
set lynow [expr [lindex [setmd box_l] 1]]
set lznow [expr [lindex [setmd box_l] 2]]
#thermostat off ;#turns off the thermostat and barostat
#turn on the thermostat for NVT simulations
if {$dpd == "yes"} {
    thermostat dpd $temp $gamma 5.0 ;#uses DPD as thermostat
}  else {
    thermostat langevin $temp $gamma ;#used Langevin thermostat
}
integrate set nvt

#setmd periodic 1 1 0 ;#pbc in xy plane
set posz 0.000
constraint plane cell -10 -10 $posz type $wall
constraint plane cell -10 -10 $lznow type $wall
  
#puts "$lxnow $lynow $lznow"
set fsurf [expr ($csurf*6.023/10000)] ;#in reduced units
set Nsurf [expr int(($lxnow*$lynow*(($lznow/2.)-2.0))*$fsurf)] ;#number of surfactants molecules
#puts "$Nsurf"
#insert random lipids
for { set ii 0} { $ii < $Nsurf } {incr ii 1} { 
    set posx [expr $lxnow*[t_random]]
    set posy [expr $lynow*[t_random]]
    set posz [expr 1+(($lznow/2.)-2.0)*[t_random] ]
    set i [expr [setmd n_part]]
    set j [expr $i+1]
    set k [expr $i+2]
    part $i pos $posx $posy $posz q 0 type $headsurf mass 1.0
    part $j pos [expr $posx+1.0] $posy $posz q 0 type $tailsurf mass 1.0
    part $k pos [expr $posx+2.0] $posy $posz q 0 type $tailsurf mass 1.0
    
    part $i bond $fene $j
    part $j bond $fene $k
    part $i bond $bend $k 
      
}
set ffluid [expr ($cfluid*6.023/10000)] ;#in reduced units
set Nfluid [expr int(($lxnow*$lynow*(($lznow/2.)-2.0))*$ffluid)] ;#number of fluid molecules

for { set ii 0} { $ii < $Nfluid } {incr ii 1} { 
    set posx [expr $lxnow*[t_random]]
    set posy [expr $lynow*[t_random]]
    if { $ii < [expr $Nfluid/2] } {
	set posz [expr (($lznow/2.)+2.0)+($lznow/2.-3.0)*[t_random] ]
	set i [expr [setmd n_part]]
	part $i pos $posx $posy $posz q 0 type $fluid1 mass 1.0
    } else {
	set posz [expr (($lznow/2.)-3.0)*[t_random] ]
	set i [expr [setmd n_part]]
	part $i pos $posx $posy $posz q 0 type $fluid2 mass 1.0
    }
}


################## membrane equilibration ##############################################
#set deg_free 3
#performs integ_steps_cap each cycle
for { set i 0 } { $i < $cap_cycles } { incr i } {
#    puts "$i $cap"
    inter forcecap $cap
    integrate $integ_steps_cap
    set cap [expr $cap*1.1]
}
inter forcecap 0
#warm up steps to reach equilibrium using NVT ensemble
#can be set zero for a test run
#integrate $integ_steps_warmNVT

################## End of NVT equilibraton steps#############################################

###### start of the main part of the simulation
##### here we evaluate the properties and save the positions for further analysis

set obs_file [open "obs-$id.dat" "w"] ;#open the observables file
puts $obs_file "\#MD simulation of a lipid bilayer"
puts $obs_file "\#Powered by Espresso! - tcl script by J. R. Bordin - December 2019"
puts $obs_file "\#---------------------------------------------------------------------"
puts $obs_file "\#All physical quantities are in reduced units"
puts $obs_file "\#Lipid model inspired in Cooke, Kremer and Deserno, PRE 011506 (2005) "
puts $obs_file "\#$n_lip lipid molecules with 3 monomers each"
puts $obs_file "\#Pressure controled with Andersen barostat"
puts $obs_file "\#pressure = $press, volume damping parameter = $gamma_v and piston mass = $piston_mass "
if {$dpd == "yes"} {
    puts $obs_file "\#Temperature controled with DPD thermostat"
} else {
    puts $obs_file "\#Temperature controled with Langevin thermostat"
}
puts $obs_file "\# Temperature = $temp and damping parameter = $gamma"
puts $obs_file "\#Time ## #Pressure_XY ## Kinectic ## Potential "

#open the snapshot file
set xyzfile [open "snap-$id.xyz" "w"]


set stress [observable new stress_tensor];# to evaluate the stress tensor. here I check the pxy parte of the tensor, that should be equal to the pressure


###################main simulation steps - production of results################################################
set nbins [expr int(2.5*$lznow)]
for {set i 0} { $i < $total_cycles} { incr i} {
    
    ###correlation for density profiles of each particle type in the system
    
    
#    set dens_headlipid [observable new density_profile type [list $headlipid] minx 0 maxx $lxnow xbins 1 miny 0 maxy $lynow ybins 1 minz 0 maxz $lznow zbins $nbins]
#    set corr_headlipid [correlation new obs1 $dens_headlipid tau_max $run_time dt 1 compress1 discard1 corr_operation componentwise_product]
#    correlation $corr_headlipid autoupdate start
    
    
#    set dens_taillipid [observable new density_profile type [list $taillipid] minx 0 maxx $lxnow xbins 1 miny 0 maxy $lynow ybins 1 minz 0 maxz $lznow zbins $nbins]
#    set corr_taillipid [correlation new obs1 $dens_taillipid tau_max $run_time dt 1 compress1 discard1 corr_operation componentwise_product]
#    correlation $corr_taillipid autoupdate start

#    set dens_fluid1 [observable new density_profile type [list $fluid1] minx 0 maxx $lxnow xbins 1 miny 0 maxy $lynow ybins 1 minz 0 maxz $lznow zbins $nbins]
#    set corr_fluid1 [correlation new obs1 $dens_fluid1 tau_max $run_time dt 1 compress1 discard1 corr_operation componentwise_product]
#    correlation $corr_fluid1 autoupdate start
#    set dens_fluid2 [observable new density_profile type [list $fluid2] minx 0 maxx $lxnow xbins 1 miny 0 maxy $lynow ybins 1 minz 0 maxz $lznow zbins $nbins]
#    set corr_fluid2 [correlation new obs1 $dens_fluid2 tau_max $run_time dt 1 compress1 discard1 corr_operation componentwise_product]
#    correlation $corr_fluid2 autoupdate start
    
    
#    if {$Nsurf > 0} {
#	set dens_headsurf [observable new density_profile type [list $headsurf] minx 0 maxx $lxnow xbins 1 miny 0 maxy $lynow ybins 1 minz 0 maxz $lznow zbins $nbins]
#	set corr_headsurf [correlation new obs1 $dens_headsurf tau_max $run_time dt 1 compress1 discard1 corr_operation componentwise_product]
#	correlation $corr_headsurf autoupdate start
	
#	set dens_tailsurf [observable new density_profile type [list $tailsurf] minx 0 maxx $lxnow xbins 1 miny 0 maxy $lynow ybins 1 minz 0 maxz $lznow zbins $nbins]
#	set corr_tailsurf [correlation new obs1 $dens_tailsurf tau_max $run_time dt 1 compress1 discard1 corr_operation componentwise_product]
#	correlation $corr_tailsurf autoupdate start
 #   }
    
    set energy [analyze energy];#evaluates the energy
    set u [lindex [lindex $energy 0] 1] ;#potential energy
    set k [lindex [lindex $energy 1] 1] ;#kinectic energy
    set pxx [expr [lindex [observable $stress print] 0]] ;#xx component of the stress tensor 
    set pyy [expr [lindex [observable $stress print] 4]] ;#yy component of the stress tensor
    set pxy [expr 0.5*($pxx+$pyy)] ;#xy component of the stress tensor
    puts $obs_file [format "%.3e %.5e %.5e %.5e " [setmd time]  $pxy  [expr $u/[setmd n_part]]  [expr $k/[setmd n_part]]   ] ;#write the quantities to the obs file
    flush $obs_file
    

    puts $xyzfile "   [setmd n_part] "
    puts $xyzfile "    Atoms   "
    
    for {set nion 0} { $nion < [setmd n_part] } { incr nion } {
	set x1  [expr [lindex [part $nion print folded_position] 0]]
	set y1  [expr [lindex [part $nion print folded_position] 1]]
	set z1  [expr [lindex [part $nion print folded_position] 2]]
	if { [lindex [part $nion ] 6] == $taillipid } {
	    puts $xyzfile " C   $x1   $y1   $z1"
	} 
	
	if { [lindex [part $nion ] 6] == $headlipid} {
	    puts $xyzfile " O  $x1   $y1   $z1"
	}
	if { [lindex [part $nion ] 6] == $tailsurf } {
	    puts $xyzfile " N   $x1   $y1   $z1"
	} 
	
	if { [lindex [part $nion ] 6] == $headsurf } {
	    puts $xyzfile " P  $x1   $y1   $z1"
	}
	if { [lindex [part $nion ] 6] == $fluid1 } {
	    puts $xyzfile " F  $x1   $y1   $z1"
	}
	if { [lindex [part $nion ] 6] == $fluid2 } {
	    puts $xyzfile " B  $x1   $y1   $z1"
	}
    }
    flush $xyzfile

    
    if { $do_config == yes} {
	#write the configurationfiles for future analyses
	set f [open "config-$id.$i" "w"]
	blockfile $f write tclvariable {box_l }
	blockfile $f write variable box_l
	blockfile $f write particles {id pos type}
	close $f
    }
    integrate $integ_steps
    ### finalize all correlations
#    for {set i 0} {$i < [correlation n_corr] } {incr i} {
#	correlation $i finalize;
#    }
    
    ### print the density profiles files
 #   set number_of_datapoints [llength [correlation $corr_headlipid print]]
    
 #   set zlist ""
 #   for {set j 0} {$j < $nbins} {incr j} {
#	lappend zlist [expr (($box_lz/$nbins)*$j)]
#    }
    
#    set fluid_den [correlation $corr_headlipid print average1]
#    set out2 [open "dens_lipid_head-$id.$i.dat" "w"]
#    foreach z $zlist c $fluid_den { puts $out2 "$z [expr $c/6.023*10000]" }
    
#    set fluid_den [correlation $corr_taillipid print average1]
#    set out3 [open "dens_lipid_tail-$id.$i.dat" "w"]
#    foreach z $zlist c $fluid_den { puts $out3 "$z [expr $c/6.023*10000]" }

#    set fluid_den [correlation $corr_fluid1 print average1]
#    set out4 [open "dens_fluid1-$id.$i.dat" "w"]
#    foreach z $zlist c $fluid_den { puts $out4 "$z [expr $c/6.023*10000]" }
#    set fluid_den [correlation $corr_fluid2 print average1]
#    set out5 [open "dens_fluid2-$id.$i.dat" "w"]
#    foreach z $zlist c $fluid_den { puts $out4 "$z [expr $c/6.023*10000]" }
    
#    if {$Nsurf > 0} {
#	set fluid_den [correlation $corr_headsurf print average1]
#	set out6 [open "dens_surf_head$id.$i.dat" "w"]
#	foreach z $zlist c $fluid_den { puts $out5 "$z [expr $c/6.023*10000]" }

#	set fluid_den [correlation $corr_tailsurf print average1]
#	set out7 [open "dens_surf_tail-$id.$i.dat" "w"]
#	foreach z $zlist c $fluid_den { puts $out6 "$z [expr $c/6.023*10000]" }
# }
#    close $out2
#    close $out3
#    close $out4
#    close $out5
#    if {$Nsurf > 0} {
#	close $out7
#	close $out6
#    }
}
############################################################### end of simulation


close $obs_file
close $xyzfile




#################end of program######################
exit


