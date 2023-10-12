###############################################################################################
########            DPD simulation of nanorods at a water-oil interface             ###########
########              Inspired in the model by....            ###########
########                    Scrit by J. R. Bordin - Apr 2023                        ########### 
########                         Powered by ESPRESSO                                ###########
###############################################################################################

#set simulation parameters
set n_pol 200 ;#número de polímeros rígidos
set n_mon 14  ;#número de monomeros no nanobastao
set tot_mon [expr int($n_pol*$n_mon)] ;# número total de monomeros no sistema
set f 0.0 ;#fração de monômeros hidrofílicos no nanobastão
set box_lx 40.
set box_ly 40.
set box_lz 80.
setmd box_l $box_lx $box_ly $box_lz;#caixa de simulação anisotrópica pra facilitar a separação
set volume [expr $box_lx*$box_ly*$box_lz];#volume total
set density 3.0 ;#densidade dos fluidos
set n_fluid [expr int(ceil($density*$volume))]
set fraction 0.5 ;#fração de água
set nwater [expr int(ceil($n_fluid*$fraction))]
set noil [expr int(ceil($n_fluid - $nwater))]
set pi [expr 4.*atan (1.0)]
cellsystem domain_decomposition
set time_step 0.02
setmd time_step $time_step;#MD time step
setmd skin 0.4;#for thermostat
setmd periodic 1 1 1 
t_random seed [ pid ]

set equil_steps 250000
set integ_steps 1000
set total_steps 250
set run_time [expr $integ_steps*$total_steps*$time_step]

set temp 1.0
set gamma 4.5 ;#friction coefficient for DPD
set cut 1.0
thermostat dpd $temp $gamma $cut
set id "npol-$n_pol-n_mon-$n_mon-temp-$temp-f-$f"

#the species needs a number to "define" his type. here I set a number to each species.
set water 1
set oil   2
set phobic   3
set philic   4
set wall 5

inter $water $water hat 131.5 $cut
inter $oil $oil hat 131.5 $cut
inter $water $oil hat 171.43 $cut

inter $water $philic hat 0.0 $cut
inter $water $phobic hat 168.9 $cut
inter $oil $philic hat 142.07 $cut
inter $oil $phobic hat 132.01 $cut

inter $philic $philic hat 131.5 $cut
inter $phobic $phobic hat 131.5 $cut
inter $philic $phobic hat 140.83 $cut

#set the bonded interactions
set harmonic 10
set kharmonic 100
set r0harmonic 0.9
inter $harmonic harmonic $kharmonic $r0harmonic

set bond_angle 40
set k_angle 500
inter $bond_angle angle_cosine $k_angle 

set n_part 0

for { set j 0 } { $j < $n_pol } {incr j 1} {
    set molid $j; #puts "Copolymer $molid"
    set posx [expr $box_lx*[t_random]]
    set posy [expr $box_ly*[t_random]]
    set posz [expr 1.0+($box_lz-2.0)*[t_random]]
    if { $f > 0.0} {
	part $n_part pos $posx $posy $posz q 0 type $philic mass 1.0
    }
    if { $f == 0.0} {
	part $n_part pos $posx $posy $posz q 0 type $phobic mass 1.0
    }
    part $n_part mol $molid
    incr n_part 1
    if { $f > 0.0} {
	for { set k 2} { $k <= [expr int($n_mon*$f)] } {incr k 1 } {
	    set posx1 [expr $posx+ $r0harmonic*($k-1)]
	    part $n_part  pos $posx1 $posy $posz q 0 type $philic mass 1.0
	    part $n_part bond $harmonic [expr $n_part -1]
	    part $n_part mol $molid
	    if {$k > 2} {
		part [expr $n_part -1] bond $bond_angle [expr $n_part -2] $n_part
	    }
	    incr n_part 1
	}
	for { set k  [expr int($n_mon*$f)+1] } { $k <= $n_mon } {incr k 1 } {
	    set posx1 [expr $posx+ ( $r0harmonic*($k-1))]
	    part $n_part  pos $posx1 $posy $posz q 0 type $phobic mass 1.0
	    part $n_part bond $harmonic [expr $n_part -1]
	    part $n_part mol $molid
	    if {$k > 2} {
		part [expr $n_part -1] bond $bond_angle [expr $n_part -2] $n_part
	    }
	    incr n_part 1
	}
    }
    if { $f == 0.0} {
	for { set k 2 } { $k <= $n_mon } {incr k 1 } {
	    set posx1 [expr $posx+ ( $r0harmonic*($k-1))]
	    part $n_part  pos $posx1 $posy $posz q 0 type $phobic mass 1.0
	    part $n_part bond $harmonic [expr $n_part -1]
	    part $n_part mol $molid
	    if {$k > 2} {
		part [expr $n_part -1] bond $bond_angle [expr $n_part -2] $n_part
	    }
	    incr n_part 1
	}
    }
}
analyze set chains 0 $n_pol $n_mon 

for { set j 0 } { $j < $nwater } {incr j 1} {
    set posx [expr $box_lx*[t_random]]
    set posy [expr $box_ly*[t_random]]
    set posz [expr 2.0+0.4*$box_lz*[t_random]]
    part $n_part pos $posx $posy $posz q 0 type $water mass 1.0
    incr n_part 1
}
for { set j 0 } { $j < $noil } {incr j 1} {
    set posx [expr $box_lx*[t_random]]
    set posy [expr $box_ly*[t_random]]
    set posz [expr $box_lz/2.0+($box_lz/2.-2.)*[t_random]]
    part $n_part pos $posx $posy $posz q 0 type $oil mass 1.0
    incr n_part 1
}


############# open the files ###########################################33
set xyzfile [open "snapshot-$id.xyz" "w"]
set agreg_file [open "agreg-$id.dat" "w"]
puts $agreg_file "\#--Time -- N_aggreg -- <n_mol> ---- err(+-) - N_max(aggreg) - N_min(aggreg)"
set obs_file [open "obs-$id.dat" "w"]
puts $obs_file "\#DPD simulation of nanorods at water-oil interface"
puts $obs_file "\#Powered by Espresso! - tcl script by J. R. Bordin - Dez 2023"
puts $obs_file "\#Universidade Federal de Pelotas -- UFPel"
puts $obs_file "\#---------------------------------------------------------------------"
puts $obs_file "\#--Time ----- (pressure ---- U/N ------- k/N )"


##################start of MD warm up steps##############################################
set deg_free 3
for {set i 0} { $i < 1 } { incr i} {
    integrate $equil_steps
}
######################### End of warm up steps#############################################

##############################################################################################################

#Start the correlations

###correlation for density profiles of each particle type in the system

set nbins [expr int(100*$box_lz)]

if {$nwater > 0} {
    set dens_water [observable new density_profile type [list $water] minx 0 maxx $box_lx xbins 1 miny 0 maxy $box_ly ybins 1 minz 0 maxz $box_lz zbins $nbins]
    set corr_water [correlation new obs1 $dens_water tau_max $run_time dt 1 compress1 discard1 corr_operation componentwise_product]
    correlation $corr_water autoupdate start
}
if {$noil > 0} {
    set dens_oil [observable new density_profile type [list $oil] minx 0 maxx $box_lx xbins 1 miny 0 maxy $box_ly ybins 1 minz 0 maxz $box_lz zbins $nbins]
    set corr_oil [correlation new obs1 $dens_oil tau_max $run_time dt 1 compress1 discard1 corr_operation componentwise_product]
    correlation $corr_oil autoupdate start
}

if {$n_pol > 0 && $f > 0.0} {
    set dens_philic [observable new density_profile type [list $philic] minx 0 maxx $box_lx xbins 1 miny 0 maxy $box_ly ybins 1 minz 0 maxz $box_lz zbins $nbins]
    set corr_philic [correlation new obs1 $dens_philic tau_max $run_time dt 1 compress1 discard1 corr_operation componentwise_product]
    correlation $corr_philic autoupdate start
}

if {$n_pol > 0 && $f < 1.0} {
    set dens_phobic [observable new density_profile type [list $phobic] minx 0 maxx $box_lx xbins 1 miny 0 maxy $box_ly ybins 1 minz 0 maxz $box_lz zbins $nbins]
    set corr_phobic [correlation new obs1 $dens_phobic tau_max $run_time dt 1 compress1 discard1 corr_operation componentwise_product]
    correlation $corr_phobic autoupdate start
}

# stress tensor for surface tension
set stress [observable new stress_tensor]



###################main simulation steps - production of results################################################

#set j 0
for {set i 0} { $i < $total_steps} { incr i} {
    integrate $integ_steps
    set energy [analyze energy];#evaluates the energy
    set u [lindex [lindex $energy 0] 1]
    set k [lindex [lindex $energy 1] 1]
    #set pxy [analyze pressure total]
    set pxx [expr [lindex [observable $stress print] 0]] ;#xx component of the stress tensor 
    set pyy [expr [lindex [observable $stress print] 4]] ;#yy component of the stress tensor
    set pzz [expr [lindex [observable $stress print] 8]] ;#yy component of the stress tensor
    set tension [expr 0.5*($pzz-(0.5*($pxx+$pyy)))*$box_lz] ;#surface tension
    puts $obs_file [format "%.3e %.5e %.5e %.5e " [setmd time]  $tension  [expr $u/$n_fluid]  [expr $k/$n_fluid]   ]
    flush $obs_file

    set agreg [analyze aggregation 1.0 0 [expr $n_pol-1] ]
    puts $agreg_file [format "%.3e %.5e %.5e %.5e %.5e %.5e" [setmd time] [lindex $agreg 9] [lindex $agreg 5]  [lindex $agreg 7] [lindex $agreg 1] [lindex $agreg 3] ]
    flush $agreg_file

    #write a .xyz file for use in jmol or vmd
    puts $xyzfile "    $tot_mon  "
    puts $xyzfile "    Atoms   "
    
    for {set nion 0} { $nion < $tot_mon } { incr nion } {
	
	set x1  [expr [lindex [part $nion print folded_position] 0]]
	set y1  [expr [lindex [part $nion print folded_position] 1]]
	set z1  [expr [lindex [part $nion print folded_position] 2]]
	if { [lindex [part $nion ] 6] == $water } {
	    puts $xyzfile " N   $x1   $y1   $z1"
	}
	if { [lindex [part $nion ] 6] == $oil } {
	    puts $xyzfile " C   $x1   $y1   $z1"
	}
	if { [lindex [part $nion ] 6] == $philic } {
	    puts $xyzfile " O   $x1   $y1   $z1"
	}
	if { [lindex [part $nion ] 6] == $phobic } {
	    puts $xyzfile " S   $x1   $y1   $z1"
	}		
	flush $xyzfile
    }
}
################### end of simulation steps - production of results################################################
close $obs_file
close $xyzfile

### finalize all correlations
for {set i 0} {$i < [correlation n_corr] } {incr i} {
    correlation $i finalize;
}
### print the density profiles files
set number_of_datapoints [llength [correlation $corr_water print]]

set zlist ""
for {set j 0} {$j < $nbins} {incr j} {
    lappend zlist [expr (($box_lz/$nbins)*$j)]
}

if {$nwater > 0} {
    set water_den [correlation $corr_water print average1]
    set out1 [open "dens_water-$id.dat" "w"]
    foreach z $zlist c $water_den { puts $out1  [format "%.3e %.3e" $z $c] }
    close $out1
}

if {$noil > 0} {
    set oil_den [correlation $corr_oil print average1]
    set out1 [open "dens_oil-$id.dat" "w"]
    foreach z $zlist c $oil_den { puts $out1 [format "%.3e %.3e" $z $c] }
    close $out1
}

if {$n_pol > 0 && $f > 0.0} {
    set philic_den [correlation $corr_philic print average1]
    set out1 [open "dens_philic-$id.dat" "w"]
    foreach z $zlist c $philic_den { puts $out1 [format "%.3e %.3e" $z $c]}
    close $out1
}

if {$n_pol > 0 && $f < 1.0} {
    set phobic_den [correlation $corr_phobic print average1]
    set out1 [open "dens_phobic-$id.dat" "w"]
    foreach z $zlist c $phobic_den { puts $out1 [format "%.3e %.3e" $z $c]}
    close $out1
}


#################end of program######################
    exit


