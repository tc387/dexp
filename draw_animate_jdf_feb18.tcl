
set nf 100
for {set i 0} {$i < $nf} {incr i} {
    echo $i
    set ind [ expr { 1000 + $i } ]
    eval topo readlammpsdata [format "test_cnf_out%04d.lammps bond" $ind ]

    echo $ind
#pbc wrap -all -compound resid -center origin


axes location Off
display projection Orthographic
menu graphics off
menu graphics on
#mol color Name
#mol representation Lines 1.000000
mol selection top
mol material Opaque
    mol addrep top
    echo AAAAAAAA
mol modstyle 0 top CPK 1.000000 0.300000 12.000000 12.000000
scale by 1.200000
mol modselect 1 top type 1
mol modstyle 1 top VDW 27.000000 72.000000
mol modstyle 0 top CPK 1.000000 0.300000 52.000000 42.000000
scale by 1.200000
    scale by 1.200000
    echo BBBBBBBB
mol modstyle 0 top CPK 1.500000 0.700000 52.000000 42.000000
mol modmaterial 0 top BrushedMetal
mol modmaterial 1 top BrushedMetal
menu color off
menu color on
color Display Background white
    display depthcue off
    menu color off

    
    display update
    render TachyonInternal [format "animate.%04d.tga" $i]
    mol delete all
}
eval "exec ffmpeg  -f image2 -i animate.%04d.tga -vb 20M animate.mpg"
