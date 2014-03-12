reset
set term pdf enhanced color size 10,5
set output "figure1_partA.pdf"
#set key 850, -0.4, 0 right top vertical Right noreverse enhanced autotitles nobox
set key bottom right
set label 1 "out of memory (8 GB)" at 360, 0.1
set label 2 "compute time exhausted" at 400, -0.7
set label 4 "compute time exhausted" at 500, -0.45
set label 3 "A)" at graph 0.02, 0.95 font "Times-Bold,20"
set arrow 1 from 395, -0.1, 0 to 395, 0.07 nohead linetype -1
set arrow 2 from 396, -0.7, 0 to 396, -0.2 nohead linetype -1
set arrow 3 from 570, -0.4, 0 to 570, -0.1 nohead linetype -1
set xtics 100
set y2tics nomirror
set ytics nomirror
set my2tics 10
set format y2 "10^{%L}"
set title "Necessary number of TDMs (T = 0.1)" font "Times,14"
set xlabel "sequence length: |s| = n" 
set ylabel "accumulated shape probability - (1 - T)" 
set y2label "|L(s)| = no. TDMs" 
set yrange [ -0.900000 : 0.300000 ] noreverse nowriteback
set y2range [ 1.00000 : 100000.0 ] noreverse nowriteback
set xrange [0:1000]
set log y2
set xtics font "Verdana,16" 
set ytics font "Verdana,16" 
set y2tics font "Verdana,16" 
plot 0 axis x1y1 with lines lt 0 lw 1 notitle, \
"figure1.data" index 0 using 1:2 with lines smooth bezier axis x1y2 lt 3 lw 2 title "|L_{sampling}(0.1, 10000, s)", \
"figure1.data" index 0 using 1:8 every ::0::113 with lines smooth bezier axis x1y2 lt 1 lw 2 title "|L_{energy}(0.1, s)|", \
"figure1.data" index 0 using 1:($3-0.9) with lines smooth bezier axis x1y1 lt 3 lw 1 notitle, \
"figure1.data" index 0 using 1:5 every ::0::113 with lines smooth bezier axis x1y2 lt 2 lw 2 title "|L_{oracle}(0.1, s)|", \
"figure1.data" index 0 using 1:($6-0.9) every ::0::113 with lines smooth bezier axis x1y1 lt 2 lw 1 notitle, \
"figure1.data" index 1 using 1:2 every ::0::24 with lines smooth bezier axis x1y2 lt 5 lw 2 title "|F(s)| lowProbFilter=0", \
"figure1.data" index 2 using 1:2 with lines smooth bezier axis x1y2 lt 4 lw 2 title "|F(s)| lowProbFilter=0.000001", \
"figure1_real.data" index 0 using 1:2 axis x1y2 with lines smooth bezier lt 7 lw 2 title "|L_{sampling}(0.1, 10000, s)| (real data)", \
"figure1_real.data" index 0 using 1:8 axis x1y2 with lines smooth bezier lt 7 lw 1 title "|L_{energy}(0.1, s)| (real data)", \
"figure1_real.data" index 0 using 1:17 axis x1y2 with lines smooth bezier lt 9 lw 2 title "|L_{sampling}(0.5, 10000, s)| (real data)", \
"figure1_real.data" index 0 using 1:20 axis x1y2 with lines smooth bezier lt 9 lw 1 title "|L_{energy}(0.5, s)| (real data)"



reset
set term pdf enhanced color size 5,5
set output "figure1_partB.pdf"
reset
set format y "10^{%L}"
set mytics 10
set key right bottom
unset label
set label 1 "out of memory (8 GB)" at 290, 100000
set label 2 "out of memory (8 GB)" at 15, 50000
set label 3 "compute time exhausted" at 500, 3e+05
set label 4 "B)" at graph 0.02, 0.95 font "Times-Bold,20"
unset arrow
set arrow 1 from 395, 20000 to 395, 60000 nohead
set arrow 2 from 125, 2000 to 125, 30000 nohead
set arrow 3 from 570, 2000 to 570, 2e+05 nohead
set logscale y
set title "Runtime comparison (T = 0.1)" font "Times,14"
set xlabel "sequence length: |s| = n" 
set ylabel "runtime in seconds"
set yrange [ 0.0100000 : 1.00000e+07 ]
set xtics font "Verdana,16" 
set ytics font "Verdana,16" 
plot "figure1.data" index 0 using 1:4 with lines smooth bezier lt 3 lw 2 title "L_{sampling}(0.1, 10000, s)", \
"figure1.data" index 0 using 1:10 every ::0::113 with lines smooth bezier lt 1 lw 2 title "L_{energy}(0.1, s)", \
"figure1.data" index 0 using 1:7 every ::0::113 with lines smooth bezier lt 2 lw 2 title "L_{oracle}(0.1, s)", \
"figure1.data" index 1 using 1:3 every ::0::24 with lines smooth bezier lt 5 lw 2 title "|F(s)| lowProbFilter=0", \
"figure1.data" index 2 using 1:3 with lines smooth bezier lt 4 lw 2 title "|F(s)| lowProbFilter=0.000001", \
"figure1.data" index 0 using 1:41 with lines smooth bezier lt 8 lw 2 title "pure sampling 10000", \
"figure1_real.data" index 0 using 1:4 with lines smooth bezier lt 7 lw 2 title "L_{sampling}(0.1, 10000, s) (real data)", \
"figure1_real.data" index 0 using 1:10 with lines smooth bezier lt 7 lw 1 title "L_{energy}(0.1, s) (real data)", \
"figure1_real.data" index 0 using 1:19 with lines smooth bezier lt 9 lw 2 title "L_{sampling}(0.5, 10000, s) (real data)", \
"figure1_real.data" index 0 using 1:22 with lines smooth bezier lt 7 lw 1 title "L_{energy}(0.5, s) (real data)"



reset
set term pdf enhanced color size 5,5
set output "figure1_partC.pdf"
reset
set format y "10^{%L}"
set format y2 "10^{%L}"
set mytics 10
set key right bottom
set label 1 "out of memory (8 GB)" at 250, 100000
set label 2 "compute time exhausted" at 590, 3e+06
set label 3 "C)" at graph 0.02, 0.95 left norotate font "Times-Bold,20"
set arrow 1 from 395, 20000 to 395, 60000 nohead linetype -1
set arrow 2 from 705, 5000 to 705, 2e+06 nohead linetype -1
set logscale y 10
set logscale y2 10
set title "Runtimes for different thresholds: 0.1 {/Symbol \243} T {/Symbol \243} 0.6" font "Times,14"
set xlabel "sequence length: |s| = n" 
set y2label "runtime in seconds"
set yrange [ 0.0100000 : 1.00000e+07 ]
set y2range [ 0.0100000 : 1.00000e+07 ]
set xtics font "Verdana,16" 
set ytics font "Verdana,16" 
set y2tics font "Verdana,16" 
plot \
"./figure1.data" index 0 using 1:($4+0.1) with lines smooth bezier lt 3 lw 2 title "|L_{sampling}(0.1 {/Symbol \243} T {/Symbol \243} 0.6, 10000, s)|", \
"./figure1.data" index 0 using 1:($13+0.1) with lines smooth bezier lt 3 lw 2 notitle, \
"./figure1.data" index 0 using 1:($19+0.1) with lines smooth bezier lt 3 lw 2 notitle, \
"./figure1.data" index 0 using 1:($25+0.1) with lines smooth bezier lt 3 lw 2 notitle, \
"./figure1.data" index 0 using 1:($31+0.1) with lines smooth bezier lt 3 lw 2 notitle, \
"./figure1.data" index 0 using 1:($37+0.1) with lines smooth bezier lt 3 lw 2 notitle, \
"./figure1.data" index 0 using 1:($10+0.1) every ::0::113 with lines smooth bezier lt 1 lw 2 title "|L_{energy}(0.1 {/Symbol \243} T {/Symbol \243} 0.6, s)|", \
"./figure1.data" index 0 using 1:($16+0.1) every ::0::199 with lines smooth bezier lt 1 lw 2 notitle, \
"./figure1.data" index 0 using 1:($22+0.1) every ::0::166 with lines smooth bezier lt 1 lw 2 notitle, \
"./figure1.data" index 0 using 1:($28+0.1) every ::0::158 with lines smooth bezier lt 1 lw 2 notitle, \
"./figure1.data" index 0 using 1:($34+0.1) every ::0::153 with lines smooth bezier lt 1 lw 2 notitle, \
"./figure1.data" index 0 using 1:($40+0.1) every ::0::113 with lines smooth bezier lt 1 lw 2 notitle, \
"./figure1.data" index 2 using 1:($3) with lines smooth bezier lt 4 lw 2 title "RNAshapes -p -F 0.000001"


