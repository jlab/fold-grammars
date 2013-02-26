Install instantbibi:
	ssh zingana
	rm ~/.bibi*
	cd /var/tmp
	mkdir StefansBIBI
	cd StefansBIBI
	hg clone ssh://hg/bibiadm/bibiserv2/main/tools/instantbibi
	cd instantbibi
	ant instant.dev -Dportbase=25000 -Dbase.dir=/var/tmp/StefansBIBI/
		
export JAVA_HOME=/media/Daten/BIBI/jdk1.7.0_10/
export PATH=/media/Daten/BIBI/jdk1.7.0_10/bin:$PATH
localhost:8080/wizard

in /media/Daten/BIBI/base:
	ant -Dxml=/home/sjanssen/Desktop/fold-grammars/Bibiserv/pkiss.bs2 -Dwithout_ws=true -Dwithout_sswap=true -Dwithout_moby=true -Dwithout_vb=true

in /tmp/pkiss_ID:
	ant deploy
	
ssh c -L 8080:129.70.160.177:8080

user: testadmin
pwd: simplepassword
