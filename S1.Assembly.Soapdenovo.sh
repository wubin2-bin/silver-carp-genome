./SOAPdenovo-63mer pregraph -s lib.cfg -K 33 -p 40 -d 1 -o bighead > pregraph.log
./SOAPdenovo-63mer contig -g bighead -R -p 40 > contig.log
./SOAPdenovo-63mer map -s lib.cfg -p 30 -g bighead -k 33 > map.log
./SOAPdenovo-63mer scaff -g bighead -F -p 40 > scaff.log
