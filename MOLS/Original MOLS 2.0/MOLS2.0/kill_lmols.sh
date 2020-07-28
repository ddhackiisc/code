#shell script to kill lmols job
for x in `cat lmols.pid`
	do
	kill -15 $x
	done
