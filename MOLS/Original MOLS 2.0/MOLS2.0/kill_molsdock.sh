#shell script to kill molsdock job
for x in `cat molsdock.pid`
	do
	kill -15 $x
	done
